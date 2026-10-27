"""σ-net v2: the v1 dilated-TCN trunk (models/tcn_noise_field.py) with a configurable
input channel set, optional FiLM conditioning on global invariants, and a configurable head family.

Heads (any subset):
  pw    per-point log σ²(s)                          [B, L]       (c) per-point + variance pooling
  rho   per-stroke lag-1 residual correlation ρ       [B]          pairs with pw -> AR(1)/OU covariance
  nig   global Evidential NIG on y = log10 rms       4 x [B]      (γ, ν, α, β), Student-t ν=2α marginal
  pyr   supervised dyadic region pyramid K=1..16     {K: [B,2,K]} (a) (mu = log10 region rms, log_var)
  cell  σ² per Haar-DWT cell (scale x position)       [B, NC]      (b) wavelet-pooled, engine-native
  intv  continuous interval head (a,b) -> (mu, lv)   [B, Q, 2]    (d)

The trunk is unchanged from v1 except: no user/calibration encoder, FiLM only when invariants are fed, and the sequence input
width is a constructor argument. All head output layers are zero-initialised (identity warm start),
like v1. Everything is static-shape at L (Core ML friendly: conv1d, GroupNorm, GELU, matmul).
"""
import math
import torch
import torch.nn as nn
import torch.nn.functional as F


class FiLM(nn.Module):
    def __init__(self, ch, cond):
        super().__init__()
        self.g = nn.Linear(cond, ch)
        self.b = nn.Linear(cond, ch)
        nn.init.zeros_(self.g.weight); nn.init.ones_(self.g.bias)
        nn.init.zeros_(self.b.weight); nn.init.zeros_(self.b.bias)

    def forward(self, x, z):
        return self.g(z).unsqueeze(2) * x + self.b(z).unsqueeze(2)


class Block(nn.Module):
    """v1 TemporalBlock: two dilated k=3 'same' convs, GroupNorm(8), optional FiLM, GELU, residual."""
    def __init__(self, ch, dilation, cond, dropout):
        super().__init__()
        pad = dilation
        self.c1 = nn.Conv1d(ch, ch, 3, dilation=dilation, padding=pad)
        self.c2 = nn.Conv1d(ch, ch, 3, dilation=dilation, padding=pad)
        self.n1 = nn.GroupNorm(8, ch)
        self.n2 = nn.GroupNorm(8, ch)
        self.f1 = FiLM(ch, cond) if cond else None
        self.f2 = FiLM(ch, cond) if cond else None
        self.drop = nn.Dropout(dropout)

    def forward(self, x, z):
        o = self.n1(self.c1(x))
        if self.f1 is not None:
            o = self.f1(o, z)
        o = self.drop(F.gelu(o))
        o = self.n2(self.c2(o))
        if self.f2 is not None:
            o = self.f2(o, z)
        o = self.drop(F.gelu(o))
        return o + x


class SigmaV2(nn.Module):
    PYR_SCALES = (1, 2, 4, 8, 16)

    def __init__(self, in_ch, n_inv, heads, hidden=128, layers=4, dropout=0.1, L=48,
                 cell_pos=None, cell_neg=None, cell_geo=None, cell_unc=False):
        super().__init__()
        self.heads = tuple(heads)
        self.L = L
        self.n_inv = n_inv
        H = hidden
        self.cond = None
        if n_inv > 0:
            self.cond = nn.Sequential(nn.Linear(n_inv, H), nn.LayerNorm(H), nn.GELU(),
                                      nn.Linear(H, H), nn.LayerNorm(H), nn.GELU())
        self.proj = nn.Conv1d(in_ch, H, 1)
        self.blocks = nn.ModuleList([Block(H, 2 ** i, H if n_inv > 0 else 0, dropout) for i in range(layers)])
        gdim = 2 * H if n_inv > 0 else H

        def zero_last(m):
            nn.init.zeros_(m.weight); nn.init.zeros_(m.bias)
            return m

        if "pw" in self.heads:
            self.pw = nn.Sequential(nn.Conv1d(H, H // 2, 1), nn.GELU(), zero_last(nn.Conv1d(H // 2, 1, 1)))
        if "rho" in self.heads:
            self.rho = nn.Sequential(nn.Linear(gdim, H), nn.GELU(), zero_last(nn.Linear(H, 1)))
        if "nig" in self.heads:
            self.nig = nn.Sequential(nn.Linear(gdim, H), nn.GELU(), zero_last(nn.Linear(H, 4)))
        if "pyr" in self.heads:
            self.pyr = nn.ModuleDict({str(k): zero_last(nn.Conv1d(H, 2, 1)) for k in self.PYR_SCALES})
        if "cell" in self.heads:
            assert cell_pos is not None
            self.register_buffer("cell_pos", cell_pos)
            self.register_buffer("cell_neg", cell_neg)
            self.register_buffer("cell_geo", cell_geo)
            G = cell_geo.shape[1]
            self.cell_unc = cell_unc
            self.cell = nn.Sequential(nn.Linear(2 * H + G + gdim, H), nn.GELU(),
                                      nn.Linear(H, H // 2), nn.GELU(), zero_last(nn.Linear(H // 2, 2 if cell_unc else 1)))
        if "intv" in self.heads:
            e = H // 4
            self.iemb = nn.Sequential(nn.Linear(2, e), nn.GELU(), nn.Linear(e, e), nn.GELU())
            self.intv = nn.Sequential(nn.Linear(H + e + gdim, H), nn.GELU(), zero_last(nn.Linear(H, 2)))

    def trunk(self, seq, inv):
        z = self.cond(inv) if self.cond is not None else None
        h = self.proj(seq)
        for b in self.blocks:
            h = b(h, z)
        g = h.mean(2)
        if z is not None:
            g = torch.cat([g, z], 1)
        return h, g

    def forward(self, seq, inv=None, intervals=None):
        h, g = self.trunk(seq, inv)
        out = {}
        if "pw" in self.heads:
            out["logv"] = self.pw(h)[:, 0]
        if "rho" in self.heads:
            out["rho"] = 0.999 * torch.tanh(self.rho(g)[:, 0])
        if "nig" in self.heads:
            r = self.nig(g)
            out["gamma"] = r[:, 0]
            out["nu"] = F.softplus(r[:, 1]) + 1e-4
            out["alpha"] = F.softplus(r[:, 2]) + 1.0 + 1e-4
            out["beta"] = F.softplus(r[:, 3]) + 1e-4
        if "pyr" in self.heads:
            out["pyr"] = {k: self.pyr[str(k)](F.adaptive_avg_pool1d(h, k)) for k in self.PYR_SCALES}
        if "cell" in self.heads:
            hp = torch.matmul(h, self.cell_pos.t()).transpose(1, 2)
            hn = torch.matmul(h, self.cell_neg.t()).transpose(1, 2)
            B, NC = hp.shape[0], hp.shape[1]
            geo = self.cell_geo.unsqueeze(0).expand(B, -1, -1)
            gg = g.unsqueeze(1).expand(-1, NC, -1)
            co = self.cell(torch.cat([hp, hn, geo, gg], 2))
            out["cell_logE"] = co[..., 0]
            if self.cell_unc:
                out["cell_logv"] = co[..., 1]
        if "intv" in self.heads and intervals is not None:
            pos = torch.linspace(0.0, 1.0, self.L, device=h.device)
            a, b = intervals[..., 0], intervals[..., 1]
            m = ((pos[None, None] >= a[..., None]) & (pos[None, None] <= b[..., None])).float()
            m = m / m.sum(-1, keepdim=True).clamp(min=1.0)
            pooled = torch.einsum("bcl,bql->bqc", h, m)
            emb = self.iemb(torch.stack([(a + b) * 0.5, torch.log((b - a).clamp(min=1e-3))], -1))
            Q = intervals.shape[1]
            out["intv"] = self.intv(torch.cat([pooled, emb, g.unsqueeze(1).expand(-1, Q, -1)], -1))
        return out

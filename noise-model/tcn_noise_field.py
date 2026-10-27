"""Dual-stream 1D-TCN with FiLM for heteroscedastic noise estimation.

Architecture:
1. Spatial Sequence Stream: 4-Layer Exponentially Dilated 1D-TCN (d = 1, 2, 4, 8) with kernel k=3 and residual skips.
2. Global Invariant Stream: MLP projecting 52 cross-disciplinary invariants z into conditioning space.
3. FiLM (Feature-wise Linear Modulation): Injects global physical/topological invariants into intermediate convolutional tokens.
4. Evidential Deep Learning (EDL) NIG Global Head: Predicts Normal-Inverse-Gamma (gamma, nu, alpha, beta) for closed-form aleatoric/epistemic uncertainty disentanglement.
5. Heavy-Tailed Pointwise Student-t Head: Predicts localized spatial noise envelope sigma_i(s).
6. Calibration Context Layer: Ingests 3-stroke onboarding summary to condition base amplitude.
"""
import math
import torch
import torch.nn as nn
import torch.nn.functional as F


class FiLMLayer(nn.Module):
    """Feature-wise Linear Modulation: h_out = gamma(z) * h + beta(z)"""
    def __init__(self, in_features, cond_features):
        super().__init__()
        self.gamma_fc = nn.Linear(cond_features, in_features)
        self.beta_fc = nn.Linear(cond_features, in_features)
        
        nn.init.ones_(self.gamma_fc.weight.data)
        nn.init.zeros_(self.gamma_fc.bias.data)
        nn.init.zeros_(self.beta_fc.weight.data)
        nn.init.zeros_(self.beta_fc.bias.data)

    def forward(self, x, z):
        gamma = self.gamma_fc(z).unsqueeze(2)
        beta = self.beta_fc(z).unsqueeze(2)
        return gamma * x + beta


class TemporalBlock(nn.Module):
    """Dilated Causal 1D Convolutional Residual Block with FiLM Modulation."""
    def __init__(self, in_channels, out_channels, kernel_size, dilation, cond_dim, dropout=0.1):
        super().__init__()
        self.padding = (kernel_size - 1) * dilation // 2
        
        self.conv1 = nn.Conv1d(in_channels, out_channels, kernel_size, dilation=dilation, padding=self.padding)
        self.norm1 = nn.GroupNorm(8, out_channels)
        self.film1 = FiLMLayer(out_channels, cond_dim)
        self.act1 = nn.GELU()
        self.drop1 = nn.Dropout(dropout)
        
        self.conv2 = nn.Conv1d(out_channels, out_channels, kernel_size, dilation=dilation, padding=self.padding)
        self.norm2 = nn.GroupNorm(8, out_channels)
        self.film2 = FiLMLayer(out_channels, cond_dim)
        self.act2 = nn.GELU()
        self.drop2 = nn.Dropout(dropout)
        
        self.downsample = nn.Conv1d(in_channels, out_channels, 1) if in_channels != out_channels else None

    def forward(self, x, z):
        res = x if self.downsample is None else self.downsample(x)
        
        out = self.conv1(x)
        if out.shape[2] != x.shape[2]:
            out = out[:, :, :x.shape[2]]
        out = self.norm1(out)
        out = self.film1(out, z)
        out = self.act1(out)
        out = self.drop1(out)
        
        out = self.conv2(out)
        if out.shape[2] != x.shape[2]:
            out = out[:, :, :x.shape[2]]
        out = self.norm2(out)
        out = self.film2(out, z)
        out = self.act2(out)
        out = self.drop2(out)
        
        return out + res


class UserContextEncoder(nn.Module):
    """Attention-pooled context encoder over a set of calibration strokes.
    
    Ingests a variable set of K calibration/previous strokes: [Batch, K, 52]
    and outputs a continuous 32-dimensional Personal Motor Fingerprint e_user in R^32.
    """
    def __init__(self, in_dim=52, embed_dim=32, hidden_dim=64):
        super().__init__()
        self.stroke_mlp = nn.Sequential(
            nn.Linear(in_dim, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.GELU(),
            nn.Linear(hidden_dim, embed_dim),
            nn.LayerNorm(embed_dim)
        )
        self.attn_score = nn.Linear(embed_dim, 1)
        self.default_prior = nn.Parameter(torch.zeros(1, embed_dim))

    def forward(self, calib_strokes=None):
        if calib_strokes is None or (torch.is_tensor(calib_strokes) and calib_strokes.shape[1] == 0):
            return self.default_prior
            
        batch_size, k_strokes, _ = calib_strokes.shape
        h = self.stroke_mlp(calib_strokes)
        
        scores = self.attn_score(h) / math.sqrt(h.shape[-1])
        weights = F.softmax(scores, dim=1)
        e_user = torch.sum(weights * h, dim=1)
        
        return e_user


class DualStreamTCNNoiseField(nn.Module):
    """Dual-stream heteroscedastic noise estimator with a user context encoder."""
    def __init__(self, in_seq_channels=4, in_invariant_dim=52, user_embed_dim=32, hidden_dim=128, num_layers=4):
        super().__init__()
        
        self.user_encoder = UserContextEncoder(in_dim=in_invariant_dim, embed_dim=user_embed_dim)
        
        total_cond_dim = in_invariant_dim + user_embed_dim
        self.cond_mlp = nn.Sequential(
            nn.Linear(total_cond_dim, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.GELU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.GELU()
        )
        
        self.seq_proj = nn.Conv1d(in_seq_channels, hidden_dim, kernel_size=1)
        
        self.tcn_blocks = nn.ModuleList([
            TemporalBlock(
                in_channels=hidden_dim,
                out_channels=hidden_dim,
                kernel_size=3,
                dilation=2**i,
                cond_dim=hidden_dim
            ) for i in range(num_layers)
        ])
        
        self.pointwise_head = nn.Sequential(
            nn.Conv1d(hidden_dim, hidden_dim // 2, kernel_size=1),
            nn.GELU(),
            nn.Conv1d(hidden_dim // 2, 2, kernel_size=1)
        )
        nn.init.zeros_(self.pointwise_head[-1].weight)
        nn.init.zeros_(self.pointwise_head[-1].bias)
        
        self.global_pool = nn.AdaptiveAvgPool1d(1)
        self.global_nig_head = nn.Sequential(
            nn.Linear(hidden_dim * 2, hidden_dim),
            nn.GELU(),
            nn.Linear(hidden_dim, 4)
        )
        nn.init.zeros_(self.global_nig_head[-1].weight)
        nn.init.zeros_(self.global_nig_head[-1].bias)

        self.pyramid_scales = [1, 2, 4, 8, 16]
        self.pyramid_heads = nn.ModuleDict({
            str(k): nn.Conv1d(hidden_dim, 2, kernel_size=1) for k in self.pyramid_scales
        })
        for k in self.pyramid_scales:
            nn.init.zeros_(self.pyramid_heads[str(k)].weight)
            nn.init.zeros_(self.pyramid_heads[str(k)].bias)

        self._iemb_dim = hidden_dim // 4
        self.interval_embed = nn.Sequential(
            nn.Linear(2, self._iemb_dim), nn.GELU(),
            nn.Linear(self._iemb_dim, self._iemb_dim), nn.GELU())
        self.interval_head = nn.Sequential(
            nn.Linear(hidden_dim + self._iemb_dim + hidden_dim, hidden_dim), nn.GELU(),
            nn.Linear(hidden_dim, 2))
        nn.init.zeros_(self.interval_head[-1].weight)
        nn.init.zeros_(self.interval_head[-1].bias)

        self.interval_infl_head = nn.Sequential(
            nn.Linear(hidden_dim + self._iemb_dim + hidden_dim, hidden_dim), nn.GELU(),
            nn.Linear(hidden_dim, 2))
        nn.init.zeros_(self.interval_infl_head[-1].weight)
        nn.init.zeros_(self.interval_infl_head[-1].bias)

    def load_state_dict(self, state_dict, strict=True):
        """Backward-compatible load: checkpoints trained before the multi-resolution pyramid head
        lack `pyramid_heads.*`. Load everything they DO have and leave the pyramid at its
        zero-init (a no-op for the existing scalar/pointwise consumers) rather than hard-failing."""
        try:
            return super().load_state_dict(state_dict, strict=strict)
        except RuntimeError:
            return super().load_state_dict(state_dict, strict=False)

    @staticmethod
    def _interval_pool(h, a, b):
        """Weighted-mean pool of h [B,C,L] over positions in [a,b] (a,b [B], in [0,1])."""
        L = h.shape[2]
        pos = torch.linspace(0.0, 1.0, L, device=h.device)
        mask = ((pos[None, :] >= a[:, None]) & (pos[None, :] <= b[:, None])).float()
        mask = mask / (mask.sum(1, keepdim=True) + 1e-6)
        return (h * mask[:, None, :]).sum(2)

    def encode_interval_context(self, X_seq, z_inv, calib_strokes=None):
        """Encode one stroke once for arbitrary batched interval queries.

        Returns the shared TCN field and its global conditioning.  This is additive to the
        training/inference `forward` contract and exists so Core ML can run the expensive trunk once
        per stroke, then invoke only the small interval heads for many `(a,b)` queries.
        """
        batch_size = X_seq.shape[0]
        e_user = self.user_encoder(calib_strokes)
        if e_user.shape[0] == 1 and batch_size > 1:
            e_user = e_user.expand(batch_size, -1)
        cond_input = torch.cat([z_inv, e_user], dim=-1)
        z_emb = self.cond_mlp(cond_input)
        h = self.seq_proj(X_seq)
        for block in self.tcn_blocks:
            h = block(h, z_emb)
        return h, z_emb, e_user

    def query_intervals(self, h, z_emb, interval):
        """Return continuous magnitude and correlation distributions from cached features."""
        if z_emb.shape[0] == 1 and interval.shape[0] > 1:
            z_emb = z_emb.expand(interval.shape[0], -1)
        a = interval[:, 0].clamp(0.0, 1.0)
        b = interval[:, 1].clamp(0.0, 1.0)
        pooled = self._interval_pool(h, a, b)
        center = (a + b) * 0.5
        log_w = torch.log((b - a).clamp(min=1e-3))
        iemb = self.interval_embed(torch.stack([center, log_w], dim=1))
        z = torch.cat([pooled, iemb, z_emb], dim=1)
        iv = self.interval_head(z)
        ivf = self.interval_infl_head(z)
        return iv[:, 0], iv[:, 1], ivf[:, 0], ivf[:, 1]

    def forward(self, X_seq, z_inv, calib_strokes=None, interval=None):
        """
        Args:
            X_seq: [Batch, 4, Length] - (x, y, signed_curvature, arc_pos_s) per
                x and y are absolute position; s is the arc-length index. Changing the channel
                set requires updating in_seq_channels and retraining.
            z_inv: [Batch, 52] - 52 global trajectory invariants for active stroke
            calib_strokes: [Batch, K, 52] - K calibration/previous stroke invariants (optional)
        Returns:
            dict containing:
                - 'pointwise_sigma': [Batch, 2, Length] - pointwise (sigma_perp, sigma_parallel);
                  see the Ch0/Ch1 semantics note on `pointwise_head` above.
                - 'log_sigma_pointwise': [Batch, 2, Length]
                - 'e_user': [Batch, 32] - 32-D Personal Motor Fingerprint vector
                - 'gamma', 'nu', 'alpha', 'beta': Evidential NIG parameters
                - 'aleatoric_var': predicted physical motor variance
                - 'epistemic_var': predicted OOD shape uncertainty
                - 'predicted_rms': global predicted RMS noise floor
        """
        batch_size = X_seq.shape[0]
        
        h, z_emb, e_user = self.encode_interval_context(X_seq, z_inv, calib_strokes)
            
        log_sigma_pt = self.pointwise_head(h)
        sigma_pt = torch.exp(log_sigma_pt)
        
        h_pool = self.global_pool(h).squeeze(2)
        h_global = torch.cat([h_pool, z_emb], dim=-1)
        
        nig_raw = self.global_nig_head(h_global)
        gamma = nig_raw[:, 0]
        nu = F.softplus(nig_raw[:, 1]) + 1e-4
        alpha = F.softplus(nig_raw[:, 2]) + 1.0 + 1e-4
        beta = F.softplus(nig_raw[:, 3]) + 1e-4
        
        aleatoric_var = beta / (alpha - 1.0 + 1e-9)
        epistemic_var = beta / (nu * (alpha - 1.0) + 1e-9)

        pyramid = {}
        for k in self.pyramid_scales:
            hk = F.adaptive_avg_pool1d(h, k)
            pyramid[k] = self.pyramid_heads[str(k)](hk)

        out = {
            "pointwise_sigma": sigma_pt,
            "log_sigma_pointwise": log_sigma_pt,
            "e_user": e_user,
            "gamma": gamma,
            "nu": nu,
            "alpha": alpha,
            "beta": beta,
            "aleatoric_var": aleatoric_var,
            "epistemic_var": epistemic_var,
            "predicted_log_rms": gamma,
            "predicted_rms": torch.pow(10.0, gamma),
            "pyramid": pyramid
        }

        if interval is not None:
            mu, log_var, infl_mu, infl_log_var = self.query_intervals(h, z_emb, interval)
            out["interval_mu"] = mu
            out["interval_log_var"] = log_var
            out["interval_infl_mu"] = infl_mu
            out["interval_infl_log_var"] = infl_log_var
        return out

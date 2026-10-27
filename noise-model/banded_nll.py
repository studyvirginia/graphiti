"""Banded (AR(1)) covariance Gaussian NLL for the noise net.

Real drawing residuals are correlated (delta_perp lag-1 autocorr ~0.87), so the standard
per-point INDEPENDENT Gaussian NLL treats a stroke's L samples as L independent observations
when the effective count is far smaller. This scores the WHOLE ordered residual vector under an
AR(1) covariance whose precision is tridiagonal -> closed form, O(L), no matrix inverse, no
architecture change. Reduces EXACTLY to the diagonal Gaussian NLL at rho=0.

Model:  r ~ N(0, Sigma),  Sigma = D^{1/2} C(rho) D^{1/2},  D = diag(exp(log_var)),
        C(rho)_{ij} = rho^{|i-j|}  (AR(1)).
NLL  =  0.5 * ( log|Sigma| + r^T Sigma^{-1} r ) + const.

For AR(1), with whitened z_i = r_i / sqrt(exp(log_var_i)):
    log|Sigma| = sum_i log_var_i + (L-1)*log(1-rho^2)
    z^T Omega z = ( sum_i z_i^2 + rho^2 * sum_{i=1}^{L-2} z_i^2 - 2*rho * sum_i z_i z_{i+1} ) / (1-rho^2)
"""
from __future__ import annotations

import math

import torch
from torch import nn


def banded_ar1_nll(residual: torch.Tensor, log_var: torch.Tensor, rho: float | torch.Tensor,
                   reduction: str = "mean") -> torch.Tensor:
    """AR(1)-banded Gaussian NLL of an ordered residual vector.

    Args:
        residual: [B, L] observed residual r = (target - mean) along the stroke's arc order.
        log_var:  [B, L] per-point ln(variance) predicted by the net.
        rho:      scalar in (-1, 1) (float or 0-dim/broadcastable tensor). rho=0 -> diagonal NLL.
        reduction: 'mean' | 'sum' | 'none' (per-sample [B] NLL).
    Returns:
        scalar (mean/sum) or [B] tensor. Excludes the additive 0.5*L*log(2*pi) constant.
    """
    assert residual.shape == log_var.shape and residual.dim() == 2, "expect [B, L]"
    L = residual.shape[1]
    if not torch.is_tensor(rho):
        rho = residual.new_tensor(float(rho))
    rho = rho.clamp(-0.999, 0.999)
    one_m_r2 = 1.0 - rho * rho

    z = residual * torch.exp(-0.5 * log_var)
    z2 = z * z
    cross = z[:, :-1] * z[:, 1:]

    if L >= 2:
        interior = z2[:, 1:-1].sum(dim=1) if L > 2 else z.new_zeros(z.shape[0])
        quad = (z2.sum(dim=1) + rho * rho * interior - 2.0 * rho * cross.sum(dim=1)) / one_m_r2
        logdet = log_var.sum(dim=1) + (L - 1) * torch.log(one_m_r2)
    else:
        quad = z2.sum(dim=1)
        logdet = log_var.sum(dim=1)

    nll = 0.5 * (logdet + quad)
    if reduction == "mean":
        return nll.mean()
    if reduction == "sum":
        return nll.sum()
    return nll


def n_eff_factor(rho: float) -> float:
    """Effective-sample deflation for an AR(1) block mean: n_eff = L*(1-rho)/(1+rho).
    Returns the multiplier (1-rho)/(1+rho) in (0,1] to down-weight a region's confidence when
    the residual is correlated."""
    return (1.0 - rho) / (1.0 + rho)


class BandedGaussianNLL(nn.Module):
    """Drop-in replacement for a diagonal Gaussian NLL on ordered residual channels.

    Wraps `banded_ar1_nll`. `rho` fixed by default; pass learn_rho=True to expose a
    single learnable scalar (still no backbone change)."""

    def __init__(self, rho: float = 0.87, learn_rho: bool = False):
        super().__init__()
        if learn_rho:
            self.rho_raw = nn.Parameter(torch.tensor(math.atanh(max(-0.999, min(0.999, rho)))))
            self._learn = True
        else:
            self.register_buffer("rho_fixed", torch.tensor(float(rho)))
            self._learn = False

    @property
    def rho(self) -> torch.Tensor:
        return torch.tanh(self.rho_raw) if self._learn else self.rho_fixed

    def forward(self, residual: torch.Tensor, log_var: torch.Tensor, reduction: str = "mean"):
        return banded_ar1_nll(residual, log_var, self.rho, reduction=reduction)

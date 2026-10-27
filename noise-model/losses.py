"""Loss Functions for the Universal Power-Law Agnostic Dual-Stream TCN.

Implements distribution-based losses:
1. Evidential Deep Learning (EDL) NIG Loss (Amini et al., 2020) for global distortion scale.
2. Heavy-Tailed Student-t Negative Log-Likelihood for pointwise spatial residuals.
3. Sobolev W^{1, 2} and W^{2, 2} smoothness regularization on variance fields.
"""
import math
import torch
import torch.nn as nn
import torch.nn.functional as F


class EvidentialNIGLoss(nn.Module):
    """Evidential loss for a Normal-Inverse-Gamma conjugate prior.
    L_NIG = -log p(y | gamma, nu, alpha, beta) + lambda_reg * |y - gamma| * (2*nu + alpha)
    """
    def __init__(self, coeff_reg=0.05):
        super().__init__()
        self.coeff_reg = coeff_reg

    def forward(self, gamma, nu, alpha, beta, y_true):
        gamma = torch.clamp(gamma, min=-10.0, max=10.0)
        nu = torch.clamp(nu, min=1e-3, max=1e4)
        alpha = torch.clamp(alpha, min=1.01, max=1e4)
        beta = torch.clamp(beta, min=1e-3, max=1e4)
        y_true = torch.clamp(y_true, min=-10.0, max=10.0)
        
        om = 2.0 * beta * (1.0 + nu)
        
        term1 = 0.5 * torch.log(math.pi / nu)
        term2 = -alpha * torch.log(om)
        diff_sq = (y_true - gamma) ** 2
        term3 = (alpha + 0.5) * torch.log(om + nu * diff_sq)
        term4 = torch.lgamma(alpha) - torch.lgamma(alpha + 0.5)
        
        loss_nll = term1 + term2 + term3 + term4
        
        error = torch.abs(y_true - gamma)
        loss_reg = error * (2.0 * nu + alpha)
        
        loss_tot = loss_nll + self.coeff_reg * loss_reg
        loss_tot = torch.nan_to_num(loss_tot, nan=10.0, posinf=50.0, neginf=-10.0)
        return torch.mean(loss_tot)


class StudentTPointwiseLoss(nn.Module):
    """Heavy-Tailed Student-t NLL Loss with Learnable Degrees of Freedom nu.
    
    Safe Parameterization:
        nu = 2.05 + softplus(log_nu_raw)
    Guarantees nu > 2.0 (finite mean and variance) while preventing numerical instability.
    """
    def __init__(self, init_nu=3.0):
        super().__init__()
        init_raw = math.log(max(math.exp(init_nu - 2.05) - 1.0, 1e-4))
        self.log_nu_raw = nn.Parameter(torch.tensor(init_raw, dtype=torch.float32))

    @property
    def nu(self):
        return 2.05 + F.softplus(torch.clamp(self.log_nu_raw, min=-5.0, max=5.0))

    def forward(self, log_sigma, delta_perp):
        nu = self.nu
        log_sigma = torch.clamp(log_sigma, min=-6.0, max=6.0)
        sigma = torch.exp(log_sigma)
        delta_perp = torch.clamp(delta_perp, min=-20.0, max=20.0)
        z = delta_perp / (sigma + 1e-4)
        
        c = torch.lgamma((nu + 1.0) / 2.0) - torch.lgamma(nu / 2.0) - 0.5 * torch.log(math.pi * nu)
        
        loss = log_sigma + 0.5 * (nu + 1.0) * torch.log(1.0 + (z ** 2) / nu) - c
        loss = torch.nan_to_num(loss, nan=5.0, posinf=20.0, neginf=-5.0)
        return torch.mean(loss)


class SobolevSmoothnessPenalty(nn.Module):
    """Penalizes high-frequency derivative fluctuations on predicted log-variance fields."""
    def __init__(self, lambda_grad=0.01):
        super().__init__()
        self.lambda_grad = lambda_grad

    def forward(self, log_sigma):
        d1 = torch.diff(log_sigma, dim=-1)
        loss_grad = torch.mean(d1 ** 2)
        return self.lambda_grad * loss_grad

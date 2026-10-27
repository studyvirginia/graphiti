#!/usr/bin/env python3
"""Export the pyramid continuous-noise model as an encode-once / field-query Core ML pair.

The encoder consumes one resampled stroke and its standardized invariant vector and
caches the shared TCN field.  The field-query model turns that field into the C++
engine's immutable point-noise field: it fuses ALL trained pyramid scales
(1,2,4,8,16) into one (mu, var) per finest-grid region and applies the exact
`_student_t_from_lognormal` moment-match, so the exported per-region sigma^2 is
byte-for-byte the `app/testbed/sigma_field.multiscale_density` value the validated
Python selection path uses.  This replaces the earlier interval-head export, whose
narrow-slice queries were strongly width-dependent (mu ~ -0.7 on a clean circle vs
the pyramid head's -2.0) and saturated the on-device selection band.

The exported checkpoint carries the pyramid heads and no interval heads.

This is an export and verification utility, not training code.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys

import coremltools as ct
import numpy as np
import torch
from torch import nn
import torch.nn.functional as F

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

NINV = 28
from models.tcn_noise_field import DualStreamTCNNoiseField

SEQUENCE_LENGTH = 48
HIDDEN_DIM = 128
REGION_COUNT = 16
SCALES = (1, 2, 4, 8, 16)
LOGVAR_CLAMP = (-8.0, 4.0)
LN10 = math.log(10.0)


class Encoder(nn.Module):
    def __init__(self, model: DualStreamTCNNoiseField):
        super().__init__()
        self.model = model

    def forward(self, sequence: torch.Tensor, invariants: torch.Tensor):
        h, z_emb, _ = self.model.encode_interval_context(sequence, invariants)
        return h, z_emb, self.model.pointwise_head(h)


class FieldQuery(nn.Module):
    """Fuse the pyramid heads on the cached TCN field into one per-region sigma^2 field.

    Mirrors `sigma_field.multiscale_density` (precision-weighted, correlated-limit fusion)
    followed by `_student_t_from_lognormal` + `_to_engine` MINUS the global ``ratio**2``
    (a stroke-dependent scalar the Swift caller applies) and MINUS the optional mean-1
    normalization (the native selector needs the ABSOLUTE engine-frame field for its
    physical-plausibility gate). Outputs are `region_s2` (net frame) and `region_log_var`
    (= ln(fused Var(log10 rms)), which `net_band_select` clips to [-8,4] for its epistemic
    factor)."""

    def __init__(self, model: DualStreamTCNNoiseField, region_count: int = REGION_COUNT,
                 scales=SCALES):
        super().__init__()
        self.model = model
        self.region_count = region_count
        self.scales = list(scales)
        for k in self.scales:
            idx = [min(j * k // region_count, k - 1) for j in range(region_count)]
            self.register_buffer(f"idx_{k}", torch.tensor(idx, dtype=torch.long), persistent=False)

    def forward(self, hidden: torch.Tensor):
        num_mu = torch.zeros(self.region_count)
        num_var = torch.zeros(self.region_count)
        den = torch.zeros(self.region_count)
        for k in self.scales:
            p = self.model.pyramid_heads[str(k)](F.adaptive_avg_pool1d(hidden, k))[0]
            mu = p[0]
            var = torch.exp(torch.clamp(p[1], LOGVAR_CLAMP[0], LOGVAR_CLAMP[1]))
            idx = getattr(self, f"idx_{k}")
            muK = mu[idx]
            varK = var[idx]
            w = 1.0 / torch.clamp(varK, min=1e-12)
            num_mu = num_mu + w * muK
            num_var = num_var + w * varK
            den = den + w
        mu_f = num_mu / den
        var_f = num_var / den
        a = 2.0 * LN10 * mu_f
        b2 = torch.clamp(4.0 * (LN10 ** 2) * var_f, 1e-12, 50.0)
        M = torch.exp(torch.clamp(a + b2 / 2.0, -60.0, 60.0))
        alpha = torch.clamp(2.0 + 1.0 / (torch.exp(b2) - 1.0), 2.01, 1e6)
        region_s2 = M * (alpha - 1.0) / alpha
        return region_s2, torch.log(var_f)


def load_model(checkpoint_path: Path) -> tuple[DualStreamTCNNoiseField, dict]:
    checkpoint = torch.load(checkpoint_path, map_location="cpu", weights_only=False)
    model = DualStreamTCNNoiseField(4, NINV, 32, HIDDEN_DIM, 4).eval()
    model.load_state_dict(checkpoint["model"], strict=False)
    return model, checkpoint


def convert(trace: torch.jit.ScriptModule, inputs: list[ct.TensorType], output_path: Path,
            output_names: list[str], precision: ct.precision):
    converted = ct.convert(
        trace,
        convert_to="mlprogram",
        minimum_deployment_target=ct.target.iOS17,
        compute_precision=precision,
        inputs=inputs,
        outputs=[ct.TensorType(name=name) for name in output_names],
    )
    converted.save(str(output_path))
    return converted


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkpoint", type=Path,
                        default=Path("models/checkpoints/fullres_c7000.pt"))
    parser.add_argument("--output-dir", type=Path,
                        default=Path("ios/GraphArt/Models/ContinuousNoise"))
    parser.add_argument("--force", action="store_true",
                        help="overwrite existing exported packages")
    args = parser.parse_args()

    model, checkpoint = load_model(args.checkpoint)
    if not any(k.startswith("pyramid_heads") for k in checkpoint["model"]):
        raise SystemExit(f"{args.checkpoint} has no pyramid heads; wrong checkpoint")
    encoder = Encoder(model).eval()
    field = FieldQuery(model).eval()
    sequence = torch.randn(1, 4, SEQUENCE_LENGTH)
    invariants = torch.randn(1, NINV)
    with torch.no_grad():
        hidden, embedding, pointwise_log_sigma = encoder(sequence, invariants)
        s2_ref, lv_ref = field(hidden)
        out = model(sequence, invariants)
        num_mu = np.zeros(REGION_COUNT); num_var = np.zeros(REGION_COUNT); den = np.zeros(REGION_COUNT)
        for k in SCALES:
            pk = out["pyramid"][k][0].numpy()
            muk = pk[0]; vk = np.exp(np.clip(pk[1], *LOGVAR_CLAMP))
            idx = np.minimum(np.arange(REGION_COUNT) * k // REGION_COUNT, k - 1)
            w = 1.0 / np.maximum(vk[idx], 1e-12)
            num_mu += w * muk[idx]; num_var += w * vk[idx]; den += w
        mu_f = num_mu / den; var_f = num_var / den
        a = 2 * LN10 * mu_f; b2 = np.clip(4 * LN10 ** 2 * var_f, 1e-12, 50.0)
        M = np.exp(np.clip(a + b2 / 2, -60, 60)); alpha = np.clip(2 + 1 / np.expm1(b2), 2.01, 1e6)
        s2_expect = M * (alpha - 1.0) / alpha
        if not np.allclose(s2_ref.numpy(), s2_expect, rtol=1e-5, atol=1e-9):
            raise RuntimeError("FieldQuery diverged from multiscale_density fusion math")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    encoder_path = args.output_dir / "ContinuousNoiseEncoder.mlpackage"
    field_path = args.output_dir / "ContinuousNoiseIntervalQuery.mlpackage"
    if not args.force and (encoder_path.exists() or field_path.exists()):
        raise SystemExit(f"Refusing to overwrite completed export in {args.output_dir} (use --force)")
    import shutil
    for p in (encoder_path, field_path):
        if p.exists():
            shutil.rmtree(p)

    traced_encoder = torch.jit.trace(encoder, (sequence, invariants), strict=False)
    traced_field = torch.jit.trace(field, (hidden,), strict=False)
    convert(traced_encoder, [
        ct.TensorType(name="sequence", shape=sequence.shape, dtype=np.float32),
        ct.TensorType(name="invariants", shape=invariants.shape, dtype=np.float32),
    ], encoder_path, ["hidden", "embedding", "log_sigma_pointwise"], ct.precision.FLOAT16)
    convert(traced_field, [
        ct.TensorType(name="hidden", shape=hidden.shape, dtype=np.float32),
    ], field_path, ["region_s2", "region_log_var"], ct.precision.FLOAT32)

    metadata = {
        "format_version": 2,
        "checkpoint": str(args.checkpoint),
        "sequence_shape": [1, 4, SEQUENCE_LENGTH],
        "invariant_shape": [1, NINV],
        "hidden_shape": [1, HIDDEN_DIM, SEQUENCE_LENGTH],
        "embedding_shape": [1, HIDDEN_DIM],
        "region_count": REGION_COUNT,
        "scales": list(SCALES),
        "encoder_outputs": ["hidden", "embedding", "log_sigma_pointwise"],
        "field_outputs": ["region_s2", "region_log_var"],
        "checkpoint_epoch": checkpoint.get("epoch"),
        "invariant_mean": np.asarray(checkpoint.get("mean", [])).tolist(),
        "invariant_std": np.asarray(checkpoint.get("std", [])).tolist(),
    }
    (args.output_dir / "continuous_noise_manifest.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(f"exported {encoder_path}")
    print(f"exported {field_path}")
    print(f"region_s2[min/med/max] = {float(s2_ref.min()):.3e} {float(s2_ref.median()):.3e} {float(s2_ref.max()):.3e}")


if __name__ == "__main__":
    main()

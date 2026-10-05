#!/usr/bin/env python3
"""Bounded geometry-to-electronic-response development campaign.

Only development chemistry is accepted by this program.  It implements the two
pre-authorized equivariant heads and chemistry-blocked LOCO evaluation using the
counterpoise response observables as the primary loss.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import random
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from ase.io import read
from torch import nn

BOHR_TO_ANGSTROM = 0.529177210903
DEBYE_PER_E_BOHR = 2.541746473


@dataclass
class Sample:
    config_id: str
    molecule_id: str
    regime: str
    positions: torch.Tensor
    numbers: torch.Tensor
    scalars: torch.Tensor
    vectors: torch.Tensor
    geometry_vectors: torch.Tensor
    design: torch.Tensor
    esp: torch.Tensor
    dipole: torch.Tensor
    q_target: torch.Tensor
    u_target: torch.Tensor
    base_q: torch.Tensor
    base_u: torch.Tensor
    train_points: torch.Tensor


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def esp_design(
    positions_angstrom: np.ndarray, points_angstrom: np.ndarray
) -> np.ndarray:
    positions = positions_angstrom / BOHR_TO_ANGSTROM
    points = points_angstrom / BOHR_TO_ANGSTROM
    displacement = points[:, None, :] - positions[None, :, :]
    distance = np.linalg.norm(displacement, axis=2)
    return np.concatenate(
        [
            1.0 / distance,
            (displacement / distance[:, :, None] ** 3).reshape(len(points), -1),
        ],
        axis=1,
    )


def geometry_basis(positions: np.ndarray, numbers: np.ndarray) -> np.ndarray:
    """Fixed E(3)-equivariant local directional basis (4 radii x 8 elements)."""
    elements = [1, 6, 7, 8, 9, 16, 17, 35]
    centers = [1.2, 2.2, 3.2, 4.4]
    n = len(numbers)
    result = np.zeros((n, len(elements) * len(centers), 3), dtype=np.float32)
    displacement = positions[None, :, :] - positions[:, None, :]
    distance = np.linalg.norm(displacement, axis=2)
    unit = displacement / np.maximum(distance[:, :, None], 1e-8)
    for element_index, element in enumerate(elements):
        mask = numbers == element
        for radial_index, center in enumerate(centers):
            weight = np.exp(-(((distance - center) / 0.75) ** 2))
            weight *= 0.5 * (np.cos(np.pi * np.minimum(distance, 6.0) / 6.0) + 1.0)
            weight[:, ~mask] = 0.0
            np.fill_diagonal(weight, 0.0)
            result[:, element_index * len(centers) + radial_index] = np.einsum(
                "ij,ijc->ic", weight, unit
            )
    return result


def load_registry(path: Path) -> dict[str, Path]:
    with (path / "observable_registry.csv").open() as handle:
        return {
            row["config_id"]: path / row["observable_file"]
            for row in csv.DictReader(handle)
        }


def load_samples(args) -> list[Sample]:
    frames = []
    for path in args.configurations:
        frames.extend(read(path, index=":"))
    frame_map = {str(a.info["config_id"]): a for a in frames}
    registry = {}
    for path in args.observables:
        overlap = set(registry) & set(load_registry(path))
        if overlap:
            raise RuntimeError(f"Duplicate response records: {sorted(overlap)}")
        registry.update(load_registry(path))
    feature = np.load(args.features)
    base_lookup = {}
    base_path = getattr(args, "base_predictions", None)
    if base_path is not None:
        base = np.load(base_path)
        for i, key in enumerate([str(x) for x in base["config_ids"]]):
            begin, end = int(base["offsets"][i]), int(base["offsets"][i + 1])
            base_lookup[key] = (
                np.asarray(base["charges_e"][begin:end]),
                np.asarray(base["dipoles_e_bohr"][begin:end]),
            )
    feature_ids = [str(x) for x in feature["config_ids"]]
    offsets = feature["offsets"]
    samples = []
    rng = np.random.default_rng(20260813)
    for config_index, config_id in enumerate(feature_ids):
        if config_id not in registry:
            continue
        atoms = frame_map[config_id]
        first, last = int(offsets[config_index]), int(offsets[config_index + 1])
        if args.encoder == "mace-polar-m":
            product = np.asarray(feature["product_1"][first:last])
            n_scalar = product.shape[1] // 4
            scalar = product[:, :n_scalar]
            # e3nn real l=1 ordering is y,z,x; store Cartesian x,y,z.
            vector = product[:, n_scalar:].reshape(len(atoms), n_scalar, 3)[
                :, :, [2, 0, 1]
            ]
        else:
            scalar = np.asarray(feature["scalars"][first:last])
            vector = np.asarray(feature["vectors_l1_cartesian"][first:last])
        values = np.load(registry[config_id])
        points = np.asarray(values["points_angstrom"])
        esp = np.asarray(values["delta_esp_hartree_per_e"])
        dipole = np.asarray(values["delta_dipole_debye"])
        design = esp_design(np.asarray(atoms.positions), points)
        q = np.asarray(values["fitted_charges_e"])
        u = np.asarray(values["fitted_dipoles_e_bohr"])
        base_q, base_u = base_lookup.get(
            config_id, (np.zeros_like(q), np.zeros_like(u))
        )
        if len(esp) > 384:
            chosen = np.sort(rng.choice(len(esp), 384, replace=False))
        else:
            chosen = np.arange(len(esp))
        samples.append(
            Sample(
                config_id=config_id,
                molecule_id=str(atoms.info["molecule_id"]),
                regime=str(atoms.info["regime"]),
                positions=torch.tensor(
                    np.asarray(atoms.positions), dtype=torch.float32
                ),
                numbers=torch.tensor(atoms.numbers, dtype=torch.long),
                scalars=torch.tensor(scalar, dtype=torch.float32),
                vectors=torch.tensor(vector, dtype=torch.float32),
                geometry_vectors=torch.tensor(
                    geometry_basis(
                        np.asarray(atoms.positions), np.asarray(atoms.numbers)
                    ),
                    dtype=torch.float32,
                ),
                design=torch.tensor(design, dtype=torch.float32),
                esp=torch.tensor(esp, dtype=torch.float32),
                dipole=torch.tensor(dipole, dtype=torch.float32),
                q_target=torch.tensor(q, dtype=torch.float32),
                u_target=torch.tensor(u, dtype=torch.float32),
                base_q=torch.tensor(base_q, dtype=torch.float32),
                base_u=torch.tensor(base_u, dtype=torch.float32),
                train_points=torch.tensor(chosen, dtype=torch.long),
            )
        )
    if len(samples) != 48 or len({s.molecule_id for s in samples}) != 14:
        raise RuntimeError(
            f"Expected exactly 48 environments / 14 chemistries, got {len(samples)} / "
            f"{len({s.molecule_id for s in samples})}"
        )
    return samples


class ResponseHead(nn.Module):
    def __init__(self, samples: list[Sample], head: str):
        super().__init__()
        scalar = torch.cat([s.scalars for s in samples])
        vector = torch.cat([s.vectors for s in samples])
        q = torch.cat([s.q_target - s.base_q for s in samples])
        u = torch.cat([s.u_target - s.base_u for s in samples])
        self.register_buffer("scalar_mean", scalar.mean(0))
        self.register_buffer("scalar_std", scalar.std(0).clamp_min(1e-3))
        self.register_buffer(
            "vector_rms", vector.square().mean((0, 2)).sqrt().clamp_min(1e-4)
        )
        self.register_buffer("q_scale", q.square().mean().sqrt().clamp_min(1e-4))
        self.register_buffer("u_scale", u.square().mean().sqrt().clamp_min(1e-4))
        residual_mu = []
        for sample in samples:
            target_mu = sample.dipole / DEBYE_PER_E_BOHR
            base_mu = torch.sum(
                sample.base_q[:, None] * (sample.positions / BOHR_TO_ANGSTROM), 0
            ) + torch.sum(sample.base_u, 0)
            residual_mu.append(target_mu - base_mu)
        self.register_buffer(
            "mu_scale", torch.cat(residual_mu).square().mean().sqrt().clamp_min(1e-3)
        )
        if head == "C":
            pooled, targets = [], []
            for sample, target in zip(samples, residual_mu):
                normalized_v = sample.vectors / self.vector_rms[None, :, None]
                # Molecular induced dipole is extensive, so the directly
                # supervised global response uses sum pooling.
                pooled.append(normalized_v.sum(0).T)
                targets.append(target)
            design = torch.cat(pooled, dim=0)
            target = torch.cat(targets, dim=0)
            alpha = 10000.0
            weight = torch.linalg.solve(
                design.T @ design
                + alpha * torch.eye(design.shape[1], device=design.device),
                design.T @ target,
            )
            self.register_buffer("dipole_ridge_weight", weight)
        self._build(scalar.shape[1], head)

    def _build(self, channels: int, head: str):
        self.head = head
        if head == "A":
            self.scalar_net = nn.Sequential(
                nn.Linear(channels, 48), nn.SiLU(), nn.Linear(48, 1)
            )
            self.vector_weight = nn.Parameter(torch.zeros(channels))
            nn.init.normal_(self.vector_weight, std=0.01)
        elif head in ("B", "C"):
            hidden, vector_channels = 64, 16
            self.local = nn.Linear(channels, hidden)
            self.global_context = nn.Linear(channels, hidden, bias=False)
            self.vector_norm = nn.Linear(channels, hidden, bias=False)
            self.q_net = nn.Sequential(
                nn.SiLU(), nn.Linear(hidden, 32), nn.SiLU(), nn.Linear(32, 1)
            )
            self.vector_projection = nn.Parameter(
                torch.empty(channels, vector_channels)
            )
            nn.init.orthogonal_(self.vector_projection)
            self.vector_gate = nn.Sequential(
                nn.SiLU(),
                nn.Linear(hidden, 48),
                nn.SiLU(),
                nn.Linear(48, vector_channels + 32),
            )
            if head == "C":
                # Start exactly from the strong frozen response prior.  This is
                # an optimization safeguard, not a post-hoc gate.
                nn.init.zeros_(self.q_net[-1].weight)
                nn.init.zeros_(self.q_net[-1].bias)
                nn.init.zeros_(self.vector_gate[-1].weight)
                nn.init.zeros_(self.vector_gate[-1].bias)
        else:
            raise ValueError(head)

    @classmethod
    def from_state_dict(cls, state: dict[str, torch.Tensor], head: str):
        """Construct an inference-only head without reading any training labels."""
        model = cls.__new__(cls)
        nn.Module.__init__(model)
        for name in (
            "scalar_mean",
            "scalar_std",
            "vector_rms",
            "q_scale",
            "u_scale",
            "mu_scale",
        ):
            model.register_buffer(name, torch.empty_like(state[name]))
        if head == "C":
            model.register_buffer(
                "dipole_ridge_weight", torch.empty_like(state["dipole_ridge_weight"])
            )
        model._build(int(state["scalar_mean"].numel()), head)
        model.load_state_dict(state)
        return model

    def forward(self, sample: Sample):
        scalar = (sample.scalars - self.scalar_mean) / self.scalar_std
        vector = sample.vectors / self.vector_rms[None, :, None]
        if self.head == "A":
            q = self.scalar_net(scalar).squeeze(-1) * self.q_scale
            u = torch.einsum("ncd,c->nd", vector, self.vector_weight) * self.u_scale
        else:
            norm = vector.square().sum(2).sqrt()
            context = scalar.mean(0, keepdim=True).expand_as(scalar)
            hidden = (
                self.local(scalar)
                + self.global_context(context)
                + self.vector_norm(norm)
            )
            q = self.q_net(hidden).squeeze(-1) * self.q_scale
            learned_vectors = torch.einsum(
                "ncd,ck->nkd", vector, self.vector_projection
            )
            candidates = torch.cat([learned_vectors, sample.geometry_vectors], dim=1)
            gate = self.vector_gate(hidden)
            u = torch.einsum("nk,nkd->nd", gate, candidates) * self.u_scale
        q = q - q.mean()  # exact correction charge conservation
        if self.head == "C":
            # A separately learned equivariant global response fixes the residual
            # molecular dipole while the local channel controls the ESP shape.
            desired = torch.einsum("ncd,c->d", vector, self.dipole_ridge_weight)
            current = torch.sum(
                q[:, None] * (sample.positions / BOHR_TO_ANGSTROM), 0
            ) + torch.sum(u, 0)
            u = u + (desired - current)[None, :] / len(q)
            q = q + sample.base_q
            u = u + sample.base_u
            q = q - q.mean()
        return q, u


def move(sample: Sample, device: torch.device) -> Sample:
    return Sample(
        **{
            field: (
                getattr(sample, field).to(device)
                if torch.is_tensor(getattr(sample, field))
                else getattr(sample, field)
            )
            for field in Sample.__dataclass_fields__
        }
    )


def predicted_observables(
    sample: Sample, q: torch.Tensor, u: torch.Tensor, training: bool
):
    vector = torch.cat([q, u.reshape(-1)])
    if training:
        design = sample.design[sample.train_points]
        target = sample.esp[sample.train_points]
    else:
        design, target = sample.design, sample.esp
    esp = design @ vector
    positions_bohr = sample.positions / BOHR_TO_ANGSTROM
    dipole = (
        torch.sum(q[:, None] * positions_bohr, 0) + torch.sum(u, 0)
    ) * DEBYE_PER_E_BOHR
    return esp, target, dipole


def loss_for(model, samples: list[Sample], training: bool):
    dipole_scale = (
        torch.cat([s.dipole for s in samples]).square().mean().sqrt().clamp_min(0.02)
    )
    terms = []
    for sample in samples:
        q, u = model(sample)
        esp, target, dipole = predicted_observables(sample, q, u, training)
        signal = target.square().mean().clamp_min(1e-10)
        esp_loss = (esp - target).square().mean() / signal
        dipole_loss = (dipole - sample.dipole).square().mean() / dipole_scale.square()
        latent = 0.5 * (
            (q - sample.q_target).square().mean() / model.q_scale.square()
            + (u - sample.u_target).square().mean() / model.u_scale.square()
        )
        dipole_weight = 2.0 if model.head == "C" else 0.5
        terms.append(esp_loss + dipole_weight * dipole_loss + 0.05 * latent)
    return torch.stack(terms).mean()


def fit(
    train: list[Sample],
    validation: list[Sample],
    head: str,
    seed: int,
    max_epochs: int = 220,
):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    train_device = [move(s, device) for s in train]
    valid_device = [move(s, device) for s in validation]
    model = ResponseHead(train_device, head).to(device)
    optimizer = torch.optim.AdamW(model.parameters(), lr=2e-3, weight_decay=2e-5)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(
        optimizer, max_epochs, eta_min=2e-5
    )
    best, best_state, best_epoch, stale = math.inf, None, 0, 0
    for epoch in range(max_epochs):
        model.train()
        optimizer.zero_grad(set_to_none=True)
        loss = loss_for(model, train_device, True)
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
        optimizer.step()
        scheduler.step()
        if epoch >= 20 and (epoch % 4 == 0 or epoch == max_epochs - 1):
            model.eval()
            with torch.no_grad():
                score = float(loss_for(model, valid_device or train_device, False))
            if score < best - 1e-5:
                best, best_epoch, stale = score, epoch, 0
                best_state = {
                    k: v.detach().cpu().clone() for k, v in model.state_dict().items()
                }
            else:
                stale += 1
            if stale >= 18:
                break
    if best_state is None:
        best_state = {
            k: v.detach().cpu().clone() for k, v in model.state_dict().items()
        }
    model.load_state_dict(best_state)
    model.eval()
    return model.cpu(), best_epoch, best


def evaluate(model, samples: list[Sample], fold: str, candidate: str):
    rows = []
    with torch.no_grad():
        for sample in samples:
            q, u = model(sample)
            esp, target, dipole = predicted_observables(sample, q, u, False)
            error = esp - target
            signal = torch.sqrt(target.square().mean()).item()
            rows.append(
                {
                    "fold": fold,
                    "candidate": candidate,
                    "config_id": sample.config_id,
                    "molecule_id": sample.molecule_id,
                    "regime": sample.regime,
                    "esp_nrmse": torch.sqrt(error.square().mean()).item()
                    / max(signal, 1e-12),
                    "esp_rmse_hartree_per_e": torch.sqrt(error.square().mean()).item(),
                    "dipole_vector_rmse_debye": torch.sqrt(
                        (dipole - sample.dipole).square().mean()
                    ).item(),
                    "dipole_error_norm_debye": torch.linalg.vector_norm(
                        dipole - sample.dipole
                    ).item(),
                    "net_charge_e": q.sum().item(),
                    "response_rms_hartree_per_e": signal,
                }
            )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--configurations", type=Path, action="append", required=True)
    parser.add_argument("--observables", type=Path, action="append", required=True)
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--encoder", choices=("mace-polar-m", "visnet"), required=True)
    parser.add_argument("--head", choices=("A", "B", "C"), required=True)
    parser.add_argument("--base-predictions", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=20260813)
    args = parser.parse_args()
    samples = load_samples(args)
    molecules = sorted({s.molecule_id for s in samples})
    rows, epochs = [], []
    candidate = f"{args.encoder}/head-{args.head}"
    for fold_index, held in enumerate(molecules):
        train = [s for s in samples if s.molecule_id != held]
        validation = [s for s in samples if s.molecule_id == held]
        model, epoch, score = fit(train, validation, args.head, args.seed + fold_index)
        rows.extend(evaluate(model, validation, held, candidate))
        epochs.append(
            {"molecule_id": held, "best_epoch": epoch, "selection_loss": score}
        )
        print(
            f"LOCO {fold_index + 1}/{len(molecules)} {held} epoch={epoch} loss={score:.5f}",
            flush=True,
        )
    args.output.mkdir(parents=True, exist_ok=True)
    table = args.output / "loco_predictions.csv"
    pd.DataFrame(rows).to_csv(table, index=False)
    pd.DataFrame(epochs).to_csv(args.output / "loco_training.csv", index=False)
    frame = pd.DataFrame(rows)
    molecule = frame.groupby("molecule_id", as_index=False).agg(
        esp_nrmse=("esp_nrmse", "mean"),
        dipole_vector_rmse_debye=(
            "dipole_vector_rmse_debye",
            lambda x: float(np.sqrt(np.mean(np.square(x)))),
        ),
    )
    summary = {
        "candidate": candidate,
        "n_chemistries": len(molecules),
        "n_configurations": len(samples),
        "molecule_mean_esp_nrmse": float(molecule.esp_nrmse.mean()),
        "molecule_rms_dipole_vector_rmse_debye": float(
            np.sqrt(np.mean(molecule.dipole_vector_rmse_debye**2))
        ),
        "prediction_sha256": sha256(table),
        "energy_or_force_supervision_used": False,
        "prospective_identities_accessed": False,
        "prospective_labels_accessed": False,
        "experimental_hydration_targets_accessed": False,
    }
    (args.output / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Check that this curated paper companion is complete and internally coherent."""

from __future__ import annotations

import csv
import hashlib
import json
import re
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
IGNORED_DIRS = {".git", ".qa", ".venv", "__pycache__", ".pytest_cache", ".ruff_cache",
                "build", "dist", "example_prediction"}


def ignored(path: Path) -> bool:
    parts = path.relative_to(ROOT).parts
    return any(part in IGNORED_DIRS or part.endswith(".egg-info") for part in parts) or \
        parts[:2] == ("third_party", "checkpoints")


def table(path: str) -> list[dict[str, str]]:
    with (ROOT / path).open(newline="") as stream:
        return list(csv.DictReader(stream))


def check_links() -> None:
    broken = []
    pattern = re.compile(r"!?\[[^\]]*\]\(([^)]+)\)")
    html_images = re.compile(r'<img\s+[^>]*src="([^"]+)"')
    for markdown in ROOT.rglob("*.md"):
        if ignored(markdown) or "provenance" in markdown.relative_to(ROOT).parts:
            continue
        content = markdown.read_text()
        for target in pattern.findall(content) + html_images.findall(content):
            target = target.split("#", 1)[0]
            if not target or ":" in target.split("/", 1)[0] or target.startswith("#"):
                continue
            if not (markdown.parent / target).exists():
                broken.append(f"{markdown.relative_to(ROOT)} -> {target}")
    assert not broken, "Broken local Markdown links:\n" + "\n".join(broken)


def check_headlines() -> None:
    rows = table("figures/figure_03/paired.csv")
    assert len(rows) == 3
    assert sum(int(r["n_molecules"]) for r in rows) == 56
    assert sum(int(r["glider_molecule_wins"]) for r in rows) == 54
    expected = [
        (0.321799354439422, 0.524476079542344, 38.6436546886518),
        (0.342174873318955, 0.571007635391217, 40.0752543204577),
        (0.426580605812233, 0.629706400723695, 32.2572225211651),
    ]
    for row, (glider, baseline, reduction) in zip(rows, expected):
        assert abs(float(row["glider_value"]) - glider) < 1e-12
        assert abs(float(row["comparator_value"]) - baseline) < 1e-12
        assert abs(float(row["relative_error_reduction_pct"]) - reduction) < 1e-10
        assert float(row["ci95_high"]) < 0
    nonwater = table("figures/figure_S15/nonwater_original72_effects.csv")
    overall = next(r for r in nonwater if r["subset"] == "Overall")
    assert int(overall["n_solutes"]) == 12
    assert abs(float(overall["esp_glider"]) - 0.450371053714298) < 1e-12
    assert abs(float(overall["esp_comparator"]) - 0.617170202219818) < 1e-12
    contact_followup = table("figures/figure_S15/nonwater_effects.csv")
    contact_overall = next(r for r in contact_followup if r["label"] == "overall")
    assert int(contact_overall["n_configurations"]) == 36
    assert abs(float(contact_overall["esp_difference"]) + 0.2884889328434251) < 1e-12
    ranks = table("figures/figure_04/energy_rank.csv")
    assert [int(r["rank"]) for r in ranks] == list(range(4, 13))
    w4 = ranks[0]
    assert abs(float(w4["glider"]) - 0.04868219341822043) < 1e-12
    assert abs(float(w4["mace_polar_l"]) - 0.11401221689887638) < 1e-12
    separation = table("figures/figure_S07/separation.csv")
    assert len(separation) == 40
    for system, glider_20, glider_100, prior_20, prior_100 in (
        ("single_water", 0.563891174, 1.540125107, 0.195369430, 0.113099723),
        ("whole_environment", 0.924588953, 2.842730048, 0.300614666, 0.239373936),
    ):
        def amplitude(method: str, distance: float) -> float:
            row = next(r for r in separation if r["system"] == system
                       and r["method"] == method and float(r["distance_A"]) == distance)
            return float(row["esp_rms_mEh_per_e"])
        assert abs(amplitude("glider", 20) - glider_20) < 1e-8
        assert abs(amplitude("glider", 100) - glider_100) < 1e-8
        assert abs(amplitude("averaged_prior", 20) - prior_20) < 1e-8
        assert abs(amplitude("averaged_prior", 100) - prior_100) < 1e-8
        assert glider_100 > glider_20 and prior_100 < prior_20
    components = table("experiments/dissociation_extended/posthoc_components.csv")
    assert len(components) == 280
    water_audit = table("experiments/water_contact_audit/panel_sensitivity.csv")
    screened = {
        row["panel"]: row for row in water_audit
        if row["analysis"] == "exclude_solute_water_contacts_below_1.5_A"
    }
    assert {panel: int(row["n_configurations"]) for panel, row in screened.items()} == {
        "I": 45, "II": 94, "III": 77
    }
    assert all(float(row["paired_95pct_ci_high"]) < 0 for row in screened.values())
    short_energies = table("experiments/water_contact_audit/cp_interaction_energies_short.csv")
    assert len(short_energies) == 8
    assert sum(float(row["cp_interaction_energy_kcal_mol"]) > 0 for row in short_energies) == 7
    pair = json.loads((ROOT / "experiments/water_contact_audit/short_water_pair_cp.json").read_text())
    assert 20.2 < float(pair["cp_pair_energy_kcal_mol"]) < 20.4


def check_separation_qm() -> None:
    base = ROOT / "experiments/dissociation_qm"
    manifest = json.loads((base / "source_manifest.json").read_text())
    for name, expected in manifest["released_file_sha256"].items():
        actual = hashlib.sha256((base / name).read_bytes()).hexdigest()
        assert actual == expected, name
    for name, expected in manifest["source_input_sha256"].items():
        actual = hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
        assert actual == expected, name

    rows = table("experiments/dissociation_qm/summary.csv")
    assert len(rows) == len(manifest["configuration_ids"]) == 20
    with np.load(base / "solute_probe_points.npz", allow_pickle=False) as probes, \
            np.load(base / "predicted_probe_potentials.npz", allow_pickle=False) as predictions:
        points = probes["dev_cyclic_carbamate"]
        assert points.shape == (512, 3)
        assert len(predictions.files) == 120
        for row in rows:
            config = (f"dissociation__dev_cyclic_carbamate__{row['system']}__"
                      f"{row['distance_A']}A")
            assert config in manifest["configuration_ids"]
            record = json.loads((base / f"{config}.json").read_text())
            arrays_file = base / record["arrays_file"]
            assert hashlib.sha256(arrays_file.read_bytes()).hexdigest() == record["arrays_sha256"]
            assert record["n_probe_points"] == 512
            assert all(attempts[-1]["converged"] for attempts in record["attempts"].values())
            with np.load(arrays_file, allow_pickle=False) as arrays:
                np.testing.assert_array_equal(arrays["points_angstrom"], points)
                components = arrays["component_esp_hartree_per_e"]
                np.testing.assert_allclose(
                    arrays["qm_response_esp_hartree_per_e"],
                    components[0] - components[1] - components[2], atol=1e-12)
                for method, key in (("glider", "glider_esp_hartree_per_e"),
                                    ("averaged_prior", "averaged_prior_esp_hartree_per_e")):
                    np.testing.assert_array_equal(arrays[key], predictions[f"{config}__{method}"])
                qm = arrays["qm_response_esp_hartree_per_e"]
                for column, values in (
                    ("qm_response_rms_mEh_per_e", qm),
                    ("glider_rms_mEh_per_e", arrays["glider_esp_hartree_per_e"]),
                    ("prior_rms_mEh_per_e", arrays["averaged_prior_esp_hartree_per_e"]),
                    ("glider_error_rms_mEh_per_e", arrays["glider_esp_hartree_per_e"] - qm),
                    ("prior_error_rms_mEh_per_e", arrays["averaged_prior_esp_hartree_per_e"] - qm),
                ):
                    value = float(np.sqrt(np.mean(values**2)) * 1000)
                    assert abs(value - float(row[column])) < 1e-10, (config, column)
                    assert abs(value - float(record[column])) < 1e-10, (config, column)
    def geometry_frames(path: Path) -> dict[str, tuple[list[str], np.ndarray]]:
        lines = iter(path.read_text().splitlines())
        frames = {}
        for count_text in lines:
            count = int(count_text)
            header = next(lines)
            match = re.search(r"\bconfig_id=(\S+)", header)
            assert match is not None
            symbols, positions = [], []
            for _ in range(count):
                fields = next(lines).split()
                symbols.append(fields[0])
                positions.append([float(x) for x in fields[1:4]])
            frames[match.group(1)] = (symbols, np.asarray(positions))
        return frames

    released_geometry = geometry_frames(base / "configurations.extxyz")
    assert sorted(released_geometry) == manifest["configuration_ids"]
    for source_name in ("dissociation", "dissociation_extended"):
        source = ROOT / "experiments" / source_name
        original_geometry = geometry_frames(source / "configurations.extxyz")
        with np.load(source / "solute_probe_points.npz", allow_pickle=False) as source_probes, \
                np.load(base / "solute_probe_points.npz", allow_pickle=False) as released_probes:
            np.testing.assert_array_equal(source_probes["dev_cyclic_carbamate"],
                                          released_probes["dev_cyclic_carbamate"])
        with np.load(source / "predicted_probe_potentials.npz", allow_pickle=False) as source_predictions, \
                np.load(base / "predicted_probe_potentials.npz", allow_pickle=False) as released_predictions:
            for config in manifest["configuration_ids"]:
                if config not in original_geometry:
                    continue
                source_symbols, source_positions = original_geometry[config]
                released_symbols, released_positions = released_geometry[config]
                assert source_symbols == released_symbols
                np.testing.assert_array_equal(source_positions, released_positions)
                for suffix in ("glider", "glider__charges_e", "glider__site_dipoles_e_bohr",
                               "averaged_prior", "averaged_prior__charges_e",
                               "averaged_prior__site_dipoles_e_bohr"):
                    key = f"{config}__{suffix}"
                    np.testing.assert_array_equal(source_predictions[key],
                                                  released_predictions[key])

    cpu_rows = table("experiments/dissociation_qm/validation/cpu_gpu_differences.csv")
    assert len(cpu_rows) == 12
    assert len({row["config_id"] for row in cpu_rows}) == 12
    for row in cpu_rows:
        config = row["config_id"]
        assert config in manifest["configuration_ids"]
        cpu_file = base / "validation/cpu" / f"{config}.npz"
        cpu_record = json.loads(cpu_file.with_suffix(".json").read_text())
        actual_hash = hashlib.sha256(cpu_file.read_bytes()).hexdigest()
        assert actual_hash == row["cpu_arrays_sha256"] == cpu_record["arrays_sha256"]
        assert all(attempts[-1]["converged"] for attempts in cpu_record["attempts"].values())
        with np.load(cpu_file, allow_pickle=False) as cpu, \
                np.load(base / f"{config}.npz", allow_pickle=False) as gpu:
            diff = (cpu["qm_response_esp_hartree_per_e"] -
                    gpu["qm_response_esp_hartree_per_e"]) * 1000
            assert abs(float(np.max(np.abs(diff))) -
                       float(row["max_pointwise_difference_mEh_per_e"])) < 1e-12
            assert abs(float(np.sqrt(np.mean(diff**2))) -
                       float(row["rms_difference_mEh_per_e"])) < 1e-12


def main():
    check_headlines()
    check_separation_qm()
    check_links()
    assert not (ROOT/".github/README.md").exists()
    for name in ['training','panel_1','panel_2','panel_3','nonwater','nonwater_contact','liquid','shell_size','heldout_water','heldout_water_pilot','distance_sweep','global_branch','dissociation','dissociation_extended','dissociation_qm','water_contact_audit']:
        assert (ROOT/'experiments'/name/'README.md').is_file(),name
    print('PASS: experiment guides, active local links, headline values and visual gallery')
if __name__=='__main__':main()

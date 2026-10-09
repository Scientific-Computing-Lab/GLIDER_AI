# Independent CPU–GPU reference check

The published 20-condition reference set was calculated with GPU4PySCF. We also recalculated 12 **distinct** conditions on CPU using PySCF and the same response definition, functional, basis, grid and SCF tolerance. There were 13 CPU runs because the one-water 20 Å condition was repeated; one copy of that condition is retained here.

[Per-condition differences](cpu_gpu_differences.csv) compare the 512-point complex-minus-fragments response potentials in mEh/e. The largest pointwise difference is **0.001120 mEh/e** and the largest per-condition RMS difference is **0.000173 mEh/e**, both at the one-water 3 Å geometry. The 100 Å one-water check agrees at the level of 0.000001 mEh/e pointwise. These comparisons support numerical consistency for the checked conditions; they do not add new solute diversity.

The [`cpu/` directory](cpu/) contains each CPU array file and JSON calculation record. The matching GPU files are one level above. Run `python scripts/verify_companion.py` from the repository root to check record hashes and rederive every difference.

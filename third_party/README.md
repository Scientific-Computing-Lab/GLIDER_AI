# External model weights

GLIDER's response head is released in [`checkpoints/`](../checkpoints/README.md). New-geometry inference also needs the **official** MACE-POLAR-1 M and L checkpoints and a compatible MACE implementation. They are not redistributed here because they come from a separate project and license.

```bash
python scripts/reproduce/fetch_mace_polar.py \
  --output-dir third_party/checkpoints --sizes M L
```

The fetcher downloads from the publisher's [`mace_polar_1` release](https://github.com/ACEsuit/mace-foundations/releases/tag/mace_polar_1) and checks each file's known SHA-256. The study used MACE source commit `91df5a2032b24ff9e23e0dc7b9407dde0da6fb31`. Review and obey the upstream license before using or redistributing those assets. The checkpoint download directory is gitignored. Frozen publication scores can be verified from the included archives without downloading these weights.

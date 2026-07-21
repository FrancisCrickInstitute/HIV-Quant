# HIV-Quant

[![License: GPL-3.0](https://img.shields.io/badge/License-GPL%203.0-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Python 3.14](https://img.shields.io/badge/Python-3.14-green.svg)](https://www.python.org/downloads/)
[![Managed by Pixi](https://img.shields.io/badge/managed%20by-pixi-yellow.svg)](https://pixi.sh)
[![Platforms](https://img.shields.io/badge/platforms-linux--64%20%7C%20win--64-lightgrey.svg)](https://pixi.sh)
![Commit activity](https://img.shields.io/github/commit-activity/y/FrancisCrickInstitute/HIV-Quant?style=plastic)

Bioimage analysis pipeline for quantifying HIV capsid/CPSF6/HA intensity in cell nuclei from multi-channel
confocal z-stacks (`.vsi` files). Nuclei are segmented in 3D from a DAPI channel, then per-nucleus intensity
statistics are measured across the remaining channels, aggregated by experimental condition, and written out as
CSVs, a summary plot, and per-slice label image overlays.

## Requirements

This project uses [pixi](https://pixi.sh) for environment management. All dependencies (numpy, pandas, scipy,
bioio, matplotlib, scikit-image, seaborn, etc.) are declared in `pixi.toml`.

```
pixi install
```

Reading `.vsi` (Olympus) files goes through `bioio-bioformats`, which relies on a JVM. On first run it will
download a JRE via `cjdk` if one isn't already available.

## Usage

```
pixi run python script_20260721_185826.py [--data-dir DATA_DIR]
```

`--data-dir` points at a folder of `.vsi` files and defaults to `./data`. Run with `--help` for the full option
list.

Filenames are expected in the form `<index>_Multichannel Z-Stack_<date>_<n>.vsi`, e.g.
`10_Multichannel Z-Stack_20260622_67.vsi`. The leading index is looked up in `CONDITION_MAPPING` (in the script)
to assign each image to an experimental condition.

## Output

Results are written to `./output/`:

- `nuclei_measurements.csv` — per-nucleus intensity metrics (mean/median/min/max/std/total) for each channel
- `summary_statistics.csv` — per-condition mean/std of each metric
- `intensity_summary.png` — per-channel boxplot with individual nuclei overlaid as a swarm plot
- `label_images/<filename_stem>/z###.png` — one PNG per z-slice per input file, showing the DAPI signal with
  segmented nuclei overlaid, for visually checking segmentation quality

## Configuration

The following are set as constants near the top of the script and should be adjusted per experiment:

- `DAPI_CHANNEL`, `HA_CHANNEL`, `CPSF6_CHANNEL`, `CAPSID_CHANNEL` — channel indices in the acquired image
- `NUCLEI_DIAMETER_PX`, `SIZE_TOLERANCE` — expected nucleus diameter (px) and acceptable size deviation, used to
  filter segmented objects
- `CONDITION_MAPPING` — maps the numeric file index parsed from each filename to an experimental condition label
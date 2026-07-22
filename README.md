# HIV-Quant

[![License: GPL-3.0](https://img.shields.io/badge/License-GPL%203.0-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Python 3.14](https://img.shields.io/badge/Python-3.14-green.svg)](https://www.python.org/downloads/)
[![Managed by Pixi](https://img.shields.io/badge/managed%20by-pixi-yellow.svg)](https://pixi.sh)
[![Platforms](https://img.shields.io/badge/platforms-linux--64%20%7C%20win--64-lightgrey.svg)](https://pixi.sh)
![Commit activity](https://img.shields.io/github/commit-activity/y/FrancisCrickInstitute/HIV-Quant?style=plastic)

A general-purpose bioimage analysis pipeline for segmenting nuclei in 3D and quantifying per-channel intensity
within them, from multi-channel confocal z-stacks (`.vsi` files). Nuclei are segmented from a DAPI channel, then
per-nucleus intensity statistics are measured across the remaining channels, aggregated by experimental condition,
and written out as CSVs, a summary plot, and per-slice label image overlays.

This repo was originally built to quantify HIV capsid/CPSF6/HA intensity in infected-cell nuclei, but the
segmentation and measurement logic is not HIV-specific — see [Adapting to other experiments](#adapting-to-other-experiments)
below to reuse it for any experiment that needs per-nucleus intensity quantification across channels.

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
pixi run python script_20260721_185826.py [OPTIONS]
```

Run with `--help` for the full option list. Available options (defaults match the original HIV experiment this
pipeline was built for):

| Option | Default | Description |
| --- | --- | --- |
| `--data-dir` | `./data` | Directory containing `.vsi` files |
| `--channel-names` | `DAPI HA CPSF6 Capsid` | Ordered list of channel names, one per channel index in the acquired image. Must include `DAPI`, which is used for nucleus segmentation |
| `--nuclei-diameter-px` | `140` | Expected nucleus diameter in pixels, used to filter segmented objects by size |
| `--size-tolerance` | `0.3` | Acceptable fractional deviation from `--nuclei-diameter-px` |
| `--condition-mapping` | built-in HIV-Quant mapping | JSON object mapping the numeric file index parsed from each filename to an experimental condition label, e.g. `'{"1": "ConditionA", "2": "ConditionB"}'` |

Filenames are expected in the form `<index>_Multichannel Z-Stack_<date>_<n>.vsi`, e.g.
`10_Multichannel Z-Stack_20260622_67.vsi`. The leading index is looked up in `--condition-mapping` to assign each
image to an experimental condition.

## Output

Results are written to `./output/`:

- `nuclei_measurements.csv` — per-nucleus intensity metrics (mean/median/min/max/std/total) for each channel
- `summary_statistics.csv` — per-condition mean/std of each metric
- `intensity_summary.png` — per-channel boxplot with individual nuclei overlaid as a swarm plot
- `label_images/<filename_stem>/z###.png` — one PNG per z-slice per input file, showing the DAPI signal with
  segmented nuclei overlaid, for visually checking segmentation quality

## Configuration

Channel names/order, expected nucleus size, and the condition mapping are all CLI options — see the table in
[Usage](#usage) — so no source changes are needed to point the pipeline at a different experiment. The
filename-parsing convention in `get_condition_from_filename` is still fixed in the script; see below if that needs
to change too.

## Adapting to other experiments

Nothing about the segmentation or measurement code is specific to HIV, capsid, CPSF6, or HA — those are just the
default `--channel-names`. To reuse the pipeline for a different multi-channel experiment:

- Pass `--channel-names` with your own stain names in acquisition order (one per channel index), e.g.
  `--channel-names DAPI GFP mCherry`. The name at each position becomes the column/plot label for that channel, and
  whichever position is named `DAPI` (required) is used for nucleus segmentation. Any number of channels is fine.
- Use `--nuclei-diameter-px`/`--size-tolerance` to match your expected nucleus size, or adjust the
  thresholding/morphology steps in `segment_nuclei_3d` if your nuclear stain behaves differently.
- Use `--condition-mapping` to match your own experimental groups, or edit `get_condition_from_filename` if
  conditions aren't identified by a leading numeric index in the filename.

Everything downstream (per-nucleus metrics, per-condition summary, plots, label image overlays) works off those
config values and needs no further changes.
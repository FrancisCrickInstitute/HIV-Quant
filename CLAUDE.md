# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project overview

HIV-Quant is a bioimage analysis pipeline for quantifying HIV capsid/CPSF6/HA intensity in cell nuclei from
multi-channel confocal z-stacks (`.vsi` files, read via `bioio`/`BioImage`). It segments nuclei in 3D from a DAPI
channel, then measures per-nucleus intensity statistics across the remaining channels, aggregates results by
experimental condition, and produces summary CSVs and a plot.

The repository is currently a single analysis script (`script_20260721_185826.py`) plus a `pixi` environment
definition — there is no package structure, test suite, or CI yet.

## Environment management (pixi)

This project uses [pixi](https://pixi.sh), not pip/conda directly. Dependencies are split across:
- `[dependencies]` in `pixi.toml` — conda-forge packages (currently just Python 3.14 and `pixi-pycharm`)
- `[pypi-dependencies]` in `pixi.toml` — PyPI packages (numpy, pandas, scipy, bioio, matplotlib, scikit-image)

Common commands:
```
pixi install              # create/update the `hiv` environment from pixi.toml/pixi.lock
pixi run python <script>  # run a script inside the environment
pixi shell                 # drop into an activated shell for the environment
```
There are no `[tasks]` defined in `pixi.toml` yet — scripts are run directly with `pixi run python ...`.

When adding a new dependency, edit `pixi.toml` (not a requirements.txt) and run `pixi install` so `pixi.lock` is
regenerated; both files should be committed together.

## Architecture of the analysis script

The pipeline (`script_20260721_185826.py`) runs as a linear sequence of stages, all currently in one file:

1. **Configuration constants** at the top of the file define channel indices (`DAPI_CHANNEL`, `HA_CHANNEL`,
   `CPSF6_CHANNEL`, `CAPSID_CHANNEL`), expected nucleus size (`NUCLEI_DIAMETER_PX`, `SIZE_TOLERANCE`), and
   `DATA_DIR`/`OUTPUT_DIR` paths. `CONDITION_MAPPING` maps a numeric index parsed from each filename to an
   experimental condition label (e.g. `D37_RR-VLPs`, `D102-VLPs`, `Uninfected`) via `get_condition_from_filename`.
2. **3D nucleus segmentation** (`segment_nuclei_3d`): normalizes the DAPI z-stack, Gaussian-smooths it, applies
   Otsu thresholding, fills holes, cleans up with binary erosion/dilation, labels connected components in 3D, and
   filters labeled objects by voxel-count bounds derived from `NUCLEI_DIAMETER_PX`/`SIZE_TOLERANCE`.
3. **Per-nucleus intensity extraction** (`extract_intensity_metrics`): for each labeled nucleus and each channel
   (`DAPI`, `HA`, `CPSF6`, `Capsid`), computes mean/median/min/max/std/total intensity over the voxels belonging to
   that nucleus.
4. **Per-file processing** (`process_vsi_file`): loads a `.vsi` file with `BioImage`, derives the DAPI stack,
   segments nuclei, and builds a per-nucleus `DataFrame` tagged with `filename` and `condition`.
5. **Pipeline entry point** (`main`): globs all `*.vsi` files under `DATA_DIR`, processes each one, concatenates
   results into `output/nuclei_measurements.csv`, computes per-condition mean/std summary statistics into
   `output/summary_statistics.csv`, and renders a 2x2 bar chart (one subplot per channel) comparing mean
   intensities across conditions.

Image arrays are expected in `(channels, z, y, x)` order; channel-to-biology mapping is controlled entirely by the
`*_CHANNEL` constants, and condition labeling is controlled entirely by `CONDITION_MAPPING` — both need to be kept
in sync with the actual acquisition/experiment setup when new data is added.

Inputs are read from `./data/*.vsi` and outputs are written to `./output/` (created automatically); neither
directory is committed to the repo.
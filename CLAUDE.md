# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project overview

HIV-Quant is a general-purpose bioimage analysis pipeline for segmenting nuclei in 3D and quantifying per-channel
intensity within them, from multi-channel confocal z-stacks (`.vsi` files, read via `bioio`/`BioImage`). It
segments nuclei in 3D from a DAPI channel, then measures per-nucleus intensity statistics across the remaining
channels, aggregates results by experimental condition, and produces summary CSVs and a plot. It was originally
built to quantify HIV capsid/CPSF6/HA intensity in infected-cell nuclei, but channel names/order, nucleus size, and
the condition mapping are all CLI options (see `parse_args` in the script), so it isn't HIV-specific.

The repository is currently a single analysis script (`quantify_nuclei_intensity.py`) plus a `pixi` environment
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

The pipeline (`quantify_nuclei_intensity.py`) runs as a linear sequence of stages, all currently in one file. Most
per-experiment settings are CLI arguments (parsed in `parse_args`, with module-level constants near the top of the
file — `CHANNEL_NAMES`, `NUCLEI_DIAMETER_PX`, `SIZE_TOLERANCE`, `CONDITION_MAPPING`, `DATA_DIR` — only supplying
their defaults) rather than hardcoded config, and are threaded through as function parameters:

1. **CLI parsing** (`parse_args`): `--data-dir`, `--channel-names` (an ordered list of channel names, one per
   channel index — must include `DAPI`, which is used for segmentation), `--nuclei-diameter-px`,
   `--size-tolerance`, and `--condition-mapping` (a JSON object string mapping the numeric file index parsed from
   each filename to a condition label, parsed by `parse_condition_mapping`).
2. **Pipeline entry point** (`main`): validates `DAPI` is present in `channel_names`, builds a `channels` list of
   `(name, index)` pairs and resolves `dapi_channel` from its position, globs all `*.vsi` files under `data_dir`,
   processes each one, concatenates results into `output/nuclei_measurements.csv`, computes per-condition mean/std
   summary statistics into `output/summary_statistics.csv`, and renders the summary plot.
3. **3D nucleus segmentation** (`segment_nuclei_3d`): Gaussian-smooths the DAPI z-stack (anisotropic sigma
   `[1.0, 2.0, 2.0]`), applies triangle thresholding (`skimage.filters.threshold_triangle` — swapped in from Otsu,
   which was under-segmenting), fills holes, cleans up with binary erosion/dilation, labels connected components in
   3D, and filters labeled objects by voxel-count bounds derived from `nuclei_diameter_px`/`size_tolerance` (via a
   vectorized `np.bincount` + lookup-table relabel, not a per-object rescan).
4. **Label image export** (`save_label_images`): percentile-normalizes the DAPI stack and saves one PNG per
   z-slice (`label_images/<filename_stem>/z###.png`) with segmented nuclei overlaid via `skimage.color.label2rgb`,
   for visually checking segmentation quality.
5. **Per-nucleus intensity extraction** (`extract_intensity_metrics`): for each labeled nucleus and each entry in
   `channels`, computes mean/median/min/max/std/total intensity over that nucleus's voxels in one vectorized pass
   per channel via `scipy.ndimage` (`mean`/`median`/`minimum`/`maximum`/`standard_deviation`/`sum_labels`, each
   given the full array of nucleus IDs at once).
6. **Per-file processing** (`process_vsi_file`): loads a `.vsi` file with `BioImage`
   (`get_image_data("CZYX", T=0)`), derives the DAPI stack from `dapi_channel`, segments nuclei, saves label
   images, extracts per-nucleus metrics, and tags the resulting `DataFrame` with `filename` and a `condition` from
   `get_condition_from_filename`.
7. **Summary plot** (`plot_intensity_summary`): for each channel, normalizes its per-nucleus mean intensity to
   that same nucleus's `DAPI_mean` (so DAPI's own panel is a ~1.0 sanity check), then renders a per-condition
   boxplot+swarmplot (seaborn) with a log-scaled y-axis in a 2x2 grid, since raw/normalized intensities span a wide
   range across channels and conditions.

Image arrays are expected in `(channels, z, y, x)` order; channel-to-biology mapping is controlled entirely by
`--channel-names` (position = channel index), and condition labeling is controlled entirely by
`--condition-mapping` — both need to be kept in sync with the actual acquisition/experiment setup when new data is
added, but neither requires editing the script.

Inputs are read from `./data/*.vsi` (or `--data-dir`) and outputs are written to `./output/` (created
automatically); neither directory is committed to the repo.
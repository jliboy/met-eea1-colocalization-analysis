# MET-EEA1 Colocalization Analysis

This repository contains a pipeline for estimating the colocalization of receptor tyrosine kinase MET with the early endosome marker EEA1 using fluorescence microscopy images. The analysis includes image segmentation, feature extraction, and colocalization quantification, with results saved to CSV and NumPy files.

## Repository structure

- `met_coloc_analysis.py` - main Python analysis script.
- `src/` - helper modules for image segmentation and filtering.
- `data_plotting.rmd`, `data_plotting.md`, `data_plotting.html` - exploratory plotting and analysis using R.
- `coloc_plot.svg` - example visualization generated from the analysis data.

## What this analysis does

`met_coloc_analysis.py` performs the following steps for each selected microscopy image:

1. Opens a directory chooser and loads image files matching `HGF*`.
2. Reads multi-channel image data from ND2 files.
3. Segments nuclei, MET signal, and EEA1 signal.
4. Counts MET and EEA1 foci and normalizes counts by nuclei number.
5. Calculates MET / EEA1 overlap as a colocalization fraction.
6. Extracts MET area and MET intensity per nucleus.
7. Infers treatment groups and timepoints from image filenames.
8. Saves results to a CSV file and masked overlays to NumPy arrays.

## Requirements

- Python 3.x
- `nd2reader`
- `numpy`
- `pandas`
- `matplotlib`
- `scikit-image`
- `tkinter` (for directory dialog)

Optional for plotting and reporting:

- R with `readr`, `dplyr`, `purrr`, `ggplot2`, and `plotrix`.

## Usage

1. Place the script in a Python environment with the required packages installed.
2. Run the script from a Python interpreter:

```bash
python met_coloc_analysis.py
```

3. When prompted, select the directory containing your ND2 image files.
4. Ensure your image filenames identify treatment and timepoint information, for example:
   - `HGF_sNRP_15min.nd2`
   - `HGF_noligand_30min.nd2`
   - `HGF_15min.nd2`

## Input data expectations

- Image files should be named with the prefix `HGF` or match the pattern used in `glob`.
- The script assumes the channels are arranged as:
  - channel 0 = EEA1
  - channel 1 = nuclei
  - channel 2 = MET

Adjust the filename pattern and channel mapping if your dataset uses a different convention.

## Output files

When `met_coloc_analysis.py` finishes, it writes:

- `met_eea1_coloc.csv` - tabulated colocalization and feature data.
- `masked_eea1.npy` - NumPy array of EEA1 mask overlays.
- `masked_met.npy` - NumPy array of MET mask overlays.

## Notes

- The current script references a local segmentation helper path; make sure the `nuclear_segmentation` helper functions are available in your Python path.
- The R notebook and exported HTML document in this repository demonstrate example plotting and statistical comparison of the generated results.

## Suggested next steps

- Inspect `data_plotting.html` to review example visualizations.
- Adapt the image selection pattern or filename parsing for your experiment naming scheme.
- Extend the analysis to export additional colocalization metrics or per-cell measurements.



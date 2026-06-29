# MTB Onkoprint

## Project description

Functions and scripts for preparing molecular tumor board (MTB) data and generating customized OncoPrint visualizations, as demonstrated in:

-   Figure 2 in the manuscript by [Klostermann *et al.* 2025](https://doi.org/10.1016/j.esmoop.2025.105497).

This repository includes data preparation scripts in Python or R, and a python-based workflow to create an OncoPrint figure with multiple annotations per gene.

## Table of Contents

-   [Installation instructions](#installation-instructions)
-   [Workflow overview](#workflow-overview)
-   [Notes on sensitive data](#notes-on-sensitive-data)
-   [Build your own Onkoprint](#build-your-own-onkoprint)
-   [License & Citation](#license)

## Installation instructions

### Project Requirements and Dependencies

The MTB-Onkoprint uses the [PyOncoPrint package](https://doi.org/10.5808/gi.22079). For this, Python is required and all packages in [requirements.txt](requirements.txt).

Data preparation can be performed using Python or R, or performed manually.

### Installation

Clone the repository and create a virtual environment with packages specified in the requirements.txt file:

``` bash
git clone https://github.com/BIH-CEI/mtb-onkoprint.git
cd mtb-onkoprint
# Unix/macOS
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements.txt

# Windows using Git Bash
py -m venv .venv
source .venv/Scripts/activate
py -m pip install -r requirements.txt
```

Activate the virtual environment:

More information on creating virtual environments from a requirements.txt file can be found on the [python](https://packaging.python.org/en/latest/guides/installing-using-pip-and-virtual-environments/#using-a-requirements-file) website.

## Workflow overview

The workflow has two main steps.

1.  **Data preparation**

    The data-preparation script reads cohort-specific MTB source data from `data/`, cleans values, standardizes gene/modality annotations, splits or detect modalities, validates marker values, and writes intermediate CSV files to `output/`.

    Main outputs include:

    ``` text
    output/<DATE>_<COHORT>_clinical.csv
    output/<DATE>_<COHORT>_ihc_heatmap.csv
    output/<DATE>_<COHORT>_complex.csv
    output/<DATE>_gene-data-single-rows.csv
    output/<DATE>_gene-data-aggregated-rows.csv
    ```

2.  **Visualization**

    The visualization script reads the intermediate CSV files, applies cohort-specific annotation settings, constructs the OnkoPrint matrix, add clinical annotations, renders the OnkoPrint, and saves the figure as PNG and SVG.

    Main outputs include:

    ``` text
    output/<DATE>_<COHORT>_oncoprint.png
    output/<DATE>_<COHORT>_oncoprint.svg
    ```

### Expected input

Input files are expected in a local `data/` directory. These files are not included in the repository, but example (synthetic, not real) data and expected outputs will be added soon (example_data.xlsx).

The data should be as follows:

-   each row contains one variable, and each column is a patient.
-   Sequence variants (gene signatures, fusions, SNVs, Deletions, and so on) need to be sorted per gene. The naming convention we used is `<gene><underscore><alteration type>`, for example `BRAF_SNV`.
-   Values of clinical metadata can be factors (for example `Sex (M/F)` with values `M` and `F`) or numeric.
-   For sequence variants, one should either specify: `assessed`, `not assessed`, or a value. `assessed` means that a method was performed for the gene to study if it was altered, and the result is "not altered". `not assessed` means that the method was not performed, meaning that we do not know whether there is an alteration. If a SNV was found, or a deletion, or a quantification of immunohistochemistry, or something else, one can either specify a numeric `value`, for example `HER2_IHC` has value 2, a specific mutation (ideally using the [HGVS Nomenclature](https://hgvs-nomenclature.org/stable/)) or just the gene name. The difference between `assessed` and `not assessed` can then be visualized in the OnkoPrint with different greyscale colors.

Code for other data input formats is currently under development.

### Using Jupyter Notebooks

The python scripts are [Jupytext](https://jupytext.org/)-compatible Python notebooks. To convert a script to an `.ipynb` notebook:

``` bash
jupytext --to ipynb 01_data_preparation.py
jupytext --to ipynb 02_visualization.py
```

For the CUP workflow:

``` bash
jupytext --to ipynb 01_data_preparation_CUP.py
jupytext --to ipynb 02_visualization_CUP.py
```

Alternatively, open the `.py` files directly in a Jupyter/Jupytext-enabled editor.

## Notes on sensitive data

The repository is designed to operate on local cohort data stored in `data/`. Do not commit patient-level source files, generated identifiable outputs, or local environment files.

## Build your own Onkoprint

We aim to write a [Tutorial](docs/Tutorial.md) with more details on preparing your own data.

## License & Citation

The license can be found here: [License](LICENSE).

If you use our code in your work, please cite the [Klostermann *et al.* 2025](https://doi.org/10.1016/j.esmoop.2025.105497) manuscript, and also consider citing:

-   cBioPortal, for designing OncoPrints, as described on [their website](https://docs.cbioportal.org/user-guide/faq/#how-do-i-cite-the-cbioportal).
-   The PyOncoPrint package by [Park & Paramasivam](https://doi.org/10.5808/gi.22079)

## 

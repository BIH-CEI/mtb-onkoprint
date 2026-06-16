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

## Installation instructions {#installation-instructions}

### Project Requirements and Dependencies

The MTB-Onkoprint uses the [PyOncoPrint package](#0). For this, Python is required and all packages in [requirements.txt](requirements.txt).

Data preparation can be performed using Python or R, or performed manually.

### Installation

Clone the repository:

``` bash
git clone https://github.com/BIH-CEI/mtb-onkoprint.git
cd mtb-onkoprint
```

*To be written*

## Workflow overview {#workflow-overview}

The workflow has two main steps.

1.  **Data preparation**

    The data-preparation scripts read cohort-specific MTB source data from `data/`, clean values, standardize gene/modality annotations, split or detect modalities, validate marker values, and write intermediate CSV files to `output/`.

    Main outputs include:

    ``` text
    output/<DATE>_<COHORT>_Entities.csv
    output/<DATE>_<COHORT>_Therapy.csv
    output/<DATE>_<COHORT>_complex.csv
    output/<DATE>_gene-data-single-rows.csv
    output/<DATE>_gene-data-aggregated-rows.csv
    ```

2.  **Visualization**

    The visualization scripts read the intermediate CSV files, apply cohort-specific annotation settings, construct the OncoPrint matrix, add clinical annotations, render the OncoPrint, and save the figure as PNG and SVG.

    Main outputs include:

    ``` text
    output/<DATE>_<COHORT>_oncoprint.png
    output/<DATE>_<COHORT>_oncoprint.svg
    ```

### Expected input

Input files are expected in a local `data/` directory. These files are not included in the repository, but example (synthetic, not real) data and expected outputs are. (example_data_CUP.xlsx)

### Using Jupyter Notebooks

The scripts are Jupytext-compatible Python notebooks. To convert a script to an `.ipynb` notebook:

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

## Notes on sensitive data {#notes-on-sensitive-data}

The repository is designed to operate on local cohort data stored in `data/`. Do not commit patient-level source files, generated identifiable outputs, or local environment files.

## Build your own Onkoprint {#build-your-own-onkoprint}

See the to-be-written [Tutorial](docs/Tutorial.md) for more details on preparing your own data.

## License & Citation {#license}

The license can be found here: [License](LICENSE).

If you use our code in your work, please cite the [Klostermann *et al.* 2025](https://doi.org/10.1016/j.esmoop.2025.105497) manuscript, and also consider citing:

-   cBioPortal, for designing OncoPrints, as described on [their website](https://docs.cbioportal.org/user-guide/faq/#how-do-i-cite-the-cbioportal).
-   The PyOncoPrint package by [Park & Paramasivam](https://doi.org/10.5808/gi.22079)

## 

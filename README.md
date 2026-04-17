<img src="./docs/assets/ocelot-logo.png" width="300px">

# OCELOT

The **O**ptimization **C**oupling **E**xpression **L**andscape **o**f **T**ranscription Factors (OCELOT) integrates gene regulatory networks (GRNs) with genome-scale metabolic models (GEMs) to support two main tasks:

1. Transcription factor knockout simulation, optionally followed by essentiality evaluation
1. Metabolic engineering, including strategy search and post hoc validation

## Installation

### Dependencies:

* MATLAB (tested with version 2025a)
* [COBRA Toolbox 3](https://opencobra.github.io/cobratoolbox/stable/index.html)
* [RAVEN Toolbox 2](https://github.com/SysBioChalmers/RAVEN/releases/tag/v2.8.3)
* Gurobi solver (tested with version 12.0.1)


### How to install:

* Clone the repository in the desired directory or download the [the latest release](https://github.com/mauricioamf/OCELOT/releases) as a ZIP file. Add folder to MATLAB path before running the code.

## How to use

The [`Protocol.m`](https://github.com/mauricioamf/OCELOT/blob/main/Protocol.m) file contains step-by-step instructions on how to run OCELOT using example data. Please refer to the documentation on the `docs/` folder for detailed instructions.

The folder `manuscript/`contains the code and data used in the manuscript for generating the presented results.

## Citation

If you use OCELOT, please cite:

> Ferreira et al. _Rational design of metabolic engineering strategies based on modulating the expression of transcription factors_. bioRxiv (2026). [DOI: pending](https://doi.org/)

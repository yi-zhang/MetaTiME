# MetaTiME: Meta-components in Tumor immune MicroEnvironment 
<p align="left"><img src="https://raw.githubusercontent.com/yi-zhang/MetaTiME/main/docs/source/_static/img/logo.png" width="290" height="240"></p>

[![Latest PyPI Version][pb]][pypi]
[![Documentation Status](https://readthedocs.org/projects/metatime/badge/?version=latest)](https://metatime.readthedocs.io/en/latest/?badge=latest)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7410180.svg)](https://doi.org/10.5281/zenodo.7410180)

[pb]: https://img.shields.io/pypi/v/metatime.svg
[pypi]: https://pypi.org/project/metatime/

MetaTiME is a framework to generate data-driven, interpretable, and reproducible gene programs by integrating millions of single cells from hundreds of tumor scRNA-seq data. Applied on large-scale tumor scRNA data with 1.7 million Tcells, MetaTiME thus utilize the meta-components to automatically annotate cell states for single-cells from tumor microenvironment. 
( * BETA version, Currently under development  :)

## Installation

Create a new virtual env and activate (optional)

`python -m venv metatime-env; 
source metatime-env/bin/activate`

Use pip to install

`pip install metatime`

Installation shall be in minutes .
## Interactive tutorials
### MetaTiME-Annotator
[Use MetaTiME to automatically annotate cell states and map signatures ![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/yi-zhang/MetaTiME/blob/main/docs/notebooks/metatime_annotator.ipynb)

## Method 
<p align="left"><img src="https://raw.githubusercontent.com/yi-zhang/MetaTiME/main/docs/source/_static/img/fig1.png" width="700" height="400"></p>


### Usage
 - [Use MetaTiME to automatically annotate cell states and map signatures](https://github.com/yi-zhang/MetaTiME/blob/main/docs/notebooks/metatime_annotator.ipynb)
 
### Dependency

- pandas
- scanpy
- anndata
- matplotlib
- adjustText
- leidenalg
- harmonypy

Dependency version tested:
- pandas==1.1.5
- scanpy==1.8.2
- anndata==0.8.0
- matplotlib==3.5.1
- adjustText==0.7.3
- leidenalg==0.8.3


### Reference
Manuscript In Revision. Repo continously being improved! More details will be updated. 

[Paper at bioRxiv](https://www.biorxiv.org/content/10.1101/2022.08.05.502989v1)

[Journal Article doi pending]

### Training Datasets

Tumor scRNAseq Data for MetaTiME @ [Zenodo](https://zenodo.org/record/7410180)

- A large collection of uniformly processed tumor single-cell RNA-seq. 

- Includes raw data and MetaTiME score for the TME cells.

### Contact


Yi Zhang, Ph.D.

yiz [AT] ds.dfci.harvard.edu

[Twitter](https://twitter.com/Wings7Spread) |  [Website](https://yi-zhang.github.io/)

Research Fellow

Department of Data Science

Dana-Farber Cancer Institute

Harvard University T.H. Chan School of Public Health




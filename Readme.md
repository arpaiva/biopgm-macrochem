# Probabilistic Macrochemical Modelling

<div align="center">
  <a href="https://www.sci.utah.edu/~arpaiva/" target="_blank">Antonio&nbsp;Paiva</a> &emsp; <b>&middot;</b> &emsp;
  <a href="https://loop.frontiersin.org/people/797425/overview" target="_blank">Giovanni&nbsp;Pilloni</a> 
</div>
<br>

Probabilistic macrochemical modeling is a general methodology that builds on 
macrochemical characterizations of microbial growth ([Heijen and 
Van&nbsp;Dijken, 1992](https://pubmed.ncbi.nlm.nih.gov/18601018/)) and
fully utilizes all of available experimental data at one for accurate and
robust estimation of parameters of interest (e.g., biomass yields).

Specifically, this repository contains the code used to generate the
validation data and perform the analysis in the article,

* Antonio R. Paiva and Giovanni Pilloni (2021). Inferring Microbial Biomass Yield and Cell Weight using Probabilistic Macrochemical Modeling. [arXiv preprint arXiv:2010.02759](https://arxiv.org/abs/2010.02759)

## Requirements
The code has been tested in Python 3.7 and requires PyStan 2.19 or newer. The 
Jupyter notebooks have been tested version 3.6.

## Overview

The repository comprises 4 main files:

* `srbsim.py` implements the class used to generate the simulated used in the 
  data. This is called in each of the Jupyter notebooks.

* `chem_model1_constCellWeight.ipynb` implements all analysis pertaining to
  Schenario 1 described in the paper.

* `chem_model2_changingCellWeights.ipynb` implements all analysis pertaining to
  Schenario 2 described in the paper.

* `chem_model3_changingCellWeights_biomassSensitivity.ipynb` tests, with
  respect to Scenario 2, the sensitivity of the results to the generic
  biomass composition formula used in the paper. This analysis is described
  in the supplementary materials, Section S1.

## License
This code is being shared under an [MIT 
license](https://github.com/arpaiva/biopgm-macrochem/blob/main/LICENSE).

# Probabilistic Macrochemical Modelling

<div align="center">
  <a href="https://www.sci.utah.edu/~arpaiva/" target="_blank">Antonio&nbsp;Paiva</a> &emsp; <b>&middot;</b> &emsp;
  <a href="https://loop.frontiersin.org/people/797425/overview" target="_blank">Giovanni&nbsp;Pilloni</a> 
</div>
<br>

Probabilistic macrochemical modeling is a general methodology that builds
on macrochemical characterizations of microbial growth ([Heijen and 
Van&nbsp;Dijken, 1992](https://pubmed.ncbi.nlm.nih.gov/18601018/)).
With it, we can build probabilistic models that infer from all available
experimental data globally consistent and accurate estimates of microbial
quantities of interest (e.g., biomass yields).

Specifically, this repository contains the code used to generate the
validation data and perform the analysis in the article,

* Antonio R. Paiva and Giovanni Pilloni (2022). Inferring Microbial Biomass Yield and Cell Weight using Probabilistic Macrochemical Modeling. _IEEE/ACM Transactions on Computational Biology and Bioinformatics_<br/>[DOI: 10.1109/TCBB.2021.3139290](https://doi.org/10.1109/TCBB.2021.3139290)  |  [arXiv preprint arXiv:2010.02759](https://arxiv.org/abs/2010.02759)

## Requirements

The code has been tested in Python 3.7 and PyStan 2.19. Note that, while newer
versions of Python 3 should work, PyStan 3 is not backwards compatible and
would require a few changes to the code.

## Overview

The repository comprises 4 main files:

* `srbsim.py` implements the class used to generate the simulated data used in the 
  paper. This is called in each of the Jupyter notebooks.

* `chem_model1_constCellWeight.ipynb` implements all analysis pertaining to
  Scenario 1 described in the paper.

* `chem_model2_changingCellWeights.ipynb` implements all analysis pertaining to
  Scenario 2 described in the paper.

* `chem_model3_changingCellWeights_biomassSensitivity.ipynb` tests, with
  respect to Scenario 2, the sensitivity of the results to the generic
  biomass composition formula used in the paper. This analysis is described
  in the paper supplementary materials, Section S1.

## License
This code is being shared under an [MIT 
license](https://github.com/arpaiva/biopgm-macrochem/blob/main/LICENSE).

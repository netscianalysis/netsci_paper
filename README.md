netsci_paper
==============================
Data and calculations used for the Netsci paper:

Stokely A, Votapka L, Hock M, Teitgen A, McCammon JA, McCullough A, Amaro R, NetSci: A Library for High Performance Biomolecular Simulation Network Analysis Computation. ChemRxiv. 2024; doi:10.26434/chemrxiv-2024-fjrpp This content is a preprint and has not been peer-reviewed.

## Installation

Make sure that Netsci is installed and tested: 
https://github.com/netscianalysis/netsci.git

Activate the Netsci Conda environment
```
conda activate netsci
```

The dependency(ies) needed for this repository must be installed manually:
```
pip install matplotlib
conda install -c conda-forge ambertools
```

## Analyses
Instructions for how to run the analyses used in the paper can be found in the READMEs
located in the various directories:

```
netsci_paper/
  toy_systems/
    README
  proteinG/
    README
  serca/
    README
```

# Citing NetSci

If you use NetSci, or any of the data within this repository, please cite the following paper:

* Stokely A, Votapka L, Hock M, Teitgen A, McCammon JA, McCullough A, Amaro R, NetSci: A Library for High Performance Biomolecular Simulation Network Analysis Computation. ChemRxiv. 2024; doi:10.26434/chemrxiv-2024-fjrpp This content is a preprint and has not been peer-reviewed.

### Copyright

Copyright (c) 2024, Lane Votapka and Andy Stokely

#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.1.

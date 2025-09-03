# A Theory of Ecological Invasions and Its Implications for Eco-Evolutionary Dynamics

## Overview

This repository contains the code related to our recent work on ecological invasions and their implications for eco-evolutionary dynamics, published on arXiv. 

## Authors

- **Zhijie Feng**
- **Emmy Blumenthal**
- **Pankaj Mehta**
- **Akshit Goyal**

## Python Script
This repository includes a Python script, `iterative_perturbation.py`**, which is imported in all of the Jupyter notebooks. It contains the key function that implements the iterative algorithm proposed in the paper for predicting the outcomes of invasions in general ecological models.
## Jupyter Notebooks

This repository includes four Jupyter Notebooks, each produces a separate result in our work:

1. **Notebook 1: `general_LV_for_figure.ipynb`**
   - produces the simulation and prediction for general Lotka-Volterra models.

2. **Notebook 2: `general_LV_for_data.ipynb`**
   - performs data analysis and prediction for experimental ecosystems fitted with general Lotka-Volterra models.

3. **Notebook 3: `general_MCRM.ipynb`**
   - produces the simulation and prediction for general MacArthur Consumer Resource models.

4. **Notebook 4: `MiCRM_and_evolution.ipynb`**
   - produces the simulation and prediction for Microbial Consumer Resource models, and also perform repeated invasion simulation to study evolved ecosystems.

5. **Notebook 4: `monod_works.ipynb`**
   - produces the simulation and prediction for Consumer Resource models with type II (Monod) growth rate.

6. **Notebook 5: `multistability_general_LV_for_figure.ipynb`**
   - produces the simulation and prediction for general Lotka-Volterra models with multiple stable fixed points.


## Data
The `data` folder in this repository is a direct copy from [DS Maynard's GitHub repository](https://github.com/dsmaynard/endpoints), and the `data` folder comes from [StefanoAllesina's GitHub repository](https://github.com/dsmaynard/endpoints](https://github.com/StefanoAllesina/lemos-costa-2023). They are included here solely to ensure compatibility with one of the Jupyter Notebooks, `general_LV_for_data.ipynb`, which relies on this data for reanalysis. 

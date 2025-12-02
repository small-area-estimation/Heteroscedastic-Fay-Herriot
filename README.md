# HFH: Heteroscedastic Fay-Herriot

This repository contains R codes for the paper:

**"Small area estimation of poverty proportions under heteroscedastic Fay-Herriot models"**

by Esteban Cabello, María Dolores Esteban, Domingo Morales, and Agustín Pérez.

---

## 📄 Overview

The HFH model is a modification of the standard Fay-Herriot model. Unlike the standard Fay-Herriot model, which assumes homoscedasticity for random effects, the HFH model explicitly models the random effect variance as a function of auxiliary variables using Box-Cox transformations. 
This approach allows for greater flexibility in capturing unobserved heterogeneity across domains.

This repository includes:

- R functions to:
  - Fit the HFH model and obtain the model parameter estimates and the EBLUPs
  - MSE estimators (analytical & bootstrap-based)
  - Example of use
  
- CSV with synthetic data

---



"Small area estimation of poverty proportions under heteroscedastic Fay-Herriot models"

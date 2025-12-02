# HFH: Heteroscedastic Fay-Herriot

This repository contains R codes for the paper:

> **Cabello, E., Esteban, M. D., Morales, D., & Pérez, A. (2025). Small area estimation of poverty proportions under heteroscedastic Fay-Herriot models. *Journal of Applied Statistics*, 1–23. [https://doi.org/10.1080/02664763.2025.2568679](https://doi.org/10.1080/02664763.2025.2568679)**
---

## 📄 Overview

The HFH model is a modification of the standard Fay-Herriot model. Unlike the standard Fay-Herriot model, which assumes homoscedasticity for random effects, the HFH model models the random effect variance as a function of auxiliary variables using Box-Cox transformations. 
This approach allows for greater flexibility in capturing unobserved heterogeneity across domains.

This repository includes:

- R functions to:
  - Fit the HFH model and obtain the model parameter estimates and the EBLUPs
  - MSE estimators (analytical & bootstrap-based)
  - Example of use
  
- CSV with synthetic data

---


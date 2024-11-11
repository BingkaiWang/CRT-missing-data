# Code for "Handling incomplete outcomes and covariates in cluster-randomized trials: doubly-robust estimation, efficiency considerations, and sensitivity analysis"

## Paper

The technical report is avaialble at https://arxiv.org/abs/2401.11278.

## Data

### Abstract 

To demonstrate our theoretical results, we performed simulation studies and re-analyzed a completed cluster-randomized experiment: The Work, Family, and Health Study (WFHS)

### Availability 
The WFHS data are publicly available at https://www.icpsr.umich.edu/web/DSDR/studies/36158.

## Code

### Abstract
The R code here can be used to reproduce all of the simulations and data applications.

### Description 
The folder `simulations` contains R code for our simulations.

 - `simulation-1-1.R` contains the code to reproduce Table 1 without within-cluster sampling. 
 - `simulation-1-2.R` contains the code to reproduce Table 1 with uniform within-cluster sampling. 
 - `simulation-2.R` contains the code to reproduce Figure 2 regarding sensitivity analysis. 
 - `simulation-1-1-tmle.R` and `simulation-1-2-tmle.R` contain the code for simulation 1 with two-stage TMLE. 
 
 The folder `data-analysis` contains R code for our data applications
 
 - `data-analysis-WFHS.R` contains the code for re-analyzing the WFHS data set.

### Reproducibility 
Tables and Figures in the main text can be reproduced using the R scripts.

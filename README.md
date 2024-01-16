# Improvements on Scalable Stochastic Bayesian Inference Methods for Multivariate Hawkes Process


## remarks

The scripts contain two parts:
- The algorithms (script files that start with "MHP")
- Observed likelihood calculation (script files that start with "calc")

## contents

- Folder 'limited_time' contains scripts for analysis in the subparagraph `Optimal subsampling ratios', Section 4.2.
- Folder 'estimation' contains scripts for analysis in the subparagraph 'Estimation accuracy', Section 4.2 and the additional methods (random walk MCMC and boundary-corrected SGLD) mentioned in subparagraph 'Sensitivity analysis for boundary approximation'
- Folder 'sensitivity/data_size' contains scripts for analysis in the subparagraph 'Sensitivity analysis for different dataset sizes', Section 4.2.
- Folder 'sensitivity/sgd_pars' contains scripts for analysis in the subparagraph 'Sensitivity analysis for the stochastic search parameters', Section 4.2.
- Folder 'sensitivity/threshold' contains scripts for analysis in the subparagraph `Sensitivity analysis for the threshold values of the boundary-corrected methods', Section 4.2.
- Folder 'application' contains scripts for analysis of the application example in Section 5.
- Folder 'dataset_intraday' contains data for the application example, and code used to generate the plots
- Folder 'goodness-of-fit' contain scripts implementing the goodness-of-fit analysis for the real-world application in Section 5.

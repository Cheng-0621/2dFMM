# Modeling Repeatedly-measured Functional Data with Fast Two-dimensional Fixed-effects Inference

## Description
We propose an innovative bivariate functional mixed model for subjects' profiles which are continuously recorded in a dense and regular grid. Such repeatedly-measure functional data symmetrically evolves along ''spatial'' and ''temporal'' directions in a dense and regular grid. Tailoring to the spatio-temporal correlation, we estimate the two-dimensional fixed effect for enhancing interpretability and obtaining pointwise confidence surfaces. The fast three-stage estimation procedure improves computational efficiency, especially when encountering large sample size and ultrahigh-dimensional functional data. Intensive simulation studies are conducted and real data application systematically demonstrates the intra-day and inter-day dynamic associations between physical activity and mental health outcomes in adolescents. Our proposed methods address current theoretical gaps and explore new ways of monitoring and early warning of health in large populations through wearable devices.


## Main Function
The main function of the two-dimensional fixed effect inference is: 

* `fmm2d`: Estimates bivariate coefficient functions and confidence surfaces.

```
fmm2d(formula, data, S, smoother = "sandwich", knots = NULL, fpca.opt = list(dataType = 'DenseWithMV', methodSelectK = 'FVE'), parallel = FALSE)
``` 

### Arguments 
*`formula`: formula two-sided formula object in lm() format, except that the response is a matrix.
*`data`:  dataframe containing variables in formula.
*`S`: number of longitudinal grids.
*`smoother`: moother sandwich smoother (sandwich) or tensor product smoother (te).
*`fpca.opt `: list of options control parameters specified by list, dataType: dense and regular (Dense); Very densely and regularly ' observed data: empirical mean and Densely recorded but irregular design, or contaminated with error: pre-smoothing for individual curves (DenseWithMV); Sparse random design (Sparse).
*`parallel`: parallel whether to run parallel computing (True only for Linux/Mac users).

### Values
* `betaHat`: The list of estimated bivariate coefficient functions.
*`betaHat.cov`: The three-way array of estimated squared standard error. 



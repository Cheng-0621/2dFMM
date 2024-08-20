# Modeling Repeatedly-measured Functional Data with Fast Two-dimensional Fixed-effects Inference


## Description
We propose an innovative bivariate functional mixed model for subjects' profiles which are continuously recorded in a dense and regular grid. Such repeatedly-measure functional data symmetrically evolves along ''spatial'' and ''temporal'' directions in a dense and regular grid. Tailoring to the spatio-temporal correlation, we estimate the two-dimensional fixed effect for enhancing interpretability and obtaining pointwise confidence surfaces. The fast three-stage estimation procedure improves computational efficiency, especially when encountering large sample sizes and ultrahigh-dimensional functional data. Intensive simulation studies are conducted and real data application systematically demonstrates the intra-day and inter-day dynamic associations between physical activity and mental health outcomes in adolescents. Our proposed methods address current theoretical gaps and explore new ways of monitoring and early warning of health in large populations through wearable devices.


## Main Function
The main function of the two-dimensional fixed effect inference is: 

* `fmm2d`: Estimates bivariate coefficient functions and confidence surfaces.

```
fmm2d(formula, data, S, smoother = "sandwich", knots = NULL,  pcb = TRUE, scb = FALSE, parallel = FALSE, silence = FALSE)
``` 

### Arguments 
* `formula`: formula two-sided formula object in lm() format, except that the response is a matrix. Noted that a covariate with $m$ categories ($m > 2$) is necessary to be transformed to $m-1$ indicator variable before fitting the model.
* `data`:  dataframe containing variables in the formula and ID column.
* `S`: number of longitudinal grids.
* `smoother`: moother sandwich smoother (sandwich) or tensor product smoother (te).
* `bootstrap.opt` bootstrapping options, B is boostrap sample size and M is for obtaining qn 
* `bandwidth.control` correlation correction and adjustment for confidence bands bandwidth control 
* `fpca.opt `: list of options control parameters specified by list, dataType: dense and regular (Dense); Very densely and regularly ' observed data: empirical mean and Densely recorded but irregular design, or contaminated with error: pre-smoothing for individual curves (DenseWithMV); Sparse random design (Sparse).
* `pcb`: whether to obtain pointwise confidence bands, default is TRUE.
* `scb`: whether to obtain simultaneous confidence bands, default is FALSE.
* `parallel`: parallel whether to run parallel computing (TRUE only for Linux/Mac users).
* `silence`: whether to show descriptions of each step.

### Values
* `betaHat`: The list of estimated bivariate coefficient functions.
* `betaHat.cov`: The three-way array of estimated squared standard error. 
* `qn`: A parameter used to construct simultaneous confidence bands.

## Examples

The bivariate functional model is simulated as follows:

$$Y_{i}(s,t) = \beta_{0}(s,t) + X_{i}(s)\beta_{1}(s,t) + \gamma_{i0}(t) + z_{i}(s)\gamma_{i1}(t) + \epsilon_{i}(s, t), \ \ s, t \in [0, 1].$$

The fixed-effect covariates are generated from $X_{i}(s) \sim N(0, 4)$ and $z_{i}(s) = 6(s-0.5)^{2} + N(0,\rho^{2})$, where $\rho$ represents how noisy the signal of repeatedly-measured visits is. For bivariate coefficient functions, the intercept function is  $\beta_{0}(s,t) = 3\sin(\pi(s+0.5)^{2})\cos(\pi t+0.5) + 1$, we take account of two different types (S1 and S2) as shown below, where S1 presents continuous non-differentiable bivariate function with local zero regions (sparse) and S2 presents smooth bivariate function. The random effects are simulated as $\gamma_{ik}(t) = a_{i1}\phi_{1}(t) + a_{i2}\phi_{2}(t)$, $k=0,1$. We use the scaled orthonormal functions

$$
\begin{bmatrix}
\phi_{1}(t)\\
\phi_{2}(t)
\end{bmatrix}
\propto
\begin{cases}
\begin{bmatrix}
1.5 - \sin(2\pi t) - \cos(2\pi t) & \sin(4\pi t)
\end{bmatrix}^{T} & \text{if } k=0 \\
\begin{bmatrix}
\cos(2\pi t) & \sin(2\pi t) 
\end{bmatrix}^{T}  & \text{if } k=1
\end{cases}
$$

to capture the subject-level fluctuations. The random coefficients are generated from $a_{i1} \sim N(0, 2\sigma_{B}^{2})$ and $a_{i2} \sim N(0, \sigma_{B}^{2})$ respectively and $\sigma_{B}^{2}$ depends on the relative importance of random effect. The measurement error $\epsilon_{ij}(s,t) \sim N(0, \sigma_{\epsilon}^{2})$, where $\sigma_{\epsilon}^{2}$ depends on signal-to-noise ratio $SNR_{\epsilon}$. Here $SNR_{B}$ is defined as the ratio of the standard deviation of fixed-effect and random-effect surfaces, while $SNR_{\epsilon}$ is the ratio of the standard deviation of all liner predictors and that of the measurement errors. We set $SNR_{B} = SNR_{\epsilon} = 1$.


<p float="left">
  <img src="https://github.com/Cheng-0621/2DFMM/blob/main/figures/3Dbeta_trueS1.jpeg" width="400" /> 
  <img src="https://github.com/Cheng-0621/2DFMM/blob/main/figures/3Dbeta_trueS2.jpeg" width="400" />
</p>


We can fit the bivariate functional mixed model. 

```  
S <- 10 #number of the longitudinal grids
T <- 100 #number of the functional grids

fit_S1 <- fmm2d(formula=Y~X, data=data, S=S, smoother="sandwich", knots=c(S-3, min(round(T/4), 35)),  parallel = TRUE)
 
fit_S2 <- fmm2d(formula=Y~X, data=data, S=S, smoother="te", knots=c(max(round(S/4),4), min(round(T/4), 35)), parallel = TRUE)
```

The estimates to bivariate coefficient functions under S1 and S2 are respectively shown as below.

<p float="left">
  <img src="https://github.com/Cheng-0621/2DFMM/blob/main/figures/3Dbeta_estS1.jpeg" width="400" /> 
  <img src="https://github.com/Cheng-0621/2DFMM/blob/main/figures/3Dbeta_estS2.jpeg" width="400" />
</p>

## Files 
* `2DFMM.R`: The main algorithm of our proposed bivariate functional mixed model.
* `simu_generate.R`: The procedures to generate simulation data.
* `demo_simu.R`: A demo script for a simulation study of cases S1 and S2.
* `example.RData`: a subset ($N=200$) of Shanghai adolescent physical activity data and their demographic and mental health outcomes. 
* `demo_app.R`: A demo script for real data application. 

## Authors
Cheng Cao, Jiguo Cao, Hao Pan, Yunting Zhang, Fan Jiang, Xinyue Li 


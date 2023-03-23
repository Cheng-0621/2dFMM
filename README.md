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
* `formula`: formula two-sided formula object in lm() format, except that the response is a matrix.
* `data`:  dataframe containing variables in formula.
* `S`: number of longitudinal grids.
* `smoother`: moother sandwich smoother (sandwich) or tensor product smoother (te).
* `fpca.opt `: list of options control parameters specified by list, dataType: dense and regular (Dense); Very densely and regularly ' observed data: empirical mean and Densely recorded but irregular design, or contaminated with error: pre-smoothing for individual curves (DenseWithMV); Sparse random design (Sparse).
* `parallel`: parallel whether to run parallel computing (True only for Linux/Mac users).

### Values
* `betaHat`: The list of estimated bivariate coefficient functions.
* `betaHat.cov`: The three-way array of estimated squared standard error. 

## Examples

The bivariate functional model is simulated as follows:

$$
    Y_{i}(s,t) = \beta_{0}(s,t) + X_{i}(s)\beta_{1}(s,t) + \gamma_{i0}(t) + z_{i}(s)\gamma_{i1}(t) + \epsilon_{i}(s, t), \ \ s, t \in [0, 1].
$$
The fixed effect covariates are generated from $X_{i}(s) \sim N(0, 4)$ and $z_{i}(s) = 6(s-0.5)^{2} + N(0,\rho^{2})$, where $\rho$ represents how noisy the signal of repeatedly-measured visits is. For bivariate coefficient functions, we take account of two different types (S1 and S2) as shown in Figure \ref{fig:true_beta}, where S1 presents continuous non-differentiable bivariate function with local zero regions and S2 presents smooth bivariate function. The random effects are simulated as $\gamma_{ik}(t) = a_{i1}\phi_{1}(t) + a_{i2}\phi_{2}(t)$, $k=0,1$. We use the scaled orthonormal functions
\[
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
\]

![alt text](3Dbeta_trueIntercept.jpeg)




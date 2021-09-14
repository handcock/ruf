---
output:
  md_document:
    variant: markdown_github
---



# Software to Implement Resource Utilization Function Estimation



The *ruf* package implements the statistical methods described in the paper
*Relating Resources to a Probabilistic Measure of Space Use: Forest Fragments and Steller's Jays*
by John M. Marzluff, J. J. Millspaugh, P. Hurvitz, and Mark S. Handcock,
*Ecology*, 2004, 85:1411-1427 [https://www.jstor.org/stable/3450181](https://www.jstor.org/stable/3450181).

It is an `R` package to estimate resource utilization functions from point measurements of use and related functions. 
There are also functions to estimate
summary measures and their uncertainties.

We recommend to install the current R version (**R4.0**) and the **Rtools4.0**.

The package can  be installed from github via:


```r
if (!requireNamespace("remotes")) install.packages("remotes")

remotes::install_github("handcock/ruf")
```
 
Once installed, it can be used by typing the call


```r
library(ruf)
```

To see help on the functions type:


```r
help(ruf.fit)
```

### Case study: The resource utilizatino of Bird D412

As an example, below is an R program to reconstruct the analysis in the paper. 
It estimates the parameters of the model and related measures of
uncertainty. This requires the above `ruf` package.


```r
#
# load the "ruf" library
#
library(ruf)
#
# attach the small test data within the library
#
data(d412)
#
# Set initial estimates at the spatial range and smoothness
#
hval <- c(0.2, 1.5)
#
# Estimate the maximum likelihood values
# with unstandardized coefficients
#
d412.fit <- ruf.fit(ruf ~ CWED + IJI + NP + MSI,
         space= ~ x + y,
         data=d412, theta=hval,
         name="Bird 412",
         standardized=FALSE)
#> Fitting completed successfully!
```
Show the details of the results


```r
summary(d412.fit)
#> 
#> Unstandardized Coefficients for name: Bird 412 
#> 
#> Matern Log-Lik = -145.6376 LS Log-Lik = -151.4624 
#> 
#> Change in Log-Lik 5.82484 p-value = 0.0029533 
#> 
#>                   MLE     s.e. LS estimate  LS s.e.
#> range        0.217498 0.041404          NA       NA
#> smoothness   1.500000       NA          NA       NA
#> (Intercept)  0.708582 1.653752    1.370407 1.650380
#> CWED         0.004491 0.046393    0.003097 0.046551
#> IJI         -0.014597 0.013885   -0.011364 0.013664
#> NP           0.368734 0.214525    0.463118 0.214353
#> MSI          0.211891 1.343349   -0.479440 1.348894
```

Estimate the maximum likelihood values
with standardized coefficients


```r
d412.fit <- ruf.fit(ruf ~ CWED + IJI + NP + MSI,
         space= ~ x + y,
         data=d412, theta=hval,
         name="Bird 412 standardized",
         standardized=T)
#> Fitting completed successfully!
summary(d412.fit)
#> 
#> Standardized Coefficients for name: Bird 412 standardized 
#> 
#> Matern Log-Lik = -145.6376 LS Log-Lik = -151.4624 
#> 
#> Change in Log-Lik 5.82484 p-value = 0.0029533 
#> 
#>                   MLE     s.e. LS estimate  LS s.e.
#> range        0.217498 0.041404          NA       NA
#> smoothness   1.500000       NA          NA       NA
#> (Intercept)  1.663544 0.106557    1.370407 1.650380
#> CWED         0.043437 0.448748    0.003097 0.046551
#> IJI         -0.234598 0.223152   -0.011364 0.013664
#> NP           0.364972 0.212336    0.463118 0.214353
#> MSI          0.035309 0.223854   -0.479440 1.348894
```

### References 

For further readings and references please consult: 
  
  + John M. Marzluff, J. J. Millspaugh, P. Hurvitz, and Mark S. Handcock, *Relating Resources to a Probabilistic Measure of Space Use: Forest Fragments and Steller's Jays* (2004) *Ecology*, 2004, 85:1411-1427 [https://www.jstor.org/stable/3450181](https://www.jstor.org/stable/3450181).

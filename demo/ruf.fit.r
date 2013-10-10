#
# Attach the library
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
#
# Show the details of the results
#
summary(d412.fit)
#
# Estimate the maximum likelihood values
# with standardized coefficients
#
d412.fit <- ruf.fit(ruf ~ CWED + IJI + NP + MSI,
         space= ~ x + y,
         data=d412, theta=hval,
         name="Bird 412 standardized",
         standardized=TRUE)
summary(d412.fit)

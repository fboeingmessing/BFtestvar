# BFtestvar
Bayes factors for testing variances. This [Shiny](https://shiny.rstudio.com) application computes the adjusted fractional Bayes factor presented in Böing-Messing, F., van Assen, M. A. L. M., Hofman, A. D., Hoijtink, H., & Mulder, J. (2017). Bayesian evaluation of constrained hypotheses on variances of multiple independent groups. *Psychological Methods*, *22*(2), 262–287. https://doi.org/10.1037/met0000116.

The app can be launched using the `runGitHub` function from the `shiny` package in R:
```r
install.packages("shiny")
library(shiny)
runGitHub("BFtestvar", "fboeingmessing")
```

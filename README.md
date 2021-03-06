# R/`survtmlerct`

> Estimators with efficiency guarantees for the RMST and survival probability in randomized trials with right censoring 
__Author:__ [Iván Díaz](http://idiaz.xyz/)

---

## Description

`survtmlerct` is an R package that computes estimates of the restricted mean survival time and survival probability at a user-given time, using targeted minimum loss based estimation, and augmented inverse probability weighting. These parameters are 
useful for evaluating efficacy of treatments on survival outcomes. The estimators implemented have a theoretical guarantee to  improve asymptotic efficiency compared to their unadjusted counterparts in randomized trials, but can also be used in observational studies.

For completeness, the package also implements several unadjusted estimators based in inverse probability weighting and the Kaplan-Meier estimator. 

The estimators implemented in the package use preliminary estimats of: (i) the probability of treatment conditional on basline covariates, (ii) the hazard of censoring conditional on baseline covariates, and (iii) the hazard of the outcome conditional on covariates. These estimators may be based on model selection or data-adaptive regression such as machine learning. The estimators are doubly robust in that they remain consistent if either (i) and (ii) are consistently estimated, or (iii) is consistently estimated. 

See [Díaz et al. (2018)](https://link.springer.com/article/10.1007/s10985-018-9428-5) for more details. 

---

## Installation

A developmental release may be installed from GitHub via
[`devtools`](https://www.rstudio.com/products/rpackages/devtools/) with:

```{r gh-installation, eval = FALSE}
devtools::install_github("idiazst/survtmlerct")
```

---

## Help Files for main functions

Help files are available for the main functions in the package and can be viewed by typing
?survtmlerct::<function_name> for one of the following functions:
```{r}
tmle_prob 	Compute TMLE of survival probability with efficiency gains
tmle_rmst 	Compute TMLE of RMST with efficiency gains
aipw_prob   	Compute AIPW of survival probability
aipw_rmst   	Compute AIPW of RMST
ipw_rmst	Compute IPW of RMST
transformData 	Transform a survival dataset from short to long form
unadjusted_prob	Compute unadjusted Kaplan-Meier estimators of survival probability
unadjusted_rmst	Compute unadjusted IPW and Kaplan-Meier estimators of RMST
```

---

## Usage

Here we demonstrate calls to `survtmlerct` to compute treatment effects using the RMST and
the dataset colon from the `survival` library. 

```{r}
library(survtmlerct)
library(survival)
library(tidyverse)
library(dummies)
library(glm2)

# load data
data(colon)

# rename variables data
data <- colon %>% dummy.data.frame(c('differ', 'extent')) %>%
   filter(rx != 'Obs') %>%
   mutate(A = rx == 'Lev+5FU', id = as.numeric(as.factor(id)),
          nanodes = is.na(nodes), nodes = ifelse(is.na(nodes), 0, nodes)) %>%
   select(-rx) %>%  group_by(id) %>% summarise_all(funs(min)) %>% select(-study) %>%
   rename(T = time, D = status)
   
# transform data to long form
dlong <- transformData(data, 30)

# fit models for the outcome, censoring, and treatment 
fitL <- glm(Lm ~ A * (m + sex + age + obstruct + perfor + adhere + nodes + D +
                        differ1 + differ2 + differ3 + differNA + extent1 + extent2 +
                        extent3 + extent4 + surg + node4 + etype),
            data = dlong, subset = Im == 1, family = binomial())
            
# in order to obtain efficiency guarantees compared to unadjusted estimators, time must be a factor in the censoring model:
fitR <- glm(Rm ~ A * (as.factor(m) + sex + age + obstruct + perfor + adhere + nodes + D +
                        differ1 + differ2 + differ3 + differNA + extent1 + extent2 +
                        extent3 + extent4 + surg + node4 + etype),
            data = dlong, subset = Jm == 1, family = binomial())
fitA <- glm(A ~ sex + age + obstruct + perfor + adhere + nodes + D +
              differ1 + differ2 + differ3 + differNA + extent1 + extent2 +
              extent3 + extent4 + surg + node4 + etype,
            data = dlong, subset = m == 1, family = binomial())

# add preliminary estimates to the data

dlong <- mutate(dlong,
                gR1 = bound01(predict(fitR, newdata = mutate(dlong, A = 1), type = 'response')),
                gR0 = bound01(predict(fitR, newdata = mutate(dlong, A = 0), type = 'response')),
                h1  = bound01(predict(fitL, newdata = mutate(dlong, A = 1), type = 'response')),
                h0  = bound01(predict(fitL, newdata = mutate(dlong, A = 0), type = 'response')),
                gA1 = bound01(predict(fitA, newdata = mutate(dlong, A = 1), type = 'response')))
                
# estimate the RMST using tmle, aipw, ipw, and unadjusted estimators
tau <- max(dlong$m)
rmst_tmle <- tmle_rmst(dlong, tau)
rmst_aipw <- aipw_rmst(dlong, tau)
rmst_ipw  <- ipw_rmst(dlong, tau)
rmst_unad <- unadjusted_rmst(dlong, tau)

# look at estimates
rmst_tmle
rmst_aipw
rmst_ipw
rmst_unad
```


You can also compute estimates of the survival probability at a given time point:

```{r}
# estimate the survival probability using tmle, aipw, ipw, and unadjusted estimators
tau <- 20
prob_tmle <- tmle_prob(dlong, tau)
prob_aipw <- aipw_prob(dlong, tau)
prob_unad <- unadjusted_prob(dlong, tau)

# look at estimates
prob_tmle
prob_aipw
prob_unad
```



## Issues

If you encounter any bugs or have any specific feature requests, please [file an
issue](https://github.com/idiazst/survtmlerct/issues).

---

## Citation

After using the `survtmlerct` R package, please cite both of the following:

    @Manual{survtmlerct,
      title = {survtmlerct: Efficiency guarantees for covariate adjustment in RCTs with survival outcomes},
      author = {Iv\'an D\'iaz},
      note = {R package version 1.0.0},
      doi = {TBA}
    }

    @article{diaz2019improved,
      title={Improved precision in the analysis of randomized trials with survival outcomes, without assuming proportional hazards},
      author={D{\'\i}az, Iv{\'a}n and Colantuoni, Elizabeth and Hanley, Daniel F and Rosenblum, Michael},
      journal={Lifetime data analysis},
      volume={25},
      number={3},
      pages={439--468},
      year={2019},
      publisher={Springer}
    }

---

## License

&copy; 2020- [Iván Díaz](http://idiaz.xyz/)

The contents of this repository are distributed under the MIT license. See
below for details:
```
The MIT License (MIT)
Copyright (c) 2020- David Benkeser
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```
© 2020 GitHub, Inc.
Terms
Privacy
Security
Status
Help
Contact GitHub
Pricing
API
Training
Blog
About

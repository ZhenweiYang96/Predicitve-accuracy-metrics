# Title: Time-dependent Predictive Accuracy Metrics in the Context of Interval Censoring and Competing Risks

**Author**: Zhenwei Yang 

**Affiliation**: Department of Biostatistics, Erasmus Medical Center

**DOI**: [under review]

****

## Content

* Predicitve-accuracy-metrics.Rproj
* [Data](#Example_ICJM1)
* [Rscript](#Rscript)
* [Output](#Output)
* [img](#poster)

****

### Data

|File name| Description|
|:----------:|:--------------|
|ICJM1.RData | the ICJM1 results|
|CV | the folder containing model results of ICJM1 based on 5-fold cross-validation|

### Rscript

* Data analysis

|File name| Description|
|:----------:|:--------------|
|evaluation metric_PASS_CV.R| the script to evaluate the ICJM1 based on 5-fold cross-validation|

* Function

|File name| Description|
|:----------:|:--------------|
|dynpred_corspe.R | the function used to generate predictions from the correctly-specified ICJM1|
|dynpred_linear.R | the function used to generate predictions from the ICJM1 with linear PSA|
|dynpred_nocovar.R | the function used to generate predictions from the ICJM1 ignoring baseline covariate PSA density|
|function.R | help and auxiliary functions|
|ICJM_correctlyspecify.R | the function used to fit the correctly-specified ICJM1|
|ICJM_linear.R | the function used to fit the ICJM1 with linear PSA|
|ICJM_nocovar.R | the function used to fit the ICJM1 ignoring baseline covariate PSA density|
|time dependent evaluation metrics.R | the function used for model evaluation based on AUC, Brier score and EPCE|
|time dependent true EPCE | the function used to calculate the EPCE in absence of censoring|

* Simulation

|File name| Description|
|:----------:|:--------------|
|Seed | including the seeds to generate the data & fit the model, the seeds for individual's personalized scheduling, and model test|
|Dataset generation.R | training and test sets generation based on three models and four scheduling scenarios|
|MM_fit.R | fit the longitudinal submodel used as an initial values in the joint model estimation|
|Model Estimation.R | model fit based on 200 simulated training sets according to the model specification and scheduling scenarios|
|Evaluation.R | evaluate the simulated models based on AUC, Brier score and EPCE|
|Evaluation_true.R | evaluate the simulated models based on reference AUC, Brier score and EPCE in absence of censoring|


* Tables&figures

|File name| Description|
|:----------:|:--------------|
|Paper plots.R | draw the plots|


### Output

* Data analysis

|File name| Description|
|:----------:|:--------------|
|AccMet_CV_1.RData | the model evaluation results (AUC, Brier score and EPCE) of the first fold version of the ICJM1 |
|AccMet_CV_2.RData | the model evaluation results (AUC, Brier score and EPCE) of the second fold version of the ICJM1 |
|AccMet_CV_3.RData | the model evaluation results (AUC, Brier score and EPCE) of the third fold version of the ICJM1 |
|AccMet_CV_4.RData | the model evaluation results (AUC, Brier score and EPCE) of the fourth fold version of the ICJM1 |
|AccMet_CV_5.RData | the model evaluation results (AUC, Brier score and EPCE) of the fifth fold version of the ICJM1 |

* Tables&figures


|File name| Description|
|:----------:|:--------------|
|Main text | all figures used in the manuscripts|
|Supplementary| all figures used in the supplementary materials|

### Package dependencies

- rajgs
- mcmcplots
- GLMMadaptive
- ggplot2
- tidyverse
- splines
- future
- mcmcse
- mvtnorm
- JMbayes2
- JMbayes
- MASS
- doParallel
- truncnorm
- Matrix
- latex2exp
- cowplot

> [!Note]
> - Please make sure to install above-mentioned packages by `install.packages()` before running the R code

### Poster

- The corresponding poster presented in ISCB 2024:

![](https://github.com/ZhenweiYang96/Predicitve-accuracy-metrics/blob/main/img/2024_ISCB_Poster.png)

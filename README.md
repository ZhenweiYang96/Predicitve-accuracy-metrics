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
|Simulation Datasets | Simulated datasets used in simulation including PASS and random schdules|

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
|time dependent true EPCE.R | the function used to calculate the EPCE in absence of censoring|

* Simulation

|File name| Description|
|:----------:|:--------------|
|Seed | including the seeds to generate the data & fit the model, the seeds for individual's personalized scheduling, and model test|
|Dataset generation.R | training and test sets generation based on three models and four scheduling scenarios|
|MM_fit.R | fit the longitudinal submodel used as an initial values in the joint model estimation|
|Model Estimation.R | model fit based on 200 simulated training sets according to the model specification and scheduling scenarios|
|Evaluation.R | evaluate the simulated models based on AUC, Brier score and EPCE|
|Evaluation_naive.R | evaluate the simulated models based on AUC, Brier score and EPCE using the naive approach (ignoring censoring)|
|Evaluation_true.R | evaluate the simulated models based on reference AUC, Brier score and EPCE in absence of censoring|


* Tables&figures

|File name| Description|
|:----------:|:--------------|
|Figure tables_main.R | draw the plots and tables in the manuscript|
|Figure tables_supplement.R | draw the plots and tables in the supplement|


### Output

* Data analysis

|File name| Description|
|:----------:|:--------------|
|AccMet_CV_1.RData | the model evaluation results (AUC, Brier score and EPCE) of the first fold version of the ICJM1 |
|AccMet_CV_2.RData | the model evaluation results (AUC, Brier score and EPCE) of the second fold version of the ICJM1 |
|AccMet_CV_3.RData | the model evaluation results (AUC, Brier score and EPCE) of the third fold version of the ICJM1 |
|AccMet_CV_4.RData | the model evaluation results (AUC, Brier score and EPCE) of the fourth fold version of the ICJM1 |
|AccMet_CV_5.RData | the model evaluation results (AUC, Brier score and EPCE) of the fifth fold version of the ICJM1 |

* Simulation Evaluation

* Simulation Models

|File name| Description|
|:----------:|:--------------|
|Correctly specified models | correctly-specified joint models |
|Fitted mixed model| mixed model part for the initials for fitting the joint models |
|Linear models | joint models assuming linear PSA|
|Linear models | joint models omiting PSA density|

* Simulation Storage

|File name| Description|
|:----------:|:--------------|
|corspe | detailed predictions for correctly-specified joint model to visualized PSA trajectories|
|linear | detailed predictions for joint model assuming linear PSA to visualized PSA trajectories|

* Tables&figures

|File name| Description|
|:----------:|:--------------|
|Main text | all figures used in the manuscripts|
|Supplementary| all figures used in the supplementary materials|

****

## Reproducibility 

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
- PCaASSim (this can be downloaded via: `devtools::install_github("https://github.com/ZhenweiYang96/PCaASSim")`)

> [!Note]
> - Please make sure to install above-mentioned packages by `install.packages()` before running the R code

### Instructions for reproducing the results

#### Section 3 - Application

- The evaluation metrics for the Canary PASS data (with cross validation) can be derived from this path: */R script/Data analysis/evaluation metric_PASS_CV.R* (In the process, there will be some saved results, e.g., the model results for each fold stored under */Data/CV/ICJM_CV_1.RData* to */Data/CV/ICJM_CV_5.RData*).

> [!Note]
> - The original Canary PASS data (*/Data/pass_id.RData* and */Data/pass.RData*) is not publically sharable and thus does not present in this repository.

#### Section 4 - Simulation

1. Data generation: run */R script/Simulation/Data Generation.R* to generate $200 \times 4$ datasets. There are four types of scnearios: the standard PASS schedule (stored under */Data/Simulation Datasets/*), three random schedules (stored under */Data/Simulation Datasets/random_3_10/*, */Data/Simulation Datasets/random_3_40/* and */Data/Simulation Datasets/random_10_20/*)

2. Preparation for model fitting: run */R script/Simulation/MM_fit.R*. 

3. Model fitting for each simulated dataset: run */R script/Simulation/Model Estimation.R*. For each dataset, three models were fitted: correctly-specified model stored under */Output/Simulation Models/Correctly specified models/*, model with linear PSA stored under */Output/Simulation Models/Linear models/*, and the model ignoring the baseline covariate PSA density stored under */Output/Simulation Models/No covariate models/*.

4. Model evaluating (two proposed approaches in the manuscript): run */R script/Simulation/Evaluation.R*. Results should be stored under */Output/Simulation Evaluation/*, and are needed for the plots.

5. Model evaluating (the naive approach in the manuscript): run */R script/Simulation/Evaluation_naive.R*. Results should be stored under */Output/Simulation Evaluation/Observed/*, and are needed for the plots.

6. Model evaluating (the reference in the manuscript): run */R script/Simulation/Evaluation_true.R*. Results should be stored under */Output/Simulation Evaluation/True linear/*, */Output/Simulation Evaluation/True models/* and */Output/Simulation Evaluation/True no cov/*, and are needed for the plots.

#### Figures and Tables in the manuscript and supplements

- run */R script/Tables & figures/Figure tables_main.R* for all the tables and figures in manuscript. This file uses the all the results from the previous steps which should be stored under */Output/*. The figures are stored under */Output/Tables & figures/Main text/*.

- run */R script/Tables & figures/Figure tables_supplement.R* for all the tables and figures in supplements. This file uses the all the results from the previous steps which should be stored under */Output/*. The figures are stored under */Output/Tables & figures/Supplementary/*.


****

## Poster

- The corresponding poster presented in ISCB 2024:

![](https://github.com/ZhenweiYang96/Predicitve-accuracy-metrics/blob/main/img/2024_ISCB_Poster.png)

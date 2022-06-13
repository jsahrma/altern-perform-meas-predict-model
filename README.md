# Replication Results for Wu & Lee (2014)

In ["Alternative Performance Measures for Prediction Models"](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0091249), Wu & Lee report results from a simulation study comparing changes in measures of model predictive performance when variables with strong explanatory power are added to a model. Their motivation is the claim that one of the most popular metrics for evaluating predictive models---the area under the receiver operating characteristics curve (AUC)---is relatively insensitive to the addition of strongly predictive variables to an existing model. They use the Gini coefficient, the Pietra index, and the scaled Brier score as comparators. Equations for calculating each of the metrics using the prevalence of the disease (or whatever is being predicted by the model), the predicted probabilities, and the true disease state of each subject are given in the Methods section.

The code in this repository is meant to replicate the results of each of Wu & Lee (2014)'s three simulation schemes---one of which is presented in the manuscript and two of which are provided in supplementary material. Each simulation hypothesizes an underlying disease mechanism, with each patient's status determined by a baseline score (confusingly referred to initially as $S$ but thereafter as $B$) and two or more markers ($M_1, M_2, \dotsc, M_{10}$), which are binary in the first and third simulation schemes and continuous in the second. Baseline scores are drawn from a standard normal distribution. Importantly, one marker is related to disease state independent of the baseline score, while the second is related to disease state in such a way that it is most predictive of the outcome when the baseline score is least predictive, i.e., when $B \approx 0$, which indicates that a patient is neither very likely or unlikely to have the disease. Wu & Lee describe these markers as having discrimination ability in the "gray zone" of a model based purely on the baseline score.

Given a simulated patient sample, Wu & Lee (2014) attempt to predict disease state using three models---one containing the baseline score only, one containing the baseline score and the independently predictive marker, and one containing the baseline score and the gray zone-discriminating marker. Wu & Lee (2014) don't explicitly state their analytic approach, but the formulation of the disease model and its popularity suggests that logistic regression was used.


## Organization of the Repository

The single R script can be found as `sim-wu-&-lee-2014.R` in the *code* subdirectory. Output consisting of CSV files with results resembling Table 1 of Wu & Lee (2014) and PNGs of histograms designed to (mostly) match Figure 2 are located in the *output* subdirectory.


## Results

| | Performance Measure | | | |
| :--- | ---: | ---: | ---: | ---: |
| | AUC | Gini | Pietra | Scaled Brier |
| Model | | | | |
| B | 0.824 | 0.644 | 0.488 | 0.309 |
| B + M1 | 0.841 | 0.684 | 0.523 | 0.347 |
| B + M2 | 0.836 | 0.671 | 0.512 | 0.335 |
| Absolute (Relative) Improvement | | | | |
| from B to B + M1 | 0.016 (2.0%) | 0.039 (6.1%) | 0.035 (7.3%) | 0.038 (12.2%) |
| from B to B + M2 | 0.012 (1.5%) | 0.027 (4.2%) | 0.025 (5.0%) | 0.026 (8.3%) |

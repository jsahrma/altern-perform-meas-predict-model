# Header -------------------------------------------------------------
#
# Replicate the simulation results from Wu & Lee (2014):
# <https://doi.org/10.1371/journal.pone.0091249>.
#
# Reference:
#
# Wu Y-C, Lee W-C (2014) Alternative Performance Measures for
# Prediction Models. PLoS ONE 9(3): e91249.
#
# John Sahrmann
# 20220512


# Setup --------------------------------------------------------------

library(arm)
library(dplyr)
library(purrr)


# Constant definitions -----------------------------------------------

nsimul <- 10000
npt <- 1000


# Function definitions -----------------------------------------------

#' Compute the area under the ROC curve.
#'
#' @param predict Vector of predicted probabilities
#' @param outcome Vector of actual outcomes
#' @return The area under the ROC curve
#' @examples
#' auc(c(.25, .3, .6, .89), c(0, 1, 1, 1))
auc <- function(predict, outcome) {
  # Reorder `predict` and `outcome` such that cases with `outcome ==
  # 1` appear first.
  temp_predict <- predict
  temp_outcome <- outcome
  predict <- predict[
    order(temp_outcome, temp_predict, decreasing = TRUE)]
  outcome <- outcome[
    order(temp_outcome, temp_predict, decreasing = TRUE)]

  n <- length(outcome)
  n1 <- sum(outcome)
  sum_scores <- 0
  for (i in 1:n1) {
    for (j in (n1 + 1):n) {
      if (predict[[i]] > predict[[j]]) {
        sum_scores <- sum_scores + 1
      } else if (predict[[i]] == predict[[j]]) {
        sum_scores <- sum_scores + 0.5
      }
    }
  }
  sum_scores / (n1 * (n - n1))
}

#' Compute the Gini index.
#'
#' @param predict Vector of predicted probabilities
#' @param outcome Vector of actual outcomes
#' @return The Gini index
#' @examples
#' gini(c(.25, .3, .6, .89), c(0, 1, 1, 1))
gini <- function(predict, outcome) {
  n <- length(outcome)
  p_bar <- mean(outcome)
  
  numer <- 0
  for (i in 1:n) {
    for (j in 1:n) {
      numer <- numer + abs(predict[[i]] - predict[[j]])
    }
  }
  denom <- 2 * n^2 * p_bar * (1 - p_bar)
  numer / denom
}

#' Compute the Pietra index.
#'
#' @param predict Vector of predicted probabilities
#' @param outcome Vector of actual outcomes
#' @return The Pietra index
#' @examples
#' pietra(c(.25, .3, .6, .89), c(0, 1, 1, 1))
pietra <- function(predict, outcome) {
  n <- length(outcome)
  p_bar <- mean(outcome)
  
  denom <- 2 * n * p_bar * (1 - p_bar)
  numer / denom
}

#' Compute the scaled Brier score.
#'
#' @param predict Vector of predicted probabilities
#' @param outcome Vector of actual outcomes
#' @return The scaled Brier score
#' @examples
#' sbrier(c(.25, .3, .6, .89), c(0, 1, 1, 1))
sbrier <- function(predict, outcome) {
  n <- length(outcome)
  p_bar <- mean(outcome)
  
  numer <- sum((predict - p_bar)^2)
  denom <- n * p_bar * (1 - p_bar)
  numer / denom
}


k1 <- function(x) exp(-x^2 / 0.5)

probab_dis_simul1 <- function(s, m1, m2) {
  arm::invlogit(-3 + (2 * s) + (1.5 * m1) + (2.2 * k1(s) * m2))
}


# Simulations --------------------------------------------------------


# Scheme 1 ------------------------

set.seed(725533)

dat <- dplyr::tibble(base_score = rnorm(npt)) %>%
  dplyr::mutate(
    marker1_85 = sample(0:1, npt, replace = TRUE, prob = c(.15, .85)),
    marker1_75 = sample(0:1, npt, replace = TRUE, prob = c(.25, .75)),
    marker1 = ifelse(base_score > 0, marker1_85, marker1_75),
    marker2_85 = sample(0:1, npt, replace = TRUE, prob = c(.15, .85)),
    marker2_75 = sample(0:1, npt, replace = TRUE, prob = c(.25, .75)),
    marker2 = ifelse(base_score > 0, marker2_85, marker2_75),
    dis = ifelse(
      runif(npt) < probab_dis_simul1(base_score, marker1, marker2),
      1L, 0L)
  ) %>%
  dplyr::select(-dplyr::matches("(_75|_85)$"))

train <- dat[1:(npt / 2), ]
valid <- dat[((npt / 2) + 1):npt, ]

model1_simul1 <- glm(
  dis ~ base_score, family = binomial, data = train)
predict1_simul1 <- predict(
  model1_simul1, newdata = valid, type = "response")

auc1_simul1 <- auc(predict1_simul1, valid$dis)
gini1_simul1 <- gini(predict1_simul1, valid$dis)
pietra1_simul1 <- pietra(predict1_simul1, valid$dis)
sbrier1_simul1 <- sbrier(predict1_simul1, valid$dis)

model2_simul1 <- glm(
  dis ~ base_score + marker1, family = binomial, data = train)
predict2_simul1 <- predict(
  model2_simul1, newdata = valid, type = "response")

auc2_simul1 <- auc(predict2_simul1, valid$dis)
gini2_simul1 <- gini(predict2_simul1, valid$dis)
pietra2_simul1 <- pietra(predict2_simul1, valid$dis)
sbrier2_simul1 <- sbrier(predict2_simul1, valid$dis)

auc2_simul1 - auc1_simul1
gini2_simul1 - gini1_simul1
pietra2_simul1 - pietra1_simul1
sbrier2_simul1 - sbrier1_simul1


run_scheme1 <- function() {
  # Simulate a data set containing the covariates and disease state.
  dat <- dplyr::tibble(base_score = rnorm(npt)) %>%
    dplyr::mutate(
      marker1_85 = sample(
        0:1, npt, replace = TRUE, prob = c(.15, .85)),
      marker1_75 = sample(
        0:1, npt, replace = TRUE, prob = c(.25, .75)),
      marker1 = ifelse(base_score > 0, marker1_85, marker1_75),
      marker2_85 = sample(
        0:1, npt, replace = TRUE, prob = c(.15, .85)),
      marker2_75 = sample(
        0:1, npt, replace = TRUE, prob = c(.25, .75)),
      marker2 = ifelse(base_score > 0, marker2_85, marker2_75),
      dis = ifelse(
        runif(npt) < probab_dis_simul1(base_score, marker1, marker2),
        1L, 0L)
    ) %>%
    dplyr::select(-dplyr::matches("(_75|_85)$"))

  # Split into training and validation sets.
  train <- dat[1:(npt / 2), ]
  valid <- dat[((npt / 2) + 1):npt, ]

  # Fit each model.
  model1 <- glm(dis ~ base_score, family = binomial, data = train)
  model2 <- glm(
    dis ~ base_score + marker1, family = binomial, data = train)
  model3 <- glm(
    dis ~ base_score + marker2, family = binomial, data = train)

  # Generate predicted probabilities of disease for the validation set
  # using each model.
  predict1 <- predict(model1, newdata = valid, type = "response")
  predict2 <- predict(model2, newdata = valid, type = "response")
  predict3 <- predict(model3, newdata = valid, type = "response")

  # Calculate the predictive performance measures under each model.
  auc1 <- auc(predict1, valid$dis)
  auc2 <- auc(predict2, valid$dis)
  auc3 <- auc(predict3, valid$dis)
  gini1 <- gini(predict1, valid$dis)
  gini2 <- gini(predict2, valid$dis)
  gini3 <- gini(predict3, valid$dis)
  pietra1 <- pietra(predict1, valid$dis)
  pietra2 <- pietra(predict2, valid$dis)
  pietra3 <- pietra(predict3, valid$dis)
  sbrier1 <- sbrier(predict1, valid$dis)
  sbrier2 <- sbrier(predict2, valid$dis)
  sbrier3 <- sbrier(predict3, valid$dis)

  list(
    dis = valid$dis,
    predict1 = predict1, predict2 = predict2, predict3 = predict3,
    auc1 = auc1, auc2 = auc2, auc3 = auc3,
    gini1 = gini1, gini2 = gini2, gini3 = gini3,
    pietra1 = pietra1, pietra2 = pietra2, pietra3 = pietra3,
    sbrier1 = sbrier1, sbrier2 = sbrier2, sbrier3 = sbrier3
  )
}


results <- purrr::rerun(10, run_scheme1())

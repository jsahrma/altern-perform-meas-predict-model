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
# 20220521


# Setup --------------------------------------------------------------

library(arm)
library(dplyr)
library(purrr)


# Constant definitions -----------------------------------------------

nsimul <- 10000                    # number of simulations
npt <- 1000                        # number of subjects per simulation


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

  numer <- sum(abs(predict - p_bar))
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


#' Reproduce Figure 3 in Wu & Lee (2014).
#'
#' Note that the histograms from the two distributions appear to be
#' stacked in the text, which doesn't seem ideal. That will not be
#' reproduced here.
#'
#' @param simul_results A list of simulation results produced by any
#'   of the `run_scheme*` functions.
#' @return The simulation results given as input (invisibly)
#' @examples
#' results <- purrr::rerun(100, run_scheme1())
#' plot_simul_results(results)
plot_simul_results <- function(simul_results) {
  # Define a helper function to produce each plot.
  plot_predict_probab_distrib <- function(predictions, outcomes) {
    hist(
      predictions[outcomes == 0], breaks = 20, freq = FALSE,
      col = "gray80", main = "", xlab = "Predicted Probability",
      ylab = "Density"
    )
    hist(
      predictions[outcomes == 1], breaks = 20, freq = FALSE,
      col = "gray30", add = TRUE
    )
    abline(v = mean(outcomes))
    abline(v = tapply(predictions, outcomes, mean), lty = "dashed")
    legend(
      "topright", legend = c("Diseased", "Non-diseased"),
      pch = 15, col = c("gray40", "gray80")
    )
    invisible(simul_results)
  }

  # Extract the outcomes and predictions from across all simulations.
  x <- unlist(purrr::map(simul_results, ~ purrr::pluck(.x, "dis")))
  p1 <- unlist(purrr::map(results, ~ purrr::pluck(.x, "predict1")))
  p2 <- unlist(purrr::map(results, ~ purrr::pluck(.x, "predict2")))
  p3 <- unlist(purrr::map(results, ~ purrr::pluck(.x, "predict3")))
  
  par(mfrow = c(3, 1))
  plot_predict_probab_distrib(p1, x)
  plot_predict_probab_distrib(p2, x)
  plot_predict_probab_distrib(p3, x)
}


# Simulations --------------------------------------------------------


# Scheme 1 ------------------------

run_scheme1 <- function() {
  # Define the (true) disease risk model for this simulation scheme.
  probab_dis <- function(s, m1, m2) {
    # Define the Gaussian kernel function used for this disease model.
    kernel <- function(x) exp(-x^2 / 0.5)
    # Return the probability of disease.
    arm::invlogit(
      -3 + (2 * s) + (1.5 * m1) + (2.2 * kernel(s) * m2))
  }

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
        runif(npt) < probab_dis(base_score, marker1, marker2), 1L, 0L)
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
    dis ~ base_score + marker1 + marker2, family = binomial, data = train)

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

## system.time({
## set.seed(781649)
## results <- purrr::rerun(nsimul, run_scheme1())
## })

set.seed(781649)
results <- purrr::rerun(100, run_scheme1())


# plot_predict_probab

pp <- unlist(purrr::map(results, ~ purrr::pluck(.x, "predict3")))
oo <- unlist(purrr::map(results, ~ purrr::pluck(.x, "dis")))

hist(pp[oo == 0], breaks = 20, freq = FALSE, col = "gray80")
hist(pp[oo == 1], breaks = 20, freq = FALSE, col = "gray30", add = TRUE)
abline(v = mean(oo))
abline(v = tapply(pp, oo, mean), lty = "dashed")
legend(
  "topright", legend = c("Diseased", "Non-diseased"), pch = 15, col = c("gray40", "gray80")
)



plot_simul_results(results)


# analyze_simul

results1 <- results

auc1 <- mean(purrr::map_dbl(results, ~ purrr::pluck(.x, "auc1")))
auc2 <- mean(purrr::map_dbl(results, ~ purrr::pluck(.x, "auc2")))
auc3 <- mean(purrr::map_dbl(results, ~ purrr::pluck(.x, "auc3")))
auc1
auc2
auc3

gini1 <- mean(purrr::map_dbl(results, ~ purrr::pluck(.x, "gini1")))
gini2 <- mean(purrr::map_dbl(results, ~ purrr::pluck(.x, "gini2")))
gini3 <- mean(purrr::map_dbl(results, ~ purrr::pluck(.x, "gini3")))
gini1
gini2
gini3

pietra1 <- mean(purrr::map_dbl(results, ~ purrr::pluck(.x, "pietra1")))
pietra2 <- mean(purrr::map_dbl(results, ~ purrr::pluck(.x, "pietra2")))
pietra3 <- mean(purrr::map_dbl(results, ~ purrr::pluck(.x, "pietra3")))
pietra1
pietra2
pietra3

sbrier1 <- mean(purrr::map_dbl(results, ~ purrr::pluck(.x, "sbrier1")))
sbrier2 <- mean(purrr::map_dbl(results, ~ purrr::pluck(.x, "sbrier2")))
sbrier3 <- mean(purrr::map_dbl(results, ~ purrr::pluck(.x, "sbrier3")))
sbrier1
sbrier2
sbrier3

d1 <- results[[1]]
mean(d1$dis); mean(d1$predict1); mean(d1$predict2); mean(d1$predict3)



x <- seq(-3, 3, length.out = 1000)
y1 <- exp(2.2 * exp(-x^2 / 0.5))
y2 <- exp(2.2 * exp(-x^2 / 2))

plot(y1 ~ x, type = "l")
lines(y2 ~ x, col = "blue")




# Define the (true) disease risk model for this simulation scheme.
probab_dis <- function(s, m1, m2) {
  # Define the Gaussian kernel function used for this disease model.
  kernel <- function(x) exp(-x^2 / 0.5)
  # Return the probability of disease.
  ## arm::invlogit(
  plogis(
    -3 + (2 * s) + (1.5 * m1) + (2.2 * kernel(s) * m2))
}

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
      runif(npt) < probab_dis(base_score, marker1, marker2), 1L, 0L)
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

auc1; auc2; auc3

list(
  dis = valid$dis,
  predict1 = predict1, predict2 = predict2, predict3 = predict3,
  auc1 = auc1, auc2 = auc2, auc3 = auc3,
  gini1 = gini1, gini2 = gini2, gini3 = gini3,
  pietra1 = pietra1, pietra2 = pietra2, pietra3 = pietra3,
  sbrier1 = sbrier1, sbrier2 = sbrier2, sbrier3 = sbrier3
)


train2 <- dat[1:500,]
valid2 <- dat[501:1000,]
model4 <- glm(dis ~ base_score + term3, data = train2, family = binomial)
pred4 <- predict(model4, newdata = valid2, type = "response")

auc(pred4, valid2$dis)
tapply(pred4, valid2$dis, mean)


# Scheme 2 ------------------------

run_scheme2 <- function() {
  # Define the (true) disease risk model for this simulation scheme.
  probab_dis <- function(s, m3, m4) {
    # Define the Gaussian kernel function used for this disease model.
    kernel <- function(x) exp(-x^2 / 0.0833)
    # Return the probability of disease.
    arm::invlogit(
      -3 + (2 * s) + (0.8 * m3) + (2.2 * kernel(s) * m4))
  }

  # Simulate a data set containing the covariates and disease state.
  dat <- dplyr::tibble(base_score = rnorm(npt)) %>%
    dplyr::mutate(
      marker3_s_gt0 = rnorm(npt, 3.65, 1),
      marker3_s_le0 = rnorm(npt, 3.55, 1),
      marker3 = ifelse(base_score > 0, marker3_s_gt0, marker3_s_le0),
      marker4_s_gt0 = rnorm(npt, 0.05, 1),
      marker4_s_le0 = rnorm(npt, -0.05, 1),
      marker4 = ifelse(base_score > 0, marker4_s_gt0, marker4_s_le0),
      dis = ifelse(
        runif(npt) < probab_dis(base_score, marker3, marker4), 1L, 0L)
    ) %>%
    dplyr::select(-dplyr::matches("(gt0|le0)$"))

  # Split into training and validation sets.
  train <- dat[1:(npt / 2), ]
  valid <- dat[((npt / 2) + 1):npt, ]

  # Fit each model.
  model1 <- glm(dis ~ base_score, family = binomial, data = train)
  model2 <- glm(
    dis ~ base_score + marker3, family = binomial, data = train)
  model3 <- glm(
    dis ~ base_score + marker4, family = binomial, data = train)

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

## system.time({
## set.seed(639320)
## results <- purrr::rerun(nsimul, run_scheme2())
## })

system.time({
set.seed(639320)
results <- purrr::rerun(1000, run_scheme2())
})


# plot_predict_probab

pp <- unlist(purrr::map(results, ~ purrr::pluck(.x, "predict1")))
oo <- unlist(purrr::map(results, ~ purrr::pluck(.x, "dis")))

hist(pp[oo == 0], breaks = 20, freq = FALSE, col = "gray80")
hist(pp[oo == 1], breaks = 20, freq = FALSE, col = "gray40", add = TRUE)
abline(v = mean(oo))
abline(v = tapply(pp, oo, mean), lty = "dashed")
legend(
  "topright", legend = c("Diseased", "Non-diseased"), pch = 15, col = c("gray40", "gray80")
)


# analyze_simul

## results1 <- results

auc1 <- mean(purrr::map_dbl(results, ~ purrr::pluck(.x, "auc1")))
auc2 <- mean(purrr::map_dbl(results, ~ purrr::pluck(.x, "auc2")))
auc3 <- mean(purrr::map_dbl(results, ~ purrr::pluck(.x, "auc3")))
auc1
auc2
auc3

gini1 <- mean(purrr::map_dbl(results, ~ purrr::pluck(.x, "gini1")))
gini2 <- mean(purrr::map_dbl(results, ~ purrr::pluck(.x, "gini2")))
gini3 <- mean(purrr::map_dbl(results, ~ purrr::pluck(.x, "gini3")))
gini1
gini2
gini3

pietra1 <- mean(purrr::map_dbl(results, ~ purrr::pluck(.x, "pietra1")))
pietra2 <- mean(purrr::map_dbl(results, ~ purrr::pluck(.x, "pietra2")))
pietra3 <- mean(purrr::map_dbl(results, ~ purrr::pluck(.x, "pietra3")))
pietra1
pietra2
pietra3

sbrier1 <- mean(purrr::map_dbl(results, ~ purrr::pluck(.x, "sbrier1")))
sbrier2 <- mean(purrr::map_dbl(results, ~ purrr::pluck(.x, "sbrier2")))
sbrier3 <- mean(purrr::map_dbl(results, ~ purrr::pluck(.x, "sbrier3")))
sbrier1
sbrier2
sbrier3

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
# 20220610


# Setup --------------------------------------------------------------

library(arm)
library(dplyr)
library(purrr)


# Constant definitions -----------------------------------------------

n_simul <- 10000                    # number of simulations
n_pt <- 1000                        # number of subjects per simulation


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


#' Reproduce Figure 3 in Wu & Lee (2014, p 3).
#'
#' Note that the histograms from the two distributions appear to be
#' stacked in the text, which doesn't seem ideal. That will not be
#' reproduced here.
#'
#' @param simul_results A list of simulation results produced by any
#'   of the `run_scheme*` functions
#' @return The simulation results given as input (invisibly)
#' @examples
#' results <- purrr::rerun(100, run_scheme1())
#' plot_simul(results)
plot_simul <- function(simul_results) {
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
  }

  # Extract the outcomes and predictions from across all simulations.
  x <- unlist(
    purrr::map(simul_results, ~ purrr::pluck(.x, "dis")))
  p1 <- unlist(
    purrr::map(simul_results, ~ purrr::pluck(.x, "predict1")))
  p2 <- unlist(
    purrr::map(simul_results, ~ purrr::pluck(.x, "predict2")))
  p3 <- unlist(
    purrr::map(simul_results, ~ purrr::pluck(.x, "predict3")))
  
  par(mfrow = c(3, 1))
  plot_predict_probab_distrib(p1, x)
  plot_predict_probab_distrib(p2, x)
  plot_predict_probab_distrib(p3, x)
  invisible(simul_results)
}


#' Reproduce Table 1 in Wu & Lee (2014, p 5).
#'
#' @param simul_results A list of simulation results produced by any
#'   of the `run_scheme*` functions
#' @return Summary of simulation results as a `data.frame`
#' @examples
#' results <- purrr::rerun(100, run_scheme1())
#' x <- summ_simul(results)
summ_simul <- function(simul_results) {
  # Define helper functions to print performance metrics and the
  # absolute and relative improvements in performance metrics.
  print_perform_meas <- function(x) {
    format(round(x, 3), nsmall = 3)
  }
  print_perform_improv <- function(new, old) {
    paste0(
      format(round(new - old, 3), nsmall = 3), " (",
      format(round((new - old) / old * 100, 1), nsmall = 1), "%)"
    )
  }

  auc1 <- mean(
    purrr::map_dbl(simul_results, ~ purrr::pluck(.x, "auc1")))
  auc2 <- mean(
    purrr::map_dbl(simul_results, ~ purrr::pluck(.x, "auc2")))
  auc3 <- mean(
    purrr::map_dbl(simul_results, ~ purrr::pluck(.x, "auc3")))
  gini1 <- mean(
    purrr::map_dbl(simul_results, ~ purrr::pluck(.x, "gini1")))
  gini2 <- mean(
    purrr::map_dbl(simul_results, ~ purrr::pluck(.x, "gini2")))
  gini3 <- mean(
    purrr::map_dbl(simul_results, ~ purrr::pluck(.x, "gini3")))
  pietra1 <- mean(
    purrr::map_dbl(simul_results, ~ purrr::pluck(.x, "pietra1")))
  pietra2 <- mean(
    purrr::map_dbl(simul_results, ~ purrr::pluck(.x, "pietra2")))
  pietra3 <- mean(
    purrr::map_dbl(simul_results, ~ purrr::pluck(.x, "pietra3")))
  sbrier1 <- mean(
    purrr::map_dbl(simul_results, ~ purrr::pluck(.x, "sbrier1")))
  sbrier2 <- mean(
    purrr::map_dbl(simul_results, ~ purrr::pluck(.x, "sbrier2")))
  sbrier3 <- mean(
    purrr::map_dbl(simul_results, ~ purrr::pluck(.x, "sbrier3")))

  col_text <- vector(mode = "list", length = 5)
  col_text[[1]] <- c(
    "", "", "Model", "B", "B + M1", "B + M2",
    "Absolute (Relative) Improvement", "from B to B + M1",
    "from B to B + M2"
  )
  col_text[[2]] <- c(
    "Performance Measure", "AUC", "",
    purrr::map_chr(
      list(auc1, auc2, auc3), print_perform_meas),
    "",
    print_perform_improv(auc2, auc1),
    print_perform_improv(auc3, auc1)
  )
  col_text[[3]] <- c(
    "", "Gini", "",
    purrr::map_chr(
      list(gini1, gini2, gini3), print_perform_meas),
    "",
    print_perform_improv(gini2, gini1),
    print_perform_improv(gini3, gini1)
  )
  col_text[[4]] <- c(
    "", "Pietra", "",
    purrr::map_chr(
      list(pietra1, pietra2, pietra3), print_perform_meas),
    "",
    print_perform_improv(pietra2, pietra1),
    print_perform_improv(pietra3, pietra1)
  )
  col_text[[5]] <- c(
    "", "Scaled Brier", "",
    purrr::map_chr(
      list(sbrier1, sbrier2, sbrier3), print_perform_meas),
    "",
    print_perform_improv(sbrier2, sbrier1),
    print_perform_improv(sbrier3, sbrier1)
  )
  table1 <- do.call(cbind, col_text)
  as.data.frame(table1)
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
  dat <- dplyr::tibble(base_score = rnorm(n_pt)) %>%
    dplyr::mutate(
      marker1_85 = sample(
        0:1, n_pt, replace = TRUE, prob = c(.15, .85)),
      marker1_75 = sample(
        0:1, n_pt, replace = TRUE, prob = c(.25, .75)),
      marker1 = ifelse(base_score > 0, marker1_85, marker1_75),
      marker2_85 = sample(
        0:1, n_pt, replace = TRUE, prob = c(.15, .85)),
      marker2_75 = sample(
        0:1, n_pt, replace = TRUE, prob = c(.25, .75)),
      marker2 = ifelse(base_score > 0, marker2_85, marker2_75),
      dis = ifelse(
        runif(n_pt) < probab_dis(base_score, marker1, marker2), 1L, 0L)
    ) %>%
    dplyr::select(-dplyr::matches("(_75|_85)$"))

  # Split into training and validation sets.
  train <- dat[1:(n_pt / 2), ]
  valid <- dat[((n_pt / 2) + 1):n_pt, ]

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

system.time({
set.seed(781649)
results_scheme1 <- purrr::rerun(n_simul, run_scheme1())
})
save(results_scheme1, file = "../data/scheme1_results.Rdata")

png(
  "../output/scheme1_predict_probab_distrib.png",
  width = 800, height = 1400, res = 200)
plot_simul(results_scheme1)
dev.off()

x <- summ_simul(results_scheme1)
write.table(
  x, "../output/table1.csv", sep = ",",
  row.names = FALSE, col.names = FALSE
)



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
  dat <- dplyr::tibble(base_score = rnorm(n_pt)) %>%
    dplyr::mutate(
      marker3_s_gt0 = rnorm(n_pt, 3.65, 1),
      marker3_s_le0 = rnorm(n_pt, 3.55, 1),
      marker3 = ifelse(base_score > 0, marker3_s_gt0, marker3_s_le0),
      marker4_s_gt0 = rnorm(n_pt, 0.05, 1),
      marker4_s_le0 = rnorm(n_pt, -0.05, 1),
      marker4 = ifelse(base_score > 0, marker4_s_gt0, marker4_s_le0),
      dis = ifelse(
        runif(n_pt) < probab_dis(base_score, marker3, marker4), 1L, 0L)
    ) %>%
    dplyr::select(-dplyr::matches("(gt0|le0)$"))

  # Split into training and validation sets.
  train <- dat[1:(n_pt / 2), ]
  valid <- dat[((n_pt / 2) + 1):n_pt, ]

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

system.time({
set.seed(639320)
results_scheme2 <- purrr::rerun(n_simul, run_scheme2())
})
save(results_scheme2, file = "../data/scheme2_results.Rdata")

png(
  "../output/scheme2_predict_probab_distrib.png",
  width = 800, height = 1400, res = 200)
plot_simul(results_scheme2)
dev.off()

x <- summ_simul(results_scheme2)
write.table(
  x, "../output/table_s1.csv", sep = ",",
  row.names = FALSE, col.names = FALSE
)


# Scheme 3 ------------------------

run_scheme3 <- function() {
  # Define the (true) disease risk model for this simulation scheme.
  probab_dis <- function(s, m5, m6, m7, m8, m9, m10) {
    # Define the Gaussian kernel function used for this disease model.
    kernel <- function(x) exp(-x^2 / 0.5)
    # Return the probability of disease.
    arm::invlogit(
      -3 + (2 * s) + (0.7 * (m5 + m6 + m7)) +
      (0.75 * kernel(s) * (m8 + m9 + m10))
    )
  }

  # Simulate a data set containing the covariates and disease state.
  dat <- dplyr::tibble(base_score = rnorm(n_pt)) %>%
    dplyr::mutate(
      marker5 = ifelse(
        base_score > 0,
        sample(0:1, n_pt, replace = TRUE, prob = c(.15, .85)),
        sample(0:1, n_pt, replace = TRUE, prob = c(.25, .75))
      ),
      marker6 = ifelse(
        base_score > 0,
        sample(0:1, n_pt, replace = TRUE, prob = c(.15, .85)),
        sample(0:1, n_pt, replace = TRUE, prob = c(.25, .75))
      ),
      marker7 = ifelse(
        base_score > 0,
        sample(0:1, n_pt, replace = TRUE, prob = c(.15, .85)),
        sample(0:1, n_pt, replace = TRUE, prob = c(.25, .75))
      ),
      marker8 = ifelse(
        base_score > 0,
        sample(0:1, n_pt, replace = TRUE, prob = c(.15, .85)),
        sample(0:1, n_pt, replace = TRUE, prob = c(.25, .75))
      ),
      marker9 = ifelse(
        base_score > 0,
        sample(0:1, n_pt, replace = TRUE, prob = c(.15, .85)),
        sample(0:1, n_pt, replace = TRUE, prob = c(.25, .75))
      ),
      marker10 = ifelse(
        base_score > 0,
        sample(0:1, n_pt, replace = TRUE, prob = c(.15, .85)),
        sample(0:1, n_pt, replace = TRUE, prob = c(.25, .75))
      ),
      dis = ifelse(
        runif(n_pt) < probab_dis(
          base_score, marker5, marker6, marker7, marker8, marker9,
          marker10
        ),
        1L, 0L)
    )

  # Split into training and validation sets.
  train <- dat[1:(n_pt / 2), ]
  valid <- dat[((n_pt / 2) + 1):n_pt, ]

  # Fit each model.
  model1 <- glm(dis ~ base_score, family = binomial, data = train)
  model2 <- glm(
    dis ~ base_score + marker5 + marker6 + marker7,
    family = binomial, data = train)
  model3 <- glm(
    dis ~ base_score + marker8 + marker9 + marker10,
    family = binomial, data = train)

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

system.time({
set.seed(866289)
results_scheme3 <- purrr::rerun(n_simul, run_scheme3())
})
save(results_scheme3, file = "../data/scheme3_results.Rdata")

png(
  "../output/scheme3_predict_probab_distrib.png",
  width = 800, height = 1400, res = 200)
plot_simul(results_scheme3)
dev.off()

x <- summ_simul(results_scheme3)
write.table(
  x, "../output/table_s2.csv", sep = ",",
  row.names = FALSE, col.names = FALSE
)

## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- eval=FALSE--------------------------------------------------------------
#  citation("SynthETIC")

## -----------------------------------------------------------------------------
library(SynthETIC)
set.seed(20200131)

## -----------------------------------------------------------------------------
set_parameters(ref_claim = 200000, time_unit = 1/4)
ref_claim <- return_parameters()[1]
time_unit <- return_parameters()[2]

## -----------------------------------------------------------------------------
years <- 10
I <- years / time_unit
E <- c(rep(12000, I)) # effective annual exposure rates
lambda <- c(rep(0.03, I))

## -----------------------------------------------------------------------------
# Number of claims occurring for each period i
# shorter equivalent code:
# n_vector <- claim_frequency()
n_vector <- claim_frequency(I = I, E = E, freq = lambda)
n_vector

# Occurrence time of each claim r, for each period i
occurrence_times <- claim_occurrence(frequency_vector = n_vector)
occurrence_times[[1]]

## -----------------------------------------------------------------------------
## input parameters
years_tmp <- 10
I_tmp <- years_tmp / time_unit
# set linearly increasing exposure, ...
E_tmp <- c(rep(12000, I)) + seq(from = 0, by = 100, length = I)
# and constant frequency per unit of exposure
lambda_tmp <- c(rep(0.03, I))

## output
# Number of claims occurring for each period i
n_vector_tmp <- claim_frequency(I = I_tmp, E = E_tmp, freq = lambda_tmp)
n_vector_tmp

# Occurrence time of each claim r, for each period i
occurrence_times_tmp <- claim_occurrence(frequency_vector = n_vector_tmp)
occurrence_times_tmp[[1]]

## -----------------------------------------------------------------------------
# simulate claim frequencies from negative binomial
# 1. using type-"r" specification (default)
claim_frequency(I = I, simfun = rnbinom, size = 100, mu = 100)
# 2. using type-"p" specification, equivalent to above
claim_frequency(I = I, simfun = pnbinom, type = "p", size = 100, mu = 100)

# simulate claim frequencies from zero-truncated Poisson
claim_frequency(I = I, simfun = actuar::rztpois, lambda = 90)
claim_frequency(I = I, simfun = actuar::pztpois, type = "p", lambda = 90)

## -----------------------------------------------------------------------------
claim_frequency(I = I, simfun = actuar::rztpois, lambda = time_unit * E_tmp * lambda_tmp)

## ----dpi=150, fig.width=7, fig.height=4, out.width=650------------------------
# sampling from non-homogeneous Poisson process
rnhpp.count <- function(no_periods) {
  rate <- 3000
  intensity <- function(x) {
    # e.g. cyclical Poisson process
    0.03 * (sin(x * pi / 2) / 4 + 1)
  }
  claim_times <- poisson::nhpp.event.times(rate, no_periods * rate * 2, intensity)
  
  as.numeric(table(cut(claim_times, breaks = 0:no_periods)))
}

n_vector_tmp <- claim_frequency(I = I, simfun = rnhpp.count)
plot(x = 1:I, y = n_vector_tmp, type = "l",
     main = "Claim frequency simulated from a cyclical Poisson process",
     xlab = "Occurrence period", ylab = "# Claims")

## -----------------------------------------------------------------------------
rate_tmp <- 3000
intensity_tmp <- function(x) {
  # e.g. cyclical Poisson process
  0.03 * (sin(x * pi / 2) / 4 + 1)
}
x_tmp <- poisson::nhpp.event.times(rate_tmp, I * rate_tmp, intensity_tmp)
event_times_tmp <- x_tmp[x_tmp <= I]

# Number of claims occurring for each period i
# by counting the number of event times in each interval (i, i + 1)
n_vector_tmp <- as.numeric(table(cut(event_times_tmp, breaks = 0:I)))
n_vector_tmp

# Occurrence time of each claim r, for each period i
occurrence_times_tmp <- to_SynthETIC(x = event_times_tmp, 
                                     frequency_vector = n_vector_tmp)
occurrence_times_tmp[[1]]

## -----------------------------------------------------------------------------
# use a power normal S^0.2 ~ N(9.5, 3), left truncated at 30
# this is the default distribution driving the claim_size() function
S_df <- function(s) {
  # truncate and rescale
  if (s < 30) {
    return(0)
  } else {
    p_trun <- pnorm(s^0.2, 9.5, 3) - pnorm(30^0.2, 9.5, 3)
    p_rescaled <- p_trun/(1 - pnorm(30^0.2, 9.5, 3))
    return(p_rescaled)
  }
}

## -----------------------------------------------------------------------------
# shorter equivalent: claim_sizes <- claim_size(frequency_vector = n_vector)
claim_sizes <- claim_size(frequency_vector = n_vector, 
                          simfun = S_df, type = "p", range = c(0, 1e24))
claim_sizes[[1]]

## ----dpi=150, fig.width=7, fig.height=4, out.width=650------------------------
## weibull
# estimate the weibull parameters to achieve the mean and cv matching that of
# the built-in test claim dataset
claim_size_mean <- mean(test_claim_dataset$claim_size)
claim_size_cv <- cv(test_claim_dataset$claim_size)
weibull_shape <- get_Weibull_parameters(target_mean = claim_size_mean, 
                                        target_cv = claim_size_cv)[1]
weibull_scale <- get_Weibull_parameters(target_mean = claim_size_mean, 
                                        target_cv = claim_size_cv)[2]
# simulate claim sizes with the estimated parameters
claim_sizes_weibull <- claim_size(frequency_vector = n_vector,
                                  simfun = rweibull, 
                                  shape = weibull_shape, scale = weibull_scale)
# plot empirical CDF
plot(ecdf(unlist(test_claim_dataset$claim_size)), xlim = c(0, 2000000),
     main = "Empirical distribution of simulated claim sizes",
     xlab = "Individual claim size")
plot(ecdf(unlist(claim_sizes_weibull)), add = TRUE, col = 2)


## inverse Gaussian
# modify actuar::rinvgauss (left truncate it @30 and right censor it @5,000,000)
rinvgauss_censored <- function(n) {
  s <- actuar::rinvgauss(n, mean = 180000, dispersion = 0.5e-5)
  while (any(s < 30 | s > 5000000)) {
    for (j in which(s < 30 | s > 5000000)) {
      # for rejected values, resample
      s[j] <- actuar::rinvgauss(1, mean = 180000, dispersion = 0.5e-5)
    }
  }
  s
}
# simulate from the modified inverse Gaussian distribution
claim_sizes_invgauss <- claim_size(frequency_vector = n_vector, simfun = rinvgauss_censored)

# plot empirical CDF
plot(ecdf(unlist(claim_sizes_invgauss)), add = TRUE, col = 3)
legend.text <- c("Power normal", "Weibull", "Inverse Gaussian")
legend("bottomright", legend.text, col = 1:3, lty = 1, bty = "n")

## ----dpi=150, fig.width=7, fig.height=4, out.width=650------------------------
# define the random generation function to simulate from the gamma GLM
sim_GLM <- function(n) {
  # simulate covariates
  age <- sample(20:70, size = n, replace = T)
  mu <- exp(27 - 0.768 * age + 0.008 * age^2)
  rgamma(n, shape = 10, scale = mu / 10)
}

claim_sizes_GLM <- claim_size(frequency_vector = n_vector, simfun = sim_GLM)
plot(ecdf(unlist(claim_sizes_GLM)), xlim = c(0, 2000000),
     main = "Empirical distribution of claim sizes simulated from GLM",
     xlab = "Individual claim size")

## ----eval=FALSE---------------------------------------------------------------
#  # install.packages("CASdatasets", repos = "http://cas.uqam.ca/pub/", type = "source")
#  library(CASdatasets)
#  data("ausautoBI8999")
#  boot <- sample(ausautoBI8999$AggClaim, size = sum(n_vector), replace = TRUE)
#  claim_sizes_bootstrap <- to_SynthETIC(boot, frequency_vector = n_vector)

## ----eval=FALSE---------------------------------------------------------------
#  sim_boot <- function(n) {
#    sample(ausautoBI8999$AggClaim, size = n, replace = TRUE)
#  }
#  claim_sizes_bootstrap <- claim_size(frequency_vector = n_vector, simfun = sim_boot)

## -----------------------------------------------------------------------------
## input
# specify the Weibull parameters as a function of claim_size and occurrence_period
notidel_param <- function(claim_size, occurrence_period) {
  # NOTE: users may add to, but not remove these two arguments (claim_size, 
  # occurrence_period) as they are part of SynthETIC's internal structure
  
  # specify the target mean and target coefficient of variation
  target_mean <- min(3, max(1, 2-(log(claim_size/(0.50 * ref_claim)))/3))/4 / time_unit
  target_cv <- 0.70
  # convert to Weibull parameters
  shape <- get_Weibull_parameters(target_mean, target_cv)[1]
  scale <- get_Weibull_parameters(target_mean, target_cv)[2]
  
  c(shape = shape, scale = scale)
}

## output
notidel <- claim_notification(n_vector, claim_sizes, 
                              rfun = rweibull, paramfun = notidel_param)

## ----dpi=150, fig.width=7, fig.height=4, out.width=650------------------------
## input
# specify the transformed gamma parameters as a function of claim_size and occurrence_period
trgamma_param <- function(claim_size, occurrence_period, rate) {
  c(shape1 = max(1, claim_size / ref_claim),
    shape2 = 1 - occurrence_period / 200,
    rate = rate)
}

## output
# simulate notification delays from the transformed gamma
notidel_trgamma <- claim_notification(n_vector, claim_sizes, 
                                      rfun = actuar::rtrgamma, 
                                      paramfun = trgamma_param, rate = 2)

# graphically compare the result with the default Weibull distribution
plot(ecdf(unlist(notidel)), xlim = c(0, 15),
     main = "Empirical distribution of simulated notification delays",
     xlab = "Notification delay (in quarters)")
plot(ecdf(unlist(notidel_trgamma)), add = TRUE, col = 2)
legend.text <- c("Weibull (default)", "Transformed gamma")
legend("bottomright", legend.text, col = 1:2, lty = 1, bty = "n")

## -----------------------------------------------------------------------------
rmixed_notidel <- function(n, claim_size) {
  # consider a mixture distribution
  # equal probability of sampling from x (Weibull) or y (transformed gamma)
  x_selected <- sample(c(T, F), size = n, replace = TRUE)
  x <- rweibull(n, shape = 2, scale = 1)
  y <- actuar::rtrgamma(n, shape1 = min(1, claim_size / ref_claim), shape2 = 0.8, rate = 2)
  result <- length(n)
  result[x_selected] <- x[x_selected]; result[!x_selected] <- y[!x_selected]
  
  return(result)
}

## -----------------------------------------------------------------------------
rmixed_params <- function(claim_size, occurrence_period) {
  # claim_size is the only "parameter" required for rmixed_notidel
  c(claim_size = claim_size)
}

## -----------------------------------------------------------------------------
notidel_mixed <- claim_notification(n_vector, claim_sizes, rfun = rmixed_notidel)

## -----------------------------------------------------------------------------
notidel_mixed <- claim_notification(n_vector, claim_sizes, 
                                    rfun = rmixed_notidel, paramfun = rmixed_params)

## -----------------------------------------------------------------------------
## input
# specify the Weibull parameters as a function of claim_size and occurrence_period
setldel_param <- function(claim_size, occurrence_period) {
  # NOTE: users may add to, but not remove these two arguments (claim_size, 
  # occurrence_period) as they are part of SynthETIC's internal structure
  
  # specify the target Weibull mean
  if (claim_size < (0.10 * ref_claim) & occurrence_period >= 21) {
    a <- min(0.85, 0.65 + 0.02 * (occurrence_period - 21))
  } else {
    a <- max(0.85, 1 - 0.0075 * occurrence_period)
  }
  mean_quarter <- a * min(25, max(1, 6 + 4*log(claim_size/(0.10 * ref_claim))))
  target_mean <- mean_quarter / 4 / time_unit
  
  # specify the target Weibull coefficient of variation
  target_cv <- 0.60

  c(shape = get_Weibull_parameters(target_mean, target_cv)[1, ],
    scale = get_Weibull_parameters(target_mean, target_cv)[2, ])
}

## output
# simulate the settlement delays from the Weibull with parameters above
setldel <- claim_closure(n_vector, claim_sizes, rfun = rweibull, paramfun = setldel_param)
setldel[[1]]

## -----------------------------------------------------------------------------
## input
# an extended parameter function for the simulation of settlement delays
setldel_param_extd <- function(claim_size, occurrence_period, notidel) {
  
  # specify the target Weibull mean
  if (claim_size < (0.10 * ref_claim) & occurrence_period >= 21) {
    a <- min(0.85, 0.65 + 0.02 * (occurrence_period - 21))
  } else {
    a <- max(0.85, 1 - 0.0075 * occurrence_period)
  }
  mean_quarter <- a * min(25, max(1, 6 + 4*log(claim_size/(0.10 * ref_claim))))
  # suppose the setldel mean is linearly related to the notidel of the claim
  target_mean <- (mean_quarter + notidel) / 4 / time_unit
  
  # specify the target Weibull coefficient of variation
  target_cv <- 0.60

  c(shape = get_Weibull_parameters(target_mean, target_cv)[1, ],
    scale = get_Weibull_parameters(target_mean, target_cv)[2, ])
}

## -----------------------------------------------------------------------------
## output
# simulate the settlement delays from the Weibull with parameters above
notidel_vect <- unlist(notidel) # convert to a vector
setldel_extd <- claim_closure(n_vector, claim_sizes, rfun = rweibull, 
                              paramfun = setldel_param_extd,
                              notidel = notidel_vect)
setldel_extd[[1]]

## -----------------------------------------------------------------------------
## input
# the default random generating function
rmixed_payment_no <- function(n, claim_size, claim_size_benchmark_1, claim_size_benchmark_2) {
  # construct the range indicators
  test_1 <- (claim_size_benchmark_1 < claim_size & claim_size <= claim_size_benchmark_2)
  test_2 <- (claim_size > claim_size_benchmark_2)

  # if claim_size <= claim_size_benchmark_1
  no_pmt <- sample(c(1, 2), size = n, replace = T, prob = c(1/2, 1/2))
  # if claim_size is between the two benchmark values
  no_pmt[test_1] <- sample(c(2, 3), size = sum(test_1), replace = T, prob = c(1/3, 2/3))
  # if claim_size > claim_size_benchmark_2
  no_pmt_mean <- pmin(8, 4 + log(claim_size/claim_size_benchmark_2))
  prob <- 1 / (no_pmt_mean - 3)
  no_pmt[test_2] <- stats::rgeom(n = sum(test_2), prob = prob[test_2]) + 4

  no_pmt
}

## -----------------------------------------------------------------------------
## output
no_payments <- claim_payment_no(n_vector, claim_sizes, rfun = rmixed_payment_no,
                                claim_size_benchmark_1 = 0.0375 * ref_claim,
                                claim_size_benchmark_2 = 0.075 * ref_claim)
no_payments[[1]]

## ---- eval=FALSE--------------------------------------------------------------
#  no_payments <- claim_payment_no(n_vector, claim_sizes)

## -----------------------------------------------------------------------------
no_payments_tmp <- claim_payment_no(n_vector, claim_sizes,
                                    claim_size_benchmark_2 = 0.1 * ref_claim)

## -----------------------------------------------------------------------------
## input
paymentNo_param <- function(claim_size) {
  no_pmt_mean <- pmax(4, pmin(8, 4 + log(claim_size / 15000)))
  c(lambda = no_pmt_mean - 3)
}

## output
no_payments_pois <- claim_payment_no(
  n_vector, claim_sizes, rfun = actuar::rztpois, paramfun = paymentNo_param)
table(unlist(no_payments_pois))

## -----------------------------------------------------------------------------
claim_dataset <- generate_claim_dataset(
  frequency_vector = n_vector,
  occurrence_list = occurrence_times,
  claim_size_list = claim_sizes,
  notification_list = notidel,
  settlement_list = setldel,
  no_payments_list = no_payments
)
str(claim_dataset)

## -----------------------------------------------------------------------------
str(test_claim_dataset)

## -----------------------------------------------------------------------------
## input
rmixed_payment_size <- function(n, claim_size) {
  # n = number of simulations, here n should be the number of partial payments
  if (n >= 4) {
    # 1) Simulate the "complement" of the proportion of total claim size 
    #    represented by the last two payments
    p_mean <- 1 - min(0.95, 0.75 + 0.04*log(claim_size/(0.10 * ref_claim)))
    p_CV <- 0.20
    p_parameters <- get_Beta_parameters(target_mean = p_mean, target_cv = p_CV)
    last_two_pmts_complement <- stats::rbeta(
      1, shape1 = p_parameters[1], shape2 = p_parameters[2])
    last_two_pmts <- 1 - last_two_pmts_complement
    
    # 2) Simulate the proportion of last_two_pmts paid in the second last payment
    q_mean <- 0.9
    q_CV <- 0.03
    q_parameters <- get_Beta_parameters(target_mean = q_mean, target_cv = q_CV)
    q <- stats::rbeta(1, shape1 = q_parameters[1], shape2 = q_parameters[2])

    # 3) Calculate the respective proportions of claim amount paid in the 
    #    last 2 payments
    p_second_last <- q * last_two_pmts
    p_last <- (1-q) * last_two_pmts

    # 4) Simulate the "unnormalised" proportions of claim amount paid 
    #    in the first (m - 2) payments
    p_unnorm_mean <- last_two_pmts_complement/(n - 2)
    p_unnorm_CV <- 0.10
    p_unnorm_parameters <- get_Beta_parameters(
      target_mean = p_unnorm_mean, target_cv = p_unnorm_CV)
    amt <- stats::rbeta(
      n - 2, shape1 = p_unnorm_parameters[1], shape2 = p_unnorm_parameters[2])

    # 5) Normalise the proportions simulated in step 4
    amt <- last_two_pmts_complement * (amt/sum(amt))
    # 6) Attach the last 2 proportions, p_second_last and p_last
    amt <- append(amt, c(p_second_last, p_last))
    # 7) Multiply by claim_size to obtain the actual payment amounts
    amt <- claim_size * amt
    
  } else if (n == 2 | n == 3) {
    p_unnorm_mean <- 1/n
    p_unnorm_CV <- 0.10
    p_unnorm_parameters <- get_Beta_parameters(
      target_mean = p_unnorm_mean, target_cv = p_unnorm_CV)
    amt <- stats::rbeta(
      n, shape1 = p_unnorm_parameters[1], shape2 = p_unnorm_parameters[2])
    # Normalise the proportions and multiply by claim_size to obtain the actual payment amounts
    amt <- claim_size * amt/sum(amt)

  } else {
    # when there is a single payment
    amt <- claim_size
  }
  return(amt)
}

## output
payment_sizes <- claim_payment_size(n_vector, claim_sizes, no_payments,
                                    rfun = rmixed_payment_size)
payment_sizes[[1]][[1]]

## ---- eval=FALSE--------------------------------------------------------------
#  payment_sizes <- claim_payment_size(n_vector, claim_sizes, no_payments)

## -----------------------------------------------------------------------------
## input
unif_payment_size <- function(n, claim_size) {
  prop <- runif(n)
  prop.normalised <- prop / sum(prop)
  
  return(claim_size * prop)
}

## output
# note that we don't need to specify a paramfun as rfun is directly a function
# of claim_size
payment_sizes_unif <- claim_payment_size(n_vector, claim_sizes, no_payments,
                                         rfun = unif_payment_size)
payment_sizes_unif[[1]][[1]]

## -----------------------------------------------------------------------------
## input
r_pmtdel <- function(n, claim_size, setldel, setldel_mean) {
  result <- c(rep(NA, n))

  # First simulate the unnormalised values of d, sampled from a Weibull distribution
  if (n >= 4) {
    # 1) Simulate the last payment delay
    unnorm_d_mean <- (1 / 4) / time_unit
    unnorm_d_cv <- 0.20
    parameters <- get_Weibull_parameters(target_mean = unnorm_d_mean, target_cv = unnorm_d_cv)
    result[n] <- stats::rweibull(1, shape = parameters[1], scale = parameters[2])

    # 2) Simulate all the other payment delays
    for (i in 1:(n - 1)) {
      unnorm_d_mean <- setldel_mean / n
      unnorm_d_cv <- 0.35
      parameters <- get_Weibull_parameters(target_mean = unnorm_d_mean, target_cv = unnorm_d_cv)
      result[i] <- stats::rweibull(1, shape = parameters[1], scale = parameters[2])
    }

  } else {
    for (i in 1:n) {
      unnorm_d_mean <- setldel_mean / n
      unnorm_d_cv <- 0.35
      parameters <- get_Weibull_parameters(target_mean = unnorm_d_mean, target_cv = unnorm_d_cv)
      result[i] <- stats::rweibull(1, shape = parameters[1], scale = parameters[2])
    }
  }

  # Normalise d such that sum(inter-partial delays) = settlement delay
  # To make sure that the pmtdels add up exactly to setldel, we treat the last one separately
  result[1:n-1] <- (setldel/sum(result)) * result[1:n-1]
  result[n] <- setldel - sum(result[1:n-1])

  return(result)
}

param_pmtdel <- function(claim_size, setldel, occurrence_period) {
  # mean settlement delay
  if (claim_size < (0.10 * ref_claim) & occurrence_period >= 21) {
    a <- min(0.85, 0.65 + 0.02 * (occurrence_period - 21))
  } else {
    a <- max(0.85, 1 - 0.0075 * occurrence_period)
  }
  mean_quarter <- a * min(25, max(1, 6 + 4*log(claim_size/(0.10 * ref_claim))))
  target_mean <- mean_quarter / 4 / time_unit

  c(claim_size = claim_size,
    setldel = setldel,
    setldel_mean = target_mean)
}

## output
payment_delays <- claim_payment_delay(
  n_vector, claim_sizes, no_payments, setldel,
  rfun = r_pmtdel, paramfun = param_pmtdel,
  occurrence_period = rep(1:I, times = n_vector))

# payment times on a continuous time scale
payment_times <- claim_payment_time(n_vector, occurrence_times, notidel, payment_delays)
# payment times in periods
payment_periods <- claim_payment_time(n_vector, occurrence_times, notidel, payment_delays,
                                      discrete = TRUE)
cbind(payment_delays[[1]][[1]], payment_times[[1]][[1]], payment_periods[[1]][[1]])

## -----------------------------------------------------------------------------
# Base inflation: a vector of quarterly rates
# In this demo we set base inflation to be at 2% p.a. constant for both past and future
# Users can choose to randominise the future rates if they wish
demo_rate <- (1 + 0.02)^(1/4) - 1
base_inflation_past <- rep(demo_rate, times = 40)
base_inflation_future <- rep(demo_rate, times = 40)
base_inflation_vector <- c(base_inflation_past, base_inflation_future)

# Superimposed inflation:
# 1) With respect to occurrence "time" (continuous scale)
SI_occurrence <- function(occurrence_time, claim_size) {
  if (occurrence_time <= 20 / 4 / time_unit) {1}
  else {1 - 0.4*max(0, 1 - claim_size/(0.25 * ref_claim))}
}
# 2) With respect to payment "time" (continuous scale)
# -> compounding by user-defined time unit
SI_payment <- function(payment_time, claim_size) {
  period_rate <- (1 + 0.30)^(time_unit) - 1
  beta <- period_rate * max(0, 1 - claim_size/ref_claim)
  (1 + beta)^payment_time
}

## -----------------------------------------------------------------------------
# shorter equivalent code:
# payment_inflated <- claim_payment_inflation(
#   n_vector, payment_sizes, payment_times, occurrence_times, claim_sizes, 
#   base_inflation_vector)
payment_inflated <- claim_payment_inflation(
  n_vector,
  payment_sizes,
  payment_times,
  occurrence_times,
  claim_sizes,
  base_inflation_vector,
  SI_occurrence,
  SI_payment
)
cbind(payment_sizes[[1]][[1]], payment_inflated[[1]][[1]])

## -----------------------------------------------------------------------------
# construct a "claims" object to store all the simulated quantities
all_claims <- claims(
  frequency_vector = n_vector,
  occurrence_list = occurrence_times,
  claim_size_list = claim_sizes,
  notification_list = notidel,
  settlement_list = setldel,
  no_payments_list = no_payments,
  payment_size_list = payment_sizes,
  payment_delay_list = payment_delays,
  payment_time_list = payment_times,
  payment_inflated_list = payment_inflated
)
transaction_dataset <- generate_transaction_dataset(
  all_claims,
  adjust = FALSE # to keep the original (potentially out-of-bound) simulated payment times
)
str(transaction_dataset)

## -----------------------------------------------------------------------------
str(test_transaction_dataset)

## ---- eval=FALSE--------------------------------------------------------------
#  claim_output(
#    frequency_vector = ,
#    payment_time_list = ,
#    payment_size_list = ,
#    aggregate_level = 1,
#    incremental = TRUE,
#    future = TRUE,
#    adjust = TRUE
#  )

## -----------------------------------------------------------------------------
# 1. Constant dollar value INCREMENTAL triangle
output <- claim_output(n_vector, payment_times, payment_sizes,
                       incremental = TRUE)

# 2. Constant dollar value CUMULATIVE triangle
output_cum <- claim_output(n_vector, payment_times, payment_sizes,
                           incremental = FALSE)

# 3. Actual (i.e. inflated) INCREMENTAL triangle
output_actual <- claim_output(n_vector, payment_times, payment_inflated,
                              incremental = TRUE)

# 4. Actual (i.e. inflated) CUMULATIVE triangle
output_actual_cum <- claim_output(n_vector, payment_times, payment_inflated,
                                  incremental = FALSE)

# Aggregate at a yearly level
claim_output(n_vector, payment_times, payment_sizes, aggregate_level = 4)

## -----------------------------------------------------------------------------
# output the past cumulative triangle
cumtri <- claim_output(n_vector, payment_times, payment_sizes,
                       aggregate_level = 4, incremental = FALSE, future = FALSE)
# calculate the age to age factors
selected <- attr(ChainLadder::ata(cumtri), "vwtd")
# complete the triangle
CL_prediction <- cumtri
J <- nrow(cumtri)
for (i in 2:J) {
  for (j in (J - i + 2):J) {
    CL_prediction[i, j] <- CL_prediction[i, j - 1] * selected[j - 1]
  }
}

CL_prediction

## ----dpi=150, fig.width=7, fig.height=6, out.width=650------------------------
plot(test_claims_object)
# compare with the "full complete picture"
plot(test_claims_object, adjust = FALSE)

## ----dpi=150, fig.width=7, fig.height=6, out.width=650------------------------
# plot by occurrence and development years
plot(test_claims_object, by_year = TRUE)

## ---- eval = FALSE------------------------------------------------------------
#  times <- 100
#  results_all <- vector("list")
#  for (i in 1:times) {
#    # Module 1: Claim occurrence
#    n_vector <- claim_frequency(I, E, lambda)
#    occurrence_times <- claim_occurrence(n_vector)
#    # Module 2: Claim size
#    claim_sizes <- claim_size(n_vector, S_df, type = "p", range = c(0, 1e24))
#    # Module 3: Claim notification
#    notidel <- claim_notification(n_vector, claim_sizes, paramfun = notidel_param)
#    # Module 4: Claim settlement
#    setldel <- claim_closure(n_vector, claim_sizes, paramfun = setldel_param)
#    # Module 5: Claim payment count
#    no_payments <- claim_payment_no(n_vector, claim_sizes, rfun = rmixed_payment_no,
#                                    claim_size_benchmark_1 = 0.0375 * ref_claim,
#                                    claim_size_benchmark_2 = 0.075 * ref_claim)
#    # Module 6: Claim payment size
#    payment_sizes <- claim_payment_size(n_vector, claim_sizes, no_payments,
#                                        rfun = rmixed_payment_size)
#    # Module 7: Claim payment time
#    payment_delays <- claim_payment_delay(n_vector, claim_sizes, no_payments, setldel,
#                                          rfun = r_pmtdel, paramfun = param_pmtdel,
#                                          occurrence_period = rep(1:I, times = n_vector))
#    payment_times <- claim_payment_time(n_vector, occurrence_times, notidel, payment_delays)
#    # Module 8: Claim inflation
#    payment_inflated <- claim_payment_inflation(
#      n_vector, payment_sizes, payment_times, occurrence_times,
#      claim_sizes, base_inflation_vector, SI_occurrence, SI_payment)
#  
#    results_all[[i]] <- generate_transaction_dataset(
#      claims(
#        frequency_vector = n_vector,
#        occurrence_list = occurrence_times,
#        claim_size_list = claim_sizes,
#        notification_list = notidel,
#        settlement_list = setldel,
#        no_payments_list = no_payments,
#        payment_size_list = payment_sizes,
#        payment_delay_list = payment_delays,
#        payment_time_list = payment_times,
#        payment_inflated_list = payment_inflated),
#      # adjust = FALSE to retain the original simulated times
#      adjust = FALSE)
#  }

## ----dpi=150, fig.width=7, fig.height=6, out.width=650------------------------
start.time <- proc.time()
times <- 10

# increase exposure to E*times to get the same results as the aggregation of
# multiple simulation runs
n_vector <- claim_frequency(I, E = E * times, lambda)
occurrence_times <- claim_occurrence(n_vector)
claim_sizes <- claim_size(n_vector)
notidel <- claim_notification(n_vector, claim_sizes, paramfun = notidel_param)
setldel <- claim_closure(n_vector, claim_sizes, paramfun = setldel_param)
no_payments <- claim_payment_no(n_vector, claim_sizes, rfun = rmixed_payment_no,
                                claim_size_benchmark_1 = 0.0375 * ref_claim,
                                claim_size_benchmark_2 = 0.075 * ref_claim)
payment_sizes <- claim_payment_size(n_vector, claim_sizes, no_payments, rmixed_payment_size)
payment_delays <- claim_payment_delay(n_vector, claim_sizes, no_payments, setldel,
                                      rfun = r_pmtdel, paramfun = param_pmtdel,
                                      occurrence_period = rep(1:I, times = n_vector))
payment_times <- claim_payment_time(n_vector, occurrence_times, notidel, payment_delays)
payment_inflated <- claim_payment_inflation(
  n_vector, payment_sizes, payment_times, occurrence_times,
  claim_sizes, base_inflation_vector, SI_occurrence, SI_payment)

all_claims <- claims(
  frequency_vector = n_vector,
  occurrence_list = occurrence_times,
  claim_size_list = claim_sizes,
  notification_list = notidel,
  settlement_list = setldel,
  no_payments_list = no_payments,
  payment_size_list = payment_sizes,
  payment_delay_list = payment_delays,
  payment_time_list = payment_times,
  payment_inflated_list = payment_inflated
)
plot(all_claims, adjust = FALSE) +
  ggplot2::labs(subtitle = paste("With", times, "simulations"))
proc.time() - start.time

## ---- eval=FALSE--------------------------------------------------------------
#  plot(claims, by_year = , inflated = , adjust = )


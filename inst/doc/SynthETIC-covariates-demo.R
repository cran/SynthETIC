## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=FALSE---------------------------------------------------------------
# citation("SynthETIC")

## -----------------------------------------------------------------------------
library(SynthETIC)
set.seed(20200131)

set_parameters(ref_claim = 200000, time_unit = 1/4)
ref_claim <- return_parameters()[1]
time_unit <- return_parameters()[2]

years <- 10
I <- years / time_unit
E <- c(rep(12000, I)) # effective annual exposure rates
lambda <- c(rep(0.03, I))

# Modelling Steps 1-2
n_vector <- claim_frequency(I = I, E = E, freq = lambda)
occurrence_times <- claim_occurrence(frequency_vector = n_vector)
claim_sizes <- claim_size(frequency_vector = n_vector)

## -----------------------------------------------------------------------------
test_covariates_obj <- SynthETIC::test_covariates_obj
print(test_covariates_obj)

## -----------------------------------------------------------------------------
claim_size_covariates <- claim_size_adj(test_covariates_obj, claim_sizes)
covariates_data_obj <- claim_size_covariates$covariates_data
head(data.frame(covariates_data_obj$data))

## -----------------------------------------------------------------------------
claim_size_w_cov <- claim_size_covariates$claim_size_adj
claim_size_w_cov[[1]]

## -----------------------------------------------------------------------------
generate_claims_dataset <- function(claim_size_list) {
    
    # SynthETIC Steps 3-5
    notidel <- claim_notification(n_vector, claim_size_list)
    setldel <- claim_closure(n_vector, claim_size_list)
    no_payments <- claim_payment_no(n_vector, claim_size_list)
    
    claim_dataset <- generate_claim_dataset(
      frequency_vector = n_vector,
      occurrence_list = occurrence_times,
      claim_size_list = claim_size_list,
      notification_list = notidel,
      settlement_list = setldel,
      no_payments_list = no_payments
    )
    
    claim_dataset
}

claim_dataset <- generate_claims_dataset(claim_size_list = claim_sizes)
claim_dataset_w_cov <- generate_claims_dataset(claim_size_list = claim_size_w_cov)

head(claim_dataset)
head(claim_dataset_w_cov)

## -----------------------------------------------------------------------------
factors_tmp <- list(
    "Vehicle Type" = c("Passenger", "Light Commerical", "Medium Goods", "Heavy Goods"),
    "Business Use" = c("Y", "N")
)

relativity_freq_tmp <- relativity_template(factors_tmp)
relativity_sev_tmp <- relativity_template(factors_tmp)

# Default Values
relativity_freq_tmp$relativity <- c(
    5, 1.5, 0.35, 0.25,
    1, 4,
    1, 0.6,
    0.35, 0.01,
    0.25, 0,
    2.5, 5
)

relativity_sev_tmp$relativity <- c(
    0.25, 0.75, 1, 3,
    1, 1,
    1, 1,
    1, 1,
    1, 1,
    1.3, 1
)

test_covariates_obj_veh <- covariates(factors_tmp)
test_covariates_obj_veh <- set.covariates_relativity(
    covariates = test_covariates_obj_veh, 
    relativity = relativity_freq_tmp, 
    freq_sev = "freq"
)
test_covariates_obj_veh <- set.covariates_relativity(
    covariates = test_covariates_obj_veh, 
    relativity = relativity_sev_tmp, 
    freq_sev = "sev"
)

claim_size_covariates_veh <- claim_size_adj(test_covariates_obj_veh, claim_sizes)

# Comparison of the same claim size except with adjustments due to covariates
data.frame(
    Claim_Size = head(round(claim_sizes[[1]]))
    ,Claim_Size_Original_Covariates = head(round(claim_size_covariates$claim_size_adj[[1]]))
    ,Claim_Size_New_Covariates = head(round(claim_size_covariates_veh$claim_size_adj[[1]]))
)

# Covariate Levels
head(claim_size_covariates$covariates_data$data)
head(claim_size_covariates_veh$covariates_data$data)


## -----------------------------------------------------------------------------
claim_sizes_known <- list(c(
    rexp(n = 100, rate = 1.5)
))

known_covariates_dataset <- data.frame(
    "Vehicle Type" = rep(rep(c("Passenger", "Light Commerical"), each = 25), times = 2),
    "Business Use" = c(rep("N", times = 50), rep("Y", times = 50))
)
colnames(known_covariates_dataset) <- c("Vehicle Type", "Business Use")

covariates_data_veh <- covariates_data(
    test_covariates_obj_veh, 
    data = known_covariates_dataset, 
    covariates_id = list(c(100:51, 1:50))
)

claim_sizes_adj_tmp <- claim_size_adj.fit(
    covariates_data = covariates_data_veh,
    claim_size = claim_sizes_known
)

head(claim_sizes_adj_tmp[[1]])



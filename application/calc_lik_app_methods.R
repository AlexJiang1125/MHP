###############################################
# Date Last modified: Feb.3, 2023 AZJ
# Program: calc_lik_app_methods.R
# Purpose: calculate likelihoods for all model fitting results for applied datasets.
###############################################

# Data inputs:
#   - the SP&500 intraday dataset
#   - file directory is: /gscratch/home/jiang14/mhp/data/sim001/data/data_app/data_app.RData
#   - fitted model .Rda files for all methods
# Data outputs
#   - evaluate marginal MHP likelihood method for all five methods
#   - needs to pre-process data
# Note:
#   - for sem/semtr: coded as "output" instead of "output_app"

# get parameter estimate out from the raw file, in matrix forms
preprocess <- function(output_app, method) {
  K <- 11
  # SGLD: get rid of burn in, get median for the rest
  if (method == "sgld") {
    mu.mean <- apply(simplify2array(output_app$mu_sgld[c(100:length(output_app$alpha_sgld))]), "median",MARGIN= c(1,2))
    alpha.mean <- apply(simplify2array(output_app$alpha_sgld[c(100:length(output_app$alpha_sgld))]), "median",MARGIN= c(1,2))
    beta.mean <- apply(simplify2array(output_app$beta_sgld[c(100:length(output_app$alpha_sgld))]), "median",MARGIN= c(1,2))
  } else if (method == "sem" | method == "semtr") {
    mu.mean <- output_app$mu_em
    alpha.mean <- output_app$alpha_em
    beta.mean <- output_app$beta_em
  } else {
    mu.mean <- output_app$eta_mu1[1, ] / output_app$eta_mu2[1, ]
    alpha.mean <- matrix(unlist(output_app$eta_alphas1), K, K, byrow = TRUE)/matrix(unlist(output_app$eta_alphas2), K, K, byrow = TRUE)
    beta.mean <- matrix(unlist(output_app$eta_betas1), K, K, byrow = TRUE)/matrix(unlist(output_app$eta_betas2), K, K, byrow = TRUE)
  }
  return(list(mu.mean = mu.mean, alpha.mean = alpha.mean, beta.mean = beta.mean))
}

# generate file name
generate_file_name <- function(method_id, seed_id) {
  filename <- paste0(RESULTS_PATH, NAME_PREFIX,
                     METHOD_NAMES[method_id], "_seed_", seed_id, "_dataset_1",
                     NAME_SUFFICES[method_id])
  return(filename)
}

# Hawkes Lik
MargLik <- function(mu,alpha,beta,out) {
  K <- 11
  n <- dim(out)[1]
  y <- out$timestamp
  res <- log(mu[out$dim[1]])
  Tdat <- 1000
  for (i in 2:n) {
    alpha_vec <- alpha[cbind(out$dim[1:(i-1)],rep(out$dim[i],(i-1)))]
    beta_vec <- beta[cbind(out$dim[1:(i-1)],rep(out$dim[i],(i-1)))]
    res <- res + log(sum(alpha_vec*beta_vec*exp(-beta_vec*(y[i] - y[1:(i-1)]))) + mu[out$dim[i]])
  }
  ns <- as.vector(table(out$dim))#c(length(out[out$dim==1,1]),length(out[out$dim==2,1]),length(out[out$dim==3,1]))
  res <- res - sum(mu)*Tdat - sum(apply(alpha, MARGIN = 1, FUN = sum)*ns)
  for (i in 1:K) {
    for (j in 1:K) {
      res <- res + alpha[i,j]*sum(exp(-beta[i,j]*(Tdat - out$timestamp[out$dim == i])))
    }
  }
  return(res)
}

# define paths
DATA_PATH <- "/gscratch/home/jiang14/mhp/data/sim001/data/data_app/data_app.RData"
RESULTS_PATH <- "/gscratch/csde/jiang14/output_app/s005/"

# name components
METHOD_NAMES <- c("sem", "semtr", "sgd", "sgdtr", "sgld")
NAME_PREFIX <- c("mhpvi")
NAME_SUFFICES <- c("_time10.RData", "_time10.RData", "newinit_time10.RData", "newinit_time10.RData", "_newinit_time10.RData")

# load in job id
taskid <- as.numeric(commandArgs(trailingOnly = TRUE)[1])
ids <- expand.grid(method_id = 1:5, seed_id = 1:16)
method_id <- ids$method_id[taskid]
seed_id <- ids$seed_id[taskid]

# loading data
load(DATA_PATH)
# loading model object
load(generate_file_name(method_id, seed_id))
if (method_id == 5) {
  output <- output_app
}
pars <- preprocess(output, METHOD_NAMES[method_id])

# calculate likelihood
mu.mean <- pars$mu.mean
alpha.mean <- pars$alpha.mean
beta.mean <- pars$beta.mean
lik <- MargLik(mu.mean, alpha.mean, beta.mean, out)

# write files
df_out <- data.frame(model = METHOD_NAMES[method_id], seed = seed_id, lik = lik)
output_filename <- paste0(RESULTS_PATH, "lik_app_", METHOD_NAMES[method_id], "_seed_", seed_id, "_10m.csv")
write.csv(df_out, file = output_filename)

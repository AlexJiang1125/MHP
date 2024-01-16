library(dplyr)
# helper function to generate the plot
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

# get parameter estimate out from the raw file, in matrix forms
preprocess <- function(output_app, method) {
  K <- 11
  # SGLD: get rid of burn in, get median for the rest
  if (method == "sgld") {
    mu.mean <- apply(simplify2array(output_app$mu_sgld[c(INIT_START_SGLD:length(output_app$alpha_sgld))]), "median",MARGIN= c(1,2))
    alpha.mean <- apply(simplify2array(output_app$alpha_sgld[c(INIT_START_SGLD:length(output_app$alpha_sgld))]), "median",MARGIN= c(1,2))
    beta.mean <- apply(simplify2array(output_app$beta_sgld[c(INIT_START_SGLD:length(output_app$alpha_sgld))]), "median",MARGIN= c(1,2))
  } else if (method == "sem" | method == "semtr") {
    mu.mean <- output_app$mu_em
    alpha.mean <- output_app$alpha_em
    beta.mean <- output_app$beta_em
  } else if (method == "mcmc" | method == "mcmc_truncated") {
    mu.mean <- apply(output_app$mus_new[1,,INIT_START:dim(output_app$mus_new)[3]], "median",MARGIN= 1)
    alpha.mean <- apply(output_app$alphas_new[,,INIT_START:dim(output_app$alphas_new)[3]], "median",MARGIN= c(1,2))
    beta.mean <- apply(output_app$betas_new[,,INIT_START:dim(output_app$betas_new)[3]], "median",MARGIN= c(1,2))
  } else {
    mu.mean <- output_app$eta_mu1[1, ] / output_app$eta_mu2[1, ]
    alpha.mean <- matrix(unlist(output_app$eta_alphas1), K, K, byrow = TRUE)/matrix(unlist(output_app$eta_alphas2), K, K, byrow = TRUE)
    beta.mean <- matrix(unlist(output_app$eta_betas1), K, K, byrow = TRUE)/matrix(unlist(output_app$eta_betas2), K, K, byrow = TRUE)
  }
  return(list(mu.mean = mu.mean, alpha.mean = alpha.mean, beta.mean = beta.mean))
}

# get parameter uncertainty estimate out from the raw file, in matrix forms
preprocess_uncertainties <- function(output_app, method, alpha) {
  K <- 11
  # SGLD: get rid of burn in, get median for the rest
  if (method == "sgld") {
    mu.upp <- apply(simplify2array(output_app$mu_sgld[c(INIT_START_SGLD:length(output_app$alpha_sgld))]),
                    "quantile", probs = 1-alpha/2, MARGIN= c(1,2))
    mu.low <- apply(simplify2array(output_app$mu_sgld[c(INIT_START_SGLD:length(output_app$alpha_sgld))]),
                    "quantile", probs = alpha/2, MARGIN= c(1,2))
    mu.cred <- mu.upp - mu.low
    alpha.upp <- apply(simplify2array(output_app$alpha_sgld[c(INIT_START_SGLD:length(output_app$alpha_sgld))]),
                       "quantile", probs = 1-alpha/2, MARGIN= c(1,2))
    alpha.low <- apply(simplify2array(output_app$alpha_sgld[c(INIT_START_SGLD:length(output_app$alpha_sgld))]),
                       "quantile", probs = alpha/2, MARGIN= c(1,2))
    alpha.cred <- alpha.upp - alpha.low
    beta.upp <- apply(simplify2array(output_app$beta_sgld[c(INIT_START_SGLD:length(output_app$alpha_sgld))]),
                      "quantile", probs = 1-alpha/2, MARGIN= c(1,2))
    beta.low <- apply(simplify2array(output_app$beta_sgld[c(INIT_START_SGLD:length(output_app$alpha_sgld))]),
                      "quantile", probs = alpha/2, MARGIN= c(1,2))
    beta.cred <- beta.upp - beta.low
  } else if (method == "sem" | method == "semtr") {
    mu.mean <- NA #output_app$mu_em
    alpha.mean <- NA #output_app$alpha_em
    beta.mean <- NA #output_app$beta_em
    mu.cred <- rep(NA, K)
    alpha.cred <- matrix(NA, K, K)
    beta.cred <- matrix(NA, K, K)
  } else if (method == "mcmc" | method == "mcmc_truncated") {
    mu.low <- apply(output_app$mus_new[1,,INIT_START:dim(output_app$mus_new)[3]],
                    "quantile", probs = alpha/2, MARGIN= 1)
    mu.upp <- apply(output_app$mus_new[1,,INIT_START:dim(output_app$mus_new)[3]],
                    "quantile", probs = 1-alpha/2, MARGIN= 1)
    mu.cred <- mu.upp - mu.low
    alpha.low <- apply(output_app$alphas_new[,,INIT_START:dim(output_app$alphas_new)[3]],
                       "quantile", probs = alpha/2, MARGIN= c(1,2))
    alpha.upp <- apply(output_app$alphas_new[,,INIT_START:dim(output_app$alphas_new)[3]],
                       "quantile", probs = 1-alpha/2, MARGIN= c(1,2))
    alpha.cred <- alpha.upp - alpha.low
    beta.low <- apply(output_app$betas_new[,,INIT_START:dim(output_app$betas_new)[3]],
                      "quantile", probs = alpha/2, MARGIN= c(1,2))
    beta.upp <- apply(output_app$betas_new[,,INIT_START:dim(output_app$betas_new)[3]],
                      "quantile", probs = 1-alpha/2, MARGIN= c(1,2))
    beta.cred <- beta.upp - beta.low
  } else {
    mu.low <- qgamma(alpha/2, shape = output_app$eta_mu1[1, ], rate = output_app$eta_mu2[1, ])
    mu.upp <- qgamma(1-alpha/2, shape = output_app$eta_mu1[1, ], rate = output_app$eta_mu2[1, ])
    mu.cred <- mu.upp - mu.low
    alpha.low <- matrix(qgamma(alpha/2, shape = unlist(output_app$eta_alphas1), rate = unlist(output_app$eta_alphas2)), K, K, byrow = TRUE)
    alpha.upp <- matrix(qgamma(1-alpha/2, shape = unlist(output_app$eta_alphas1), rate = unlist(output_app$eta_alphas2)), K, K, byrow = TRUE)
    alpha.cred <- alpha.upp - alpha.low
    beta.low <- matrix(qgamma(alpha/2, shape = unlist(output_app$eta_betas1), rate = unlist(output_app$eta_betas2)), K, K, byrow = TRUE)
    beta.upp <- matrix(qgamma(1-alpha/2, shape = unlist(output_app$eta_betas1), rate = unlist(output_app$eta_betas2)), K, K, byrow = TRUE)
    beta.cred <- beta.upp - beta.low
  }
  return(list(mu.cred = mu.cred, alpha.cred = alpha.cred, beta.cred = beta.cred))
}

# generate file name
generate_file_name <- function(method_id, seed_id) {
  if (method_id <= 5) {
    filename <- paste0(RESULTS_PATH, NAME_PREFIX,
                       METHOD_NAMES[method_id], "_seed_", seed_id, "_dataset_1",
                       NAME_SUFFICES[method_id])
  } else if (method_id == 6) {
    filename <- paste0(RESULTS_PATH, "mcmc_", seed_id, ".RData")
  } else {
    filename <- paste0(RESULTS_PATH, "mcmc_", seed_id, "_truncated.RData")
  }
  return(filename)
}

slice <- function(x,pos) {
  K <- 11
  res <- rep(0, K^2)
  for (i in 1:K) {
    for (j in 1:K) {
      res[(i-1)*K + j] <- x[[i]][[j]][pos]
    }
  }
  return(matrix(res, nrow = K, byrow = TRUE))
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

compensator <- function(mu,alpha,beta,dim,out) {
  y_dim <- out[out$dim == dim, "timestamp"]
  n_dim <- length(y_dim)
  z_out <- rep(0, n_dim)
  for (i in 1:n_dim) {
    z_out[i] <- mu[dim]*y_dim[i]
    if (i > 1) {
      out_filtered <- out[out$timestamp < y_dim[i], ]
      n_dim_filtered <- dim(out_filtered)[1]
      alpha_vec <- alpha[cbind(out_filtered$dim,rep(dim,n_dim_filtered))]
      beta_vec <- beta[cbind(out_filtered$dim,rep(dim,n_dim_filtered))]
      z_out[i] <- z_out[i] + sum(alpha_vec*(1 - exp(-beta_vec*(y_dim[i] - out_filtered$timestamp))))
    }
  }
  return(z_out)
}

compensator_mcmc <- function(results,dim,out,method,iter_id) {
  tempfun <- function(mu, alpha, beta) {
    z_out <- rep(0,n_dim)
    for (i in 1:n_dim) {
      z_out[i] <- mu[dim]*y_dim[i]
      if (i > 1) {
        out_filtered <- out[out$timestamp < y_dim[i], ]
        n_dim_filtered <- dim(out_filtered)[1]
        alpha_vec <- alpha[cbind(out_filtered$dim,rep(dim,n_dim_filtered))]
        beta_vec <- beta[cbind(out_filtered$dim,rep(dim,n_dim_filtered))]
        z_out[i] <- z_out[i] + sum(alpha_vec*(1 - exp(-beta_vec*(y_dim[i] - out_filtered$timestamp))))
      }
    }
    return(z_out)
  }
  
  range(out$timestamp)
  
  out_test <- out
  
  y_dim <- out_test[out_test$dim == dim, "timestamp"]
  n_dim <- length(y_dim)
  iters <- dim(results$alphas_new)[3]
  mu <- lapply(seq(dim(results$mus_new)[3]), function(x) results$mus_new[ , , x])
  alpha <- lapply(seq(dim(results$alphas_new)[3]), function(x) results$alphas_new[ , , x])
  beta <- lapply(seq(dim(results$betas_new)[3]), function(x) results$betas_new[ , , x])
  z_out_list <- list()
  
  id_start <- (iter_id-1)*40 + 1
  id_end <- iter_id*40
  
  for (iter in id_start:id_end) {
    mu <- results$mus_new[,,iter]
    alpha <- results$alphas_new[,,iter]
    beta <- results$betas_new[,,iter]
    z_out_list[[iter - id_start + 1]] <- tempfun(mu, alpha, beta)
  }
  return(z_out_list)
}


compensator_ld <- function(results,dim,out,method,iter_id) {
  tempfun <- function(mu, alpha, beta) {
    z_out <- rep(0,n_dim)
    for (i in 1:n_dim) {
      z_out[i] <- mu[dim]*y_dim[i]
      if (i > 1) {
        out_filtered <- out[out$timestamp < y_dim[i], ]
        n_dim_filtered <- dim(out_filtered)[1]
        alpha_vec <- alpha[cbind(out_filtered$dim,rep(dim,n_dim_filtered))]
        beta_vec <- beta[cbind(out_filtered$dim,rep(dim,n_dim_filtered))]
        z_out[i] <- z_out[i] + sum(alpha_vec*(1 - exp(-beta_vec*(y_dim[i] - out_filtered$timestamp))))
      }
    }
    return(z_out)
  }
  
  out_test <- out
  
  y_dim <- out_test[out_test$dim == dim, "timestamp"]
  n_dim <- length(y_dim)
  z_out_list <- list()
  
  id_start <- (iter_id-1)*40 + 1
  id_end <- min(iter_id*40, length(results$alpha_sgld))
  
  for (iter in id_start:id_end) {
    mu <- results$mu_sgld[[iter]]
    alpha <- results$alpha_sgld[[iter]]
    beta <- results$beta_sgld[[iter]]
    z_out_list[[iter - id_start + 1]] <- tempfun(mu, alpha, beta)
  }
  return(z_out_list)
}


repackage <- function(x) {
  K <- 11
  pos <- length(x[[1]][[1]])
  res <- rep(0, K^2)
  for (i in 1:K) {
    for (j in 1:K) {
      res[(i-1)*K + j] <- x[[j]][[i]][pos]
    }
  }
  return(res)
}

compensator_vi <- function(output,dim,out) {
  tempfun <- function(mu, alpha, beta) {
    z_out <- rep(0,n_dim)
    for (i in 1:n_dim) {
      z_out[i] <- mu[dim]*y_dim[i]
      if (i > 1) {
        out_filtered <- out[out$timestamp < y_dim[i], ]
        n_dim_filtered <- dim(out_filtered)[1]
        alpha_vec <- alpha[cbind(out_filtered$dim,rep(dim,n_dim_filtered))]
        beta_vec <- beta[cbind(out_filtered$dim,rep(dim,n_dim_filtered))]
        z_out[i] <- z_out[i] + sum(alpha_vec*(1 - exp(-beta_vec*(y_dim[i] - out_filtered$timestamp))))
      }
    }
    return(z_out)
  }
  
  y_dim <- out[out$dim == dim, "timestamp"]
  n_dim <- length(y_dim)
  B <- 1
  z_out_list <- list()
  for (b in 1:B) {
    mu <- rgamma(K, shape = output$eta_mu1[1,], rate = output$eta_mu2[1,])
    alpha <- matrix(rgamma(K^2, shape = repackage(output$eta_alphas1),
                           rate = repackage(output$eta_alphas2)), nrow = K)
    beta <- matrix(rgamma(K^2, shape = repackage(output$eta_betas1),
                          rate = repackage(output$eta_betas2)), nrow = K)
    z_out_list[[b]] <- tempfun(mu, alpha, beta)
  }
  return(z_out_list)
}

# load in data
taskid <- as.numeric(commandArgs(trailingOnly = TRUE)[1])
K <- 11
pars <- expand.grid(dim = 1:K)
mydim <- pars$dim[taskid]

set.seed(123)
METHOD_NAMES <- c("sem", "semtr", "sgd", "sgdtr", "sgld", "mcmc", "mcmc_truncated")
METHOD_NAMES_FORPLOT <- c("SEM", "SEM-c", "SGVI", "SGVI-c", "SGLD", "MCMC", "MCMC-c")
NAME_PREFIX <- c("mhpvi")
NAME_SUFFICES <- c("_time30.RData", "_time30.RData", "newinit_time30.RData", "newinit_time30.RData", "_newinit_time30.RData")

# graphics
PLOT_HEIGHT <- 10
PLOT_WIDTH <- 22


setwd("/mnt/beegfs/homes/jiang14/mhp/goodnessoffit")
# parameters for the dataset
K <- 11 # number of dimensions
return_names <- c("S5COND", "S5CONS", "S5ENRS", "S5FINL",
                  "S5HLTH", "S5INDU", "S5MATR", "S5INFT", "S5TELS", "S5UTIL", "S5RLST")
return_names_p <- sub('..', '', return_names)

load("data_app.RData")

plist <- list()
load("mhpvisem_seed_12_dataset_1_time30.RData")
res <- compensator(alpha = output$alpha_em, 
                   beta = output$beta_em, mu = output$mu_em, 
                   dim = mydim, out = out)
fileoutname <- paste0("sgem_dim_", mydim, "_", 1, ".RData")
save(res, file = file.path("output/sgem", fileoutname))

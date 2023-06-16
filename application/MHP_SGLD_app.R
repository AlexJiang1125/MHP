## calculate likelihood for model output_apps
rm(list = ls())
library(dplyr)
library(ggplot2)
## data re-formatting

MargLik_grad <- function(mu,alpha,beta,outdat) {
  ## denominator: a vector of n elements
  n <- dim(outdat)[1]
  y <- outdat$timestamp
  Tdat <- T*split_ratio
  denoms <- rep(0, n)
  denoms[1] <- mu[outdat$dim[1]]
  Marg_partial_alpha <- matrix(0,K,K)
  Marg_partial_beta <- matrix(0,K,K)
  Marg_partial_mu <- rep(0,K)
  for (i in 2:n) {
    alpha_vec <- alpha[cbind(outdat$dim[1:(i-1)],rep(outdat$dim[i],(i-1)))]
    beta_vec <- beta[cbind(outdat$dim[1:(i-1)],rep(outdat$dim[i],(i-1)))]
    denoms[i] <- sum(alpha_vec*beta_vec*exp(-beta_vec*(y[i] - y[1:(i-1)]))) + mu[outdat$dim[i]]
  }
  for (i in 1:K) {
    t_compensator <- outdat[outdat$dim == i, "timestamp"]
    for (j in 1:K) {
      temp_alpha <- rep(0,n)
      temp_beta <- rep(0,n)
      for (k in 2:n) {
        if (outdat$dim[k] == j) {
          ttemp <- outdat[outdat$dim == i,"timestamp"]
          ttemp_t <- ttemp[ttemp < outdat$timestamp[k]]
          temp_alpha[k] <- sum(exp(-beta[i,j]*(outdat$timestamp[k] - ttemp_t))*beta[i,j])
          temp_beta[k] <- sum(alpha[i,j]*(1-beta[i,j]*(outdat$timestamp[k] - ttemp_t))*exp(-beta[i,j]*(outdat$timestamp[k] - ttemp_t)))
        }
      }
      Marg_partial_alpha[i,j] <- sum(temp_alpha/denoms) - sum(outdat$dim == i) + sum(exp(-beta[i,j]*(Tdat - t_compensator)))
      Marg_partial_beta[i,j] <- sum(temp_beta/denoms) - sum(alpha[i,j]*(Tdat - t_compensator)*exp(-beta[i,j]*(Tdat - t_compensator)))
    }
  }
  for (i in 1:K) {
    Marg_partial_mu[i] <- sum((outdat$dim == i)/denoms) - Tdat
  }
  res_list = list(
    alpha_grad = Marg_partial_alpha,
    beta_grad = Marg_partial_beta,
    mu_grad = Marg_partial_mu
  )
  return(res_list)
}


MargLik <- function(mu,alpha,beta,out) {
  n <- dim(out)[1]
  y <- out$timestamp
  res <- log(mu[out$dim[1]])
  Tdat <- T*split_ratio
  for (i in 2:n) {
    alpha_vec <- alpha[cbind(out$dim[1:(i-1)],rep(out$dim[i],(i-1)))]
    beta_vec <- beta[cbind(out$dim[1:(i-1)],rep(out$dim[i],(i-1)))]
    res <- res + log(sum(alpha_vec*beta_vec*exp(-beta_vec*(y[i] - y[1:(i-1)]))) + mu[out$dim[i]])
  }
  ns <- as.vector(table(out$dim))
  res <- res - sum(mu)*Tdat - sum(apply(alpha, MARGIN = 1, FUN = sum)*ns)
  for (i in 1:K) {
    for (j in 1:K) {
      res <- res + alpha[i,j]*sum(exp(-beta[i,j]*(Tdat - out$timestamp[out$dim == i])))
    }
  }
  return(res)
}


All_grad <- function(mu, alpha, beta, out) {
  temp_grad <- MargLik_grad(mu,alpha,beta,out)
  mu_grad <- -temp_grad$mu_grad + b_mu - (a_mu - 1)/mu
  alpha_grad <- -temp_grad$alpha_grad + 1*f_alpha - (e_alpha - 1)/alpha
  #alpha_grad <- -temp_grad$alpha_grad + (pis_true + 100*(1-pis_true))*f_alpha - (e_alpha - 1)/alpha
  beta_grad <- -temp_grad$beta_grad + s_beta - (r_beta - 1)/beta
  res_list = list(
    alpha_grad = alpha_grad,
    beta_grad = beta_grad,
    mu_grad = mu_grad
  )
  return(res_list)
}

# gradient for parameters on a log scale

All_grad_log <- function(mu, alpha, beta, out) {
  temp_grad <- MargLik_grad(mu,alpha,beta,out)
  mu_grad <- -mu*temp_grad$mu_grad + b_mu*mu - (a_mu - 1) - 1
  alpha_grad <- -temp_grad$alpha_grad*alpha + alpha*f_alpha - (e_alpha - 1) - 1
  #alpha_grad <- -temp_grad$alpha_grad + (pis_true + 100*(1-pis_true))*f_alpha - (e_alpha - 1)/alpha
  beta_grad <- -temp_grad$beta_grad*beta + s_beta*beta - (r_beta - 1) - 1
  res_list = list(
    alpha_grad = alpha_grad,
    beta_grad = beta_grad,
    mu_grad = mu_grad
  )
  return(res_list)
}


# scaled by the factor
All_sgrad_log <- function(mu, alpha, beta, outdat, fct) {
  temp_grad <- MargLik_grad(mu,alpha,beta,outdat)
  mu_grad <- -mu*temp_grad$mu_grad*fct + b_mu*mu - (a_mu - 1) #- 1
  alpha_grad <- -temp_grad$alpha_grad*alpha*fct + alpha*f_alpha - (e_alpha - 1) #- 1
  #alpha_grad <- -temp_grad$alpha_grad + (pis_true + 100*(1-pis_true))*f_alpha - (e_alpha - 1)/alpha
  beta_grad <- -temp_grad$beta_grad*beta*fct + s_beta*beta - (r_beta - 1) #- 1
  res_list = list(
    alpha_grad = alpha_grad,
    beta_grad = beta_grad,
    mu_grad = mu_grad
  )
  return(res_list)
}

#load(paste0("/Users/alexziyujiang/Documents/data/data/data_",1,".RData"))
#load(paste0("/Users/alexziyujiang/Documents/data/EM_inits/inits_",1,".RData"))
type <- "cluster"
modes <- c("short","medium","long")
mode <- modes[2]
K <- 11

if (type == "cluster") {
  taskid <- as.numeric(commandArgs(trailingOnly = TRUE)[1])
} else {
  taskid <- 1
}

for (dataset_seed in 1: 1) {
  if (type == "cluster") {
    if (mode == "medium") {
      load(
        paste0(
          "/gscratch/home/jiang14/mhp/data/sim001/data/data_app/",
          "data_app",
          ".RData"
        )
      )
    } else if (mode == "short") {
      load(
        paste0(
          "/gscratch/home/jiang14/mhp/data/sim001/data/data_short/data_",
          dataset_seed,
          ".RData"
        )
      )
    } else {
      load("/gscratch/home/jiang14/mhp/data/sim001/data/data_app/data_app.RData")
    }
  } else {
    load("/Users/alexziyujiang/Documents/GitHub/HawkesSP500/dataset/data_app.RData")
  }
  files <- c("100",
             "090",
             "080",
             "070",
             "060",
             "050",
             "040",
             "030",
             "020",
             "010",
             "005",
             "001")
  split_ratios <- as.numeric(files) / 100
  split_id <- 12
  split_ratio <- split_ratios[split_id]
  NBATCH <- 1 / split_ratio # number of batches
  T <- 1000
  stepsize <- min(0.1/(T*split_ratios[split_id]), 0.002)
  mus <- list()
  alphas <- list()
  betas <- list()
  # sgd pars
  sgd_par1 <- 1
  sgd_par2 <- 0.51
  print(paste0("learning parameters: ", sgd_par1, ",", sgd_par2))
  out_all <- out
  N <- length(out$timestamp)
  Nall <- N
  y <- (out$timestamp)
  N_each <- N
  print(taskid)
  set.seed(taskid)
  if (type == "cluster") {
    load(
      paste0(
        "/gscratch/home/jiang14/mhp/data/sim001/SGLD_app_inits/inits_",
        taskid,
        ".RData"
      )
    )
  } else {
    #load("/Users/alexziyujiang/Dropbox/data/mhp_app/output_new/mhpvisem_seed_5_dataset_1_time30.RData")
    load("/Users/alexziyujiang/Dropbox/data/mhp_app/EM_app_inits/inits_1.RData")
  }
  alpha <- inits[[1]]
  beta <- inits[[2]]
  mu <- inits[[3]]
  rm(output)
  # hyperpars
  a_mu <- 2
  b_mu <- 4
  e_alpha <- 2
  f_alpha <- 4
  r_beta <- 2
  s_beta <- 2 / 4
  count <- 1
  offset <- 1
  # initializing time
  time <- 0
  time_old <- time
  #mode <- "log"
  clock_start <- Sys.time()
  count_warmup <- 0
  if (mode == "long") {
    T <- 2000
  } else if (mode == "short") {
    T <- 500
  } else {
    T <- 1000
  }
  while (time < 10) {
    if (count <= offset) {
      rho <- stepsize
    } else {
      rho <-
        stepsize * (count - offset - 1 + sgd_par1) ^ (-sgd_par2) # step size parameter
    }
    # sample sub dataset
    time_start <- runif(1, 0, T - T * split_ratio)
    time_end <- time_start + T * split_ratio
    out_temp <-
      out[out$timestamp < time_end &
            out$timestamp >= time_start ,]
    if (dim(out_temp)[1] < 10) {
      while (dim(out_temp)[1] < 10) {
        time_start <- runif(1,0,T-T*split_ratio)
        time_end <- time_start + T*split_ratio
        out_temp <- out[out$timestamp < time_end & out$timestamp >= time_start ,]
      }
    }
    out_temp$timestamp <- out_temp$timestamp - time_start
    #print(mean(alpha))
    newgrad <-
      All_sgrad_log(mu, alpha, beta, out_temp, 1 / split_ratio)
    beta <- exp(log(beta) - newgrad$beta_grad * rho / 2 + sqrt(rho)*matrix(rnorm(K^2, mean = 0, sd = 1), K, K))
    newgrad <-
      All_sgrad_log(mu, alpha, beta, out_temp, 1 / split_ratio)
    mu <- exp(log(mu) - newgrad$mu_grad * rho / 2 + sqrt(rho)*rnorm(K, mean = 0, sd = 1))
    newgrad <-
      All_sgrad_log(mu, alpha, beta, out_temp, 1 / split_ratio)
    alpha <- exp(log(alpha) - newgrad$alpha_grad * rho / 2 + sqrt(rho)*matrix(rnorm(K^2, mean = 0, sd = 1), K, K))
    
    mus[[count]] <- mu
    alphas[[count]] <- alpha
    betas[[count]] <- beta
    print(mean(alpha))
    count <- count + 1
    time_old <- time
    time <- difftime(Sys.time(), clock_start, units = 'mins')
    if (time_old < 0.5 && time >= 0.5) {
      count_warmup <- count
    }
  }
  #temp = simplify2array(mus[c(2:count)])
  #if (length(dim(temp)) == 2) {
  #   mu.mean <- apply(temp, "median",MARGIN= 1)
  #} else {
  #   mu.mean <- apply(temp, "median",MARGIN= 1)
  #}
  #alpha.mean <- apply(simplify2array(alphas[c(2:count)]), "median",MARGIN= c(1,2))
  #beta.mean <- apply(simplify2array(betas[c(2:count)]), "median",MARGIN= c(1,2))
  
  if (type == "cluster") {
    output_app <- list(alpha_sgld = alphas,
                       beta_sgld = betas,
                       mu_sgld = mus,
                       warmup = count_warmup)
    if (mode == "long") {
      save(
        output_app,
        file = paste0(
          "/gscratch/csde/jiang14/output_app_long/s",
          files[split_id],
          "/mhpvisgld_seed_",
          taskid,
          "_dataset_" ,
          dataset_seed,
          "_newinit_time10.RData"
        ))
    } else if (mode == "short") {
      save(
        output_app,
        file = paste0(
          "/gscratch/csde/jiang14/output_app_short/s",
          files[split_id],
          "/mhpvisgld_seed_",
          taskid,
          "_dataset_" ,
          dataset_seed,
          "_newinit_time30.RData"
        ))
    } else {
      save(
        output_app,
        file = paste0(
          "/gscratch/csde/jiang14/output_app/s",
          files[split_id],
          "/mhpvisgld_seed_",
          taskid,
          "_dataset_" ,
          dataset_seed,
          "_newinit_time10.RData"
        ))
    }
  }
}

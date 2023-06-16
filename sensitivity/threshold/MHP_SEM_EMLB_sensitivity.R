## variational algorithm
## updating branching structures
# new version: does not store old parameters
# alpha: KxK kernels, B reps
rm(list = ls())
library(dplyr)
library(ggplot2)
library(tidyr)
modes <- c("short","medium","long")
mode <- modes[2]
type <- "cluster"
if (type == "local") {
  setwd("/Users/alexziyujiang/Dropbox/Alex_RP/scripts")
  source("MHP.R")
  taskid <- 1
} else {
  #source("/gscratch/home/jiang14/mhp/R/sim001/MHP.R")
  taskid <- as.numeric(commandArgs(trailingOnly=TRUE)[1])
}
a_mu <- 2
b_mu <- 4
e_alpha <- 2
f_alpha <- 4
r_beta <- 2
s_beta <- 2/4
type <- "cluster"
dat_tmp <- ceiling(taskid/16)
dat_start <- (dat_tmp-1)*10 + 1
dat_end <- dat_tmp*10
print(paste0(dat_start, " ", dat_end))

taskid <- (taskid - 1) %% 16 + 1
print(taskid)

for (dataset_seed in dat_tmp: dat_tmp) {
  if (type == "cluster") {
    load(
      paste0(
        "/gscratch/home/jiang14/mhp/data/sim001/data/data/data_",
        dataset_seed,
        ".RData"
      )
    )
  } else {
    taskid <- 1
    load(
      paste0(
        "/Users/alexziyujiang/Documents/data/data/data_",
        dataset_seed,
        ".RData"
      )
    )
  }
  #load(paste0("/gscratch/home/jiang14/mhp/data/sim001/data/data_",dataset_seed,".RData"))
  print(taskid)
  set.seed(taskid)
  time_EMLB <- 0
  sgd_par1 <- 1
  sgd_par2 <- 0.5
  print(paste0("learning parameters: ",sgd_par1,",",sgd_par2))
  out_all <- out
  N <- length(out$timestamp)
  y <- (out$timestamp)
  N_each <- N
  ## split the dataset into
  files <- c("100","090","080","070","060","050",
             "040","030","020","010","005","001")
  split_ratios <- as.numeric(files)/100
  split_id <- 11
  split_ratio <- split_ratios[split_id]
  NBATCH <- 1/split_ratio # number of batches
  #pis_true <- matrix(c(1,1,0,0,1,1,1,1,1),3,3, byrow = TRUE)
  ### EM parameters
  type <- "cluster"
  if (type == "cluster") {
    load(
      paste0(
        "/gscratch/home/jiang14/mhp/data/sim001/MCMC_inits/inits_",
        taskid,
        ".RData"
      )
    )
  } else {
    load(
      paste0(
        "/Users/alexziyujiang/Documents/data/MCMC_inits/inits_",
        taskid,
        ".RData"
      )
    )
  }
  if (mode == "long") {
    T <- 2000
  } else if (mode == "short") {
    T <- 500
  } else {
    T <- 1000
  }
  alpha_em <- inits[[1]]
  beta_em <- inits[[2]]
  mu_em <- inits[[3]]
  alpha_em_new <- alpha_em
  beta_em_new <- beta_em
  mu_em_new <- mu_em
  s_mu1 <- rep(0,3)
  s_mu2 <- rep(0,3)
  s_alpha1 <- matrix(0,3,3)
  s_alpha2 <- matrix(0,3,3)
  s_beta1 <- matrix(0,3,3)
  s_beta2 <- matrix(0,3,3)
  Nall <- length(out_all$id)
  pos <- 1
  count <- 1
  offset <- 1
  maxiters <- 1000
  elbo <- 0
  time <- 0
  clock_start <- Sys.time()
  time <- 0
  flag <- 0
  stepsize <- 0.02
  thres_id <- 3
  thresholds <- c(0.17, 0.6, 0.75, 1.15)
  threshold_names <- c("t05", "t01", "t005", "t001")
  threshold <- thresholds[thres_id]
  # step 0
  time_start <- runif(1,0,T-T*split_ratio)
  time_end <- time_start + T*split_ratio
  out <- out_all[out_all$timestamp < time_end & out_all$timestamp >= time_start ,]
  if (dim(out)[1] < 10) {
    while (dim(out)[1] < 10) {
      time_start <- runif(1,0,T-T*split_ratio)
      time_end <- time_start + T*split_ratio
      out <- out_all[out_all$timestamp < time_end & out_all$timestamp >= time_start ,]
    }
  }
  out$timestamp <- out$timestamp - min(out$timestamp)
  y <- out$timestamp
  N <- length(out$timestamp)
  ## E step
  B <- matrix(0, N, N)
  B[1,1] <- 1
  for (ii in 2:N) {
    ptemp <- rep(0,ii)
    ptemp[ii] <- log(mu_em[out$dim[ii]])
    ptemp[1:(ii-1)] <- log(alpha_em[out$dim[1:(ii-1)],out$dim[ii]]) + log(beta_em[out$dim[1:(ii-1)],out$dim[ii]]) - beta_em[out$dim[1:(ii-1)], out$dim[ii]]*(y[ii]-y[1:(ii-1)])
    B[ii,1:ii] <- exp(ptemp - max(ptemp))/sum(exp(ptemp - max(ptemp)))
  }
  for (j in 1:K) {
    s_mu1[j] <- sum(diag(B)*(out$dim == j))/split_ratio
    s_mu2[j] <- T
  }
  for (iter_i in 1:K) {
    for (iter_j in 1:K) {
      ntrans <- 0
      if (iter_i == iter_j) {
        ntrans <- sum(B[out$dim == iter_j, out$dim == iter_i]) - sum(diag(B)*(out$dim == iter_i))
      } else {
        ntrans <- sum(B[out$dim == iter_j, out$dim == iter_i])
      }
      ts <- out$timestamp[out$dim == iter_i]
      s_alpha1[iter_i,iter_j] <- ntrans/split_ratio
      s_alpha2[iter_i,iter_j] <- (length(ts) - sum(exp(-beta_em[iter_i, iter_j]*(T * split_ratio - ts))))/split_ratio
      dists <- 0
      indices <- expand.grid(which(out$dim == iter_i), which(out$dim == iter_j))
      indices <- indices[indices[,1] < indices[,2],]
      dists <- sum(B[cbind(indices$Var2, indices$Var1)]*(y[indices$Var2] - y[indices$Var1]))
      s_beta1[iter_i,iter_j] <- ntrans/split_ratio
      s_beta2[iter_i,iter_j] <- dists/split_ratio
    }
  }
  for (j in 1:K) {
    mu_em[j] <- (s_mu1[j] + (a_mu - 1))/(s_mu2[j] + b_mu)
  }
  #mu_em <- mu_em*(1-rho) + mu_em_new*rho
  for (ii in 1:K) {
    for (jj in 1:K) {
      alpha_em[ii,jj] <- (s_alpha1[ii,jj] + (e_alpha - 1))/(f_alpha + s_alpha2[ii,jj])
      beta_em[ii,jj] <- (s_beta1[ii,jj] + (r_beta - 1))/(s_beta + s_beta2[ii,jj])
    }
  }
  while (time < 30) {
    count <- count + 1
    #
    if (count <= offset) {
      rho <- 1
    } else {
      rho <- 1*(count - offset - 1 + sgd_par1)^(-sgd_par2) # step size parameter
    }
    ## scan through batch
    # set statistics for global
    for (sb in 1:1) {
      time_start <- runif(1,0,T-T*split_ratio)
      time_end <- time_start + T*split_ratio
      out <- out_all[out_all$timestamp < time_end & out_all$timestamp >= time_start ,]
      if (dim(out)[1] < 10) {
        while (dim(out)[1] < 10) {
          time_start <- runif(1,0,T-T*split_ratio)
          time_end <- time_start + T*split_ratio
          out <- out_all[out_all$timestamp < time_end & out_all$timestamp >= time_start ,]
        }
      }
      out$timestamp <- out$timestamp - min(out$timestamp)
      y <- out$timestamp
      N <- length(out$timestamp)
      truncate_ratio <- 0.95
      out_truncated <- out[out$timestamp < T*split_ratio - threshold, ]
      out_tail <- out[out$timestamp >= T*split_ratio - threshold, ]
      n_truncated <- as.vector(table(out_truncated$dim))
      ys_truncated <- list()
      for (id_k in 1:K) {
        ys_truncated[[id_k]] <- out_tail[out_tail$dim == id_k,]
      }
      ## E step
      B <- matrix(0, N, N)
      B[1,1] <- 1
      for (ii in 2:N) {
        ptemp <- rep(0,ii)
        ptemp[ii] <- log(mu_em[out$dim[ii]])
        ptemp[1:(ii-1)] <- log(alpha_em[out$dim[1:(ii-1)],out$dim[ii]]) + log(beta_em[out$dim[1:(ii-1)],out$dim[ii]]) - beta_em[out$dim[1:(ii-1)], out$dim[ii]]*(y[ii]-y[1:(ii-1)])
        B[ii,1:ii] <- exp(ptemp - max(ptemp))/sum(exp(ptemp - max(ptemp)))
      }
      for (j in 1:K) {
        s_mu1[j] <- (1-rho)*s_mu1[j] + rho*sum(diag(B)*(out$dim == j))/split_ratio
        s_mu2[j] <- T
      }
      for (iter_i in 1:K) {
        for (iter_j in 1:K) {
          ntrans <- 0
          if (iter_i == iter_j) {
            ntrans <- sum(B[out$dim == iter_j, out$dim == iter_i]) - sum(diag(B)*(out$dim == iter_i))
          } else {
            ntrans <- sum(B[out$dim == iter_j, out$dim == iter_i])
          }
          ts <- out$timestamp[out$dim == iter_i]
          s_alpha1[iter_i,iter_j] <- s_alpha1[iter_i,iter_j]*(1-rho) + rho*ntrans/split_ratio
          s_alpha2[iter_i,iter_j] <-  s_alpha2[iter_i,iter_j]*(1-rho) + (length(ts) - sum(exp(-beta_em[iter_i, iter_j]*(T * split_ratio - ts))))/split_ratio*rho
          dists <- 0
          indices <- expand.grid(which(out$dim == iter_i), which(out$dim == iter_j))
          indices <- indices[indices[,1] < indices[,2],]
          dists <- sum(B[cbind(indices$Var2, indices$Var1)]*(y[indices$Var2] - y[indices$Var1])) + alpha_em[iter_i,iter_j]*sum(T*split_ratio - ys_truncated[[iter_i]]$timestamp)
          s_beta1[iter_i,iter_j] <- s_beta1[iter_i,iter_j]*(1-rho) + rho*ntrans/split_ratio
          s_beta2[iter_i,iter_j] <- s_beta2[iter_i,iter_j]*(1-rho) + rho*dists/split_ratio
        }
      }
      ## M step
      for (j in 1:K) {
        mu_em[j] <- (s_mu1[j] + (a_mu - 1))/(s_mu2[j] + b_mu)
      }
      #mu_em <- mu_em*(1-rho) + mu_em_new*rho
      for (ii in 1:K) {
        for (jj in 1:K) {
          alpha_em[ii,jj] <- (s_alpha1[ii,jj] + (e_alpha - 1))/(f_alpha + s_alpha2[ii,jj])
          beta_em[ii,jj] <- (s_beta1[ii,jj] + (r_beta - 1))/(s_beta + s_beta2[ii,jj])
        }
      }
    }
    # calculate the time used to calculate elbo
    time_old <- time
    time <- difftime(Sys.time(), clock_start, units = 'mins') - time_EMLB
  }
  # update ELBO
  elbo_start <- Sys.time()
  # out_EMLB <-
  #   out_all
  #out_EMLB$timestamp <- out_EMLB$timestamp - 0.45 * T
  # y_EMLB <- out_EMLB$timestamp
  # N_EMLB <- length(out_EMLB$timestamp)
  # Bnew <- matrix(0, N_EMLB, N_EMLB)
  # Bnew[1,1] <- 1
  # for (ii in 2:N_EMLB) {
  #   ptemp <- rep(0,ii)
  #   ptemp[ii] <- log(mu_em[out$dim[ii]])
  #   ptemp[1:(ii-1)] <- log(alpha_em[out$dim[1:(ii-1)],out$dim[ii]]) + log(beta_em[out$dim[1:(ii-1)],out$dim[ii]]) - beta_em[out$dim[1:(ii-1)], out$dim[ii]]*(y_EMLB[ii]-y_EMLB[1:(ii-1)])
  #   Bnew[ii,1:ii] <- exp(ptemp - max(ptemp))/sum(exp(ptemp - max(ptemp)))
  # }
  # Bnew_EMLB <- Bnew
  # emlb <- EMLB(alpha_em = alpha_em, beta_em = beta_em, mu_em = mu_em, Bnew_EMLB, out_EMLB, 1)
  time_EMLB <- time_EMLB + difftime(Sys.time(), elbo_start, units = 'mins')
  time <- difftime(Sys.time(), clock_start, units = 'mins') - time_EMLB
  output <- list(alpha_em = alpha_em,
                 beta_em = beta_em,
                 mu_em = mu_em, threshold = threshold)
  if (type == "cluster") {
    if (mode == "long") {
      save(
        output,
        file = paste0(
          "/gscratch/home/jiang14/mhp/output_long/s",
          files[split_id],
          "/mhpvisemtr_seed_",
          taskid,
          "_dataset_" ,
          dataset_seed,
          "_newinit_time30.RData"
        ))
    } else if (mode == "short") {
      save(
        output,
        file = paste0(
          "/gscratch/home/jiang14/mhp/output_short/sensitivity",
          files[split_id],
          "/mhpvisemtr_seed_",
          taskid,
          "_dataset_" ,
          dataset_seed,
          "_newinit_time30.RData"
        ))
    } else {
      save(
        output,
        file = paste0(
          "/gscratch/home/jiang14/mhp/output/sensitivity/",
          threshold_names[thres_id],
          "/mhpvisemtr_seed_",
          taskid,
          "_dataset_" ,
          dataset_seed,
          "_newinit_time30.RData"
        ))
    }
  } else {
    # save(
    #   output,
    #   file = paste0(
    #     "/Users/alexziyujiang/Documents/data/mhp/SGD/s",
    #     files[split_id],
    #     "/mhpsemtr_seed_",
    #     taskid,
    #     "_dataset_" ,
    #     dataset_seed,
    #     "_time05.RData"
    #   ))
  }
}

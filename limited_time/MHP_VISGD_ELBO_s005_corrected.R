## variational algorithm
## updating branching structures
# new version: does not store old parameters
# alpha: KxK kernels, B reps
rm(list = ls())
library(dplyr)
library(ggplot2)
library(tidyr)
type <- "cluster"
if (type == "local") {
  setwd("/Users/alexziyujiang/Dropbox/Alex_RP/scripts")
  source("MHP.R")
  taskid <- 1
} else {
  #source("/gscratch/home/jiang14/mhp/R/sim001/MHP.R")
  taskid_org <- as.numeric(commandArgs(trailingOnly = TRUE)[1])
  ids<-expand.grid(seed = 1:16, dataset = 1:50)
  dataset_seed <-ids$seed[taskid_org]
  taskid <- ids$dataset[taskid_org]
}
type <- "cluster"
slice <- function(x,pos) {
  return(matrix(c(x[[1]][[1]][pos], x[[1]][[2]][pos], x[[1]][[3]][pos],
                  x[[2]][[1]][pos], x[[2]][[2]][pos], x[[2]][[3]][pos],
                  x[[3]][[1]][pos], x[[3]][[2]][pos], x[[3]][[3]][pos]), nrow = 3, byrow = TRUE))
}
ELBO <-
  function(eta_mu1,
           eta_mu2,
           eta_alphas1,
           eta_alphas2,
           eta_betas1,
           eta_betas2,
           B_dat,
           out_dat,
           nb) {
    T <- 1000
    a_mu <- 1
    b_mu <- 2
    e_alpha <- 1
    f_alpha <- 2
    r_beta <- 1
    s_beta <- 1 / 4
    N <- length(out_dat[, 1])
    res <- 0
    #pos <- length(eta_betas1[[1]][[1]])
    for (i in 1:K) {
      for (j in 1:K) {
        ntrans <- 0
        if (i == j) {
          ntrans <-
            sum(B_dat[out_dat$dim == j, out_dat$dim == i]) - sum(diag(B_dat) * (out_dat$dim == i))
        } else {
          ntrans <- sum(B_dat[out_dat$dim == j, out_dat$dim == i])
        }
        res <-
          res + (digamma(eta_alphas1[[i]][[j]][pos]) - log(eta_alphas2[[i]][[j]][pos])) *
          (a_mu - eta_alphas1[[i]][[j]][pos] + ntrans)
        res <-
          res + (digamma(eta_betas1[[i]][[j]][pos]) - log(eta_betas2[[i]][[j]][pos])) *
          (e_alpha - eta_betas1[[i]][[j]][pos] + ntrans)
        res <-
          res - exp(log(eta_alphas1[[i]][[j]][pos]) - log(eta_alphas2[[i]][[j]][pos])) *
          (sum(out_dat$dim == i) + f_alpha - eta_alphas2[[i]][[j]][pos])
        res <-
          res + lgamma(eta_alphas1[[i]][[j]][pos]) - eta_alphas1[[i]][[j]][pos] *
          log(eta_alphas2[[i]][[j]][pos])
        res <-
          res + lgamma(eta_betas1[[i]][[j]][pos]) - eta_betas1[[i]][[j]][pos] * log(eta_betas2[[i]][[j]][pos])
        dists <- 0
        for (ki in 2:N) {
          for (kj in 1:(ki - 1)) {
            if (out_dat$dim[ki] == j & out_dat$dim[kj] == i) {
              dists <-
                dists + (out_dat$timestamp[ki] - out_dat$timestamp[kj]) * B_dat[ki, kj]
            }
          }
        }
        res <-
          res - exp(log(eta_betas1[[i]][[j]][pos]) - log(eta_betas2[[i]][[j]][pos])) *
          (dists + s_beta - eta_betas2[[i]][[j]][pos])
      }
    }
    for (i in 1:K) {
      res <-
        res + (digamma(eta_mu1[pos, i]) - log(eta_mu1[pos, i])) * (sum(diag(B_dat) *
                                                                         (out_dat$dim == i)) + a_mu - eta_mu1[pos, i])
      res <-
        res + lgamma(eta_mu1[pos, i]) - eta_mu1[pos, i] * log(eta_mu2[pos, i])
      res <-
        res - eta_mu1[pos, i] / eta_mu2[pos, i] * (T / 10 + b_mu - eta_mu2[pos, i])
    }
    res <- res - sum(B_dat[B_dat != 0] * log(B_dat[B_dat != 0]))
    return(res)
  }

type <- "cluster"
for (split_id in 7:11) {
for (dataset_seed in dataset_seed:dataset_seed) {
  if (type == "cluster") {
    load(
      paste0(
        "/gscratch/home/jiang14/mhp/data/sim001/data/data_long/data_",
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
  print(taskid)
  set.seed(taskid)
  time_ELBO <- 0
  sgd_par1 <- 1
  sgd_par2 <- 0.5
  print(paste0("learning parameters: ", sgd_par1, ",", sgd_par2))
  out_all <- out
  N <- length(out$timestamp)
  y <- (out$timestamp)
  N_each <- N
  ## split the dataset into
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
  #split_id <- 11
  split_ratio <- split_ratios[split_id]
  NBATCH <- 1 / split_ratio # number of batches
  #pis_true <- matrix(c(1, 1, 0, 0, 1, 1, 1, 1, 1), 3, 3, byrow = TRUE)
  if (type == "cluster") {
    load(
      paste0(
        "/gscratch/home/jiang14/mhp/data/sim001/VI_inits/inits_",
        taskid,
        ".RData"
      )
    )
  } else {
    load(
      paste0(
        "/Users/alexziyujiang/Documents/data/VI_inits/inits_",
        taskid,
        ".RData"
      )
    )
  }
  # prior pars
  a_mu <- 2
  b_mu <- 4
  e_alpha <- 2
  f_alpha <- 4
  r_beta <- 2
  s_beta <- 2/4
  T <- 2000
  eta_alphas1 <- inits$eta_alphas1
  eta_alphas2 <- inits$eta_alphas2
  eta_betas1 <- inits$eta_betas1
  eta_betas2 <- inits$eta_betas2
  eta_mu1 <- inits$eta_mu1
  eta_mu2 <- inits$eta_mu2
  Nall <- length(out_all$id)
  count <- 1
  offset <- 1
  maxiters <- 1000
  elbo <- 0
  time <- 0
  times <- 0
  clock_start <- Sys.time()
  time <- 0
  flag <- 0
  stepsize <- 0.02
  while (time < 5) {
    eta_mu1 <- rbind(eta_mu1, rep(0, K))
    eta_mu2 <- rbind(eta_mu2, rep(0, K))
    count <- count + 1
    #
    if (count <= offset) {
      rho <- 1
    } else {
      rho <-
        1 * (count - offset - 1 + sgd_par1) ^ (-sgd_par2) # step size parameter
    }
    ## scan through batch
    # set statistics for global
    for (sb in 1:1) {
      new_mu1 <- new_mu2 <- matrix(0, nrow = 1, ncol = K)
      new_alphas1 <- new_alphas2 <- list()
      new_betas1 <- new_betas2 <- list()
      for (i in 1:K) {
        temp_all <- list()
        for (j in 1:K) {
          temp_all[[j]] <- 0
        }
        new_alphas1[[i]]  <-
          new_alphas2[[i]] <- new_betas1[[i]] <- new_betas2[[i]] <- temp_all
      }
      # preprocessing data
      # randomly sample
      time_start <- runif(1, 0, T - T * split_ratio)
      time_end <- time_start + T * split_ratio
      out <-
        out_all[out_all$timestamp < time_end &
                  out_all$timestamp >= time_start , ]
      out$timestamp <- out$timestamp - time_start
      y <- out$timestamp
      N <- length(out$timestamp)
      out_truncated <- out[out$timestamp < T*split_ratio - 0.25, ]
      out_tail <- out[out$timestamp >= T*split_ratio - 0.25, ]
      n_truncated <- as.vector(table(out_truncated$dim))
      ys_truncated <- list()
      for (id_k in 1:K) {
        ys_truncated[[id_k]] <- out_tail[out_tail$dim == id_k,]
      }
      
      ## update B
      B <- matrix(0, N, N)
      B[1, 1] <- 1
      sliced_eta_alphas1 <- slice(eta_alphas1,1)
      sliced_eta_alphas2 <- slice(eta_alphas2,1)
      sliced_eta_betas1 <- slice(eta_betas1,1)
      sliced_eta_betas2 <- slice(eta_betas2,1)
      for (ii in 2:N) {
        ptemp <- rep(0, ii)
        ptemp[ii] <-
          digamma(eta_mu1[1, out$dim[ii]]) - log(eta_mu2[1, out$dim[ii]])
        # for (jj in 1:(ii - 1)) {
        #   ptemp[jj] <-
        #     digamma(eta_alphas1[[out$dim[jj]]][[out$dim[ii]]][1]) -
        #     log(eta_alphas2[[out$dim[jj]]][[out$dim[ii]]][1]) +
        #     digamma(eta_betas1[[out$dim[jj]]][[out$dim[ii]]][1]) -
        #     log(eta_betas2[[out$dim[jj]]][[out$dim[ii]]][1]) -
        #     eta_betas1[[out$dim[jj]]][[out$dim[ii]]][1] / eta_betas2[[out$dim[jj]]][[out$dim[ii]]][1] * (y[ii] - y[jj])
        # }
        ptemp[1:(ii-1)] <- digamma(sliced_eta_alphas1[cbind(out$dim[1:(ii-1)],rep(out$dim[ii],(ii-1)))])-
          log(sliced_eta_alphas2[cbind(out$dim[1:(ii-1)],rep(out$dim[ii],(ii-1)))]) +
          digamma(sliced_eta_betas1[cbind(out$dim[1:(ii-1)],rep(out$dim[ii],(ii-1)))]) -
          log(sliced_eta_betas2[cbind(out$dim[1:(ii-1)],rep(out$dim[ii],(ii-1)))]) -
          sliced_eta_betas1[cbind(out$dim[1:(ii-1)],rep(out$dim[ii],(ii-1)))]/sliced_eta_betas2[cbind(out$dim[1:(ii-1)],rep(out$dim[ii],(ii-1)))]*(y[ii]-y[1:(ii-1)])
        B[ii, 1:ii] <-
          exp(ptemp - max(ptemp)) / sum(exp(ptemp - max(ptemp)))
      }
      for (j in 1:K) {
        # update sufficient statistics
        new_mu1[j] <- (sum(diag(B) * (out$dim == j)))
        new_mu2[j] <- (T * split_ratio)
        # update parameters
        eta_mu1[1, j] <-
          eta_mu1[1, j] * (1 - rho) + rho * (a_mu + new_mu1[j] / split_ratio)
        eta_mu2[1, j] <-
          eta_mu2[1, j] * (1 - rho) + rho * (b_mu + new_mu2[j] / split_ratio)
      }
      for (ii in 1:K) {
        for (jj in 1:K) {
          ntrans <- 0
          if (ii == jj) {
            ntrans <-
              sum(B[out$dim == jj, out$dim == ii]) - sum(diag(B) * (out$dim == ii))
          } else {
            ntrans <- sum(B[out$dim == jj, out$dim == ii])
          }
          ts <- out$timestamp[out$dim == ii]
          new_alphas1[[ii]][[jj]][1] <- ntrans
          new_betas1[[ii]][[jj]][1] <- ntrans
          eta_alphas1[[ii]][[jj]][1] <-
            eta_alphas1[[ii]][[jj]][1] * (1 - rho) + (e_alpha + new_alphas1[[ii]][[jj]][1] /
                                                              split_ratio) * rho
          eta_betas1[[ii]][[jj]][1] <-
            eta_betas1[[ii]][[jj]][1] * (1 - rho) + (r_beta + new_betas1[[ii]][[jj]][1] /
                                                             split_ratio) * rho
          new_alphas2[[ii]][[jj]][1] <-
            sum(1 - exp(-log(
              1 - (ts - T * split_ratio) / (eta_betas2[[ii]][[jj]][1])
            ) * eta_betas1[[ii]][[jj]][1]))
          eta_alphas2[[ii]][[jj]][1] <-
            eta_alphas2[[ii]][[jj]][1] * (1 - rho) + (f_alpha + new_alphas2[[ii]][[jj]][1] / split_ratio) *
            rho
          dists <- 0
          indices <- expand.grid(which(out$dim == ii), which(out$dim == jj))
          indices <- indices[indices[,1] < indices[,2],]
          dists <- sum(B[cbind(indices$Var2, indices$Var1)]*(y[indices$Var2] - y[indices$Var1]))
          # dists <- 0
          # for (ki in 2:N) {
          #   for (kj in 1:(ki - 1)) {
          #     if (out$dim[ki] == jj & out$dim[kj] == ii) {
          #       dists <- dists + (y[ki] - y[kj]) * B[ki, kj]
          #     }
          #   }
          # }
          new_betas2[[ii]][[jj]][1] <- dists + eta_alphas1[[ii]][[jj]][1]/eta_alphas2[[ii]][[jj]][1]*sum(T*split_ratio - ys_truncated[[ii]]$timestamp)
          eta_betas2[[ii]][[jj]][1] <-
            eta_betas2[[ii]][[jj]][1] * (1 - rho) + (s_beta + new_betas2[[ii]][[jj]][1] /
                                                             split_ratio) * rho
        }
      }
      ## update B
    }
    # calculate the time used to calculate elbo
    time_old <- time
    time <- difftime(Sys.time(), clock_start, units = 'mins') - time_ELBO
    if (time_old < 1 & time > 1) {
      {
        elbo_start <- Sys.time()
        time_ELBO <-
          time_ELBO + difftime(Sys.time(), elbo_start, units = 'mins')
      }
      output <- list(
        eta_mu1 = eta_mu1,
        eta_mu2 = eta_mu2,
        eta_alphas1 = eta_alphas1,
        eta_alphas2 = eta_alphas2,
        eta_betas1 = eta_betas1,
        eta_betas2 = eta_betas2,
        times = times
      )
      if (type == "cluster") {
        save(
          output,
          file = paste0(
            "/gscratch/home/jiang14/mhp/output_long/s",
            files[split_id],
            "/mhpvisgdtr_seed_",
            taskid,
            "_dataset_" ,
            dataset_seed,
            "_time01.RData"
          ))
      }
    }
    if (time_old < 3 & time > 3) {
      {
        elbo_start <- Sys.time()
        time_ELBO <-
          time_ELBO + difftime(Sys.time(), elbo_start, units = 'mins')
      }
      output <- list(
        eta_mu1 = eta_mu1,
        eta_mu2 = eta_mu2,
        eta_alphas1 = eta_alphas1,
        eta_alphas2 = eta_alphas2,
        eta_betas1 = eta_betas1,
        eta_betas2 = eta_betas2,
        times = times
      )
      if (type == "cluster") {
        save(
          output,
          file = paste0(
            "/gscratch/home/jiang14/mhp/output_long/s",
            files[split_id],
            "/mhpvisgdtr_seed_",
            taskid,
            "_dataset_" ,
            dataset_seed,
            "_time03.RData"
          ))
      }
    }
  }
  {
    elbo_start <- Sys.time()
    time_ELBO <-
      time_ELBO + difftime(Sys.time(), elbo_start, units = 'mins')
  }
  # update ELBO
  output <- list(
    eta_mu1 = eta_mu1,
    eta_mu2 = eta_mu2,
    eta_alphas1 = eta_alphas1,
    eta_alphas2 = eta_alphas2,
    eta_betas1 = eta_betas1,
    eta_betas2 = eta_betas2,
    times = times
  )
  if (type == "cluster") {
    save(
      output,
      file = paste0(
        "/gscratch/home/jiang14/mhp/output_long/s",
        files[split_id],
        "/mhpvisgdtr_seed_",
        taskid,
        "_dataset_" ,
        dataset_seed,
        "_time05.RData"
      )
    )
  }
}
}

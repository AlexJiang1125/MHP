### MCMC (latent variable version)
# version: use original likelihood


rm(list = ls())
library(dplyr)

CompleteLik_diff_jk <- function(alpha, beta, alpha_new, beta_new, trans, dist, ndim, out, dim_from) {
  sum1_delta <- trans*(log(alpha_new) + log(beta_new) - log(alpha) - log(beta))
  sum2_delta <- -dist*(beta_new - beta) - ndim*(alpha_new - alpha)
  sum3_delta <- alpha_new*sum(exp(-beta_new*(T_all - out[out$dim==dim_from,"timestamp"]))) - alpha*sum(exp(-beta*(T_all - out[out$dim==dim_from,"timestamp"])))
  return(sum1_delta + sum2_delta + sum3_delta)
}


# hyperparameters

a_mu <- 2
b_mu <- 4
e_alpha <- 2
f_alpha <- 4
r_beta <- 2
s_beta <- 2 / 4
T_all <- 1000

# load in data
taskid <- as.numeric(commandArgs(trailingOnly = TRUE)[1])

pars <- expand.grid(
  dat_id = 1:50,
  MH = c(TRUE, FALSE)
)

dat_id <- pars$dat_id[taskid]
MH <- pars$MH[taskid]

load(
  paste0(
    "/mnt/beegfs/homes/jiang14/mhp/data/sim001/data/data/data_",
    dat_id,
    ".RData"
  )
)


# get MLE estimate as starting vals

mu_MLE <- NA
beta_MLE <- alpha_MLE <- NA
out <- out[order(out$timestamp),]
N <- length(out$timestamp)
y <- out$timestamp

tbl <- table(out$dim, out$parentdim)
for (ii in 1:K) {
  mu_MLE[ii] <- sum(out$dim == ii & out$parentdim == 0)/T_all
}
for (ii in 1:K) {
  for (jj in 1:K) {
    pos <- (ii-1)*3 + jj
    out_sub <- out[out$dim == jj & out$parentdim == ii,]
    res <- 0
    if (length(out_sub$id) == 0) {
      beta_MLE[pos] <- 0
    } else {
      for (kk in 1:length(out_sub$id)) {
        pid <- out_sub$parent[kk]
        res[kk] <- out_sub$timestamp[kk] - out[out$id == pid,]$timestamp
      }
      beta_MLE[pos] <- 1/mean(res)
    }
    alpha_MLE[pos] <- tbl[jj,ii+1]/sum(out$dim == ii)
  }
}
alpha_beta_mtp <- matrix(1, 3, 3)
alpha_MLE <- matrix(alpha_MLE,3, byrow = T)
beta_MLE <- matrix(beta_MLE,3, byrow = T)

# storing the chains

mus <- list()
alphas <- list()
betas <- list()

# latent variables

# a change in data structure
Bvec <- rep(0, N)
# assign values to true latent structure
Bvec[1] <- 1
for (i in 2:N) {
  if (out$parentdim[i] == 0) {
    Bvec[i] <- i
  } else {
    parentid <- out$parent[i]
    parentrowid <- which(out$id == parentid)
    Bvec[i] <- parentrowid
  }
}

col_max <- 200
B <- matrix(0, N, col_max)
for (i in 1:N) {
  if (Bvec[i] == i) {
    B[i,1] <- 1
  } else {
    x <- min(i - Bvec[i] + 1, col_max)
    B[i,x] <- 1
  }
}

B_true <- B

# get transition pairs
ids <- list()
pos <- 1
for (j in 1:3) {
  for (k in 1:3) {
    id_j <- which(out$dim == j)
    id_k <- which(out$dim == k)
    id_jk <- expand.grid(id_j, id_k)
    id_jk <- id_jk[id_jk$Var2 - id_jk$Var1 < col_max, ]
    id_jk <- id_jk[id_jk$Var2 - id_jk$Var1 > 0, ]
    ids[[pos]] <- id_jk
    pos <- pos + 1
  }
}


# not updating anything else

alpha <- alpha_MLE
beta <- beta_MLE
mu <- mu_MLE

alphas_old  <- list()#output$alpha_sgld[c(1:count)]
betas_old  <- list()#output$beta_sgld[c(1:count)]
mus_old  <- list()#output$beta_sgld[c(1:count)]


count <- 1

while (count <= 15000) {
  for (id_i in 1:3) {
    ndim <- sum(out$dim == id_i)
    ts <- out$timestamp[out$dim == id_i]
    
    for (id_j in 1:3) {
      pos <- (id_i - 1) * 3 + id_j
      ntrans <- sum(B[cbind(ids[[pos]]$Var2, ids[[pos]]$Var2 - ids[[pos]]$Var1 + 1)])
      
      dists_B <- B[cbind(ids[[pos]]$Var2, ids[[pos]]$Var2 - ids[[pos]]$Var1 + 1)]
      dists_y <- y[ids[[pos]]$Var2] - y[ids[[pos]]$Var1]
      dists <- sum(dists_B * dists_y)
      
      alpha_s2 <- length(ts) - sum(exp(-beta[id_i, id_j]*(T_all - ts)))
      
      alpha[id_i, id_j] <- rgamma(1, e_alpha + ntrans, f_alpha + alpha_s2)
      if (MH == FALSE) {
        beta[id_i, id_j] <- rgamma(1, r_beta + ntrans, s_beta + dists)
      } else {
        beta_new <- beta
        beta_new[id_i, id_j] <- exp(rnorm(1, mean = log(beta[id_i, id_j]), sd = 0.05))
        loglik_ratio <- CompleteLik_diff_jk(alpha = alpha[id_i, id_j], beta = beta[id_i, id_j], 
                                            alpha_new = alpha[id_i, id_j], beta_new = beta_new[id_i, id_j], 
                                            trans = ntrans, dist = dists, ndim = ndim, out = out,
                                            dim_from = id_i)
        if (log(runif(1, 0, 1)) < loglik_ratio) {
          beta <- beta_new
        }
      }
    }
    mu[id_i] <- rgamma(1, shape = sum(B[out$dim == id_i,1]) + a_mu, rate = T_all + b_mu)
  }
  
  #ids_B <- sample(2:N, size = 500)
  
  for (ii in 2:N) {
    ptemp <- rep(0,min(col_max, ii))
    ptemp[1] <- log(mu[out$dim[ii]])
    id_start <- max(1, ii - col_max + 1)
    id_end <- ii - 1
    ptemp[2:min(col_max, ii)] <-
      log(alpha[out$dim[id_end:id_start], out$dim[ii]]) + log(beta[out$dim[id_end:id_start], out$dim[ii]]) - beta[out$dim[id_end:id_start], out$dim[ii]] *
      (y[ii] - y[id_end:id_start])
    res_temp <- exp(ptemp - max(ptemp))/sum(exp(ptemp - max(ptemp)))
    sampled <- sample(1:min(col_max,ii), size = 1, prob = res_temp)
    B[ii,] <- 0
    B[ii,sampled] <- 1#sample(1:min(col_max,ii), size = 1, prob = res_temp)
  }
  #print(beta)
  alphas_old[[count]] <- alpha
  betas_old[[count]] <- beta
  mus_old[[count]] <- mu
  
  count <- count + 1
  
  if (count == 10 || count %% 1000 == 0) {
    count_start <- count / 2
    alphas_new <- simplify2array(alphas_old[c((count_start+1):count)])
    betas_new <- simplify2array(betas_old[c((count_start+1):count)])
    mus_new <- simplify2array(mus_old[c((count_start+1):count)])
    results <- list(
      alphas_new = alphas_new,
      betas_new = betas_new,
      mus_new = mus_new
      #B = B,
      #B_true = B_true
    )
    save(
      results,
      file = paste0("/mnt/beegfs/homes/jiang14/mhp/output_mcmc/mcmc_dataset_", dat_id, "_MH_", MH, ".RData")
    )
  }
}


alphas_new <- simplify2array(alphas_old[c(5001:15000)])
betas_new <- simplify2array(betas_old[c(5001:15000)])
mus_new <- simplify2array(mus_old[c(5001:15000)])
results <- list(
  alphas_new = alphas_new,
  betas_new = betas_new,
  mus_new = mus_new
  #B = B,
  #B_true = B_true
)
save(
  results,
  file = paste0("/mnt/beegfs/homes/jiang14/mhp/output_mcmc/mcmc_dataset_", dat_id, "_MH_", MH, ".RData")
)
# update beta


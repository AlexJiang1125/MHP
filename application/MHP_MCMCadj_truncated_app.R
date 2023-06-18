### MCMC (latent variable version)
# version: use original likelihood


rm(list = ls())
library(dplyr)

CompleteLik_diff_jk <- function(alpha, beta, alpha_new, beta_new, trans, dist, ndim, out) {
  sum1_delta <- trans*(log(alpha_new) + log(beta_new) - log(alpha) - log(beta))
  sum2_delta <- -dist*(beta_new - beta) - ndim*(alpha_new - alpha)
  sum3_delta <- alpha_new*sum(exp(-beta_new*(T - out[out$dim==j,"timestamp"]))) - alpha*sum(exp(-beta*(T - out[out$dim==j,"timestamp"])))
  return(sum1_delta + sum2_delta + sum3_delta)
}


# hyperparameters

a_mu <- 2
b_mu <- 4
e_alpha <- 2
f_alpha <- 4
r_beta <- 2
s_beta <- 2 / 4


# load in data
taskid <- as.numeric(commandArgs(trailingOnly = TRUE)[1])
dat_id <- taskid

mode <- "cluster"

if (mode == "local") {
  load("/Users/alexziyujiang/Dropbox/data/playground/data_app.RData")
} else {
  load("/mnt/beegfs/homes/jiang14/mhp/data/sim001/data/data_app/data_app.RData")
}



# get MLE estimate as starting vals

if (mode == "local") {
  load("/Users/alexziyujiang/Dropbox/data/playground/EM_app_inits/inits_1.RData")
} else {
  load(paste0("/mnt/beegfs/homes/jiang14/mhp/data/sim001/EM_app_inits/inits_", taskid ,".RData"))
}

alpha <- inits[[1]]
beta <- inits[[2]]
mu <- inits[[3]]


out <- out[order(out$timestamp),]
N <- length(out$timestamp)
y <- out$timestamp


# storing the chains

mus <- list()
alphas <- list()
betas <- list()

# latent variables

# a change in data structure
Bvec <- 1:N#rep(0, N)
# assign values to true latent structure
# Bvec[1] <- 1
# for (i in 2:N) {
#   if (out$parentdim[i] == 0) {
#     Bvec[i] <- i
#   } else {
#     parentid <- out$parent[i]
#     parentrowid <- which(out$id == parentid)
#     Bvec[i] <- parentrowid
#   }
# }

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
K <- 11
# get transition pairs
ids <- list()
pos <- 1
for (j in 1:K) {
  for (k in 1:K) {
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

alphas_old  <- list()#output$alpha_sgld[c(1:count)]
betas_old  <- list()#output$beta_sgld[c(1:count)]
mus_old  <- list()#output$beta_sgld[c(1:count)]

T <- 1000
count <- 1

out_tail <- out[out$timestamp >= T - 0.25, ]
ys_truncated <- list()
for (id_k in 1:K) {
  ys_truncated[[id_k]] <- out_tail[out_tail$dim == id_k,]
}

while (count <= 20000) {
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
  for (id_i in 1:K) {
    ndim <- sum(out$dim == id_i)
    for (id_j in 1:K) {
      pos <- (id_i - 1) * K + id_j
      ntrans <- sum(B[cbind(ids[[pos]]$Var2, ids[[pos]]$Var2 - ids[[pos]]$Var1 + 1)])
      
      dists_B <- B[cbind(ids[[pos]]$Var2, ids[[pos]]$Var2 - ids[[pos]]$Var1 + 1)]
      dists_y <- y[ids[[pos]]$Var2] - y[ids[[pos]]$Var1]
      dists <- sum(dists_B * dists_y)
      
      alpha_compensator_remainder <- sum(exp(-beta[id_i, id_j]*(T - out$timestamp[out$dim == id_i])))
      truncation_correction <- alpha[id_i, id_j]*sum(T - ys_truncated[[id_i]]$timestamp)
      alpha[id_i, id_j] <- rgamma(1, e_alpha + ntrans, f_alpha + ndim - alpha_compensator_remainder)
      beta[id_i, id_j] <- rgamma(1, r_beta + ntrans, s_beta + dists + truncation_correction)
    }
    mu[id_i] <- rgamma(1, shape = sum(B[out$dim == id_i,1]) + a_mu, rate = T + b_mu)
  }
  
  #ids_B <- sample(2:N, size = 500)
  
  print(mean(alpha))
  
  
  alphas_old[[count]] <- alpha
  betas_old[[count]] <- beta
  mus_old[[count]] <- mu
  
  count <- count + 1
}


alphas_new <- simplify2array(alphas_old[c(5001:20000)])
betas_new <- simplify2array(betas_old[c(5001: 20000)])
mus_new <- simplify2array(mus_old[c(5001: 20000)])


results <- list(
  alphas_new = alphas_new,
  betas_new = betas_new,
  mus_new = mus_new
  #B = B,
  #B_true = B_true
)
save(
  results,
  file = paste0("/mnt/beegfs/homes/jiang14/mhp/output_app/mcmc_", dat_id,"_truncated.RData")
)
# update beta


###############################################
# Date Last modified: Feb.5, 2023 AZJ
# Program: gengraphs_app_methods.R
# Purpose: generate graphs for model fitting results
###############################################

# Data inputs:
#   - the SP&INIT_START intraday dataset
#   - fitted model .Rda files for all methods
# Data outputs
#   - model output graphs
# Note:
#   - run locally

library(dplyr)
library(ggplot2)
library(ggrepel)
library(vegan)
library(cowplot)


rm(list = ls())

theme_set(theme_bw())

INIT_START <- 1000
INIT_START_SGLD <- 1000
UPP_LIMIT_ALPHA <- 0.3
UPP_LIMIT_ALPHA_UNC <- 0.18
UPP_LIMIT_BETA <- 120
UPP_LIMIT_BETA_UNC <- 100
UPP_LIMIT_MU <- 4
UPP_LIMIT_MU_UNC <- 1.8
BAR_HEIGHT <- 5
METHOD_ORDER <- c(7,6,4,3,2,1,5)

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

MODEL <- "local"

if (MODEL == "local") {
  DATA_PATH <- "/Users/alexziyujiang/Documents/GitHub/HawkesSPINIT_START/dataset_intraday/data_app.RData"
  RESULTS_PATH <- "/Users/alexziyujiang/Dropbox/data/mhp_app/output_new/"
  LIK_PATH <- "likelihood/"
} else {
  DATA_PATH <- "/gscratch/home/jiang14/mhp/data/sim001/data/data_app/data_app.RData"
  RESULTS_PATH <- "/gscratch/csde/jiang14/output_app/s005/"
}

setwd(RESULTS_PATH)

# name components

METHOD_NAMES <- c("sem", "semtr", "sgd", "sgdtr", "sgld", "mcmc", "mcmc_truncated")
METHOD_NAMES_FORPLOT <- c("SEM", "SEM-c", "SGVI", "SGVI-c", "SGLD", "MCMC", "MCMC-c")
NAME_PREFIX <- c("mhpvi")
NAME_SUFFICES <- c("_time30.RData", "_time30.RData", "newinit_time30.RData", "newinit_time30.RData", "_newinit_time30.RData")

# graphics
PLOT_HEIGHT <- 10
PLOT_WIDTH <- 22

# parameters for the dataset
K <- 11 # number of dimensions
return_names <- c("S5COND", "S5CONS", "S5ENRS", "S5FINL",
                  "S5HLTH", "S5INDU", "S5MATR", "S5INFT", "S5TELS", "S5UTIL", "S5RLST")
return_names_p <- sub('..', '', return_names)

# plot heatmap for the alpha matrix
# input: alpha matrix, a vector for the seed that c
# output: the graph
# name: "app_alpha_heatmap.pdf"

REODERED_DIMS <- c(1, 10, 7, 4, 8,6, 5,2,11, 9, 3)

plot_mu_heatmap <- function(mu.mat, method_id, upp_limit, bar_height) {
  mus <- mu.mat %>% data.frame()
  if (nrow(mus) == 1) {
    mus <- mus %>% tidyr::pivot_longer(cols = c(1:K))
    colnames(mus) <- c("dimensions", "value")
    mus$dimensions <- 1:K
  } else {
    mus$dimensions <- 1:K
    colnames(mus) <- c("value", "dimensions")
  }
  mus$dimensions <- factor(mus$dimensions, levels = REODERED_DIMS)
  mus$y <- factor(1, levels = "")
  mus <- mus |> mutate(dimensions = plyr::mapvalues(dimensions, 1:K, return_names_p))

  p <- ggplot(mus, aes(x = dimensions, y=y, fill=(value))) +
    geom_tile()+ scale_fill_gradient(low = "white", high = "black", name = "", limits=c(0, upp_limit))+ ggtitle(METHOD_NAMES_FORPLOT[method_id]) + ylab("")+
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.x = element_text(size = 4.5)
    ) + theme(legend.position="right") +
    guides(fill = guide_colourbar(barheight = bar_height))
  return(p)
}

plot_alpha_heatmap <- function(alpha.mat, method_id, upp_limit,bar_height) {
  alphas <- alpha.mat %>% as.data.frame()
  colnames(alphas) <- 1:K
  alphas <- alphas %>% tidyr::pivot_longer(cols = c(1:K)) %>% mutate(send = rep(1:K, each = K))
  colnames(alphas) <- c("receiver", "value", "triggerer")
  alphas$receiver <- factor(alphas$receiver, levels = REODERED_DIMS)
  alphas$triggerer <- factor(alphas$triggerer, levels = REODERED_DIMS)
  alphas <- alphas |> mutate(receiver = plyr::mapvalues(receiver, 1:K, return_names_p),
                             triggerer = plyr::mapvalues(triggerer, 1:K, return_names_p))
  p <- ggplot(alphas, aes(x = receiver, y = triggerer, fill=(value))) +
    geom_tile()+ scale_fill_gradient(low = "white", high = "black", name = "", limits=c(0, upp_limit)) +
    ggtitle(METHOD_NAMES_FORPLOT[method_id]) + theme(legend.position="right") +
    theme(axis.text = element_text(size = 4.5)) + theme(legend.position = "right") +
    guides(fill = guide_colorbar(barheight = bar_height, barwidth = 2, label.theme = element_text(size = 12)))
  return(p)
}

plot_alpha_heatmap_forall <- function(seed_ids,bar_height) {
  # preprocess alpha
  plot_list <- list()
  # loop over all five methods
  for (i in 1:length(METHOD_NAMES)) {
    load(generate_file_name(i, seed_ids[i]))
    if (i == 5) {
      output <- output_app
    } else if (i > 5) {
      output <- results
    }
    pars <- preprocess(output, METHOD_NAMES[i])
    # calculate likelihood
    alpha.mean <- pars$alpha.mean
    plot_list[[i]] <- plot_alpha_heatmap(alpha.mat = alpha.mean,
                                         method_id = i,
                                         upp_limit = UPP_LIMIT_ALPHA,
                                         bar_height = bar_height)
  }
  return(plot_list)
}

# calculate a 7x7 table of estimation results, on the log scale
table_alpha_dist <- function(seed_ids) {
  # preprocess alpha
  mat_list <- list()
  # loop over all five methods
  for (i in 1:length(METHOD_NAMES)) {
    load(generate_file_name(i, seed_ids[i]))
    if (i == 5) {
      output <- output_app
    } else if (i > 5) {
      output <- results
    }
    pars <- preprocess(output, METHOD_NAMES[i])
    # calculate likelihood
    alpha.mean <- pars$alpha.mean
    mat_list[[i]] <- alpha.mean
  }
  dist_mat <- matrix(0, nrow = 7, ncol = 7)
  for (i in 1:7) {
    for (j in 1:7) {
      dist_mat[i,j] <- mean((log(mat_list[[i]]) - log(mat_list[[j]]))^2)
    }
  }
  return(dist_mat)
}

table_beta_dist <- function(seed_ids) {
  # preprocess alpha
  mat_list <- list()
  # loop over all five methods
  for (i in 1:length(METHOD_NAMES)) {
    load(generate_file_name(i, seed_ids[i]))
    if (i == 5) {
      output <- output_app
    } else if (i > 5) {
      output <- results
    }
    pars <- preprocess(output, METHOD_NAMES[i])
    # calculate likelihood
    beta.mean <- pars$beta.mean
    mat_list[[i]] <- beta.mean
  }
  dist_mat <- matrix(0, nrow = 7, ncol = 7)
  for (i in 1:7) {
    for (j in 1:7) {
      dist_mat[i,j] <- mean((log(mat_list[[i]]) - log(mat_list[[j]]))^2)
    }
  }
  return(dist_mat)
}

table_mu_dist <- function(seed_ids) {
  # preprocess alpha
  mat_list <- list()
  # loop over all five methods
  for (i in 1:length(METHOD_NAMES)) {
    load(generate_file_name(i, seed_ids[i]))
    if (i == 5) {
      output <- output_app
    } else if (i > 5) {
      output <- results
    }
    pars <- preprocess(output, METHOD_NAMES[i])
    # calculate likelihood
    mu.mean <- pars$mu.mean
    mat_list[[i]] <- mu.mean
  }
  dist_mat <- matrix(0, nrow = 7, ncol = 7)
  for (i in 1:7) {
    for (j in 1:7) {
      dist_mat[i,j] <- mean((log(mat_list[[i]]) - log(mat_list[[j]]))^2)
    }
  }
  return(dist_mat)
}



# plot the heatmaps of uncertainty interval lengths
plot_alpha_heatmap_forall_uncertainties <- function(seed_ids, alpha, bar_height) {
  # preprocess alpha
  plot_list <- list()
  # loop over all five methods
  for (i in 1:length(METHOD_NAMES)) {
    load(generate_file_name(i, seed_ids[i]))
    if (i == 5) {
      output <- output_app
    } else if (i > 5) {
      output <- results
    }
    pars <- preprocess_uncertainties(output, METHOD_NAMES[i], alpha)
    # calculate likelihood
    alpha.cred <- pars$alpha.cred
    plot_list[[i]] <- plot_alpha_heatmap(alpha.mat = alpha.cred,
                                         method_id = i,
                                         upp_limit = UPP_LIMIT_ALPHA_UNC,
                                         bar_height = bar_height)
  }
  return(plot_list)
}

plot_alpha_mds <- function(mds.coor, method_id) {
  mds.coor <- as.data.frame(mds.coor)
  mds.coor$sectors <- return_names_p
  colnames(mds.coor) <- c("X1", "X2", "sectors")
  p <- ggplot(mds.coor, aes(x = X1, y = X2, label = sectors)) + geom_text(hjust=0, vjust=0) + geom_point()
  p <- p + xlim(c(min(mds.coor$X1) - 0.1, max(mds.coor$X1) + 0.1)) +
    ylim(c(min(mds.coor$X2) - 0.05, max(mds.coor$X2) + 0.05))
  p <- p + geom_hline(yintercept = 0, lty = 2, col = 'gray50') +
    geom_vline(xintercept = 0, lty = 2, col ="gray50")
  p <- p + ggtitle(METHOD_NAMES_FORPLOT[method_id]) + xlab("") + ylab("")
  return(p)
}

plot_alpha_mds_forall  <- function(seed_ids) {
  # preprocess alpha
  plot_list <- list()
  # list of alpha matrices (for storage)
  mds_list <- list()

  # loop over all five methods
  for (i in 1:length(METHOD_NAMES)) {
    load(generate_file_name(i, seed_ids[i]))
    if (i == 5) {
      output <- output_app
    } else if (i > 5) {
      output <- results
    }
    pars <- preprocess(output, METHOD_NAMES[i])
    # calculate likelihood
    alpha.mean <- pars$alpha.mean
    alpha.mean.dist <- exp(-alpha.mean)
    diag(alpha.mean.dist) <- 0
    mds.coor <- cmdscale(alpha.mean.dist)
    rownames(mds.coor) <- return_names_p
    mds.coor <- data.frame(mds.coor)
    if (i == 1) {
      mds_list[[i]] <- mds.coor
    } else {
      mds_list[[i]] <- vegan::procrustes(mds_list[[1]], mds.coor, scale = TRUE)$Yrot
    }
  }
  for (i in 1:length(METHOD_NAMES)) {
    plot_list[[i]] <- plot_alpha_mds(mds_list[[i]], i)
  }
  return(plot_list)
}

plot_beta_heatmap <- function(beta.mean, method_id, upp_limit, bar_height) {
  betas <- beta.mean %>% as.data.frame()
  colnames(betas) <- 1:K
  betas <- betas %>% tidyr::pivot_longer(cols = c(1:K)) %>% mutate(send = rep(1:K, each = K),)
  colnames(betas) <- c("receiver", "value", "triggerer")
  betas$receiver <- factor(betas$receiver, levels = REODERED_DIMS)
  betas$triggerer <- factor(betas$triggerer, levels = REODERED_DIMS)
  betas <- betas |> mutate(receiver = plyr::mapvalues(receiver, 1:K, return_names_p),
                           triggerer = plyr::mapvalues(triggerer, 1:K, return_names_p))
  p <- ggplot(betas, aes(receiver, triggerer, fill=(value))) +
    geom_tile()+ scale_fill_gradient(low = "white", high = "black", name = "", limits=c(0, upp_limit)) + ggtitle(METHOD_NAMES_FORPLOT[method_id]) +
    theme(axis.text = element_text(size = 4.5)) + theme(legend.position = "right") +
    guides(fill = guide_colorbar(barheight = bar_height, barwidth = 2, label.theme = element_text(size = 12)))
  return(p)
}

plot_beta_heatmap_forall <- function(seed_ids, bar_height) {
  # preprocess beta
  plot_list <- list()
  # loop over all five methods
  for (i in 1:length(METHOD_NAMES)) {
    load(generate_file_name(i, seed_ids[i]))
    if (i == 5) {
      output <- output_app
    } else if (i > 5) {
      output <- results
    }
    pars <- preprocess(output, METHOD_NAMES[i])
    # calculate likelihood
    beta.mean <- pars$beta.mean
    plot_list[[i]] <- plot_beta_heatmap(beta.mean = beta.mean,
                                        method_id = i,
                                        upp_limit = UPP_LIMIT_BETA,
                                        bar_height = bar_height)
  }
  return(plot_list)
}

plot_beta_heatmap_forall_uncertainties <- function(seed_ids, alpha, bar_height) {
  # preprocess beta
  plot_list <- list()
  # loop over all five methods
  for (i in 1:length(METHOD_NAMES)) {
    load(generate_file_name(i, seed_ids[i]))
    if (i == 5) {
      output <- output_app
    } else if (i > 5) {
      output <- results
    }
    pars <- preprocess_uncertainties(output, METHOD_NAMES[i], alpha = alpha)
    # calculate likelihood
    beta.cred <- pars$beta.cred
    plot_list[[i]] <- plot_beta_heatmap(beta.mean = beta.cred,
                                        method_id = i,
                                        upp_limit = UPP_LIMIT_BETA_UNC,
                                        bar_height = bar_height)
  }
  return(plot_list)
}

plot_beta_mds <- function(mds.coor, method_id) {
  mds.coor <- as.data.frame(mds.coor)
  mds.coor$sectors <- return_names_p
  colnames(mds.coor) <- c("X1", "X2", "sectors")
  p <- ggplot(mds.coor, aes(x = X1, y = X2, label = sectors)) + geom_text(hjust=0, vjust=0) + geom_point()
  p <- p + xlim(c(min(mds.coor$X1) - 0.1, max(mds.coor$X1) + 0.1)) +
    ylim(c(min(mds.coor$X2) - 0.05, max(mds.coor$X2) + 0.05))
  p <- p + geom_hline(yintercept = 0, lty = 2, col = 'gray50') +
    geom_vline(xintercept = 0, lty = 2, col ="gray50")
  p <- p + ggtitle(METHOD_NAMES_FORPLOT[method_id]) + xlab("") + ylab("")
  return(p)
}

plot_beta_mds_forall  <- function(seed_ids) {
  # preprocess beta
  plot_list <- list()
  # list of beta matrices (for storage)
  mds_list <- list()

  # loop over all five methods
  for (i in 1:length(METHOD_NAMES)) {
    load(generate_file_name(i, seed_ids[i]))
    if (i == 5) {
      output <- output_app
    } else if (i > 5) {
      output <- results
    }
    pars <- preprocess(output, METHOD_NAMES[i])
    # calculate likelihood
    beta.mean <- pars$beta.mean
    beta.mean.dist <- exp(-beta.mean)
    diag(beta.mean.dist) <- 0
    mds.coor <- cmdscale(beta.mean.dist)
    rownames(mds.coor) <- return_names_p
    mds.coor <- data.frame(mds.coor)
    if (i == 1) {
      mds_list[[i]] <- mds.coor
    } else {
      mds_list[[i]] <- vegan::procrustes(mds_list[[1]], mds.coor, scale = TRUE)$Yrot
    }
  }
  for (i in 1:length(METHOD_NAMES)) {
    plot_list[[i]] <- plot_beta_mds(mds_list[[i]], i)
  }
  return(plot_list)
}

plot_mu_heatmap_forall <- function(seed_ids, bar_height) {
  # preprocess mu
  plot_list <- list()
  # loop over all five methods
  for (i in 1:length(METHOD_NAMES)) {
    load(generate_file_name(i, seed_ids[i]))
    if (i == 5) {
      output <- output_app
    } else if (i > 5) {
      output <- results
    }
    pars <- preprocess(output, METHOD_NAMES[i])
    # calculate likelihood
    mu.mean <- pars$mu.mean
    plot_list[[i]] <- plot_mu_heatmap(mu.mat = mu.mean,
                                      method_id = i,
                                      upp_limit = UPP_LIMIT_MU,
                                      bar_height = bar_height)
  }
  return(plot_list)
}

plot_mu_heatmap_forall_uncertainties <- function(seed_ids, alpha, bar_height) {
  # preprocess mu
  plot_list <- list()
  # loop over all five methods
  for (i in 1:length(METHOD_NAMES)) {
    load(generate_file_name(i, seed_ids[i]))
    if (i == 5) {
      output <- output_app
    } else if (i > 5) {
      output <- results
    }
    pars <- preprocess_uncertainties(output, METHOD_NAMES[i], alpha = alpha)
    # calculate likelihood
    mu.cred <- pars$mu.cred
    plot_list[[i]] <- plot_mu_heatmap(mu.mat = mu.cred,
                                      method_id = i,
                                      upp_limit = UPP_LIMIT_MU_UNC,
                                      bar_height = bar_height)
  }
  return(plot_list)
}

# read in likelihood data
# determine seed ids
dfs <- list()
pos <- 1
for (method_id in 1:5) {
  for (seed_id in 1:16) {
    dfs[[pos]] <- read.csv(paste0(LIK_PATH, "lik_app_", METHOD_NAMES[method_id], "_seed_", seed_id, ".csv"))[,-1]
    pos <- pos + 1
  }
}
df <- dplyr::bind_rows(dfs, .id = "column_label")[,-1]
df <- df %>%  as_tibble() %>% group_by(model) %>% summarise(max = max(lik), seed_id = which.max(lik))

# check if names match
METHOD_NAMES
df$seed_id

seed_ids <- c(df$seed_id,1,1)

p <- plot_alpha_heatmap_forall(seed_ids, bar_height = 18*3)
p_unc <- plot_alpha_heatmap_forall_uncertainties(seed_ids, alpha = 0.05, bar_height = 20*2)
p2 <- plot_alpha_mds_forall(seed_ids)

ggsave(gridExtra::grid.arrange(gridExtra::arrangeGrob(p[[1]] + theme(legend.position="none"), p2[[1]],
                              p[[2]] + theme(legend.position="none"), p2[[2]],
                              p[[3]] + theme(legend.position="none"), p2[[3]],
                              p[[4]] + theme(legend.position="none"), p2[[4]],
                              p[[7]] + theme(legend.position="none"), p2[[7]],
                              p[[6]] + theme(legend.position="none"), p2[[6]],
                              p[[5]] + theme(legend.position="none"), p2[[5]], ncol = 4), g_legend(p[[1]]), widths = c(8,1), ncol = 2),
       file = "figures/app_alpha_heatmap.pdf", width = 19, height = 18)

ggsave(gridExtra::grid.arrange(gridExtra::arrangeGrob(p_unc[[3]] + theme(legend.position="none"),
                                                      p_unc[[4]] + theme(legend.position="none"),
                                                      p_unc[[5]] + theme(legend.position="none"),
                                                      p_unc[[6]] + theme(legend.position="none"),
                                                      p_unc[[7]] + theme(legend.position="none"), ncol = 3),g_legend(p_unc[[3]]), widths = c(8,1), ncol = 2),
       file = "figures/app_alpha_unc_heatmap.pdf", width = 16, height = 10)

p <- plot_beta_heatmap_forall(seed_ids, bar_height = 18*3)
p_unc <- plot_beta_heatmap_forall_uncertainties(seed_ids, alpha = 0.05, bar_height = 20*2)

p2 <- plot_beta_mds_forall(seed_ids)

ggsave(gridExtra::grid.arrange(gridExtra::arrangeGrob(p[[1]] + theme(legend.position="none"),# p2[[1]],
                                                      p[[2]] + theme(legend.position="none"),# p2[[2]],
                                                      p[[3]] + theme(legend.position="none"),# p2[[3]],
                                                      p[[4]] + theme(legend.position="none"),# p2[[4]],
                                                      p[[7]] + theme(legend.position="none"),# p2[[7]],
                                                      p[[6]] + theme(legend.position="none"),# p2[[6]],
                                                      p[[5]] + theme(legend.position="none"), ncol = 2), g_legend(p[[1]]), widths = c(8,1), ncol = 2),
       file = "figures/app_beta_heatmap.pdf", width = 10, height = 18)

ggsave(gridExtra::grid.arrange(gridExtra::arrangeGrob(p_unc[[3]] + theme(legend.position="none"),
                                                      p_unc[[4]] + theme(legend.position="none"),
                                                      p_unc[[5]] + theme(legend.position="none"),
                                                      p_unc[[6]] + theme(legend.position="none"),
                                                      p_unc[[7]] + theme(legend.position="none"), ncol = 3),g_legend(p_unc[[3]]), widths = c(8,1), ncol = 2),
       file = "figures/app_beta_unc_heatmap.pdf", width = 16, height = 10)

p <- plot_mu_heatmap_forall(seed_ids, bar_height = 18)
p_unc <- plot_mu_heatmap_forall_uncertainties(seed_ids, alpha = 0.05, bar_height = 14)

ggsave(gridExtra::grid.arrange(gridExtra::arrangeGrob(p[[1]] + theme(legend.position="none"),
                                                      p[[2]] + theme(legend.position="none"),
                                                      p[[3]] + theme(legend.position="none"),
                                                      p[[4]] + theme(legend.position="none"),
                                                      p[[6]] + theme(legend.position="none"),
                                                      p[[7]] + theme(legend.position="none"),
                                                      p[[5]] + theme(legend.position="none"),
                                                      ncol = 2), g_legend(p[[1]]), widths = c(8,1), ncol = 2),
       file = "figures/app_mu_heatmap.pdf", width = PLOT_HEIGHT, height = 4.6)


# ggsave(gridExtra::arrangeGrob(p[[1]],p[[2]],p[[3]],p[[4]],p[[6]], p[[7]], p[[5]], ncol = 2),
#        file = "figures/app_mu_heatmap.pdf", width = PLOT_HEIGHT, height = 5)

ggsave(gridExtra::grid.arrange(gridExtra::arrangeGrob(p_unc[[3]] + theme(legend.position="none"),
                              p_unc[[4]] + theme(legend.position="none"),
                              p_unc[[6]] + theme(legend.position="none"),
                              p_unc[[7]] + theme(legend.position="none"),
                              p_unc[[5]] + theme(legend.position="none"),
                              ncol = 2), g_legend(p_unc[[3]]), widths = c(8,1), ncol = 2),
       file = "figures/app_mu_unc_heatmap.pdf", width = PLOT_HEIGHT , height = 3.5)



# generate latex tables

alpha_dist <- round(table_alpha_dist(seed_ids)[METHOD_ORDER,METHOD_ORDER],3)
colnames(alpha_dist) <- METHOD_NAMES_FORPLOT[METHOD_ORDER]
rownames(alpha_dist) <- METHOD_NAMES_FORPLOT[METHOD_ORDER]
xtable::xtable(alpha_dist, digits = 3)

beta_dist <- round(table_beta_dist(seed_ids)[METHOD_ORDER,METHOD_ORDER],3)
colnames(beta_dist) <- METHOD_NAMES_FORPLOT[METHOD_ORDER]
rownames(beta_dist) <- METHOD_NAMES_FORPLOT[METHOD_ORDER]
xtable::xtable(beta_dist, digits = 3)

mu_dist <- round(table_mu_dist(seed_ids)[METHOD_ORDER,METHOD_ORDER],3)
colnames(mu_dist) <- METHOD_NAMES_FORPLOT[METHOD_ORDER]
rownames(mu_dist) <- METHOD_NAMES_FORPLOT[METHOD_ORDER]
xtable::xtable(mu_dist, digits = 3)





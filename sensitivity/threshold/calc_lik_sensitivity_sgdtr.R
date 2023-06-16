## calculate likelihood for model outputs
rm(list = ls())
## data re-formatting

repackage <- function(x) {
  pos <- length(x[[1]][[1]])
  return(c(x[[1]][[1]][pos], x[[2]][[1]][pos], x[[3]][[1]][pos], x[[1]][[2]][pos], x[[2]][[2]][pos], x[[3]][[2]][pos], x[[1]][[3]][pos],
           x[[2]][[3]][pos], x[[3]][[3]][pos]))
}

repackage_pos <- function(x, pos) {
  return(c(x[[1]][[1]][pos], x[[2]][[1]][pos], x[[3]][[1]][pos], x[[1]][[2]][pos], x[[2]][[2]][pos], x[[3]][[2]][pos], x[[1]][[3]][pos],
           x[[2]][[3]][pos], x[[3]][[3]][pos]))
}

MargLik <- function(mu,alpha,beta,out) {
  n <- dim(out)[1]
  y <- out$timestamp
  res <- log(mu[out$dim[1]])
  Tdat <- T*1
  for (i in 2:n) {
    alpha_vec <- alpha[cbind(out$dim[1:(i-1)],rep(out$dim[i],(i-1)))]
    beta_vec <- beta[cbind(out$dim[1:(i-1)],rep(out$dim[i],(i-1)))]
    res <- res + log(sum(alpha_vec*beta_vec*exp(-beta_vec*(y[i] - y[1:(i-1)]))) + mu[out$dim[i]])
  }
  ns <- c(length(out[out$dim==1,1]),length(out[out$dim==2,1]),length(out[out$dim==3,1]))
  res <- res - sum(mu)*Tdat - sum(apply(alpha, MARGIN = 1, FUN = sum)*ns)
  for (i in 1:3) {
    for (j in 1:3) {
      res <- res + alpha[i,j]*sum(exp(-beta[i,j]*(Tdat - out$timestamp[out$dim == i])))
    }
  }
  return(res)
}

type <- "cluster"

if (type == "local") {
  id_i <- 1
  setwd("/Users/alexziyujiang/Documents/data/mhp/SGD/delay1learning51")
} else if (type == "cluster") {
  taskid <- as.numeric(commandArgs(trailingOnly = TRUE)[1])
  print(taskid)
}
seeds <- 1:16
datasets <- 1:50
foldernames <- c("s001", "s005", "s010","s020","s030","s040","s050","s060","s070","s080","s090",
                 "s100")
filenametypes <- c("", "", "", "")
times <- c("01", "03","05")
count <- 1
mliks <- NA
sub_ratios <- NA
datasets_list <- NA
seeds_list <- NA
times_list <- NA
sens <- NA
T <- 1000
ids <- expand.grid(
  id_i = 1:16,
  id_k = 1:4)


id_k <- ids$id_k[taskid]
id_i <- ids$id_i[taskid]
id_l <- 2

sens_foldernames <- c("t001", "t005", "t01", "t05")

for (id_j in 1:50) {
  print(id_j)
  if (type == "local") {
    load(paste0("/Users/alexziyujiang/Documents/data/data_short/data_",datasets[id_j],".RData"))
  } else {
    load(
      paste0(
        "/gscratch/home/jiang14/mhp/data/sim001/data/data/data_",
        datasets[id_j],
        ".RData"
      )
    )
  }
  filename <- paste0("/gscratch/csde/jiang14/output/sensitivity/", sens_foldernames[id_k],"/mhpvisgdtr_seed_",id_i,"_dataset_",id_j,"newinit_time30.RData")
  print(filename)
  if (file.exists(filename)) {
    print(paste(id_i, id_j, id_k, id_l, sep = " "))
    if (type == "local") {
      load(paste0("/mhpvisgd_seed_",seeds[id_i],"_dataset_",datasets[id_j],"_time",times[id_k],".RData"))
    } else {
      load(filename)
    }
    mu.mean <- output$eta_mu1[1, ] / output$eta_mu2[1, ]
    alpha.mean <- matrix(repackage(output$eta_alphas1) / repackage(output$eta_alphas2),nrow=3)
    beta.mean <- matrix(repackage(output$eta_betas1) / repackage(output$eta_betas2),nrow=3)
    mliks[count] <- MargLik(mu.mean, alpha.mean, beta.mean, out)
    print(mliks[count])
    datasets_list[count] <- id_j
    seeds_list[count] <- id_i
    sens[count] <- id_k
    count <- count + 1
  }
}

df_out <- data.frame(data = datasets_list,
                     seed = seeds_list,
                     mlik = mliks, sen = sens)

if (type == "local") {
  write.csv(df_out, file = paste0("mhpvisgd_seed_",seeds[id_i],".csv"))
} else if (type == "cluster") {
  write.csv(df_out, file = paste0("/gscratch/csde/jiang14/output/sensitivity/mhpvisgdtr_seed_",id_i,id_k,id_l,".csv"))
}


# 这里用的是老师上课讲的东西
rm(list = ls())
# Generate fixed effect and observed data-----------------------------------------------------------
library(doParallel)
library(foreach)
library(parallel)
set.seed(1)
n <- 100
Time <- 10
beta1 <- 1
beta2 <- 1
pi1 <- 0.6
sigma1 <- 2
sigma2 <- 10
U_i <- rbinom(n = n,size = 1,prob = 1 - pi1) + 1
X1 <- matrix(rnorm(n = n*Time),nrow = n,ncol = Time) # fixed effect
X2 <- matrix(rnorm(n = n*Time),nrow = n,ncol = Time)
z1 <- rnorm(n = n,mean = 0,sd = sigma1) # random effect
z2 <- rnorm(n = n,mean = 0,sd = sigma2)
g <- function(x){
  return(exp(x)/(1+exp(x)))
}

Y_ij <- array(dim = c(n,Time))
for(i in 1:n) {
  for(j in 1:Time) {
    # generate the Y_ij AKA observed data
    # the group item indicates the cluster for each subject
    if(U_i[i] == 1) {
      Y_ij[i,j] <- g(beta1 * X1[i,j] + z1[i])
    }
    else {
      Y_ij[i,j] <- g(beta2 * X2[i,j] + z2[i])
    }
    Y_ij[i,j] <- rbinom(n = 1,size = 1,prob = Y_ij[i,j])
  }
}
#



logLikelihood <- function(Z, # 这里的Z1，Z2都是1*n的列向量。
                          U,
                          Omega,
                          m = 1) {# EM进行到第几步
  wi1 <- which(U == 1)
  wi2 <- which(U == 2)
  #if(Omega[5] > 1) Omega[5] <- 1
  l1 <- sapply(wi1,function(i){
    TT <- sapply(1:Time,function(j){
      eta <- Omega[1] * X1[i,j] + Z[i]
      if(eta >= 600) eta <- 600
      TTT <- Y_ij[i,j] * (eta - log(1 + exp(eta))) + (1 - Y_ij[i,j]) * (- log(1 + exp(eta)))
      return(TTT)
    })
    lll <- (log(Omega[5]) - log(sqrt(2*pi) * Omega[3]) - Z[i]^2/(2*Omega[3]^2) + sum(TT))
    return(lll)
  })
  l1 <- sum(l1)
  
  l2 <- sapply(wi2,function(i){
    TT <- sapply(1:Time,function(j){
      eta <- Omega[2] * X2[i,j] + Z[i]
      
      if(eta >= 600) eta <- 600
      TTT <- Y_ij[i,j] * (eta - log(1 + exp(eta))) + (1 - Y_ij[i,j]) * (- log(1 + exp(eta)))
      return(TTT)
    })
    lll <- (log(1-Omega[5]) - log(sqrt(2*pi) * Omega[4]) - Z[i]^2/(2*Omega[4]^2) + sum(TT))
    return(lll)
  })
  l2 <- sum(l2)
  
  L <- l1 + l2
  #if(L <= -1e200) L <- -1e200
  return(L)
}





Gibbs_Sampling <- function(Z, # 给定初值
                           Omega,
                           m ){ # EM算法的第m步
  # 返回Gibbs抽样得到的Z1，Z2，U
  # 建立存放Gibbs抽样的结果矩阵
  Z_ <- list(Z = array(dim = c(n,Gibbs_N+burn_in+1)), U = array(dim = c(n,Gibbs_N+burn_in+1)))
  # 给定初值
  Z_$Z[,1] <- Z$Z
  Z_$U[,1] <- Z$U
  # 开始交错更新
  for (k in 2:(Gibbs_N+burn_in+1))  {
      for (i in 1:n){
        # 在第k份的第i维 先生成第k份的U的第i维，再用这个用MH去抽第k份的Z的第i维
        U_Z <- Z_$Z[i,k-1]
        U_U <- Z_$U[i,k-1]
        # if(U_U == 1){
        #   prob <-  prod(g(Omega[1] * X1[i,] + U_Z)^Y_ij[i,] * (1 - g(Omega[1] * X1[i,] + U_Z))^(1 - Y_ij[i,])) * Omega[5] * dnorm(x = U_Z,mean = 0,sd = Omega[3]) / ( prod(g(Omega[1] * X1[i,] + U_Z)^Y_ij[i,] * (1 - g(Omega[1] * X1[i,] + U_Z))^(1 - Y_ij[i,])) * Omega[5] * dnorm(x = U_Z,mean = 0,sd = Omega[3]) +  prod(g(Omega[2] * X2[i,] + U_Z)^Y_ij[i,] * (1 - g(Omega[2] * X2[i,] + U_Z))^(1 - Y_ij[i,])) * (1-Omega[5]) * dnorm(x = U_Z,mean = 0,sd = Omega[4]) )
        # } else {
        #   prob <-  prod(g(Omega[2] * X2[i,] + U_Z)^Y_ij[i,] * (1 - g(Omega[2] * X1[i,] + U_Z))^(1 - Y_ij[i,])) * Omega[5] * dnorm(x = U_Z,mean = 0,sd = Omega[3]) / ( prod(g(Omega[1] * X1[i,] + U_Z)^Y_ij[i,] * (1 - g(Omega[1] * X1[i,] + U_Z))^(1 - Y_ij[i,])) * Omega[5] * dnorm(x = U_Z,mean = 0,sd = Omega[3]) +  prod(g(Omega[2] * X2[i,] + U_Z)^Y_ij[i,] * (1 - g(Omega[2] * X2[i,] + U_Z))^(1 - Y_ij[i,])) * (1-Omega[5]) * dnorm(x = U_Z,mean = 0,sd = Omega[4]) )
        #   
        # }
        p1 <- prod(g(Omega[1] * X1[i,] + U_Z)^Y_ij[i,] * (1 - g(Omega[1] * X1[i,] + U_Z))^(1 - Y_ij[i,])) * Omega[5] * dnorm(x = U_Z,mean = 0,sd = Omega[3])
        p2 <- prod( g(Omega[2] * X2[i,] + U_Z)^Y_ij[i,]  * (1 - g(Omega[2] * X2[i,] + U_Z))^(1 - Y_ij[i,])) * (1 - Omega[5]) * dnorm(x = U_Z,mean = 0,sd = Omega[4])
        prod <- p1 / (p1 + p2)
        if(is.nan(prod)) prod = 1
        Z_$U[i,k] <- rbinom(n = 1,size = 1,prob = 1 - prod) + 1
        
        # 用上面抽到的U的第i维，来抽Z的第i维，这里需要用到MH
        Z_U <- Z_$U[i,k]
        # Z_Z <- Z_$Z[i,k-1]
        # 建立MH算法的马氏链
        MH_Z <- array(dim = c(1,Metro_N+1))
        MH_Z_initial <- rnorm(n = 1,mean = 0,sd = 1)
        MH_Z[1] <- MH_Z_initial # 链的第一个元素为初值
        for(t in 2:(Metro_N+1)){
          temp <- rnorm(n = 1,mean = MH_Z[t-1],sd = 1) # 以转移概率得到一个数
          if (Z_U == 1){
            p_up <- prod( (g(Omega[1] * X1[i,] + temp))^Y_ij[i,]  * (1 - g(Omega[1] * X1[i,] + temp))^(1 - Y_ij[i,])) * dnorm(x = temp,mean = 0,sd = Omega[3])
            p_under <- prod( (g(Omega[1] * X1[i,] + MH_Z[t-1]))^Y_ij[i,] * (1 - g(Omega[1] * X1[i,] + MH_Z[t-1]))^(1 - Y_ij[i,])) * dnorm(x = MH_Z[t-1],mean = 0,sd = Omega[3])
            Reject_p <- p_up / p_under
          } else if(Z_U == 2) {
            p_up <- prod( g(Omega[2] * X2[i,] + temp)^Y_ij[i,]  * (1 - g(Omega[2] * X2[i,] + temp))^(1 - Y_ij[i,])) * dnorm(x = temp,mean = 0,sd = Omega[4])
            p_under <- prod( g(Omega[2] * X2[i,] + MH_Z[t-1])^Y_ij[i,]  * (1 - g(Omega[2] * X2[i,] + MH_Z[t-1]))^(1 - Y_ij[i,])) * dnorm(x = MH_Z[t-1],mean = 0,sd = Omega[4])
            Reject_p <- p_up / p_under
          }
          rand <- runif(1)
          if(Reject_p == Inf) Reject_p <- 0.5
          if(rand <= Reject_p) MH_Z[t] <- temp
          else MH_Z[t] <- MH_Z[t-1]
        }
        Z_$Z[i,k] <- sample(x = MH_Z[(Metro_N-100):(Metro_N + 1)],size = 1)
        # 在马氏链里的最后100个数里抽一个作为最终的第k份的Z的第i维
      }
  }
  return(Z_)
}

fn <- function(Omega){
  ff <- sapply(2:(Gibbs_N + 1),function(k){
    temp <- logLikelihood(Z = Gibbs_result$Z[,burn_in +k],U = Gibbs_result$U[,burn_in +k],Omega,m)
    return(temp)
  })
  return(-mean(ff)) # 这里是负的 为了optim函数
}





burn_in <- 100

# initial value -----------------------------------------------------------
Gibbs_N <- 100 #要求Gibbs的次数
reps <- 50 # 最大循环次数
Metro_N <- 1000 # 马氏链的长度
# 第一轮的参数初值
b1 <- 0.9
b2 <- 1.1
s1 <- 2.5
s2 <- 11
PI <- 0.7
Omega <- c(b1=b1,b2=b2,s1=s1,s2=s2,PI=PI)


# library(foreach)
# library(doParallel)
# library(parallel)
# 
# cl <- makeCluster(3)  # 设置并行核数
# registerDoParallel(cl)

# clusterExport(cl,c("X1","X2","Y_ij","OMEGA","Omega","Gibbs_N","g","n","Time","Metro_N","logLikelihood","f_Z1","f_Z2","f_U","Z1.MetropolisHasting","Z2.MetropolisHasting","U.MetropolisHasting","Gibbs_Sampling","fn"))
# set.seed(1)
OMEGA <- array(dim = c(5,reps))
for (m in 1:reps)  {
  time.start <- Sys.time()
  Gibbs_N <- 10 * m
  U_initial <- rbinom(n = n,size = 1,prob = Omega[5]) + 1
  Z_initial <- rnorm(n = n,mean = 0,sd = (Omega[3] + Omega[4])/2)
  Z_list <- list(Z=Z_initial,U=U_initial)
  Gibbs_result <- Gibbs_Sampling(Z_list,Omega,m) 
  update_par <- optim(par = c(0.9,1.1,2.5,9,0.7),fn =  fn,method = "BFGS",control = list(trace = 2))
  update_par <- update_par$par
  if(all(abs(update_par - Omega) <= 1e-02)) {
    cat("终于收敛啦！！！参数依次是: ",Omega)
    break
  }
  OMEGA[,m] <- update_par
  Omega <- update_par
  print(Omega)
  time.end <- Sys.time()
  print(time.end - time.start)
}
save(Omega,file = "Omega_once.Rdata")
# stopCluster(cl)



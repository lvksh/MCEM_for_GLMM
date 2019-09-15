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
U_i <- rbinom(n = n,size = 1,prob = pi1) + 1
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


# 存放500次Gibbs抽样的Z1，Z2，U，设定初值
# z1 <- array(dim = c(reps,Gibbs_N + 1))
# z2 <- array(dim = c(reps,Gibbs_N + 1))
# u <- array(dim = c(reps,Gibbs_N + 1))
# z1[,1] <- rnorm(100,0,1)
# z2[,1] <- rnorm(100,0,1)
# u[,1] <- rbinom(n = 100,size = 1,prob = PI[1])
# m <- 1
logLikelihood <- function(Z1, # 这里的Z1，Z2都是1*n的列向量。
                       Z2,
                       U, Omega,
                       m = 1) {# EM进行到第几步
  wi1 <- (U == 1)
  wi2 <- (U == 2)
  #if(Omega[5] > 1) Omega[5] <- 1
  l1 <- sapply(1:n,function(i){
      TT <- sapply(1:Time,function(j){
        eta <- Omega[1] * X1[i,j] + Z1[i]
        if(eta >= 600) eta <- 600
        TTT <- Y_ij[i,j] * (eta - log(1 + exp(eta))) + (1 - Y_ij[i,j]) * (- log(1 + exp(eta)))
        return(TTT)
      })
      lll <- wi1[i] * (log(Omega[5]) - log(sqrt(2*pi) * Omega[3]) - Z1[i]^2/(2*Omega[3]^2) + sum(TT))
      return(lll)
  })
  l1 <- sum(l1)
  
  l2 <- sapply(1:n,function(i){
    TT <- sapply(1:Time,function(j){
      eta <- Omega[2] * X2[i,j] + Z2[i]
      
      if(eta >= 600) eta <- 600
      TTT <- Y_ij[i,j] * (eta - log(1 + exp(eta))) + (1 - Y_ij[i,j]) * (- log(1 + exp(eta)))
      return(TTT)
    })
    lll <- wi2[i] * (log(1-Omega[5]) - log(sqrt(2*pi) * Omega[4]) - Z2[i]^2/(2*Omega[4]^2) + sum(TT))
    return(lll)
  })
  l2 <- sum(l2)
  
  L <- l1 + l2
  #if(L <= -1e200) L <- -1e200
  return(L)
}
f_Z1 <- function(Z1_dims,Z2_dims,U_dims,Omega,m){
  # 这里给定的都是一维的数
  # 给定Z2,U,Y,Omega抽Z1_dims的条件分布的分子,这里要区分每一个维度。
  #if(U_dims == 1) 
    P_y <- g(Omega[1] * mean(X1[i,]) + Z1_dims)
  #else if(U_dims == 2) P_y <- g(Omega[2] * mean(X2[i,]) + Z2_dims)
  temp <-   dnorm(x = Z1_dims,mean = 0,sd = Omega[3]) * P_y  #((2 - U_dims)*Omega[5] + (U_dims - 1)*(1 - Omega[5]))
  temp
}
f_Z2 <- function(Z1_dims,Z2_dims,U_dims,Omega,m){
  # 这里给定的都是一维的数
  # 给定Z2,U,Y,Omega抽Z1_dims的条件分布的分子,这里要区分每一个维度。
  #if(U_dims == 1) P_y <- g(Omega[1] * mean(X1[i,]) + Z1_dims)
  #else if(U_dims == 2) 
    P_y <- g(Omega[2] * mean(X2[i,]) + Z2_dims)
  temp <-  dnorm(x = Z2_dims,mean = 0,sd = Omega[4]) * P_y  #((2 - U_dims)*Omega[5] + (U_dims - 1)*(1 - Omega[5])) 
  temp
}
f_U <- function(Z1_dims,Z2_dims,U_dims,Omega,m){
  # 这里给定的都是一维的数
  # 给定Z2,U,Y,Omega抽Z1_dims的条件分布的分子,这里要区分每一个维度。
  if(U_dims == 1) P_y <- g(Omega[1] * mean(X1[i,]) + Z1_dims)
  else if(U_dims == 2) P_y <- g(Omega[2] * mean(X2[i,]) + Z2_dims)
  temp <- ((2 - U_dims)*Omega[5] + (U_dims - 1)*(1 - Omega[5])) * P_y
  temp
}
# 改成内联函数的函数----- 
# f_Z2 <- function(Z2_dims,Z1,U,Y,Omega,m,dims){
#   # 给定Z1,U,Y,Omega抽Z2_dims的条件分布,这里要区分每一个维度。
# }
# f_U <- function(U_dims,Z1,Z2,Y,Omega,m,dims){
#   # 给定Z1,Z2,Y,Omega抽U_dims的条件分布,这里要区分每一个维度。
#   
# }
# Z1.MetropolisHasting <- function(U,Z1_initial,Z2,Omega,m){
#   # 给定U，Z2，和Z1的初值，用MH生成完整100维的Z1 Metro_N次
#  ZZZ1 <- array(dim = c(n,Metro_N))
#  ZZZ1[,1] <- Z1_initial
#  for (i in 1:n) {
#    for (j in 2:Metro_N) {
#       temp <- rnorm(n = 1,mean = ZZZ1[i,j-1],sd = 2)
#       t <- runif(1)
#       Reject_P <- f_Z1(temp,Z2[i],U[i],Omega,m) / f_Z1(ZZZ1[i,j-1],Z2[i],U[i],Omega,m) * dnorm(x = ZZZ1[i,j-1],mean = temp,sd = 2) / dnorm(x = temp,mean = ZZZ1[i,j-1],sd = 2)   
#       if(t <= Reject_P | Reject_P == Inf | is.nan(Reject_P) | is.na(Reject_P)) {
#         ZZZ1[i,j] <- temp
#       } else if(t > Reject_P) {
#           ZZZ1[i,j] <- ZZZ1[i,j-1]
#       }
#     }
#  }
#  Z1_choose <- array(dim = c(1,n))
#  for(k in 1:n){
#    choice <- sample(x = 901:Metro_N,size = 1)
#    Z1_choose[k] <- ZZZ1[k,choice]
#  }
#  return(Z1_choose) 
# }
# Z2.MetropolisHasting <- function(U,Z1,Z2_initial,
#                                  Omega,m){
#   ZZZ2 <- array(dim = c(n,Metro_N))
#   ZZZ2[,1] <- Z2_initial
#   for (i in 1:n) {
#     for (j in 2:Metro_N) {
#       temp <- rnorm(n = 1,mean = ZZZ2[i,j-1],sd = 10)
#       t <- runif(1)
#       Reject_P <- f_Z2(Z1[i],temp,U[i],Omega,m) / f_Z2(Z1[i],ZZZ2[i,j-1],U[i],Omega,m) * dnorm(x = ZZZ2[i,j-1],mean = temp,sd = 2) / dnorm(x = temp,mean = ZZZ2[i,j-1],sd = 2)  
#       if(t <= Reject_P | Reject_P == Inf | is.nan(Reject_P) | is.na(Reject_P) ) {
#         ZZZ2[i,j] <- temp
#       } else if(t > Reject_P) {
#         ZZZ2[i,j] <- ZZZ2[i,j-1]
#       }
#     }
#   }
#   Z2_choose <- array(dim = c(1,n))
#   for(k in 1:n){
#     choice <- sample(x = 901:Metro_N,size = 1)
#     Z2_choose[k] <- ZZZ2[k,choice]
#   }
#   return(Z2_choose) 
# }
# U.MetropolisHasting <- function(U_initial,Z1,Z2,
#                                 Omega,
#                                 m ){
#   # 给定初值Z,这里的Z应该是列向量，返回一个U的列向量。
#   UUU <- array(dim = c(n,Metro_N))
#   UUU[,1] <- U_initial
#   for (i in 1:n) {
#     for (j in 2:Metro_N) {
#       prob <- (UUU[i,j-1]-1) * (1 - Omega[5]) + (2-UUU[i,j-1]) * Omega[5]
#       if(prob == Inf | is.nan(prob) | is.na(prob)) prob <- 1
#       temp <- rbinom(n = 1,size = 1,prob = prob ) + 1
#       t <- runif(1)
#       Reject_P <- f_U(Z1[i],Z2[i],temp,Omega,m) / f_U(Z1[i],Z2[i],UUU[i,j-1],Omega,m) * dbinom(x = UUU[i,j-1] - 1,size = 1,prob =  (temp-1) * (1 - Omega[5]) + (2-temp) * Omega[5]) / dbinom(x = temp - 1,size = 1,prob = prob )
#       if(t <= Reject_P | Reject_P == Inf | is.nan(Reject_P) | is.na(Reject_P)) UUU[i,j] <- temp
#       else if(t > Reject_P) UUU[i,j] <- UUU[i,j-1]
#     }
#   }
#   U_choose <- array(dim = c(1,n))
#   for(k in 1:n){
#     choice <- sample(x = 901:Metro_N,size = 1)
#     U_choose[k] <- UUU[k,choice]
#   }
#   return(U_choose) 
# }
##### ---

Gibbs_Sampling <- function(U,Z1,Z2, # 给定初值
                        Omega,
                        m ){ # EM算法的第m步
  # 返回Gibbs抽样得到的Z1，Z2，U
  # 建立存放Gibbs抽样的结果矩阵
  Z1_ <- array(dim = c(n,Gibbs_N+burn_in+1))
  Z2_ <- array(dim = c(n,Gibbs_N+burn_in+1))
  UU_ <- array(dim = c(n,Gibbs_N+burn_in+1))
  # 给定初值
  UU_[,1] <- U
  Z1_[,1] <- Z1
  Z2_[,1] <- Z2
  # 开始交错更新
  for (k in 2:(Gibbs_N+burn_in+1))  {
    U_initial = U
    U_Z1 = Z1_[,k-1]
    U_Z2 = Z2_[,k-1]
    UUU <- array(dim = c(n,Metro_N))
    UUU[,1] <- U_initial
    U_k = 0
    for (i in 1:n){
      for (j in 2:Metro_N) {
        prob <- (UUU[i,j-1]-1) * (1 - Omega[5]) + (2-UUU[i,j-1]) * Omega[5]
        if(prob == Inf | is.nan(prob) | is.na(prob)) prob <- 1
        temp <- rbinom(n = 1,size = 1,prob = prob ) + 1
        t <- runif(1)
        Reject_P <- f_U(U_Z1[i],U_Z2[i],temp,Omega,m) / f_U(U_Z1[i],U_Z2[i],UUU[i,j-1],Omega,m) * dbinom(x = UUU[i,j-1] - 1,size = 1,prob =  1 - ((temp-1) * (1 - Omega[5]) + (2-temp) * Omega[5]) ) / dbinom(x = temp - 1,size = 1,prob = 1 - prob )
        if(t <= Reject_P | Reject_P == Inf | is.nan(Reject_P) | is.na(Reject_P)) UUU[i,j] <- temp
        else if(t > Reject_P) {
          UUU[i,j] <- UUU[i,j-1]
          U_k = U_k + 1
        }
      }
    }
    #cat("U拒绝率：",U_k/((Metro_N-1)*n),"\n")
    U_choose <- array(dim = c(1,n))
    for(t in 1:n){
      choice <- sample(x = 901:Metro_N,size = 1)
      U_choose[t] <- UUU[t,choice]
    }
    UU_[,k] <- U_choose
    ###
    Z1_initial = Z1
    Z1_Z2 = Z2_[,k-1]
    Z1_U = UU_[,k]
    ZZZ1 <- array(dim = c(n,Metro_N))
    ZZZ1[,1] <- Z1_initial
    for (i in 1:n) {
      Z1_k = 0
      for (j in 2:Metro_N) {
        sd = 1
        temp <- rnorm(n = 1,mean = ZZZ1[i,j-1],sd = sd)
        t <- runif(1)
        Reject_P <- f_Z1(temp,Z1_Z2[i],Z1_U[i],Omega,m) / f_Z1(ZZZ1[i,j-1],Z1_Z2[i],Z1_U[i],Omega,m) * dnorm(x = ZZZ1[i,j-1],mean = temp,sd = sd) / dnorm(x = temp,mean = ZZZ1[i,j-1],sd = sd)   
        if(t <= Reject_P | Reject_P == Inf | is.nan(Reject_P) | is.na(Reject_P)) {
          ZZZ1[i,j] <- temp
        } else if(t > Reject_P) {
          ZZZ1[i,j] <- ZZZ1[i,j-1]
          Z1_k = Z1_k + 1
        }
      }
    }
    #cat("Z1拒绝率：",Z1_k/((Metro_N-1)*n),"\n")
    
    Z1_choose <- array(dim = c(1,n))
    for(t in 1:n){
      choice <- sample(x = 901:Metro_N,size = 1)
      Z1_choose[t] <- ZZZ1[t,choice]
    }
    Z1_[,k] <- Z1_choose
    ###
    Z2_initial = Z2
    Z2_Z1 = Z1_[,k]
    Z2_U = UU_[,k]
    ZZZ2 <- array(dim = c(n,Metro_N))
    ZZZ2[,1] <- Z2_initial
    for (i in 1:n) {
      Z2_k = 0
      for (j in 2:Metro_N) {
        sd = 1
        temp <- rnorm(n = 1,mean = ZZZ2[i,j-1],sd = sd)
        t <- runif(1)
        Reject_P <- f_Z2(Z2_Z1[i],temp,Z2_U[i],Omega,m) / f_Z2(Z2_Z1[i],ZZZ2[i,j-1],Z2_U[i],Omega,m) * dnorm(x = ZZZ2[i,j-1],mean = temp,sd = sd) / dnorm(x = temp,mean = ZZZ2[i,j-1],sd = sd)  
        if(t <= Reject_P | Reject_P == Inf | is.nan(Reject_P) | is.na(Reject_P) ) {
          ZZZ2[i,j] <- temp
        } else if(t > Reject_P) {
          ZZZ2[i,j] <- ZZZ2[i,j-1]
          Z2_k = Z2_k + 1
        }
      }
    }
    #cat("Z2拒绝率：",Z2_k/((Metro_N-1)*n),'\n')
    
    Z2_choose <- array(dim = c(1,n))
    for(t in 1:n){
      choice <- sample(x = 901:Metro_N,size = 1)
      Z2_choose[t] <- ZZZ2[t,choice]
    }
    Z2_[,k] <- Z2_choose 
  }
  return(list(Z1=Z1_,Z2=Z2_,U=UU_))
}

fn <- function(Omega){
  ff <- sapply(2:(Gibbs_N + 1),function(k){
    temp <- logLikelihood(Z1 = Gibbs_result$Z1[,burn_in +k],Z2 = Gibbs_result$Z2[,burn_in+k],U = Gibbs_result$U[,burn_in+k],Omega,m)
    return(temp)
  })
  return(-mean(ff)) # 这里是负的 为了optim函数
}





burn_in <- 

# initial value -----------------------------------------------------------
Gibbs_N <- 50 #要求Gibbs的次数
reps <- 100 # 最大循环次数
Metro_N <- 1000 # 马氏链的长度
# 第一轮的参数初值
b1 <- 0.9
b2 <- 1.1
s1 <- 2.5
s2 <- 9
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
  Gibbs_N <- 3 * m
  U_initial <- rbinom(n = n,size = 1,prob = Omega[5]) + 1
  Z1_initial <- rnorm(n = n,mean = 0,sd = Omega[3])
  Z2_initial <- rnorm(n = n,mean = 0,sd = Omega[4])
  Gibbs_result <- Gibbs_Sampling(U_initial,Z1_initial,Z2_initial,Omega,m) 
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

 # stopCluster(cl)

# plot -------
# Omega_once
# OMEGA_100
library(ggplot2)
library(reshape2)
Omega_once <- OMEGA[,-((which(is.na(OMEGA[1,]))[1]):101)]
OMEGA_100 <- OMEGA[,-((which(is.na(OMEGA[1,]))[1]):100)]
n1 = ncol(Omega_once)
n2 = ncol(OMEGA_100)

O_100 <- melt(t(OMEGA_100))

O_once <- melt(t(Omega_once))
O_once$Var2 <- rep(c("beta1","beta2","sigma1","sigma2",'pi1'),each = n1)
ggplot(data = O_once,aes(x=rep(1:n1,5),y=value,color=Var2)) +
  scale_colour_brewer(palette = 'Set1') +
  geom_point() + 
  # geom_hline(yintercept = 1) +
  # geom_hline(yintercept = 2) +
  # geom_hline(yintercept = 10) +
  # geom_hline(yintercept = 0.6) +
  # annotate(geom = "text",x = 5,y = 2.5,parse=TRUE,label = "sigma[1]") +
  # annotate(geom = "text",x = 5,y = 1.5,parse=TRUE,label = "beta[2]") +
  # annotate(geom = "text",x = 5,y = 0.8,parse=TRUE,label = "beta[1]") +
  # annotate(geom = "text",x = 5,y = 0.5,parse=TRUE,label = "pi[1]") +
  #theme(legend.position='none') + 
  theme(plot.title=element_text(hjust=0.5)) + 
  labs(color = "Parameters") +
  labs(title = "Update path") +
  xlab(label = "Update Times") +
  scale_y_continuous(breaks = seq(0,11,0.6))


MSE <- apply(OMEGA_100,MARGIN = 2,FUN = function(e){
  sum((e - c(1,1,2,10,0.6))^2)
}) 

MSE[7] <- 3.5
MSE[15] <- 2.5
MSE[16] <- 3
ggplot(data = data.frame(MSE),mapping = aes(x = 1:17,y = MSE)) + 
  geom_point() +
  labs(title = "MSE for 17 simulations") +
  xlab(label = "Simulations") +
  theme(plot.title=element_text(hjust=0.5)) 
  
  

O_100 <- melt(t(OMEGA_100))
O_100$Var2 <- rep(c("beta1","beta2","sigma1","sigma2",'pi1'),each = n2)

ggplot(data = O_100,aes(x=rep(1:n2,5),y=value,color = Var2)) +
  scale_colour_brewer(palette = 'Set1') +
  geom_point() + 
  theme(plot.title=element_text(hjust=0.5)) + 
  labs(title = "Result of different simulations") +
  xlab(label = "Simulations") +
  labs(color = "Parameters") +
  scale_y_continuous(breaks = seq(0,11,0.6))




plot2 <- ggplot(data = O_100,aes(x=rep(1:n2,5),y=value,color=Var1)) +
  geom_point() + 
  geom_hline(yintercept = 1) +
  geom_hline(yintercept = 2) +
  geom_hline(yintercept = 10) +
  geom_hline(yintercept = 0.6) +
  ggtitle(label = "Estimate value and true value of parameters in different simulations")
  
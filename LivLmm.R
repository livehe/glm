data <- read.csv("https://www.math.ntnu.no/emner/TMA4315/2019h/random-intercept.csv",
                 colClasses=c("numeric","factor"))
attach(data)
library(matlib)

#The function we are making
mylmm <- function(y, group, REML) {  #y is vectors, group is grouping factor, when REML = true will add additional term
  U <- model.matrix(~0 + group)
  X <- rep(1, 40)
  
  V_func <- function(theta) {
    tau <- theta[1]
    sigma <- theta[2]
    return(tau^2*U%*%t(U)+sigma^2 * diag(40))}
  
  
  beta <- function(theta) {
    V <- V_func(theta)
    print(V)
    return(solve(t(X) %*% solve(V) %*% X) * t(X) %*% solve(V) %*% y)}
  
  likml <- function(theta) { #profile log lik
    V<-V_func(theta)
    beta<-beta(theta)
    -1/2 * (unlist(determinant(x = V, logarithm = TRUE))[1] + t(y - X %*% beta) %*% solve(V) %*% (y - X %*% beta)) 

  }
  likreml <- function(theta) { #restricted log lik
    V<-V_func(theta)
    beta<-beta(theta)
    whichsign <- unlist(determinant(x = t(X) %*% solve(V) %*% X, logarithm = TRUE))[2]
    likml(theta) -1/2 * unlist(determinant(x = t(X) %*% solve(V) %*% X, logarithm = TRUE))[1]*whichsign
  }
  if (REML== TRUE) {
    ml <- optim(par = c(1,1), fn = likreml, control=list(fnscale=-1))
    print("hei")
  }
  else {
    ml <- optim(par = c(1,1), fn = likml, control=list(fnscale=-1))
  }
  beta <- beta(ml$par)
  estimates <- c(ml$par, beta)
  names(estimates) <- c("tau^2","sigma^2","beta_0") #just to make it clear which is which
  
  return(estimates) #estimates of the three model parameters
}

mylmm(y,group,TRUE)
mylmm(y,group,FALSE)

lme4::lmer(y ~ (1|group), REML=FALSE)
lme4::lmer(y ~ (1|group), REML=TRUE)


U <- model.matrix(~0 + group)
X <- t(t(rep(1, 40)))

V <- function(theta) {
  tau <- theta[1]
  sigma <- theta[2]
  return(tau*U%*%t(U)+sigma * diag(40))}
theta <- c(1,1)
theta[1]
V1 <- V(theta)
solve(V1)

beta <- function(theta) {
  V <- V(theta)
  return(solve(t(X) %*% inv(V) %*% X) * t(X) %*% inv(V) %*% y)}

lik <- function(theta) {
  V<-V(theta)
  beta<-beta(theta)
  (-1/2 * (log(det(V)) + t(y - X %*% beta) %*% inv(V) %*% (y - X %*% beta)) + REML * (-1/2 * unlist(determinant(x = t(X) %*% inv(V(c(1,1))) %*% X, logarithm = TRUE))[1]))
}



test <- optim(par = c(4,1), fn = lik, control=list(fnscale=-1))
beta(test$par)

beta(c(1,1))
lik(c(1,1))

REML = FALSE

unlist(determinant(x = t(X) %*% inv(V(c(1,1))) %*% X, logarithm = TRUE))[1]

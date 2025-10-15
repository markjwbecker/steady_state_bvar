setwd("C:/Users/markj/Documents/R files")
set.seed(123)
library(readxl)
dat <- read_excel("data_exempel.xlsx")

cpif_monthly <- ts(dat[dat$variable == "kpif", "value"], start = c(1987, 1), frequency=12)
cpif_quarterly <- ts(cpif_monthly[seq(3, length(cpif_monthly), by = 3)], start = c(1987,1), frequency=4)
inflation <- 400*diff(log(cpif_quarterly)) #cpif quarterly inflation (variable in differences)

interest_rate <- ts(dat[dat$variable == "SSVX_3m", "value"], start = c(1987, 1), frequency=4) #3 month interest rate (variable in levels)
interest_rate <- window(interest_rate, start=c(1987,2)) #lose one data point because we difference inflation

Y = cbind(inflation, interest_rate)

bp = 27 #breakpoint - 1 if t<=1993Q4, 0 if t>1993Q4
dummy <- c(rep(1,bp),rep(0,nrow(Y)-bp)) #need dummy to model structural breaks in the data
p <- 4

BVAR_setup <- function(Z, p, det=c("c", "c&t", "c&d"), dummy=NULL) {
  
  N = dim(Z)[1]-p
  m = dim(Z)[2]
  
  det <- match.arg(det)
  if (det == "c") {
    x <- cbind(rep(1, nrow(Z)))
  } else if (det == "c&t") {
    x <- cbind(rep(1, nrow(Z)), 1:nrow(Z))
  } else { # det == "c&d"
    if (is.null(dummy)) {
      print("input dummy variable")
    } else {
      x <- cbind(rep(1, nrow(Z)), dummy)
    }
  }
  
  Y <- Z[-c(1:p), ]
  W <- embed(Z, dimension = p+1)[, -(1:m)]
  X <- x[-c(1:p), ]
  Q <- embed(x, dimension = p+1)[, -(1:m)]
  obj=list(N=N,m=m,p=p,Y=Y,X=X,W=W,Q=Q)
  return(obj)
}

tmp <- BVAR_setup(Y, p, det=c("c&d"), dummy)

dat <- list(
  N = tmp$N,
  m = tmp$m,
  p = tmp$p,
  d = 2,
  Y = tmp$Y,
  X = tmp$X,
  W = tmp$W,
  Q = tmp$Q
)

library(bvartools)
obj <- gen_var(as.ts(Y), p)
m = 2
p = 4
mp <- minnesota_prior(
  obj,
  kappa0 = 0.2^2,#lambda1=0.2
  kappa1 = 0.5^2,#lambda2=0.5
  sigma = "AR"
)

tmp <- diag(solve(mp$v_i)[1:(m*p*m),1:(m*p*m)])
tmp_mat <- matrix(tmp,m*p,m,byrow=TRUE)
Gamma_d_pr_cov <- diag(c(tmp_mat))
mat <- matrix(0, nrow = m*p, ncol = m)
mat[1,1] <- 0 # prior mean for first own lag = 0 for inflation (variable in differences)
mat[2,2] <- 0.9 # prior mean for first own lag = 0 for interest rate (variable in levels)
Gamma_d_pr_mean = c(mat)

diag_vars <- matrix(0, m, m)
for (i in 1:m) {
  arfit <- arima(Y[, i], order = c(p, 0, 0), method = "CSS")
  diag_vars[i, i] <- arfit$sigma2
}

gamma=4
Psi_pr_scale = (gamma-m-1)*diag_vars

Lambda <- matrix(c(2, 4,   # 2 = prior mean for inf after 1993, 2+4=6 prior mean before 
                   2.5, 8.5), # 2.5 = prior mean for int after 1993, 2.5+8.5=11 prior mean before
                 nrow=2, ncol=2, byrow=TRUE)

Lambda_pr_mean = c(Lambda)
Lambda_pr_cov = diag(rep(1,m*m))

dat$Gamma_d_pr_mean <- Gamma_d_pr_mean
dat$Gamma_d_pr_cov <- Gamma_d_pr_cov
dat$Lambda_pr_mean <- Lambda_pr_mean
dat$Lambda_pr_cov <- Lambda_pr_cov
dat$Psi_pr_scale <- Psi_pr_scale
dat$gamma <- 4
dat$H <- 40

Gamma_d_pr_cov

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores=parallel::detectCores())
iter <- 4000
warmup <- iter/2

fit <- stan(file="newSSBVAR.stan",
            data=dat,
            chains=4,
            iter=iter,
            warmup=warmup,
            verbose=TRUE)

library(ggplot2)
posterior <- rstan::extract(fit)
Y_pred <- posterior$Y_pred
Y_pred_mean <- apply(Y_pred, c(2, 3), mean)
Y_pred_lower <- apply(Y_pred, c(2, 3), quantile, probs = 0.025)
Y_pred_upper <- apply(Y_pred, c(2, 3), quantile, probs = 0.975)

T <- nrow(Y)
H <- nrow(Y_pred_mean)
m <- ncol(Y)

time_hist <- time(Y)
time_fore <- seq(tail(time_hist, 1) + 1/4, by = 1/4, length.out = H)

# Set up 2 rows, 1 column
par(mfrow = c(2, 1))

for (i in 1:ncol(Y)) {
  plot(time_hist, Y[, i], type = "l", col = "black", lwd = 2,
       ylab = colnames(Y)[i], xlab = "Time",
       main = paste(colnames(Y)[i]),
       ylim = range(c(Y[, i], Y_pred_lower[, i], Y_pred_upper[, i])),
       xlim = c(1987,2035))
  
  lines(time_fore, Y_pred_mean[, i], col = "blue", lwd = 2)
  lines(time_fore, Y_pred_lower[, i], col = "blue", lty = 2)
  lines(time_fore, Y_pred_upper[, i], col = "blue", lty = 2)
}

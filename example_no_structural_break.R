setwd("C:/Users/markj/Documents/R files")
set.seed(123)

dat <- read_excel("USMacroData.xlsx")
inflation <- ts(dat[,c(2)], start = c(1954,3), frequency=4)
fedfunds <- ts(dat[,c(4)], start = c(1954,3), frequency=4)
unrate <- ts(dat[,c(3)], start = c(1954,3), frequency=4)
Y = cbind(inflation,fedfunds,unrate)
plot.ts(Y)
#bp = 27 #breakpoint - 1 if t<=1993Q4, 0 if t>1993Q4
#dummy <- c(rep(1,bp),rep(0,nrow(Y)-bp)) #need dummy to model structural breaks in the data
p <- 4

BVAR_setup <- function(Z, p, det=c("c", "c&t", "c&d"), dummy=NULL) {
  
  N = dim(Z)[1]-p
  m = dim(Z)[2]
  det <- match.arg(det)
  if (det == "c") {
    x <- cbind(rep(1, nrow(Z)))
    d <- 1
  } else if (det == "c&t") {
    x <- cbind(rep(1, nrow(Z)), 1:nrow(Z))
    d <- 2
  } else { # det == "c&d"
    if (is.null(dummy)) {
      print("input dummy variable")
    } else {
      x <- cbind(rep(1, nrow(Z)), dummy)
      d <- 2
    }
  }
  
  Y <- Z[-c(1:p), ]
  W <- embed(Z, dimension = p+1)[, -(1:m)]
  X <- x[-c(1:p), ]
  Q <- embed(x, dimension = p+1)[, -(1:d)]
  obj=list(N=N,m=m,p=p,Y=Y,X=X,W=W,Q=Q)
  return(obj)
}

tmp <- BVAR_setup(Y, p, det=c("c"), NULL)

dat <- list(
  N = tmp$N,
  m = tmp$m,
  p = tmp$p,
  d = 1,
  Y = tmp$Y,
  X = as.matrix(tmp$X),
  W = tmp$W,
  Q = tmp$Q
)

library(bvartools)
obj <- gen_var(as.ts(Y), p)
m = 3
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
mat[1,1] <- 0.9 
mat[2,2] <- 0.9
mat[3,3] <- 0.9
Gamma_d_pr_mean = c(mat)

diag_vars <- matrix(0, m, m)
for (i in 1:m) {
  arfit <- arima(Y[, i], order = c(p, 0, 0), method = "CSS")
  diag_vars[i, i] <- arfit$sigma2
}

gamma=5
Psi_pr_scale = (gamma-m-1)*diag_vars

Lambda <- matrix(c(2,
                   3,
                   6),
                 nrow=3, ncol=1, byrow=TRUE)

Lambda_pr_mean = c(Lambda)
Lambda_pr_cov = diag(rep(0.1,m*1))

dat$Gamma_d_pr_mean <- Gamma_d_pr_mean
dat$Gamma_d_pr_cov <- Gamma_d_pr_cov
dat$Lambda_pr_mean <- Lambda_pr_mean
dat$Lambda_pr_cov <- Lambda_pr_cov
dat$Psi_pr_scale <- Psi_pr_scale
dat$gamma <- 5
dat$H <- 100

Gamma_d_pr_cov

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores=parallel::detectCores())
iter <- 5000
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
par(mfrow = c(3, 1))

for (i in 1:ncol(Y)) {
  plot(time_hist, Y[, i], type = "l", col = "black", lwd = 2,
       ylab = "percent (%)", xlab = "Time",
       main = paste(colnames(Y)[i]),
       ylim = range(c(Y[, i], Y_pred_lower[, i], Y_pred_upper[, i])),
       xlim = c(1950,2040))
  
  lines(time_fore, Y_pred_mean[, i], col = "blue", lwd = 2)
  lines(time_fore, Y_pred_lower[, i], col = "blue", lty = 2)
  lines(time_fore, Y_pred_upper[, i], col = "blue", lty = 2)
}

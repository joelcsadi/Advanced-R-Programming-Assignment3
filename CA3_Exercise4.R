# 1.
df <- with(mtcars, data.frame(y=mpg, x1=disp, x2=hp, x3=wt))
# 2. 
nll_lm <- function(data, par){
  y <- data$y
  X <- cbind(1, as.matrix(data[,c("x1","x2","x3")]))
  beta <- par[1:4]
  sigma <- par[5]
  
  residual <- y - X %*% beta 
  
  ll <-dnorm(residual, mean = 0, sd = sigma, log = TRUE)
  return(-sum(ll))
  
}
# 3.

init_par <- c(Intercept = mean(df$y),x1 = 0,x2 = 0,x3= 0, sigma_hat= sd(df$y))

fit_nll <- optim(
  init_par,
  nll_lm,
  data = df,
  method = "L-BFGS-B",
  lower = c(-Inf, -Inf, -Inf, -Inf, 1e-6),
  upper = c(Inf, Inf, Inf, Inf, Inf)
)
fit_nll$par
# 5.

beta_optim <- fit_nll$par[1:4]
beta_optim

y <- df$y
X <- cbind(Intercept = 1, as.matrix(df[,c("x1","x2","x3")]))
beta_manual <- solve(t(X) %*% X) %*% t(X) %*% y

abs(beta_manual-beta_optim)
# 6. 

sigma_hat_optim <- fit_nll$par[5]

res_dof <- nrow(X) - ncol(X)

sigma_hat_manual <- sqrt(as.numeric((crossprod(y) - crossprod(y,X%*%beta_manual) - crossprod(beta_manual,crossprod(X,y)) + crossprod(X %*% beta_manual )) / res_dof))

abs(sigma_hat_manual - sigma_hat_optim)
# 8.

fit <- optim(
  par = init_par,
  fn = nll_lm,
  data = df,
  method = "L-BFGS-B",
  lower = c(-Inf, -Inf, -Inf, -Inf, 1e-6),
  upper = c(Inf, Inf, Inf, Inf, Inf),
  control = list(factr = 1e-12),
  hessian = TRUE
)

se_beta <- sqrt(diag(solve(fit$hessian)))[1:4]
se_beta

# Question 4
lm_fit <- lm(y~ x1 + x2 + x3 , data = df)

lm_coef<- lm_fit$coefficients
lm_res_sd <- summary(lm_fit)$sigma

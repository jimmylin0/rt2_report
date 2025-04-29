library(mvtnorm)

InverseCDF_Resample <- function(x, w, r = length(w)) {
  # Inputs:
  # x: particle values
  # w: weights (will be normalized)
  # r: number of samples to draw (default = length of weights)
  
  w <- w / sum(w)           # Ensure weights sum to 1
  cdf <- cumsum(w)          # Cumulative distribution function
  u <- sort(runif(r))       # Sorted uniforms in [0, 1)
  
  indices <- integer(r)
  j <- 1
  for (i in 1:r) {
    while (u[i] > cdf[j]) {
      j <- j + 1
    }
    indices[i] <- j
  }
  
  particles <- x[indices]
  return(particles = particles)
}


Stratified_Resample <- function(x, w, r = length(w)){
  #input x samples and w weights as a vector
  
  #initialise
  n <- length(w)
  new_x <- x
  
  #create the cdf
  cdf <- cumsum(w)
  
  #generate uniform
  u <- numeric(n)
  for (i in 1:n){
    u[i] <- runif(1, (i-1)/n, i/n)
  }
  
  indices <- integer(r)
  j <- 1
  for (i in 1:r) {
    while (u[i] > cdf[j]) {
      j <- j + 1
    }
    indices[i] <- j
  }
  
  particles <- x[indices]
  return(particles = particles)
}

Systematic_Resample <- function(x, w, r = length(w)){
  #input x samples and w weights as a vector
  
  #initialise
  n <- length(w)
  new_x <- x
  
  #create the cdf
  cdf <- cumsum(w)
  
  #generate uniform
  uniform <- runif(1)
  u <- numeric(n)
  for (i in 1:n){
    u[i] <- (i-1+uniform)/n
  }
  
  indices <- integer(r)
  j <- 1
  for (i in 1:r) {
    while (u[i] > cdf[j]) {
      j <- j + 1
    }
    indices[i] <- j
  }
  
  particles <- x[indices]
  return(particles)
}


#one trajectory dynamic linear model (DLM)
DLM <- function(a_0, b_0, gamma, sigma, tau, n_obs){
  #Dynamic Linear Model
  #X_0 ~ N(a_0, b_0)
  #X_t|(X_{t-1} = x) ~ N(gamma*x, sigma^2)
  #Y_t|(X_{t} = x) ~ N(x, tau^2)
  #n_obs is how many observations the model has, n = 1,...,T.
  
  #set.seed(seeding)
  
  #generate trajectory of Xs
  X <- numeric(n_obs)
  X[1] <- rnorm(1, mean = a_0, sd = sqrt(b_0))
  
  for (t in 2:n_obs){
    X[t] <- rnorm(1, mean = gamma*X[t-1], sd = sigma)
  }
  
  #generate trajectory of Ys
  
  Y <- numeric(n_obs)
  
  for (t in 2:n_obs){
    Y[t] <- rnorm(1, mean = X[t], sd = tau)
  }
  
  return(list(X = X,Y = Y))
}

test <- DLM(a_0 = 0, b_0 = 1, gamma = 1, sigma = 2, tau = 5, n_obs = 50)

plot(test$Y, type = "l", col = "red")
lines(pf_results2$x[2,], type = "l")
lines(test$X, type = "l", col = "blue", lty = "dotted")

#Sequential Importance Sampling
SISfilter <- function(M, y, simX0, simXt, loglik, a_0, b_0, gamma, sigma, tau){
  n_obs <- length(y)
  
  #trajectory
  x_all <- matrix(nrow = M, ncol = n_obs, data = 0)
  w_all <- matrix(nrow = M, ncol = n_obs, data = 0)
  #weights
  
  #prior
  for (t in 1:n_obs){
    #prior simulation
    if (t == 1){
      for (i in 1:M){
        x_all[i,t] <- simX0(a_0, b_0)
        w_all[i,t] <- 1/M
      }
    }
    
    else{
      logwstar <- numeric(M)
      for (i in 1:M){
        x_all[i,t] <- simXt(x_all[i,t-1], gamma, sigma)
        logwstar[i] <- log(w_all[i,t-1]) + loglik(y[t], x_all[i,t], tau)
      }
      
      wstar <- exp(logwstar)
      w_all[,t] <- wstar / sum(wstar)
    }
  }
  
  return(list(x = x_all, w = w_all))
}

pf_results <- SISfilter(M = 10, y = test$Y, simX0 = simX0_DLM, simXt = simXfwd_DLM, loglik = LL_DLM, 
                       a_0 = 0, b_0 = 1, gamma = 1, sigma = 2, tau = 5)

pf_results$x

##############################
simX0_DLM <- function(a0, b0){
  return(rnorm(1, a0, b0))
}

simXfwd_DLM <- function(x, gamma, sigma){
  return(rnorm(1, x*gamma, sigma))
}

LL_DLM <- function(y, x, tau){
  return(dnorm(y, x, tau, log = TRUE))
}

##############################
#GPT ggplot

library(ggplot2)
library(reshape2)

df <- as.data.frame(t(pf_results$x))
df$time <- 1:5

# Melt the data to long format
df_long <- melt(df, id.vars = "time", variable.name = "particle", value.name = "position")

# Add normalized importance weights (repeating for each time step)
df_long$weight <- 0
x <- 1
for (i in 1:10){
  for (j in 1:5){
    df_long[x,4] <- pf_results$w[i,j]
    x <- x + 1
  }
}

x <- 1
for (i in 1:nrow(df_long)){
  df_long[i,2] <- x
  if (i %% 5 == 0){
    x <- x + 1
  }
}

observed_df <- data.frame(observed = test$Y, time = 1:5)

# Plot using ggplot2
ggplot(df_long, aes(x = time, y = position, size = weight, color = particle)) +
  geom_point() +
  geom_line(size = 0.1) +
  geom_line(data = observed_df, aes(x = time, y = observed), color = "red", size = 1.2, lty = "longdash") +
  theme_minimal() +
  scale_size(range = c(1, 10)) +  # Adjust size range as needed
  labs(x = "Time Step",
       y = "Particle Location",
       size = "Weight",
       color = "Particle") +
  theme(legend.position = "right")
##############################

#ESS Calculation
ESS <- numeric(5)
for (i in 1:50){
  ESS[i] <- 1 / sum((pf_results2$w_preresample[,i])^2)
}

ess.df <- data.frame(ESS = ESS, time = 1:50)
ggplot(ess.df, aes(x = time, y = ESS)) +
  geom_line(size = 0.8) +
  geom_point(size = 2.5) +
  theme_minimal() +
  labs(x = "Time Step",
       y = "Effective Sample Size") +
  ylim(0, 10)

InverseCDF_Resample <- function(x, w){
  #input x samples and w weights as a vector
  
  #initialise
  n <- length(w)
  new_x <- x
  
  #create the cdf
  c <- numeric(n)
  c[1] <- w[1]
  for (i in 2:n){
    c[i] <- w[i] + c[i-1]
  }
  
  #generate uniform
  u <- runif(n, 0, 1)
  
  #final step
  for (i in 1:n){
    j <- which.min(c[c > u[i]])
    new_x[j] <- x[i]
  }
  
  return(new_x)
}

Stratified_Resample <- function(x, w){
  #input x samples and w weights as a vector
  
  #initialise
  n <- length(w)
  new_x <- x
  
  #create the cdf
  c <- numeric(n)
  c[1] <- w[1]
  for (i in 2:n){
    c[i] <- w[i] + c[i-1]
  }
  
  #generate uniform
  u <- numeric(n)
  for (i in 1:n){
    u[i] <- runif(1, (i-1)/n, i/n)
  }
  
  #final step
  for (i in 1:n){
    j <- which(u[i] <= c)[1]
    new_x[i] <- x[j]
  }
  
  return(new_x)
}

Systematic_Resample <- function(x, w){
  #input x samples and w weights as a vector
  
  #initialise
  n <- length(w)
  new_x <- x
  
  #create the cdf
  c <- numeric(n)
  c[1] <- w[1]
  for (i in 2:n){
    c[i] <- w[i] + c[i-1]
  }
  
  #generate uniform
  uniform <- runif(1)
  u <- numeric(n)
  for (i in 1:n){
    u[i] <- (i-1+uniform)/n
  }
  
  #final step
  for (i in 1:n){
    j <- which(u[i] <= c)[1]
    new_x[i] <- x[j]
  }
  
  return(new_x)
}

Residual_Resample <- function(x, w){
  #input x samples and w weights as a vector
  
  #initialise
  n <- length(w)
  new_x <- numeric(n)
  k <- 1
  
  #deterministic phase
  floor_nw <- floor(n*w)
  for (i in 1:n){
    if (floor_nw[i] > 0){
      for (j in 1:floor_nw[i]){
        new_x[k] <- x[i]
        k <- k + 1
      }
    }
  }
  
  #stochastic phase
  #initialise stochastic weightings
  r = n - sum(floor_nw)
  if (r > 0){
    epsilon <- n*w - floor(n*w)
    epsilon <- epsilon / (sum(epsilon)) #normalisation
    multinomial_x <- InverseCDF_Resample(x, epsilon)
    new_x[k:n] <- multinomial_x[1:r]
  }
  
  return(new_x)
}

#Bootstrap Filter
BSfilter <- function(M, y, simX0, simXt, loglik, a_0, b_0, gamma, sigma, tau, resample){
  #resample is the function used for the resampling algorithm
  n_obs <- length(y)
  
  #trajectory
  x_preresample <- matrix(nrow = M, ncol = n_obs, data = 0)
  w_preresample <- matrix(nrow = M, ncol = n_obs, data = 0)
  x_all <- matrix(nrow = M, ncol = n_obs, data = 0)
  w_all <- matrix(nrow = M, ncol = n_obs, data = 0)
  #weights
  
  #prior
  for (t in 1:n_obs){
    #prior simulation
    if (t == 1){
      for (i in 1:M){
        x_all[i,t] <- simX0(a_0, b_0)
        x_preresample[i,t] <- x_all[i,t]
        w_all[i,t] <- 1/M
        w_preresample[i,t] <- w_all[i,t]
      }
    }
    
    else{
      logwstar <- numeric(M)
      for (i in 1:M){
        x_all[i,t] <- simXt(x_all[i,t-1], gamma, sigma)
        x_preresample[i,t] <- x_all[i,t]
        logwstar[i] <- log(w_all[i,t-1]) + loglik(y[t], x_all[i,t], tau)
      }
      
      wstar <- exp(logwstar - max(logwstar))
      
      w_all[,t] <- wstar / sum(wstar)
      w_preresample[,t] <- w_all[,t]
      
      #resampling
      x <- as.vector(x_all[,t])
      w <- as.vector(w_all[,t])
      new_x <- as.matrix(resample(x, w))
      x_all[,t] <- t(new_x)
      w_all[,t] <- 1/M
    }
  }
  
  return(list(x = x_all, w = w_all, w_preresample = w_preresample))
}

pf_results2 <- BSfilter(M = 10, y = test$Y, simX0 = simX0_DLM, simXt = simXfwd_DLM, loglik = LL_DLM, 
                        a_0 = 0, b_0 = 1, gamma = 1, sigma = 2, tau = 5, resample = Systematic_Resample)

##########################################
#Numerical Experiment

set.seed(1234)
DLM_ts <- replicate(20, DLM(a_0 = 0, b_0 = 1, gamma = 1, sigma = 2, tau = 0.5, n_obs = 50))
#So DLM_ts has [,1] to [,20] simulations each with $X latent and $Y observed

set.seed(12345)
ESS_SIS <- matrix(0, nrow = 20, ncol = 50)
ESS_Multinomial <- matrix(0, nrow = 20, ncol = 50)
ESS_Residual <- matrix(0, nrow = 20, ncol = 50)
ESS_Stratified <- matrix(0, nrow = 20, ncol = 50)
ESS_Systematic <- matrix(0, nrow = 20, ncol = 50)
for (i in 1:20){
  #SIS approach
  pf_results <- SISfilter(M = 50, y = DLM_ts[,i]$Y, simX0 = simX0_DLM, simXt = simXfwd_DLM, loglik = LL_DLM, 
                          a_0 = 0, b_0 = 1, gamma = 1, sigma = 2, tau = 0.5)
  #Multinomial Resampling
  pf_results2 <- BSfilter(M = 50, y = DLM_ts[,i]$Y, simX0 = simX0_DLM, simXt = simXfwd_DLM, loglik = LL_DLM, 
                          a_0 = 0, b_0 = 1, gamma = 1, sigma = 2, tau = 0.5, resample = InverseCDF_Resample)
  #Residual Resampling
  pf_results3 <- BSfilter(M = 50, y = DLM_ts[,i]$Y, simX0 = simX0_DLM, simXt = simXfwd_DLM, loglik = LL_DLM, 
                         a_0 = 0, b_0 = 1, gamma = 1, sigma = 2, tau = 0.5, resample = Residual_Resample)
  #Stratified Resampling
  pf_results4 <- BSfilter(M = 50, y = DLM_ts[,i]$Y, simX0 = simX0_DLM, simXt = simXfwd_DLM, loglik = LL_DLM, 
                          a_0 = 0, b_0 = 1, gamma = 1, sigma = 2, tau = 0.5, resample = Stratified_Resample)
  #Systematic Resampling
  pf_results5 <- BSfilter(M = 50, y = DLM_ts[,i]$Y, simX0 = simX0_DLM, simXt = simXfwd_DLM, loglik = LL_DLM, 
                          a_0 = 0, b_0 = 1, gamma = 1, sigma = 2, tau = 0.5, resample = Systematic_Resample)
  
  for (j in 1:50){
    #SIS
    ESS_SIS[i,j] <- 1 / sum((pf_results$w[,j])^2)
    
    #Multinomial
    ESS_Multinomial[i,j] <- 1 / sum((pf_results2$w_preresample[,j])^2)
    
    #Residual
    ESS_Residual[i,j] <- 1 / sum((pf_results3$w_preresample[,j])^2)
    
    #Stratified
    ESS_Stratified[i,j] <- 1 / sum((pf_results4$w_preresample[,j])^2)
    
    #Systematic
    ESS_Systematic[i,j] <- 1 / sum((pf_results5$w_preresample[,j])^2)
  }
}

mean_SIS <- colMeans(ESS_SIS)
mean_Multinomial <- colMeans(ESS_Multinomial)
mean_Residual <- colMeans(ESS_Residual)
mean_Stratified <- colMeans(ESS_Stratified)
mean_Systematic <- colMeans(ESS_Systematic)

ESS_means_df <- data.frame(
  Time = 1:50,
  No_Resample = mean_SIS,
  Multinomial = mean_Multinomial,
  Residual = mean_Residual,
  Stratified = mean_Stratified,
  Systematic = mean_Systematic
)

ESS_long <- ESS_means_df %>%
  pivot_longer(cols = -Time, names_to = "Method", values_to = "ESS")

ggplot(ESS_long, aes(x = Time, y = ESS, color = Method)) +
  geom_line(size = 1) +
  labs(x = "Time",
       y = "Effective Sample Size (ESS)",
       color = "Method") +
  theme_minimal()

##########################################
#Stochastic Volatility Particle MCMC

Stochastic_Volatility <- function(gamma, sigma, n_obs){
  #Dynamic Linear Model
  #X_0 ~ N(0, sigma^2)
  #X_t|(X_{t-1} = x) ~ N(gamma*x, sigma^2)
  #Y_t|(X_{t} = x) ~ N(0, exp(x_t))
  #n_obs is how many observations the model has, n = 1,...,T.
  
  #generate trajectory of Xs
  X <- numeric(n_obs)
  X[1] <- rnorm(1, mean = 0, sd = sigma)
  
  for (t in 2:n_obs){
    X[t] <- rnorm(1, mean = gamma*X[t-1], sd = sigma)
  }
  
  #generate trajectory of Ys
  
  Y <- numeric(n_obs)
  
  for (t in 2:n_obs){
    Y[t] <- rnorm(1, mean = 0, sd = exp(X[t]))
  }
  
  return(list(X = X,Y = Y))
}

test <- Stochastic_Volatility(gamma = 0.95, sigma = exp(-1), n_obs = 400)

plot(test$Y, type = "h", col = "red", xlab = "Time", ylab = "X and Y values")
lines(test$X, type = "l", col = "blue",lwd = 1.5)

#BS Filter for SV

simX0_SV <- function(sigma){
  return(rnorm(1, mean = 0, sd = exp(sigma)))
}

simXfwd_SV <- function(x, gamma, sigma){
  return(rnorm(1, x*gamma, exp(sigma)))
}

LL_SV <- function(y, x){
  ll <- dnorm(y, 0, sd = exp(x), log = TRUE)
  if (!is.finite(ll)){
    ll <- -1e10
  }
  return(ll)
}

BSfilter_SV <- function(M, y, simX0, simXt, loglik, gamma,  sigma, resample){
  #resample is the function used for the resampling algorithm
  n_obs <- length(y)
  
  #trajectory
  x_preresample <- matrix(nrow = M, ncol = n_obs, data = 0)
  w_preresample <- matrix(nrow = M, ncol = n_obs, data = 0)
  x_all <- matrix(nrow = M, ncol = n_obs, data = 0)
  w_all <- matrix(nrow = M, ncol = n_obs, data = 0)
  ancestry <- matrix(nrow = M, ncol = n_obs)
  #weights
  
  #no ancestry in first 
  ancestry[,1] <- NA
  
  #prior
  for (t in 1:n_obs){
    #prior simulation
    if (t == 1){
      for (i in 1:M){
        x_all[i,t] <- simX0(sigma)
        x_preresample[i,t] <- x_all[i,t]
        w_all[i,t] <- 1/M
        w_preresample[i,t] <- w_all[i,t]
      }
    }
    
    else{
      logwstar <- numeric(M)
      for (i in 1:M){
        x_all[i,t] <- simXt(x_all[i,t-1], gamma, sigma)
        x_preresample[i,t] <- x_all[i,t]
        logwstar[i] <- log(w_all[i,t-1]) + loglik(y[t], x_all[i,t])
      }
      
      wstar <- exp(logwstar - max(logwstar))
      
      w_preresample[,t] <- wstar*exp(max(logwstar))
      w_all[,t] <- wstar / sum(wstar)
      
      
      #resampling
      x <- as.vector(x_all[,t])
      w <- as.vector(w_all[,t])
      resampled <- resample(x, w)
      x_all[, t] <- resampled$particles
      w_all[, t] <- rep(1 / M, M)
      ancestry[, t] <- resampled$ancestry
    }
  }
  
  return(list(x = x_all,
              w = w_all,
              w_preresample = w_preresample,
              ancestry = ancestry))
}

Residual_Resample <- function(x, w){
  #input x samples and w weights as a vector
  
  #initialise
  
  n <- length(w)
  new_x <- numeric(n)
  ancestry <- integer(n)
  k <- 1
  
  #deterministic phase
  floor_nw <- floor(n*w)
  for (i in 1:n){
    if (floor_nw[i] > 0){
      for (j in 1:floor_nw[i]){
        new_x[k] <- x[i]
        ancestry[k] <- i
        k <- k + 1
      }
    }
  }
  
  #stochastic phase
  #initialise stochastic weightings
  r = n - sum(floor_nw)
  if (r > 0){
    epsilon <- n*w - floor(n*w)
    epsilon <- epsilon / (sum(epsilon)) #normalisation
    multinomial_x <- InverseCDF_Resample(x, epsilon, r)
    new_x[k:n] <- multinomial_x$particles
    ancestry[k:n] <- multinomial_x$indices
  }
  
  return(list(particles = new_x, ancestry = ancestry))
}

InverseCDF_Resample <- function(x, w, r = length(x)) {
  # x: particle values
  # w: normalized weights (should sum to 1)
  # r: number of samples to draw
  
  n <- length(w)
  cdf <- cumsum(w)
  u <- sort(runif(r))
  
  indices <- integer(r)
  j <- 1
  for (i in 1:r) {
    while (u[i] > cdf[j]) {
      j <- j + 1
    }
    indices[i] <- j
  }
  
  x_resampled <- x[indices]
  return(list(particles = x_resampled, indices = indices))
}
# 
py <- numeric(50)
for (i in 1:50){
  print(i)
  out <- BSfilter_SV(
    M = 500,
    y = test$Y,
    simX0 = simX0_SV,
    simXt = simXfwd_SV,
    loglik = LL_SV,
    gamma = 0.95,
    sigma = exp(-1),
    resample = Residual_Resample
  )
  py[i] <- sum(log(colMeans(out$w_preresample)))
}

prior <- function(theta){
  return(dmvnorm(theta, mean = c(0.5, 0), sigma = 1*diag(2), log = TRUE))
}

pmmh <- function(prior, pf, theta0, n_its, M, y, simX0, simXt, loglik, resample, h){
  #initialisation
  acc <- 0
  theta <- matrix(nrow = length(theta0), ncol = n_its + 1)
  t <- length(y)
  x <- matrix(nrow = t, ncol = n_its + 1)
  
  #first step
  theta[,1] <- theta0
  u <- pf(M, y, simX0, simXt, loglik, theta0[1], theta0[2], resample) #generate auxiliary
  k0 <- sample(1:M, 1, prob = u$w_preresample[,t]) #sample a path x_{1:T}
  
  #trace the backwards path for x_{1:T} 
  counter <- t
  a <- k0
  while (counter > 0){
    x[counter,1] <- u$x[a,counter]
    a <- u$ancestry[a,counter]
    counter <- counter - 1
  }
  
  #p(y_{1:T}|theta)
  log_p_hat <- sum(log(colMeans(u$w_preresample)))
  log_pi_old <- log_p_hat + prior(theta0)
  
  #iterating for new thetas
  for (i in 2:(n_its+1)){
    print(i)
    theta_old <- theta[,i-1]
    #generate new proposal
    theta_new <- theta_old + t(rmvnorm(1, mean = c(0,0), sigma = h))
    u_new <- pf(M, y, simX0, simXt, loglik, theta_new[1], theta_new[2], resample)

    #calculate probability p(y_{1:T}|theta)
    log_p_hat <- sum(log(colMeans(u_new$w_preresample)))
    log_pi_new <- log_p_hat + prior(t(theta_new))
    
    #acceptance
    log_alpha <- log_pi_new - log_pi_old
    U <- runif(1)
    #accept
    if (U <= exp(log_alpha)){
      theta[,i] <- theta_new
      log_pi_old <- log_pi_new
      
      #sample a new path
      k <- sample(1:M, 1, prob = u_new$w_preresample[,t])
      
      #trace the path
      counter <- t
      a <- k
      while (counter > 0){
        x[counter,i] <- u_new$x[a,counter]
        a <- u_new$ancestry[a,counter]
        counter <- counter - 1
      }
      
      #update u
      u <- u_new
      acc <- acc + 1
    }
    #reject
    else{
      theta[,i] <- theta_old
      x[,i] <- x[,i-1]
    }
  }
  return(list(theta = theta, x = x, acc = acc))
}

# pm_result <- pmmh(
#   prior = prior,
#   pf = BSfilter_SV,
#   theta0 = c(0.5,0),
#   n_its = 1000,
#   M = 500,
#   y = test$Y,
#   simX0 = simX0_SV,
#   simXt = simXfwd_SV,
#   loglik = LL_SV,
#   resample = Residual_Resample,
#   h = 0.005*diag(2)
# )

pm_result2 <- pmmh(
  prior = prior,
  pf = BSfilter_SV,
  theta0 = c(0.5,0),
  n_its = 10000,
  M = 500,
  y = test$Y,
  simX0 = simX0_SV,
  simXt = simXfwd_SV,
  loglik = LL_SV,
  resample = Residual_Resample,
  h = 0.005*diag(2)
)


plot(pm_result2$theta[1,1:500], type = "l", ylab = "gamma", xlab = "Iteration")
abline(h = 0.95)
plot(pm_result2$theta[2,1:500], type = "l", ylab = "log(sigma)", xlab = "Iteration")
abline(h = -1)

plot(pm_result2$theta[1,100:10001], pm_result2$theta[2,100:10001], pch = 20, col = rgb(0, 0, 1, 0.02),
     xlab = "gamma",
     ylab = "log(sigma)")
abline(v = 0.95, lty = 2, col = "red", lwd = 1.5)
abline(h = -1, lty = 2, col = "red", lwd = 1.5)

hist(pm_result2$theta[1,100:10001], xlab = "gamma", main = "")
hist(pm_result2$theta[2,100:10001], xlab = "log(sigma)", main = "")


x_post <- pm_result2$x[,100:10001]

x_mean <- rowMeans(x_post)
x_lower <- apply(x_post, 1, quantile, probs = 0.025)
x_upper <- apply(x_post, 1, quantile, probs = 0.975)



plot(test$X, type = "l", col = "red", lwd = "1", xlab = "Time", ylab = "X")

lines(x_lower, col = "blue", lty = 2)
lines(x_upper, col = "blue", lty = 2)

plot(x_mean, type = "n", xlab = "Time", ylab = "X", ylim = c(-3, 3))
polygon(c(1:length(x_mean), rev(1:length(x_mean))),
        c(x_lower, rev(x_upper)),
        col = rgb(0, 0, 1, 0.2), border = NA)
lines(x_mean, col = "blue", lwd = 1)
lines(test$X, col = "red", lwd = 1.3)


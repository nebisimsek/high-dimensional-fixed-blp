#Cleaning
rm(list = ls())

#install.packages('extRemes')
#install.packages('matrixcalc')
#library(matrixcalc)
#install.packages('gmm')
#install.packages('glmnet')

library(glmnet)
library(gmm)
library(MASS)
library(AER)

start.time <- Sys.time()

set.seed(42)
#Number of markets
t <- 5
#Number of products
j <- 20
#Number of observations
obs <- j*t

#Number of characteristics
k <- 50

#Number of instruments
l <- 50

#Number of relevant instruments
r <- 5

#Assign the mean values of X
mu.x <- rep(0, k)

#Assign the mean values of Z
mu.z <- rep(0, (l+1))

#Variance covariance matrix for X vlues
sigma.x <- matrix(rep(0, k*k), ncol = k, nrow = k)
for (i.x in 1:k) {
  for (l.x in 1:k) {
    sigma.x[i.x, l.x] <- (0.8^abs(i.x-l.x))/8
  }
}

#Relevance of instruments
rel <- 0.2

#Variance covariance matrix for Z values
sigma.z <- diag(l+1)
for (i.z in 1:(1+r)) {
  for (l.z in 1:(1+r)) {
    sigma.z[i.z, l.z] <- (rel^abs(i.z-l.z))
  }
}
# sigma.z <- matrix(rep(0, l*l), ncol = l, nrow = l)
# for (i.z in 1:l) {
#   for (l.z in 1:l) {
#     sigma.z[i.z, l.z] <- (0.5^abs(i.z-l.z))
#   }
#}


nsim = 2500

#slope values for x characteristics
coeff.x <- as.vector(c(4, 4, 4, 4, rep(0, k-4)))
#At this point we assume errors are extreme value distributed

#price coefficient
coeff.p <- as.vector(c(1, 1, 0, 0, rep(0, k-4)))

beta.p <- -5

#Set the insensity of endogeneity
endo <- 0.5

# #instrument coefficient
# coeff.z <- as.vector(c(rep(1, 4), rep(0, k-4)))

#market subscript is
t.subs <- rep(c(1:t), each = j)
#Product subscript js
j.subs <- rep(1:j, t)

#Assign names for the columns
xnam <- paste("X", 1:k, sep = "")
znam <- paste("Z", 1:l, sep = "")

#Assign names for the releavant X and Z values
rel.xnam <- paste("X", 1:4, sep = "")
rel.znam <- paste("Z", 1:r, sep = "")

#OLS regression formula
fm.olsreg = paste0("delta ~ ", 
                   paste(xnam, collapse = " + "), " + ", 
                   paste("p", collapse = " + "))

#IV regression formula
fm.ivreg = paste0("delta ~ ", 
                  paste(xnam, collapse = " + "), " + ", 
                  paste("p", collapse = " + "), " | ",
                  paste(xnam, collapse = " + "), " + ",
                  paste(znam, collapse = " + "))

fm.oracle = paste0("delta ~ ", 
                   paste(rel.xnam, collapse = " + "), " + ", 
                   paste("p", collapse = " + "), " | ",
                   paste(rel.xnam, collapse = " + "), " + ",
                   paste(rel.znam, collapse = " + "))

#Storing Units
#Storing price values for all cases
p.values <- matrix(NaN, nrow = nsim, ncol = 6)
colnames(p.values) <- c('oracle', 'full_model','triple-lasso',
                        'only_first', 'wo_second', 'wo_third')

#Confidence interval for beta
coverage <- matrix(NaN, nrow = nsim, ncol = 5)
colnames(coverage) <- c('oracle', 'full_model','triple-lasso',
                        'only_first', 'wo_second', 'wo_third')

#Storing number of X characteristics included in the data
#And check whether relevant X char. included or not
x.count <- matrix(NaN, nrow = nsim, ncol = 4)
colnames(x.count) <- c('triple-lasso_count_rel_x',
                       'triple-lasso_count-nonrel_x',
                       'first_stage_count_rel_x',
                       'first_stage_count_nonrel_x')

z.count <- matrix(NaN, nrow = nsim, ncol = 2)
z.selection <- matrix(NaN, nrow = nsim, ncol = l)

lambda.set <- c(seq(0.02, 3, by = 0.02), seq(3.5,5,0.5), seq(6,20,2))


for (sim in 1:nsim) {
  #Creating x characteristics for each of j products in t markets
  X.jt <- kronecker(rep(1, t),mvrnorm(j, mu = mu.x, Sigma = sigma.x))
  
  #The base values for price and instrumental variable
  Z.jt <- mvrnorm(obs, mu = mu.z, Sigma = sigma.z)
  
  #Assigning the endogenous part
  xi.jt <- as.matrix(rnorm(obs, 0, endo))
  
  #Uniform error for price 
  nu.jt <- runif(obs)
  
  #real p data correlated with xi
  p.jt <- xi.jt + Z.jt[,1] + X.jt%*%coeff.p + nu.jt
  
  #instrumental variable correlated with p
  Z.jt <- Z.jt[, 2:(l + 1)]
  
  #putting everything together
  data <- as.matrix(cbind(X.jt, p.jt, xi.jt, Z.jt, t.subs, j.subs))
  #remove separate data to save space
  #rm(X.jt, p.jt, xi.jt, z.jt, t.subs, j.subs)
  #Assign names to related columns
  colnames(data)[1:as.integer(k+l+4)] <- c(xnam, 'p', 'xi', znam, 't', 'j')
  
  ###Obtaining shares by assuming EVT1 distribution##
  
  #Prepare delta
  delta <- X.jt %*% coeff.x + p.jt*beta.p + xi.jt
  
  #Prepare the numerator part
  numer <- exp(delta)
  
  #Prepare the denominator part
  temp <- cumsum(numer)
  temp1 <- temp[(1:t)*j]
  temp1[2:t] <- diff(temp1)
  denom <- 1 + temp1
  
  #assign shares
  s.jt <- as.vector(numer) / as.vector(kronecker(denom, rep(1, j)))
  
  #Assign the outside shares
  s.j0 <- kronecker(1/denom, rep(1, j))
  
  #attach share to data
  data <- cbind(data, s = s.jt)
  
  #Attach the share of outside good
  data <- cbind(data, s.j0)
  
  #Obtain delta due to log-linear shares ?? 
  data <- cbind(data, delta = log(s.jt) - log(s.j0))
  
  #OLS regression
  #try(ols <- lm(data = as.data.frame(data),
  #         formula = fm.olsreg))
  
  
  ##First step(Regressing delta on X values to find relevant ones)##
  #Specify the model
  x.model <- cv.glmnet(data[,1:k], data[, 'delta'], lambda = lambda.set)
  #Obtain the best  model
  x.best <- coef(x.model, s = 'lambda.min')
  #Get the non-zero coefficients 
  x.non.zero <- rownames(coef(x.model, s = 'lambda.min'))[coef(x.model, s = 'lambda.min')[,1]!= 0][-1] ### returns nonzero coefs
  
  ####store coefficients in matrix of course####
  #coef(b.model)
  
  ##Second step(Regressing delta on p values to find correlated ones)
  p.model <- cv.glmnet(data[,1:k], data[, 'p'], lambda = lambda.set)
  
  #Obtain the best model
  p.best <- coef(p.model, s = 'lambda.min')
  
  #Get the non-zero coefficients 
  p.non.zero <- rownames(coef(p.model, s = 'lambda.min'))[coef(p.model, s = 'lambda.min')[,1]!= 0][-1] 
  
  ##Third step(Regressing p on z values to find relevant instruments)
  
  z.model <- glmnet(data[,(k+3):(k+l+2)], data[, 'p'], lambda = 0.1)
  
  #Obtain the best model
  #z.best <- coef(z.model, s = z.model$lambda.min/100)
  #z.best <- coef(z.model, s = 'lambda.min')
  #Get the non-zero coefficients 
  
  z.non.zero <- rownames(coef(z.model))[coef(z.model)[,1]!= 0][-1] 
  
  
  #a <- glmnet(data[,(k+3):(k+l+2)], data[, 'p'], lambda = 0.1)
  #Get the union of both non-zero coefficients
  post.est.var <- union(x.non.zero, p.non.zero)
  
  #third step includes instrumental variable
  
  #prepare the triple_Lasso regression
  fm.post.est <- paste0("delta ~ ", 
                        paste(post.est.var, collapse = " + "), " + ", 
                        paste("p", collapse = " + "), " | ",
                        paste(post.est.var, collapse = " + "), " + ",
                        paste(z.non.zero, collapse = " + "))
  
  #Preparing regression for only first stage implementation
  fm.only.first <- paste0("delta ~ ", 
                          paste(x.non.zero, collapse = " + "), " + ", 
                          paste("p", collapse = " + "), " | ",
                          paste(x.non.zero, collapse = " + "), " + ",
                          paste(znam, collapse = " + "))
  
  #Preparing regression without selecting char. corr. with price 
  fm.wo.price <- paste0("delta ~ ", 
                        paste(x.non.zero, collapse = " + "), " + ", 
                        paste("p", collapse = " + "), " | ",
                        paste(x.non.zero, collapse = " + "), " + ",
                        paste(z.non.zero, collapse = " + "))
  
  #Preparing regression without selecting relevant instruments
  fm.wo.last <- paste0("delta ~ ", 
                       paste(post.est.var, collapse = " + "), " + ", 
                       paste("p", collapse = " + "), " | ",
                       paste(post.est.var, collapse = " + "), " + ",
                       paste(znam, collapse = " + "))
  
  #Oracle Benchmark regression
  try(oracle <- ivreg(data = as.data.frame(data),
                      formula = fm.oracle))
  
  #IV regression with all characteristics
  try(mo.ivreg <- ivreg(data = as.data.frame(data),
                        formula = fm.ivreg))
  
  #Triple Lasso regression
  try(post.iv.reg <- ivreg(data = as.data.frame(data),
                           formula = fm.post.est))
  
  #Regression with only first stage
  try(only.first.reg <- ivreg(data = as.data.frame(data),
                              formula = fm.only.first))
  
  #Regression without selecting char.corr. with price
  try(wo.price.reg <- ivreg(data = as.data.frame(data),
                            formula = fm.wo.price))
  
  #Regression without selecting relevant instruments
  try(wo.last.reg <- ivreg(data = as.data.frame(data),
                           formula = fm.wo.last))
  
  #Get the standard errors
  se.oracle <- sqrt(diag(vcov(oracle)))['p']
  se.mo.ivreg <- sqrt(diag(vcov(mo.ivreg)))['p']
  se.post.iv.reg <- sqrt(diag(vcov(post.iv.reg)))['p']
  se.only.first.reg <- sqrt(diag(vcov(only.first.reg)))['p']
  se.wo.price.reg <- sqrt(diag(vcov(wo.price.reg)))['p']
  se.wo.last.reg <- sqrt(diag(vcov(wo.last.reg)))['p']
  
  
  
  p.values.oracle <- oracle$coefficients['p']
  #Store p values with high dimensionality
  p.values.mo.ivreg <- mo.ivreg$coefficients['p']
  #Store p values for Triple-Lasso
  p.values.post.iv.reg <- post.iv.reg$coefficients['p']
  #Store p values for only first stage
  p.values.only.first.reg <- only.first.reg$coefficients['p']
  #Store p values for without second stage
  p.values.wo.price.reg <- wo.price.reg$coefficients['p']
  #Store p values for without third stage
  p.values.wo.last.reg <- wo.last.reg$coefficients['p']
  
  
  #Get the confidence interval
  CI.lower.oracle <- as.numeric(p.values.oracle - 1.96 * se.oracle)
  CI.upper.oracle <- as.numeric(p.values.oracle + 1.96 * se.oracle)
  coverage[sim,1] <- ifelse(CI.lower.oracle <= beta.p & beta.p <= CI.upper.oracle, 1, 0)
  
  #Get the confidence interval
  CI.lower.mo.ivreg <- as.numeric(p.values.mo.ivreg - 1.96 * se.mo.ivreg)
  CI.upper.mo.ivreg <- as.numeric(p.values.mo.ivreg + 1.96 * se.mo.ivreg)
  coverage[sim,2] <- ifelse(CI.lower.mo.ivreg <= beta.p & beta.p <= CI.upper.mo.ivreg, 1, 0)
  
  #Get the confidence interval
  CI.lower.post.iv.reg <- as.numeric(p.values.post.iv.reg - 1.96 * se.post.iv.reg)
  CI.upper.post.iv.reg <- as.numeric(p.values.post.iv.reg + 1.96 * se.post.iv.reg)
  coverage[sim,3] <- ifelse(CI.lower.post.iv.reg <= beta.p & beta.p <= CI.upper.post.iv.reg, 1, 0)
  
  #Get the confidence interval
  CI.lower.only.first.reg <- as.numeric(p.values.only.first.reg - 1.96 * se.only.first.reg)
  CI.upper.only.first.reg <- as.numeric(p.values.only.first.reg + 1.96 * se.only.first.reg)
  coverage[sim,4] <- ifelse(CI.lower.only.first.reg <= beta.p & beta.p <= CI.upper.only.first.reg, 1, 0)
  
  #Get the confidence interval
  CI.lower.wo.last.reg <- as.numeric(p.values.wo.last.reg - 1.96 * se.wo.last.reg)
  CI.upper.wo.last.reg <- as.numeric(p.values.wo.last.reg + 1.96 * se.wo.last.reg)
  coverage[sim,5] <- ifelse(CI.lower.wo.last.reg <= beta.p & beta.p <= CI.upper.wo.last.reg, 1, 0)
  
  
  
  
  #Store p values for the oracle benchmark
  p.values[sim, 1] <- oracle$coefficients['p']
  #Store p values with high dimensionality
  p.values[sim, 2] <- mo.ivreg$coefficients['p']
  #Store p values for Triple-Lasso
  p.values[sim, 3] <- post.iv.reg$coefficients['p']
  #Store p values for only first stage
  p.values[sim, 4] <- only.first.reg$coefficients['p']
  #Store p values for without second stage
  p.values[sim, 5] <- wo.price.reg$coefficients['p']
  #Store p values for without third stage
  p.values[sim, 6] <- wo.last.reg$coefficients['p']
  
  
  #Store whether relevant characteristics included or not
  #A vector return TRUE if the rel. char in post.est.var
  is.in.post.est.var <- rel.xnam %in% post.est.var
  #Store the proportion of included relevant characteristics
  x.count[sim, 1] <- sum(is.in.post.est.var, na.rm = TRUE)/length(rel.xnam)
  
  #Store how many non-relevant characteristics included 
  x.count[sim, 2] <- length(setdiff(post.est.var, rel.xnam))
  #----------------------------------------------
  #A vector return TRUE if the rel. char in x.non.zero
  is.in.x.non.zero <- rel.xnam %in% x.non.zero
  #Store the proportion of included relevant characteristics
  x.count[sim, 3] <- sum(is.in.x.non.zero, na.rm = TRUE)/length(rel.xnam)
  
  #Store how many non-relevant characteristics included 
  x.count[sim, 4] <- length(setdiff(x.non.zero, rel.xnam))
  
  #--------------------------------------------------
  
  #A vector return TRUE if the rel. char in x.non.zero
  is.in.z.non.zero <- rel.znam %in% z.non.zero
  #Store the proportion of included relevant characteristics
  z.count[sim, 1] <- sum(is.in.z.non.zero, na.rm = TRUE)/length(rel.znam)
  
  #Store how many non-relevant characteristics included 
  z.count[sim, 2] <- length(setdiff(z.non.zero, rel.znam))
  
  #--------------------------------------------------
  #To see which Zs are selected
  z.selection[sim, ] <- znam %in%z.non.zero
  
}

end.time <- Sys.time()
running.time <- start.time - end.time

summary(p.values)
summary(x.count)
summary(z.count)
summary(z.selection)
sum(z.count[,1] == 0)

summary(x.count[,2] >= 17)
#plot(density(p.values[,'wo_third']))

dens1 <- density(p.values[,'oracle'])
dens2 <- density(p.values[,'full_model'])
dens3 <- density(p.values[,'triple-lasso'])
dens4 <- density(p.values[,'only_first'])
dens5 <- density(p.values[,'wo_third'])


plot(density(p.values[,'oracle']), xlim = range(p.values[,'oracle'], 
                                                p.values[,'full_model'],
                                                p.values[,'triple-lasso'],
                                                p.values[,'only_first'],
                                                p.values[,'wo_third']),
     ylim = range(dens1$y,dens2$y,dens3$y, dens4$y, dens5$y), 
     col = 'ivory4', lwd = 2, lty = 5)
lines(density(p.values[,'full_model']), col = 'red', lwd = 2, lty = 3)
lines(density(p.values[,'triple-lasso']), col = 'blue', lwd = 2)
lines(density(p.values[,'only_first']), col = 'lightgreen', lwd = 2, lty = 5)
#lines(density(p.values[,'wo_second']), col = 'darkmagenta', lwd = 2)
lines(density(p.values[,'wo_third']), col = 'goldenrod1', lwd = 2,lty = 1)



legend("topright",legend = c("Oracle", "Unpenalized", "Triple-Lasso",
                             "Single Selection", "Wo Instrument Selection"),
       col = c('ivory4', 'red', 'blue', 'lightgreen',' goldenrod1' ),
       lty = c(2, 3, 1, 2, 1), lwd = 2, cex = 0.8)



#Forcing R to add at least some variable to model
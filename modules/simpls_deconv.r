#!/usr/bin/env R
#title      : simpls.r
#decription : Deconvolution by SIMPLS or Linear Regression method
#author     : Huamei Li
#date       : 27/11/2018
#type       : module

#####################################################################################
# Check package installed or not

required_pkgs = c('MASS', 'pls', 'transport')

for(p in required_pkgs){
  if(!require(p,character.only = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}

#####################################################################################
# Get relative fraction of each cell type

getfrac <- function(model, start = 1) {
	coeffs = coef(model)
	coeffs[coeffs < 0] = 0
	coeffs = as.vector(coeffs[start : length(coeffs)] / sum(coeffs[start : length(coeffs)]))
	return(coeffs)
}

######################################################################################
# Select model for deconvolution

select_models <- function(X_norm, Y, method = 'SIMPLS') {
	ncols  <- dim(X)[2]
	Y_norm <- as.matrix((Y - mean(Y)) / sd(Y))
	
	if (method == 'SIMPLS') {
		model <- mvr(Y_norm ~ X_norm, ncols, scale = FALSE, validation = 'CV', center = FALSE, method = 'simpls')
	} else {
		model <- lm(Y_norm ~ X_norm)
	}
  
	if (method == 'SIMPLS') {
		Y_predict <- model$fitted.values[, , ncols]
		rsquared  <- R2(model)$val[,,ncols]
		rmsep     <- RMSEP(model)$val[,,ncols][2]
	} else {
	    Y_predict <- model$fitted.values
		rsquared  <- summary(model)$r.squared
		rmsep     <- sqrt(mean(model$residuals^2))
	}

	start <- ifelse(method == 'SIMPLS', 1, 2)
	dist  <- wasserstein1d(Y_predict, Y_norm)
	return(list(Y_predict = Y_predict, Y = Y_norm, model = model, rsquared = rsquared, rmse = rmsep, start = start, dist = dist))
}

#####################################################################################
# Calculate P-value of estimate results


null_distribution <- function(Y, X, method = 'SIMPLS', iter_num = 1000) {	
	null_dists <- c()
	for (i in 1 : iter_num) {
		random_Y   <- runif(length(Y), min(Y), max(Y))
		rand_dist  <- select_models(X, random_Y, method = method)$dist
		null_dists <- c(null_dists, rand_dist)
	}
	null_dists <- sort(null_dists, decreasing = TRUE)
	return(null_dists)
}

calc_pval <- function(Y, X_norm, dist_obs, method = 'SIMPLS', iter_num = 1000) {
	null_dists <- null_distribution(Y, X_norm, method, iter_num)
	pval <- 1 - (which.min(abs(null_dists - dist_obs)) / length(null_dists))
	return(pval)
}

######################################################################################
# Deconvolution for mixed samples based on SIMPLS algorithm

simpls_deconv <- function(Y, X, method = 'SIMPLS', est.pval = FALSE, iter_num = 1000){
	X <- as.matrix(X) + 1; Y <- as.matrix(Y) + 1
	
	bc <- boxcox(Y ~ X, lambda = seq(0, 1.0, 0.1), plotit = FALSE)
	lambda <- bc$x[which.max(bc$y)]
	if (lambda == 0) { Y = log(Y); X = log(X); print(paste('[WARN] Lambda = 0, Log transfer for', method)) }
	if (lambda == 0.5) { Y = sqrt(Y); X = sqrt(X); print(paste('[WARN] Lambda = 0.5, Sqrt transfer for', method)) }
	
	X_norm <- as.matrix((X - mean(X)) / sd(as.vector(X)))
	model  <- select_models(X_norm, Y, method = method)
	coeffs <- getfrac(model$model, start = model$start)
	
	if (est.pval) {
		pval <- calc_pval(model$Y, X_norm, model$dist, iter_num = iter_num)
	} else {
		pval <- 9999
	}
	return(list('coeffs' = coeffs, 'R2' = model$rsquared, 'rmse' = model$rmse, 'pval' = pval))
}

#####################################################################################################


#!/usr/bin/env R
#title      : simpls.r
#decription : Deconvolution by SIMPLS or PCR method
#author     : Huamei Li
#date       : 27/11/2018
#type       : module

#####################################################################################
# Check package installed or not

required_pkgs = c('MASS', 'pls')

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

#####################################################################################
# Calculate P-value of estimate results

calc_pval <- function(X, L) {
	k   <- length(unique(L))
	n_j <- data.frame(table(L))$Freq
	s_j <- aggregate(X, by = list(L), FUN = sd)[2]
	x_j <- aggregate(X, by = list(L), FUN = mean)[2]

	m_j <- n_j - 1
	D <- (sum(m_j * (s_j / x_j))) / sum(m_j)
	D_AD <- (sum(m_j * (s_j / x_j - D)^2 )) / ( D^2 * (0.5 + D^2) )
	p_value <-  pchisq(D_AD, k - 1, lower.tail = TRUE)
	return(list(D_AD = D_AD, p_value = p_value))
}

######################################################################################
# Deconvolution for mixed samples based on SIMPLS algorithm

simpls_deconv <- function(Y, X, method = 'SIMPLS', est.pval = FALSE){
	X <- as.matrix(X) + 1; Y <- as.matrix(Y) + 1
	nrows <- dim(X)[1]; ncols <- dim(X)[2]; start <- 1
	
	bc <- boxcox(Y ~ X, lambda = seq(0, 1.0, 0.1), plotit = FALSE)
	lambda <- bc$x[which.max(bc$y)]
	if (lambda == 0) {Y = log(Y); X = log(X); print(paste('[WARN] Lambda = 0, Log transfer for', method))}
	if (lambda == 0.5) {Y = sqrt(Y); X = sqrt(X); print(paste('[WARN] Lambda = 0.5, Sqrt transfer for', method))}
	
	if ((lambda != 0) && (lambda != 0.5)) {
		if (est.pval) {
			X_norm <- as.matrix((X - min(X)) / (max(X) - min(X)))
			Y_norm <- as.matrix((Y - min(Y)) / (max(Y) - min(Y)))
		} else {
			X_norm <- as.matrix((X - mean(X)) / sd(as.vector(X)))
			Y_norm <- as.matrix((Y - mean(Y)) / sd(Y))
		}
	} else {
		X_norm <- X; Y_norm <- Y
	}

	if (method == 'SIMPLS') {
		model <- mvr(Y_norm ~ X_norm, ncols, scale = FALSE, validation = 'CV', center = FALSE, method = 'simpls')
	} else if (method == 'PCR') {
		model <- pcr(Y_norm ~ X_norm, ncols, scale = FALSE, validation = 'CV', center = FALSE)
	} else {
		model <- lm(Y_norm ~ X_norm); start = 2
	}
	
	if (method != 'LR') {
		Y_predict <- model$fitted.values[, , ncols]
		rsquared  <- R2(model)$val[,,ncols]
		rmsep     <- RMSEP(model)$val[,,ncols][2]
	} else {
		Y_predict <- model$fitted.values
		rsquared  <- summary(model)$r.squared
		rmsep     <- sqrt(mean(model$residuals^2))
	}
	coeffs <- getfrac(model, start = start)
	if (est.pval) {
		pval <- calc_pval(c(Y_norm, Y_predict), c(rep('a', nrows), rep('b', nrows)))$p_value
	} else {
		pval <- 9999
	}
	return(list('coeffs' = coeffs, 'R2' = rsquared, 'rmse' = rmsep, 'pval' = pval))
}

#####################################################################################################


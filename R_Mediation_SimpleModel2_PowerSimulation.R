########################################################################################################################################
### P o w e r   A n a l y s i s   f o r   t h e   S i m p l e   M e d i a t i o n   M o d e l   u s i n g   C o r r e l a t i o n s  ###
### Author: Thomas Ledermann                                                                                                         ###
### Created: Januray 3, 2025                                                                                                         ###
### Last update: Januray 22, 2025                                                                                                    ###
########################################################################################################################################

sampleSize <- 250
alphaLevel <- .05		# significance level
nsim <- 10000		# number of iterations
setSeed <- 123

rXM <- .14	# correlation between X and M
rXY <- .426795	# correlation between X and Y
rMY <- .426795	# correlation between M and Y

vX <- 1		# VAR(X), can be set to 1
vM <- 1		# VAR(M), can be set to 1
vY <- 1		# VAR(Y), can be set to 1

# install and load package
if(!require("lavaan")) install.packages("lavaan")
if(!require("paramtest")) install.packages("paramtest")
if(!require("simsem")) install.packages("simsem")
if(!require("dplyr")) install.packages("dplyr")

library(lavaan)
library(paramtest)
library(simsem)
library(dplyr)

# Covariance matrix
rMat <- matrix(NA, 3, 3)
rMat[lower.tri(rMat, diag = TRUE)] <- c(1, rXM, rXY, 1, rMY, 1)
rMat[upper.tri(rMat)] <- t(rMat)[upper.tri(rMat)]
vNames <- c('x', 'm', 'y')                    # Variable names
dimnames(rMat) <- list(vNames, vNames)
sds <- sqrt(c(vX, vM, vY))
covMat <- sds %*% t(sds) * rMat
covMat

## Estimate the parameters using the covariance matrix
med <- '	
	m ~ a * x
	y ~ b * m + c * x

	x ~~ vX * x
	m ~~ vM * m
	y ~~ vY * y 

	ab := a * b
	tot := ab + c
	ab_c := ab - c
'
fit <- sem(med, sample.cov = covMat, sample.nobs = 100000)
summary(fit, standardized = TRUE, rsquare = TRUE)

ests <- parameterEstimates(fit)
ests

# extract results from the lavaan output
ea <- ests[ests$label == 'a', 'est']
eb <- ests[ests$label == 'b', 'est']
ec <- ests[ests$label == 'c', 'est']
varX <- ests[ests$label == 'vX', 'est']
varM <- ests[ests$label == 'vM', 'est']
varY <- ests[ests$label == 'vY', 'est']

popEst <- cbind.data.frame(ea, eb, ec, varX, varM, varY)
popEst

## MC Simulation
models <- popEst %>% rowwise() %>% do({
	genModel <- paste0('
		Y ~ ', .$eb, '*M
		Y ~ ', .$ec, '*X
		M ~ ', .$ea, '*X
		X ~~ ', .$varX, '*X
		M ~~ ', .$varM, '*M
		Y ~~ ', .$varY, '*Y
	')
	fitModel <-'	
		M ~ ea * X
		Y ~ ec * X + eb * M

		X ~~ varX * X
		M ~~ varM * M
		Y ~~ varY * Y

		eab := ea * eb
		total := eab + ec
		eab_c := eab - ec
	'
	data.frame(ea = .$ea, eb = .$eb, ec = .$ec,
		     gen = genModel, fit = fitModel, stringsAsFactors = FALSE)
})
models

powerSim <- models %>% do({
	pSim <- sim(nRep = nsim, model = .$fit[1], n = sampleSize, generate = .$gen[1], lavaanfun = "sem")
	data_frame(ea = .$ea, eb = .$eb, ec = .$ec,
	pSimEst = list(pSim))
})

# Point estimates of the population model
powerSim

# Power estimates and other statistics
round(summaryParam(powerSim$pSimEst[[1]], detail = TRUE, alpha = 0.05), 3)

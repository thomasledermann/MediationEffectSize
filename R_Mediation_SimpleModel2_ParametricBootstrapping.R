library(lavaan)
library(MASS)

# covariance matrix
m <- matrix(c(
	0.9999900, 0.1399986, 0.4445956,
	0.1399986, 0.9999900, 0.4445956,
	0.4445956, 0.4445956, 1.0851651
	),
nrow = 3, byrow = TRUE)
rownames(m) <- colnames(m) <- c("X", "M", "Y")
m

# Model specification
modmed <- '
	M ~ a * X
	Y ~ b * M + c * X
	X ~~ X
	M ~~ M
	Y ~~ Y
	
	ab := a * b
	tot := ab + c
	ab_c := ab - c 
'

# Fit the model to get parameter estimates
n <- 50
fit <- sem(modmed, sample.cov = m, sample.nobs = n)
summary(fit, standardized = TRUE, rsq = TRUE)

out <- parameterEstimates(fit, standardized = TRUE, rsq = TRUE)	# $est[out$op == "r2"]
out
out[, "est"]
nrow(est)

### Parametric Bootstrapping 
paraBoot <- function(covMat, n = 50, nBoot = 10000) {
	medMod <- '
		M ~ a * X
		Y ~ b * M + c * X
		X ~~ vX * X
		M ~~ vM * M
		Y ~~ vY * Y
	
		ab := a * b
		tot := ab + c
		ab_c := ab - c 
	'
	fitM <- sem(medMod, sample.cov = covMat, sample.nobs = n)
	outpM <- parameterEstimates(fitM, standardized = TRUE, rsq = TRUE)

	# r_y(m.x)
	rspbMod <- '
		lx =~ NA*X
		lm =~ NA*M
		ly =~ NA*Y

		lx ~~ 1*lx
		lm ~~ 1*lm
		ly ~~ 1*ly

		X ~~ 0*X
		M ~~ 0*M
		Y ~~ 0*Y

		M ~ lx

		lm ~~ rspb*ly
		lx ~~ 0*lm
	'

	rspcMod <- '
		lx =~ NA*X
		lm =~ NA*M
		ly =~ NA*Y

		lx ~~ 1*lx
		lm ~~ 1*lm
		ly ~~ 1*ly

		X ~~ 0*X
		M ~~ 0*M
		Y ~~ 0*Y

		X ~ lm

		lx ~~ rspc*ly
		lm ~~ 0*lx
	'
	bResults <- matrix(NA, nrow = nBoot, ncol = nrow(outpM) + 10 + 2)
	colnames(bResults) <- c(unlist(outpM["label"])[1:9], "mR2", "yR2", "as", "bs", "cs", "vXs", "vMs", "vYs", "abs", "tots", "ab_cs", "rXY", "rspb", "rspc")

	for (i in 1:nBoot) {
		datR <- MASS::mvrnorm(n, mu = rep(0, nrow(covMat)), Sigma = covMat)
		colnames(datR) <- rownames(covMat)
		fitM <- sem(medMod, data = datR)
		fitRb <- sem(rspbMod, data = datR)
		fitRc <- sem(rspcMod, data = datR)
		if(fitM@Fit@converged == TRUE & fitRb@Fit@converged == TRUE & fitRc@Fit@converged == TRUE) { # check for convergence
			bResults[i, 1:11] <- parameterEstimates(fitM, rsq = TRUE)[, "est"]
			bResults[i, 12:20] <- parameterEstimates(fitM, standardized = TRUE)[, "std.all"]
			bResults[i, "rXY"] <- cor(datR[, "X"], datR[, "Y"])
			outpRb <- parameterEstimates(fitRb, standardized = TRUE)
			bResults[i, "rspb"] <- outpRb[outpRb$label == "rspb", "std.all"]
			outpRc <- parameterEstimates(fitRc, standardized = TRUE)
			bResults[i, "rspc"] <- outpRc[outpRc$label == "rspc", "std.all"]
		}
	}
  return(bResults)
}
bEst50 <- as.data.frame(paraBoot(covMat = m, n = 50, nBoot = 10000))
bEst50$upsilon <- bEst50$bs^2 - (bEst50$yR2 - bEst50$rXY^2)
bEst50$rspb2 <- bEst50$rspb^2
bEst50$rspc2 <- bEst50$rspc^2
bEst50$f2b <- bEst50$rspb^2/(1 - bEst50$yR2)
bEst50$f2c <- bEst50$rspc^2/(1 - bEst50$yR2)

bEst100 <- as.data.frame(paraBoot(covMat = m, n = 100, nBoot = 10000))
bEst100$upsilon <- bEst100$bs^2 - (bEst100$yR2 - bEst100$rXY^2)

bEst150 <- as.data.frame(paraBoot(covMat = m, n = 150, nBoot = 10000))
bEst150$upsilon <- bEst150$bs^2 - (bEst150$yR2 - bEst150$rXY^2)
bEst150$rspb2 <- bEst150$rspb^2
bEst150$rspc2 <- bEst150$rspc^2
bEst150$f2b <- bEst150$rspb^2/(1 - bEst150$yR2)
bEst150$f2c <- bEst150$rspc^2/(1 - bEst150$yR2)

bEst200 <- as.data.frame(paraBoot(covMat = m, n = 200, nBoot = 10000))
bEst200$upsilon <- bEst200$bs^2 - (bEst200$yR2 - bEst200$rXY^2)

bEst250 <- as.data.frame(paraBoot(covMat = m, n = 250, nBoot = 10000))
bEst250$upsilon <- bEst250$bs^2 - (bEst250$yR2 - bEst250$rXY^2)
bEst250$rspb2 <- bEst250$rspb^2
bEst250$rspc2 <- bEst250$rspc^2
bEst250$f2b <- bEst250$rspb^2/(1 - bEst250$yR2)
bEst250$f2c <- bEst250$rspc^2/(1 - bEst250$yR2)

# Remove non-converged bootstrap samples - complete cases
bEst50 <- bEst50[complete.cases(bEst50), ]
bEst100 <- bEst100[complete.cases(bEst100), ]
bEst150 <- bEst150[complete.cases(bEst150), ]
bEst200 <- bEst200[complete.cases(bEst200), ]
bEst250 <- bEst250[complete.cases(bEst250), ]

# Calculate 95% percentile confidence intervals 
ci50 <- apply(bEst50, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
ci100 <- apply(bEst100, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
ci150 <- apply(bEst150, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
ci200 <- apply(bEst200, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
ci250 <- apply(bEst250, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
round(ci50, 3)
round(ci100, 3)
round(ci150, 3)
round(ci200, 3)
round(ci250, 3)


### Sensitivity Analysis: Parametric Bootstrapping Residual Correlation M-Y
paraBootSA <- function(covMat, n = 50, nBoot = 10000) {
	modSA <- '
		M ~ a * X
		Y ~ c * X

		M ~~ r * Y

		X ~~ X
		M ~~ M
		Y ~~ Y 
	'

	fit0 <- sem(modSA, sample.cov = covMat, sample.nobs = n)
	outp <- parameterEstimates(fit0, standardized = TRUE)
	bResults <- matrix(NA, nrow = nBoot, ncol = nrow(outp))
	colnames(bResults) <- c(unlist(outp["label"])[1:6])

	for (i in 1:nBoot) {
		datR <- MASS::mvrnorm(n, mu = rep(0, nrow(covMat)), Sigma = covMat)
		colnames(datR) <- rownames(covMat)
		fitBoot <- sem(modSA, data = datR)
		if(fitBoot@Fit@converged == TRUE) { # check for convergence
			bResults[i, 1:6] <- parameterEstimates(fitBoot, standardized = TRUE)[, "std.all"]
		}
	}
  return(bResults)	# (bResults)
}
bEst50 <- paraBootSA(covMat = m, n = 50, nBoot = 10000)

bEst100 <- paraBootSA(covMat = m, n = 100, nBoot = 10000)

bEst150 <- paraBootSA(covMat = m, n = 150, nBoot = 10000)

bEst200 <- paraBootSA(covMat = m, n = 200, nBoot = 10000)

bEst250 <- paraBootSA(covMat = m, n = 250, nBoot = 10000)

# Remove non-converged bootstrap samples - complete cases

bEst50 <- bEst50[complete.cases(bEst50), ]
bEst100 <- bEst100[complete.cases(bEst100), ]
bEst150 <- bEst150[complete.cases(bEst150), ]
bEst200 <- bEst200[complete.cases(bEst200), ]
bEst250 <- bEst250[complete.cases(bEst250), ]

# Calculate 95% percentile confidence intervals 
ci50 <- apply(bEst50, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
ci100 <- apply(bEst100, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
ci150 <- apply(bEst150, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
ci200 <- apply(bEst200, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
ci250 <- apply(bEst250, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
round(ci50, 3)
round(ci100, 3)
round(ci150, 3)
round(ci200, 3)
round(ci250, 3)




#### Calculate the mean estimates and standard errors
means <- apply(bEst, 2, mean, na.rm = TRUE)
sds <- apply(bEst, 2, sd, na.rm = TRUE)
print(cbind(means, sds))

#### Example of how to get p-values from the bootstrap results:
p_value_ab <- mean(abs(bEst[, "ab"]) > abs(coef(fit1)["a"]*coef(fit1)["b"]))
print(paste("P-value for indirect effect (ab):", p_value_ab))

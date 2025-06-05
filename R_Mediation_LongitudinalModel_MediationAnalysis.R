#################################################################
### L o n g i t u d i n a l   M e d i a t i o n   M o d e l   ###
### Author: Thomas Ledermann                                  ###
### Created: April 4, 2024                                    ###
### Last update: June 6, 2025                                 ###
#################################################################

### This illustration uses data from Midlife in the United States (MIDUS) study.
### Datsets used are:	MIDUS2 Biomarker 29282-0001-Data.sav
### 				MIDUS_W2 04652-0001-Data.sav
###				MIDUS_W3 M3_P1_SURVEY_N3294_20181219.sav
### The data can be downladed from: https://midus.wisc.edu/?page_id=89
### Variables used: "B1PAGE_M2", "B4QCT_EA", "B1SESTEE", "B1SNEGPA", "C1SNEGPA"

# emotional abuse 
library(car)
datall$B4QCT_EA
table(datall$B4QCT_EA)
datall$eAbuse <- as.numeric(recode(datall$B4QCT_EA,"10:25 = 1; else = -1")) # 1 = yes, -1 = no
table(datall$eAbuse)

dat <- datall[c("B1PAGE_M2", "eAbuse", "B1SESTEE", "B1SNEGPA", "C1SNEGPA")]
names(dat) <- c("age", "X", "M1", "M2", "Y")

dat$age <- as.numeric(dat$age)
dat$X <- as.numeric(dat$X)
dat$M1 <- as.numeric(dat$M1)
dat$M2 <- as.numeric(dat$M2)
dat$Y <- as.numeric(dat$Y)
round(cor(dat),3)
nrow(dat)	# 396

## descriptive statistics
descriptive.table(vars = d(age, X, M1, M2, Y), data = dat, func.names = c("Mean", "St. Deviation", "Min", "Max", "Valid N", "Skew", "Kurtosis"))

### Correlations
## Point-biserial correlations
library(psych)
psych::biserial(dat$M1, dat$X)	# -.327
psych::biserial(dat$M2, dat$X)	# .273

# Bootstrapping point-biseral r 
boot <- function(data, Xd, M1c, M2c, nBoot = 5000) {
	n <- nrow(data)
	bRes <- matrix(NA, nBoot, 2)
	for (i in 1:nBoot) {
		bcases <- sample(1:n, n, replace = TRUE)
		bdata <- data[bcases, ]
		bRes[i, 1] <- psych::biserial(bdata[[M1c]], bdata[[Xd]])
		bRes[i, 2] <- psych::biserial(bdata[[M2c]], bdata[[Xd]])
	}
	colnames(bRes) <- c(paste("rpb.", Xd, M1c, sep = ""), paste("rpb.", Xd, M2c, sep = ""))
	return(bRes)
}
bres <- boot(dat, "X", "M1", "M2", nBoot = 5000)
head(bres)
rpbCI <- round(apply(bres, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE), 3)
rpbCI


### squared semi-partial correlations
# full model
lmfull <- summary(lm(Y ~ X + M1 + M2, dat))
R2full <- lmfull$r.squared
R2full

# with M1 and M2
lmM1M2 <- summary(lm(Y ~ M1 + M2, dat))
R2M1M2 <- lmM1M2$r.squared
R2M1M2

# with X and M1
lmXM1 <- summary(lm(Y ~ X + M1, dat))
R2XM1 <- lmXM1$r.squared
R2XM1

# with X and M2
lmXM2 <- summary(lm(Y ~ X + M2, dat))
R2XM2 <- lmXM2$r.squared
R2XM2

# b1
r2_spb1 <- R2full - R2XM2
r2_spb1

# b2
r2_spb2 <- R2full - R2XM1
r2_spb2

# c'
r2_spc <- R2full - R2M1M2
r2_spc

# alternative method
lmX.M12 <- lm(X ~ M1 + M2, dat)
lmM1.XM2 <- lm(M1 ~ X + M2, dat)
lmM2.XM1 <- lm(M2 ~ X + M1, dat)

rm1y <- cor(residuals(lmM1.XM2), dat$Y)
rm1y^2									# .1532865
rm2y <- cor(residuals(lmM2.XM1), dat$Y)
rm2y^2									# .1269534
rxy <- cor(residuals(lmX.M12), dat$Y)
rxy^2

# alternative method 2
lmr <- lm(M1 ~ X + M2, dat)
r_spb1 <- cor(dat$Y, residuals(lmr))
r_spb1
r_spb1^2

lmr <- lm(M2 ~ X + M1, dat)
r_spb2 <- cor(dat$Y, residuals(lmr))
r_spb2
r_spb2^2

lmr <- lm(X ~ M1 + M2, dat)
summary(lmr)
r_spc <- cor(dat$Y, residuals(lmr))
r_spc
r_spc^2

## Cohen's f^2
f2_b1 <- r2_spb1/(1 - R2full)
f2_b1

f2_b2 <- r2_spb2/(1 - R2full)
f2_b2

f2_c <- r2_spc/(1 - R2full)
f2_c

# Bootstrapping point-biseral r 
boot <- function(data, x, m1, m2, y, nBoot = 5000) {
	n <- nrow(data)
	bRes <- matrix(NA, nBoot, 6)
	for (i in 1:nBoot) {
		bcases <- sample(1:n, n, replace = TRUE)
		bdata <- data[bcases, ]
  
		lmfull <- summary(lm(as.formula(paste(y, "~", x, "+", m1, "+", m2)), data = bdata))
		R2full <- lmfull$r.squared

		# with M1 and M2
		lmM1M2 <- summary(lm(as.formula(paste(y, "~", m1, "+", m2)), data = bdata))
		R2M1M2 <- lmM1M2$r.squared

		# with X and M1
		lmXM1 <- summary(lm(as.formula(paste(y, "~", x, "+", m1)), data = bdata))
		R2XM1 <- lmXM1$r.squared

		# with X and M2
		lmXM2 <- summary(lm(as.formula(paste(y, "~", x, "+", m2)), data = bdata))
		R2XM2 <- lmXM2$r.squared

		# b1
		r2_spb1 <- R2full - R2XM2

		# b2
		r2_spb2 <- R2full - R2XM1

		# c'
		r2_spc <- R2full - R2M1M2

		bRes[i, 1] <- R2full - R2XM2
		bRes[i, 2] <- R2full - R2XM1
		bRes[i, 3] <- R2full - R2M1M2
		bRes[i, 4] <- r2_spb1/(1 - R2full)
  		bRes[i, 5] <- r2_spb2/(1 - R2full)
  		bRes[i, 6] <- r2_spc/(1 - R2full)
	}
	colnames(bRes) <- c("r2_spb1", "r2_spb2", "r2_spc", "f2_b1", "f2_b2", "f2_c")
	return(bRes)
}
bres <- boot(dat, "X", "M1", "M2", "Y", nBoot = 5000)
head(bres)
rpbCI <- round(apply(bres, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE), 3)
rpbCI

### Mediation Analyses
library(lavaan)
mod2 <- '
	M1 ~ a1 * X
	M2 ~ a2 * X
	Y ~ b1 * M1 + b2 * M2 + c * X

	M1 ~~ cov*M2
		
	a1b1 := a1 * b1
	a2b2 := a2 * b2
	IEtot := a1 * b1 + a2 * b2
	total := a1 * b1 + a2 * b2 + c
	a1b1_a2b2 := a1 * b1 - a2 * b2
	a1b1_c := a1 * b1 - c
	a2b2_c := a2 * b2 - c
	IEtot_c := IEtot - c
	tot_c := total - c
	tot_a1b1 := total - a1b1
'
fit <- sem(mod2, data = dat, meanstructure = TRUE)
summary(fit, standardized = TRUE, rsquare = TRUE)
parameterEstimates(fit)
parameterEstimates(fit)[c("label", "est")]

# partial standardized point estimates
coef(fit, type = "user")["a1"]/sd(dat$M1, na.rm = TRUE)
coef(fit, type = "user")["a2"]/sd(dat$M2, na.rm = TRUE)
coef(fit, type = "user")["c"]/sd(dat$Y, na.rm = TRUE)
coef(fit, type = "user")["a1b1"]/sd(dat$Y, na.rm = TRUE)
coef(fit, type = "user")["a2b2"]/sd(dat$Y, na.rm = TRUE)
coef(fit, type = "user")["IEtot"]/sd(dat$Y, na.rm = TRUE)
coef(fit, type = "user")["total"]/sd(dat$Y, na.rm = TRUE)

# bootstrapping
fitb <- sem(mod2, data = dat, meanstructure = TRUE, se = "bootstrap", bootstrap = 5000)
bfit <- parameterEstimates(fitb, boot.ci.type = "perc", level = 0.95, ci = TRUE)
bfit


### Bootstrapping partial standardized effects
nBoot <- 5000
bRes <- matrix(NA, nBoot, 9)	# 9 effects

# Perform the bootstrapping
set.seed(123)  # For reproducibility
for (i in 1:nBoot) {
	bcases <- sample(1:nrow(dat), nrow(dat), replace = TRUE)
	bdata <- dat[bcases, ]
  
	fit <- sem(mod2, data = bdata, meanstructure = TRUE)
  
	bRes[i, 1] <- coef(fit, type = "user")["a1"] / sd(bdata$M1, na.rm = TRUE)
	bRes[i, 2] <- coef(fit, type = "user")["a2"] / sd(bdata$M2, na.rm = TRUE)
	bRes[i, 3] <- coef(fit, type = "user")["b1"] * sd(bdata$M1, na.rm = TRUE) / sd(bdata$Y, na.rm = TRUE)
	bRes[i, 4] <- coef(fit, type = "user")["b2"] * sd(bdata$M2, na.rm = TRUE) / sd(bdata$Y, na.rm = TRUE)
	bRes[i, 5] <- coef(fit, type = "user")["c"]  / sd(bdata$Y, na.rm = TRUE)
	bRes[i, 6] <- coef(fit, type = "user")["a1b1"] / sd(bdata$Y, na.rm = TRUE)
	bRes[i, 7] <- coef(fit, type = "user")["a2b2"] / sd(bdata$Y, na.rm = TRUE)
	bRes[i, 8] <- coef(fit, type = "user")["IEtot"] / sd(bdata$Y, na.rm = TRUE)
	bRes[i, 9] <- coef(fit, type = "user")["total"] / sd(bdata$Y, na.rm = TRUE)
}

colnames(bRes) <- c("a1", "a2", "b1", "b2", "c", "a1b1", "a2b2", "IEtot", "total")

# Calculate bootstrap confidence intervals
round(apply(bRes, 2, quantile, probs = c(0.025, 0.975)), 3)

### Sensitivity Analysis
library(lavaan)
mod2sen <- '
	M1 ~ a1 * X
	M2 ~ a2 * X
	Y ~ b2 * M2 + c * X

	M1 ~~ rM1Y*Y

	M1 ~~ cov*M2
'

fit <- sem(mod2sen, data = dat, meanstructure = TRUE)
summary(fit, standardized = TRUE, rsquare = TRUE)
parameterEstimates(fit, standardized = TRUE)
parameterEstimates(fit, standardized = TRUE)[c("label", "est", "std.all")]
parameterEstimates(fit, standardized = TRUE)[label = "rM1Y", "std.all"]
coef(fit, type = "user")["rM1Y"]
coef(fit)

# bootstrapping
fitb <- sem(mod2sen, data = dat, meanstructure = TRUE, se = "bootstrap", bootstrap = 5000)
bfit <- parameterEstimates(fitb, boot.ci.type = "perc", level = 0.95, ci = TRUE)
bfit

install.packages("semhelpinghands")
library(semhelpinghands)
standardizedSolution_boot_ci(fitb)


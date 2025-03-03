######################################################################################################
### S i m p l e   M e d i a t i o n   M o d e l   u s i n g   P o p u l a t i o n   E f f e c t s  ###
### Author: Thomas Ledermann                                                                       ###
### Created: January 3, 2025                                                                       ###
### Last update: January 22, 2025                                                                  ###
######################################################################################################

### a = 0.14, b = c' = 0.39
covDat <- matrix(c(1.00, 0.14, 0.39, 0.14, 1.00, 0.39, 0.39, 0.39, 1.00 ), nrow = 3, byrow = TRUE)
covDat

rownames(covDat) <- colnames(covDat) <- c("X", "M", "Y")
covDat

library(lavaan)
medmod0 <- ' 
	M ~ 0.14 * X
	Y ~ 0.39 * M + 0.39 * X

	X ~~ X
	M ~~ M
	Y ~~ Y 
'
fit0 <- sem(medmod0, sample.cov = covDat, sample.nobs = 100000)
summary(fit0, standardized = TRUE, rsq = TRUE)

# Get the model-implied covariance matrix - lavInspect(fit0, "implied") is an alternative
iCov <- fitted(fit0)
iCov
fitted(fit0)

iCov <- as.matrix(as.data.frame(iCov))
colnames(iCov) <- rownames(iCov)
iCov

# model implied correlation matrix
iCor <- cov2cor(iCov)
iCor

medmod <- '
	M ~ a * X
	Y ~ b * M + c * X

	X ~~ X
	M ~~ M
	Y ~~ Y 

	ab := a * b
	tot := ab + c
	ab_c := ab - c
'
fit1 <- sem(medmod, sample.cov = iCov, sample.nobs = 100000)
summary(fit1, standardized = TRUE, rsq = TRUE)

# Sensitivity analysis - rRexXM
modSA <- '
	M ~ a * X
	Y ~ 0 * M + c * X

	M ~~ r * Y

	X ~~ X
	M ~~ M
	Y ~~ Y 
'
fitSA <- sem(modSA, sample.cov = iCov, sample.nobs = 100000)
summary(fitSA, standardized = TRUE)
parameterEstimates(fitSA, standardized = TRUE)

# N = 50
fit1.50 <- sem(medmod, sample.cov = iCov, sample.nobs = 50)
summary(fit1.50, standardized = TRUE, rsq = TRUE)
pe1.50 <- parameterEstimates(fit1.50, standardized = TRUE, rsq = TRUE)
pe1.50
bs <- pe1.50[pe1.50$label == "b", "std.all"]
bs
as <- pe1.50[pe1.50$label == "a", "std.all"]
as
yR2.50 <- pe1.50[pe1.50$lhs == "Y" & pe1.50$op == "r2", "est"]
yR2.50

rxy <- iCor["X", "Y"]
rxy

upsilon50 <- bs^2	- (yR2.50 - rxy^2)	# b_ym.x^2 - (R2_y.mx - r_yx^2)
upsilon50

fitSA.50 <- sem(modSA, sample.cov = iCov, sample.nobs = 50)
summary(fitSA, standardized = TRUE)

# N = 100
fit1.100 <- sem(medmod, sample.cov = iCov, sample.nobs = 100)
summary(fit1.100, standardized = TRUE, rsq = TRUE)
pe1.100 <- parameterEstimates(fit1.100, standardized = TRUE, rsq = TRUE)
pe1.100
yR2.100 <- pe1.100[pe1.100$lhs == "Y" & pe1.100$op == "r2", "est"]
yR2.100

# N = 150
fit1.150 <- sem(medmod, sample.cov = iCov, sample.nobs = 150)
summary(fit1.150, standardized = TRUE, rsq = TRUE)
pe1.150 <- parameterEstimates(fit1.150, standardized = TRUE, rsq = TRUE)
pe1.150
yR2.150 <- pe1.150[pe1.150$lhs == "Y" & pe1.150$op == "r2", "est"]
yR2.150

# N = 200
fit1.200 <- sem(medmod, sample.cov = iCov, sample.nobs = 200)
summary(fit1.200, standardized = TRUE, rsq = TRUE)
pe1.200 <- parameterEstimates(fit1.200, standardized = TRUE, rsq = TRUE)
pe1.200
yR2.200 <- pe1.200[pe1.200$lhs == "Y" & pe1.200$op == "r2", "est"]
yR2.200

# N = 250
fit1.250 <- sem(medmod, sample.cov = iCov, sample.nobs = 250)
summary(fit1.250, standardized = TRUE, rsq = TRUE)
pe1.250 <- parameterEstimates(fit1.250, standardized = TRUE, rsq = TRUE)
pe1.250
yR2.250 <- pe1.250[pe1.250$lhs == "Y" & pe1.250$op == "r2", "est"]
yR2.250

fit1.300 <- sem(medmod, sample.cov = iCov, sample.nobs = 300)
summary(fit1.300, standardized = TRUE, rsq = TRUE)
pe1.300 <- parameterEstimates(fit1.300, standardized = TRUE, rsq = TRUE)
pe1.300
yR2.300 <- pe1.300[pe1.300$lhs == "Y" & pe1.300$op == "r2", "est"]
yR2.300

### semi-partial correlation
# r_y(m.x)
rspb <- '
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

	lm ~~ rsp*ly
	lx ~~ 0*lm
'
## N = 50
fitrspb50 <- sem(rspb, sample.cov = iCov, sample.nobs = 50)
summary(fitrspb50, standardized = TRUE)
peb50 <- parameterEstimates(fitrspb50, standardized = TRUE)
peb50
rspb50 <- peb50[peb50$label == "rsp", "std.all"]
rspb50
rspb2.50 <- rspb50^2
rspb2.50

# f2
f2.50 <- rspb2.50/(1 - yR2.50)
f2.50	# 0.01702596

## N = 100
fitrspb100 <- sem(rspb, sample.cov = iCov, sample.nobs = 100)
summary(fitrspb100, standardized = TRUE)
peb100 <- parameterEstimates(fitrspb100, standardized = TRUE)
peb100
rspb100 <- peb100[peb100$label == "rsp", "std.all"]
rspb100
rspb2.100 <- rspb100^2
rspb2.100

# f2
f2.100 <- rspb2.100/(1 - yR2.100)
f2.100	# 0.01702596

## N = 150
fitrspb150 <- sem(rspb, sample.cov = iCov, sample.nobs = 150)
summary(fitrspb150, standardized = TRUE)
peb150 <- parameterEstimates(fitrspb150, standardized = TRUE)
peb150
rspb150 <- peb150[peb150$label == "rsp", "std.all"]
rspb150
rspb2.150 <- rspb150^2
rspb2.150

# f2
f2.150 <- rspb2.150/(1 - yR2.150)
f2.150	# 0.01702596


fitrspb200 <- sem(rspb, sample.cov = iCov, sample.nobs = 200)
summary(fitrspb200, standardized = TRUE)
peb200 <- parameterEstimates(fitrspb200, standardized = TRUE)
peb200
rspb200 <- peb200[peb200$label == "rsp", "std.all"]
rspb200
rspb2.200 <- rspb200^2
rspb2.200

# f2
f2.200 <- rspb2.200/(1 - yR2.200)
f2.200	# 0.01702596


fitrspb250 <- sem(rspb, sample.cov = iCov, sample.nobs = 250)
summary(fitrspb250, standardized = TRUE)
peb250 <- parameterEstimates(fitrspb250, standardized = TRUE)
peb250
rspb250 <- peb250[peb250$label == "rsp", "std.all"]
rspb250
rspb2.250 <- rspb250^2
rspb2.250

# f2
f2.250 <- rspb2.250/(1 - yR2.250)
f2.250	# 0.01702596


fitrspb300 <- sem(rspb, sample.cov = iCov, sample.nobs = 300)
summary(fitrspb300, standardized = TRUE)
peb300 <- parameterEstimates(fitrspb300, standardized = TRUE)
peb300
rspb300 <- peb300[peb300$label == "rsp", "std.all"]
rspb300
rspb2.300 <- rspb300^2
rspb2.300

# f2
f2.300 <- rspb2.300/(1 - yR2.300)
f2.300	# 0.01702596


# r_y(x.m)
rspc <- '
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

	lx ~~ rsp*ly
	lm ~~ 0*lx
'
fitrspc50 <- sem(rspc, sample.cov = iCov, sample.nobs = 50)
summary(fitrspc50, standardized = TRUE)
pec50 <- parameterEstimates(fitrspc50, standardized = TRUE)
pec50
rspc50 <- pec50[pec50$label == "rsp", "std.all"]
rspc50
rspc2.50 <- rspc50^2
rspc2.50

# f2
f2.50 <- rspc2.50/(1 - yR2.50)
f2.50	# 0.01702596


fitrspc100 <- sem(rspc, sample.cov = iCov, sample.nobs = 100)
summary(fitrspc100, standardized = TRUE)
pec100 <- parameterEstimates(fitrspc100, standardized = TRUE)
pec100
rspc100 <- pec100[pec100$label == "rsp", "std.all"]
rspc100
rspc2.100 <- rspc100^2
rspc2.100

# f2
f2.100 <- rspc2.100/(1 - yR2.100)
f2.100	# 0.01702596


fitrspc150 <- sem(rspc, sample.cov = iCov, sample.nobs = 150)
summary(fitrspc150, standardized = TRUE)
pec150 <- parameterEstimates(fitrspc150, standardized = TRUE)
pec150
rspc150 <- pec150[pec150$label == "rsp", "std.all"]
rspc150
rspc2.150 <- rspc150^2
rspc2.150

# f2
f2.150 <- rspc2.150/(1 - yR2.150)
f2.150	# 0.01702596


fitrspc200 <- sem(rspc, sample.cov = iCov, sample.nobs = 200)
summary(fitrspc200, standardized = TRUE)
pec200 <- parameterEstimates(fitrspc200, standardized = TRUE)
pec200
rspc200 <- pec200[pec200$label == "rsp", "std.all"]
rspc200
rspc2.200 <- rspc200^2
rspc2.200

# f2
f2.200 <- rspc2.200/(1 - yR2.200)
f2.200	# 0.01702596


fitrspc250 <- sem(rspc, sample.cov = iCov, sample.nobs = 250)
summary(fitrspc250, standardized = TRUE)
parameterEstimates(fitrspc250, standardized = TRUE)[label == "rsp"]
pec250 <- parameterEstimates(fitrspc250, standardized = TRUE)
pec250
rspc250 <- pec250[pec250$label == "rsp", "std.all"]
rspc250
rspc2.250 <- rspc250^2
rspc2.250

# f2
f2.250 <- rspc2.250/(1 - yR2.250)
f2.250	# 0.01702596

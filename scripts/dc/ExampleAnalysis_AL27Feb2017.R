###########
# WECARE example code for:
# 1) Burden test via logistic regression and conditional logistic
# 2) SKAT test (Burden, SKAT, SKAT-O)
# 3) ordinal regression (proportional odds)
# 4) multinomial trend test 
# updated 11/07/16
#################################

# note that I am assuming that the order of individuals is the same in all files

#----------------------------------------------------------------#
#                            Read data                           #
#----------------------------------------------------------------#

#setwd("C:/Google Drive/WECARE/DC/email_of_18Nov2016") # AL27Feb2017
setwd("/scratch/medgen/scripts/wecare_stat_01.17/scripts/dc") # AL27Feb2017

# Input individual-level genotype information
d.g <- read.table("genotypes_1.txt", header=T, sep="\t")
d.g <- t(d.g)

# input information on variants
d.v <- read.table("variants_1.txt", header=T, sep="\t")

# input individual-level phenotype information
d.1 <- read.table("phenotypes.txt", header=T, sep="\t")
d.1 <- data.frame(d.1, d.g)

d.2 <- read.table("WECARE.Exome.DemographicVariables.txt", header=T, sep="\t")

dim(d.1)
dim(d.2)

# Merge phenotype and genotype tables
d  <- merge(x=d.1, y=d.2, by.x="gwas_id", by.y="labid.x", all.x=T, all.y=F)
dim(d)

# Remove individuals with missed data
d[is.na(d$ID.x),]  # 4 individuals not in d.2 dataset????
d <- d[!is.na(d$ID.x),]
dim(d)

# Cleanup
rm(d.g, d.1, d.2) # AL27Feb2017

#----------------------------------------------------------------#
#                  Create varibles for analysis                  #
#----------------------------------------------------------------#
# Why Q, W and Z ? 

# we will adjust by the first 3 components
Q <- as.matrix(d[,c("Eigen_1", "Eigen_2","Eigen_3")])  

# covariates to adjust for... note that these should be jointly associated with
# the genetic variants and outcome to technically warrant adjustment
W <- as.matrix(d[,c("stage", "hormone.x", "family_history")]) 

# Covariates to compare to matching ?
Z <- d[,c("age_dx", "registry.x")] 

# Matching information
pairID <- d$setno.x

# OUtcome
Y <- d$cc

# Keep only individuals with complete covariates and outcome
# Impute missed values for covariates??

keep <- apply(cbind(Y, Q, W, pairID), 1, FUN=function(v) { ifelse(sum(is.na(v))>0, F, T) })  
sum(!keep) # 1

Q <- Q[keep,]
W <- W[keep,]
Z <- Z[keep,]
Y <- Y[keep]
pairID <- pairID[keep]

# Prepare ordered factor for trend analysis
Y.m <- as.factor(ifelse(d$cc==1, ifelse(d$family_history==1, 2, 1), 0))
Y.m <- ordered(Y.m, levels=c(0,1,2), labels=c("UBC", "CBC-", "CBC+"))
Y.m <- Y.m[keep]
Y.m[1:20]

#----------------------------------------------------------------#
#                   Analysis for ATM or CHEK2                    #
#----------------------------------------------------------------#

# --- Select variants for the gene --- #

v.list <- d.v$VarID[d.v$SYMBOL=="ATM"]
G <- as.matrix(d[, (names(d) %in% v.list)])
G <- G[keep,]

M <- ncol(G) # used later to generate uniform weights

# --- Impute missing values for variants --- #

# (we have a large number of individuals missing variants)
# note that there are better ways to impute like using haplotype infor etc. 
# However, because these are rare, there may be very little difference. 
# Also note that SKAT has its own default...

# Imputing variant-wise means
G <- apply(G, 2, function(v) { ifelse(is.na(v), mean(v, na.rm=T), v) })

# Check the imputation
#H <- matrix(c(1:20), nrow=4)
#NA -> H[1,5]
#NA -> H[2,3]
#NA -> H[4,4]
#NA -> H[3,1]
#H
#K <- apply(H, 2, function(v) { ifelse(is.na(v), mean(v, na.rm=T), v) })
#K
#rm(H,K)

# --- Prepare weights for variants --- # 

# Calculate allelic frequencies
p <- apply(G, 2, mean, na.rm=T)
p <- p/2 # Change to DC example: to account for two alleles per each genotype

# Check allelic frequencies calculations
#cnt <- apply(G, 2, sum, na.rm=T)
#notisna <- apply(G, 2, function(v){sum(!is.na(v))})
#cnt
#notisna
#af <- cnt/(2*notisna)
#p-af

# Calculate weights (se in SKAT default)
w=dbeta(p, 1, 25)

# --- Burden type analysis using logistic regression with a LRT --- #
# Unweighted, adjusted by eigenvectors (W) and covariates (Q)

# X is a sum of vatiants multiplyed by weights
# (each variant count is multiplyed by the individual variant's weight)

X <- G%*%rep(1,M) # uniform weights = 1 ( = no weighting)
reg.null <- glm(Y ~ W+Q, family=binomial(link=logit)) # default link for bionomial is logit
reg <- glm(Y ~ X + W+Q, family=binomial(link=logit))
X.coef <- summary(reg)$coef["X", ]
chi.stat.LRT <- 2*(logLik(reg) - logLik(reg.null))
P.LRT = 1-pchisq(chi.stat.LRT, df=1)
X.coef
P.LRT # 'log Lik.' 0.1985253 (df=8) : what does this mean??

# Check the matrix multiplication for weights
#a <- c(1:5)
#b <- rep(1,5)
#a
#b
#K %*% a
#K %*% b

# --- Burden type analysis using logistic regression with a LRT and weight --- #
# Weighted by default SKAT weight
# Adjusted only for age_dx and registry to directly compare to matching
# It seems that refage could be a better adjustement for matching 
# (although it is not essential) ... Why not adjusted by rs.time?

X <- G%*%w
dat <- data.frame(Y, X, Z)
reg.null <- glm(Y ~ age_dx + registry.x, family=binomial, data=dat)
reg <- glm(Y ~ X + age_dx + registry.x, family=binomial, data=dat)
X.coef <- summary(reg)$coef["X", ]
chi.stat.LRT = 2*(logLik(reg) - logLik(reg.null))
P.LRT = 1-pchisq(chi.stat.LRT, df=1)
X.coef
P.LRT

########### MATCHED analysis #################
##### Burden type analysis using conditional logistic regression with a LRT and weight
library(survival)
X <- G%*%w
reg <- clogit(Y ~ X + strata(pairID))
X.coef <- summary(reg)$coef["X", ]
P.LRT = anova(reg)["X","Pr(>|Chi|)"]
X.coef
P.LRT
# we need to compare the X.coef and P.LRT between matched and unmatched analysis



########### below are unmatched analyses

##### Burden type analysis using SKAT
library(SKAT)
skat.reg.null <- SKAT_Null_Model(Y ~ W+Q, out_type="D")
skat.reg <- SKAT(Z=G, obj=skat.reg.null, kernel="linear.weighted", method="davies", weights=w, is_dosage=T, r.corr=1)
P.SKAT <- skat.reg$p.value


##### SKAT variance component type test
library(SKAT)
skat.reg.null <- SKAT_Null_Model(Y ~ W+Q, out_type="D")
skat.reg <- SKAT(Z=G, obj=skat.reg.null, kernel="linear.weighted", method="davies", weights=w, is_dosage=T)
P.SKAT <- skat.reg$p.value

# SKAT-O
skat.reg <- SKAT(Z=G, obj=skat.reg.null, kernel="linear.weighted", method="optimal.adj", weights=w, is_dosage=T)
P.SKAT <- skat.reg$p.value

##### Weighted burden with proportional odds regression for CBC.FH+, CBC.FH-, UBC
library(MASS)
X <- G%*%w
m.reg.null <- polr(Y.m ~ W + Q, Hess=T)
m.reg <- polr(Y.m ~ X + W + Q, Hess=T)
chi.stat.LRT = 2*(logLik(m.reg) - logLik(m.reg.null))
P.LRT = 1-pchisq(chi.stat.LRT, df=1)
anova(m.reg.null, m.reg) # alternative test...

##### Weighted burden with multinomial regression for CBC.FH+, CBC.FH-, UBC
library(nnet)
X <- G%*%w
m.reg.null <- multinom(Y.m ~ W + Q, Hess=T)
m.reg <- multinom(Y.m ~ X + W + Q, Hess=T)
chi.stat.LRT = 2*(logLik(m.reg) - logLik(m.reg.null))
P.LRT = 1-pchisq(chi.stat.LRT, df=2)   # note that this test both contrasts (Wald tests can be used to test each contrast)
anova(m.reg.null, m.reg) # alternative test...



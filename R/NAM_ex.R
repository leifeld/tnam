#written by @caramnix on 10.12.21 

#also NEED to run  Estimation.R 

#Bayes TNAM simple example 
library(Matrix) 
library(matrixcalc) 
library(mvtnorm) 
library(rARPACK)
library(tmvtnorm)
library(tnam) # need this and then I've been overwriting the old functions with the updated functions 

############################################## DATA
data("knecht")
#prepare the dependent variable y
delinquency <- as.data.frame(delinquency)
rownames(delinquency) <- letters

delinquency1<- delinquency[,1]
delinquency1<- as.data.frame(delinquency1)
rownames(delinquency1) <- letters

# replace structural zeros (denoted as 10) and add row labels
friendship[[3]][friendship[[3]] == 10] <- NA
friendship[[4]][friendship[[4]] == 10] <- NA
for (i in 1:length(friendship)) {
  rownames(friendship[[i]]) <- letters
}
#get rid of NA's--> 0 
t2<- friendship[2]$t2
t2[is.na(t2)] <- 0 
friendship[2]$t2 <- t2
friendship1<- friendship[1]
#change to double-- 
friendship1$t1<- matrix(as.numeric(friendship[1]$t1), ncol=26, nrow=26)

# prepare the covariates sex and religion
sex <- demographics$sex
names(sex) <- letters
sex <- list(t1 = sex, t2 = sex, t3 = sex, t4 = sex)
sex1<- sex[1]
#religion <- demographics$religion
#names(religion) <- letters
#religion <- list(t1 = religion, t2 = religion, t3 = religion, 
#                 t4 = religion)
#NOTE: identical 

ethnicity <- demographics$ethnicity
names(ethnicity) <- letters
ethnicity <- list(t1 = ethnicity, t2 = ethnicity, t3 = ethnicity, t4 = ethnicity)
ethnicity1<- ethnicity[1]

############################################## FORMULA w/ simulated data from NAM_ex.R -- Dino
#10.12.21

set.seed(1)
W1=matrix(0,50,50); W1[sample(1:2500,250)]=1; diag(W1)=0; W1=W1/rowSums(W1)
W2=matrix(0,50,50); W2[sample(1:2500,250)]=1; diag(W2)=0; W2=W2/rowSums(W2)
X=matrix(c(rep(1,50),rnorm(50*3)),50,4)                   #Draw some covariates
rho=c(.2,.2); beta=rep(1,4); sigma=1                      #Set the model parameters
eps=rnorm(50,0,sigma)                                     #Draw the disturbances
y=c(solve(diag(50)-rho[1]*W1-rho[2]*W2)%*%(X%*%beta+eps)) #Draw the response vector y
mu.prior=c(0,0); Sigma.prior=50*diag(2)                   #Set non-informative priors

X[,1]

nam_formula <- y ~  
  covariate(X[,1], coefname = "intercept") +
  covariate(X[,2], coefname = "X2") +
  covariate(X[,3], coefname = "X3") +
  covariate(X[,4], coefname = "X4") + 
  W(W1) + 
  W(W2)

class(list(W1))
class(friendship1)
test1<- tnamdata_mod(nam_formula) ##data is in correct format now for betnam 

fitted<- tnam(formula= nam_formula) 

#look @ distributions of coeffeicents

class(fitted["beta"]$beta) ["covariate.intercept"]

intecerpt_samples<- vector()
X1_samples <- vector() 
X2_samples <- vector() 
X3_samples <- vector() 

for (i in 1:1000) {
  intecerpt_samples<- append(intecerpt_samples, fitted[["beta"]][i])
  X1_samples<- append(X1_samples, fitted[["beta"]][i+1000]) 
  X2_samples<- append(X2_samples, fitted[["beta"]][i+2000]) 
  X3_samples<- append(X3_samples, fitted[["beta"]][i+3000]) 
}


attach(mtcars)
par(mfrow=c(2,2))
#mtext("Bayesian Estimated TNAM", outer = TRUE, cex = 1.5)
hist(intecerpt_samples)
hist(X1_samples)
hist(X2_samples)
hist(X3_samples)
mtext("Bayesian Estimated TNAM",                
      side = 3,
      line = -1.5,
      outer = TRUE, 
      cex=1.5)


sigma_samples<- vector()
rho1_samples <- vector() 
rho2_samples <- vector() 

for (i in 1:1000) {
  sigma_samples<- append(sigma_samples, fitted[["sigma"]][i])
  rho1_samples<- append(rho1_samples, fitted[["rho"]][i]) 
  rho2_samples<- append(rho2_samples, fitted[["beta"]][i+1000]) 
}

par(mfrow=c(2,2))
hist(sigma_samples)
hist(rho1_samples)
hist(rho2_samples)
mtext("Bayesian Estimated TNAM",                   # Add main title
      side = 3,
      line = -1.5,
      outer = TRUE, 
      cex=1.5)

fitted[["rho"]]

##########################################

#now test with knecht data instead of simulated data
ethnicity$t5 <- rep(1,26) 
names(ethnicity$t5) <- letters
intercept_term<- ethnicity[5] # hacky way to do this-- sorry

formula4<- delinquency1 ~ #just one t for now...
  covariate(intercept_term, coefname = "intercept")+ 
  covariate(sex1, coefname = "sex") + 
  covariate(ethnicity1, coefname = "ethnicity") +
  W(friendship1) + 
  structsim(delinquency1, friendship1) 


dat4<- tnamdata_mod(formula4)

mu.prior4=c(0,0); Sigma.prior4=26*diag(2) #-- attempt to fix X.tilde error 

fitted<- tnam(formula= formula4, mu.prior= mu.prior4, Sigma.prior= Sigma.prior4)





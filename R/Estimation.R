## code from #https://github.com/DittrichD/BayesNAM, @DittrichD

#1. Dittrich D, Leenders RTAJ, Mulder J. Network Autocorrelation Modeling: Bayesian Techniques for Estimating and 
#Testing Multiple Network Autocorrelations. Sociological Methodology. 2020;50(1):168-214. doi:10.1177/0081175020913899

# commented by @caramnix 

# mh.rho.1.beta0 Metropolis Hastings Function for when W = 1 -- see mh.rho.R.beta0 for more thorough comments 
mh.rho.1.beta0=function(mu.cand,Sigma.cand,lb,ub,y,Wy,ev,lndet.curr,s2.curr,X.beta.tilde.curr,g,rss,val.curr,mu.prior,sigma.prior,Ay.curr) {
  val.cand=c(rtmvnorm(1,mean=mu.cand,sigma=Sigma.cand,lower=c(lb,-Inf),upper=c(ub,Inf))) # R+1-variate density 
  rho.cand=val.cand[1]
  beta0.cand=val.cand[2]
  Ay.cand=y-rho.cand*Wy
  lndet.cand=sum(log(1-rho.cand*ev))
  if(log(runif(1))<
     lndet.cand-lndet.curr
     -(1/(2*s2.curr))*(sum(Ay.cand**2)-2*beta0.cand*sum(Ay.cand)-2*sum(Ay.cand*X.beta.tilde.curr)+g*beta0.cand**2+2*beta0.cand*sum(X.beta.tilde.curr)+sum(X.beta.tilde.curr**2)-rss)
     -(1/(2*sigma.prior**2))*(((rho.cand-mu.prior)**2)-((val.curr[1]-mu.prior)**2))
     +dmvnorm(val.curr,mean=mu.cand,sigma=Sigma.cand,log=T)-dmvnorm(val.cand,mean=mu.cand,sigma=Sigma.cand,log=T)
  )
  {val=val.cand;lndet=lndet.cand;Ay=Ay.cand}
  else
  {val=val.curr;lndet=lndet.curr;Ay=Ay.curr}
  parlist=list(val,lndet,Ay)
  names(parlist)=c("value","lndet","Ay")
  return(parlist)
}
################################################################################
# Bayesian Estimation when len (W) =1 --> see nam.Bayes() for commented code 
nam.Bayes.1=function(y,X,W.list,mu.prior,Sigma.prior,N=100,burnin=0){
  W=as.matrix(W.list)
  sigma.prior=Sigma.prior
  if(is.null(y)>0){stop("The response vector 'y' must not be empty.")}
  if(is.vector(y)<1){stop("The response vector 'y' must be a vector.")}
  if(is.numeric(y)<1){stop("The response vector 'y' must be numeric.")}
  if(is.numeric(X)<1){stop("The covariance matrix 'X' must be numeric.")}
  if(is.matrix(X)<1){stop("The covariance matrix 'X' must be a matrix.")}
  if(is.numeric(W)<1){stop("The connectivity matrix must be numeric.")}
  if(abs(length(y)-nrow(X))>.5){stop("Number of observations in 'X' must match length of 'y'.")}
  g=length(y) 
  if(length(unique(g,nrow(W),ncol(W)))>1){stop("Order of the connectivity matrix must match length of 'y'.")}
  if(sum(is.vector(mu.prior)+length(mu.prior))<2){stop("The prior mean 'mu.prior' must be a scalar.")}
  if(is.numeric(mu.prior)<1){stop("The prior mean 'mu.prior' must be numeric.")}
  if(sum(is.vector(sigma.prior)+length(sigma.prior)+as.numeric(sigma.prior>0))<3){stop("The prior standard deviation 'sigma.prior' must be a positive scalar.")}
  if((N%%1==0)<1|N<1){stop("The number of desired posterior draws 'N' must be a positive integer.")}
  if((burnin%%1==0)<1|burnin<0){stop("The 'burnin' must be a non-negative integer.")}
  
  g=length(y)
  Id=diag(g)
  ones=rep(1,g)
  k=ncol(X)
  X.tilde=X[,-1]
  XtXi=solve(t(X)%*%X)
  XtXiXt=XtXi%*%t(X)
  M=Id-X%*%XtXiXt
  Wy=c(W%*%y) #this has to be subsetting Wi*yi
  yWones=sum(Wy)
  yWWy=sum(Wy**2)
  ev=Re(eigen(W)$values)
  tau2=sum(ev**2)
  lb=1/min(ev)
  ub=max(ev)
  Sigmai12h=yWones
  Sigmai21h=Sigmai12h
  Sigmai22h=g
  
  # Set starting values
  vals.start=lm(y~X.tilde) 
  rho.curr=0
  beta0.curr=vals.start$coefficients[1]
  beta.tilde.curr=vals.start$coefficients[-1]
  s2.curr=sigma(vals.start)**2
  
  A.curr=Id
  lndet.curr=0
  Ay.curr=y
  print("here 4")
  print(X.tilde)
  print(beta.tilde.curr)
  print("here 5")
  
  X.beta.tilde.curr=c(X.tilde%*%beta.tilde.curr)
  rss=sum(Ay.curr**2)-2*beta0.curr*sum(Ay.curr)-2*sum(Ay.curr*X.beta.tilde.curr)+g*beta0.curr**2+2*beta0.curr*sum(X.beta.tilde.curr)+sum(X.beta.tilde.curr**2)
  
  # Create output vectors and matrices
  rho.samples=NULL
  beta.samples=matrix(NA,nrow=burnin+N,ncol=k)  
  sigma.samples=NULL
  fitted.vals=list()  
  res=list()  
  rho.samples[1]=rho.curr
  beta.samples[1,]=c(beta0.curr,beta.tilde.curr)
  sigma.samples[1]=sqrt(s2.curr)
  fitted.vals[[1]]=c(solve(A.curr)%*%(X.beta.tilde.curr+rnorm(g,mean=0,sd=sqrt(s2.curr))))
  res[[1]]=y-fitted.vals[[1]]
  
  #MCMC step 
  for (i in 2:(burnin+N)){
    #loop is here 
    # Draw rho
    Sigmai11h=yWWy+s2.curr*tau2+s2.curr/sigma.prior**2 #page 204 
    Sigmaih=rbind(c(Sigmai11h,Sigmai12h),c(Sigmai21h,Sigmai22h))
    Sigmah=(1/(Sigmai11h*Sigmai22h-Sigmai12h**2))*rbind(c(Sigmai22h,-Sigmai12h),c(-Sigmai21h,Sigmai11h))
    Sigma.cand=s2.curr*Sigmah
    if(is.positive.definite(Sigma.cand)<1){Sigma.cand=nearPD(Sigma.cand)$mat}
    r.tilde=y-c(X.tilde%*%beta.tilde.curr)
    SS.tilde=sum(Wy*r.tilde)
    mu1=sum(r.tilde*yWones-SS.tilde-mu.prior*s2.curr/sigma.prior**2)/(yWones**2-g*Sigmai11h)
    if(mu1>ub|mu1<lb){
      a=(lb-mu1)/sqrt(Sigma.cand[1,1]);b=(ub-mu1)/sqrt(Sigma.cand[1,1]);z=pnorm(b)-pnorm(a)
      mu1=mu1+(sqrt(Sigma.cand[1,1])/z)*(dnorm(a)-dnorm(b))
    }
    mu2=(SS.tilde+mu.prior*s2.curr/sigma.prior**2-mu1*Sigmai11h)/yWones
    mu.cand=c(mu1,mu2)
    draw=mh.rho.1.beta0(mu.cand,Sigma.cand,lb,ub,y,Wy,ev,lndet.curr,s2.curr,X.beta.tilde.curr,g,rss,c(rho.curr,beta0.curr),mu.prior,sigma.prior,Ay.curr)
    rho.curr=draw$value[1]
    beta0.curr=draw$value[2]
    A.curr=Id-rho.curr*W
    lndet.curr=draw$lndet
    Ay.curr=draw$Ay
    rho.samples[i]=rho.curr
    beta.samples[i,1]=beta0.curr
    # Draw betatilde
    mu.beta=c(XtXiXt%*%Ay.curr)
    Sigma.beta=s2.curr*XtXi
    Sigma.beta11=Sigma.beta[1:1]
    Sigma.beta12=Sigma.beta[1,2:k]
    Sigma.beta22=Sigma.beta[2:k,2:k]
    mu.beta.tilde=mu.beta[-1]+Sigma.beta12*(beta0.curr-mu.beta[1])/Sigma.beta11
    Sigma.beta.tilde=Sigma.beta22-outer(Sigma.beta12,Sigma.beta12)/Sigma.beta11
    beta.tilde.curr=c(rmvnorm(1,mean=mu.beta.tilde,sigma=Sigma.beta.tilde))
    X.beta.tilde.curr=c(X.tilde%*%beta.tilde.curr)
    rss=sum(Ay.curr**2)-2*beta0.curr*sum(Ay.curr)-2*sum(Ay.curr*X.beta.tilde.curr)+g*beta0.curr**2+2*beta0.curr*sum(X.beta.tilde.curr)+sum(X.beta.tilde.curr**2)
    beta.samples[i,-1]=beta.tilde.curr
    # Draw Sigma2
    s2.curr=1/rgamma(1,shape=g/2,rate=rss/2)
    sigma.samples[i]=sqrt(s2.curr)
    fitted.vals[[i]]=c(solve(A.curr)%*%(c(X%*%beta.samples[i,])+rnorm(g,mean=0,sd=sqrt(s2.curr))))
    res[[i]]=y-fitted.vals[[i]]
  }
  
  if(burnin>0){
    rho.samples=rho.samples[-(1:burnin)]
    beta.samples=beta.samples[-(1:burnin),]
    sigma.samples=sigma.samples[-(1:burnin)]
    fitted.vals=fitted.vals[-(1:burnin)]
    res=res[-(1:burnin)]
  }
  if(is.null(colnames(X))>0){colnames(beta.samples)=paste("X",1:k,sep="")}else{colnames(beta.samples)=colnames(X)}
  theta.samples=cbind(rho.samples,beta.samples,sigma.samples)
  colnames(theta.samples)=c("rho1",colnames(beta.samples),"sigma")
  
  o=list()
  o$y=y
  o$X=X
  o$W.list=W.list
  o$mu.prior=mu.prior
  o$Sigma.prior=sigma.prior
  o$N=N
  o$burnin=burnin
  o$cvm=cov(theta.samples)
  o$rho=rho.samples
  o$beta=beta.samples
  o$sigma=sigma.samples
  o$rho.mean=mean(rho.samples)
  names(o$rho.mean)="rho1"
  o$beta.mean=apply(beta.samples,2,mean)
  o$sigma.mean=mean(sigma.samples)
  names(o$sigma.mean)="sigma"
  o$rho.se=sd(rho.samples)
  o$beta.se=apply(beta.samples,2,sd)
  o$sigma.se=sd(sigma.samples)
  o$fitted.values=fitted.vals
  o$residuals=res
  o$call <- match.call()
  class(o) <- c("nam.Bayes")
  o
}
################################################################################
# functions in.space.2.rsd,in.space.2, and in.space.R check if (rho, Beta) in theta_rho x R -- page 205, step 3 first bullet point. 

# CASE 1 when len(W)= 2 and is row standardized.
in.space.2.rsd=function(rho,W.list,R){ 
  inout=0
  if(sum(rho>=0)>1.5){
    if(sum(rho)<1){inout=1}
    else{inout=0}
  }
  else{
    mf=c(1,rho[2]/rho[1]); W=W.list[[1]]+mf[2]*W.list[[2]]
    if(rho[1]>0){
      rho1b=1/Re(eigs(W,1,which="LR",opts=list(retvec=FALSE))$values) 
      if(identical(rho1b,complex(0))>0){rho1b=1/max(Re(eigen(W,only.values=TRUE)$values))}
    }
    else{
      rho1b=1/Re(eigs(W,1,which="SR",opts=list(retvec=FALSE))$values) 
      if(identical(rho1b,complex(0))>0){rho1b=1/min(Re(eigen(W,only.values=TRUE)$values))}
    }
    if(sum(abs(rho)<=abs(mf*rho1b))==2){inout=1}
  }
  return(inout)
}

# CASE 2 when len(W)= 2 and is not row standardized.
in.space.2=function(rho,W.list,R){ 
  inout=0
  mf=c(1,rho[2]/rho[1]); W=W.list[[1]]+mf[2]*W.list[[2]]
  if(rho[1]>0){
    rho1b=1/Re(eigs(W,1,which="LR",opts=list(retvec=FALSE))$values) 
    if(identical(rho1b,complex(0))>0){rho1b=1/max(Re(eigen(W,only.values=TRUE)$values))}
  }
  else{
    rho1b=1/Re(eigs(W,1,which="SR",opts=list(retvec=FALSE))$values) 
    if(identical(rho1b,complex(0))>0){rho1b=1/min(Re(eigen(W,only.values=TRUE)$values))}
  }
  if(sum(abs(rho)<=abs(mf*rho1b))==2){inout=1}
  return(inout)
}

# CASE 3 when len(W) > 2 
in.space.R=function(rho,W.list,R){ 
  inout=0
  alpha=NULL; prod.cos=NULL; mf=NULL; ss=NULL
  alpha[1]=atan2(rho[2],rho[1])
  prod.cos[1:2]=rep(1,2)
  mf[1:2]=c(1,tan(alpha[1]))
  ss[1]=rho[1]**2; ss[2]=ss[1]+rho[2]**2
  W=mf[1]*W.list[[1]]+mf[2]*W.list[[2]]
  for(i in 3:R){
    ss[i]=ss[i-1]+rho[i]**2
    if(rho[i]>0){alpha[i-1]=acos(sqrt(ss[i-1]/ss[i]))}
    else{alpha[i-1]=-acos(sqrt(ss[i-1]/ss[i]))}
    prod.cos[i]=prod.cos[i-1]*cos(alpha[i-2])
    mf[i]=tan(alpha[i-1])/prod.cos[i]
    W=W+mf[i]*W.list[[i]]
  }
  if(rho[1]>0){
    rho1b=1/Re(eigs(W,1,which="LR",opts=list(retvec=FALSE))$values) 
    if(identical(rho1b,complex(0))>0){rho1b=1/max(Re(eigen(W,only.values=TRUE)$values))}
  } 
  else{
    rho1b=1/Re(eigs(W,1,which="SR",opts=list(retvec=FALSE))$values) 
    if(identical(rho1b,complex(0))>0){rho1b=1/min(Re(eigen(W,only.values=TRUE)$values))}
  }
  if(sum(abs(rho)<abs(mf*rho1b))==R){inout=1}
  return(inout)
}

# metropolis hastings function, rho0, beta1.  returns (rho,beta) pair, ln(det) and Ay ----- page 205 of Dittrich (2020)
# also see Equation 10, page 180 of Dittrich (2020) 
# "The conditional posterior in equation (10) does not have a well-known form and cannot be directly sampled from.
# Instead, we use the Metropolis- Hastings algorithm to generate draws from the conditional posterior for rho, Beta1 "

# "The density in equation (10) can be approximated by an (R+1)-variate normal candidate-generating density for rho, Beta1 that is tailored
# to the conditional posterior for rho, Beta1. 

# R = length(W.list) 
# g = length(y)

mh.rho.R.beta0=function(mu.cand,Sigma.cand,R,W.list,Id,y,Wy,lndet.curr,s2.curr,X.beta.tilde.curr,g,rss,val.curr,mu.prior,Sigma.prior,Ay.curr,FUN=in.space) {
  inout=0
  while(inout<1){
    val.cand=c(rmvnorm(1,mean=mu.cand,sigma=Sigma.cand)) # R+1-variate density 
    rho.cand=val.cand[-(R+1)]
    beta0.cand=val.cand[R+1]
    if(FUN(rho.cand,W.list,R)>0){inout=1} # check to see if in parameter space theta_rho x R (reals), pg 204. FUN =in.space is inputted into function and is one of (in.space.2.rsd, in.space.2, in.space.2.R) defined above. 
    # keep sampling until in parameter space, then you've found a rho and beta0 candidate 
  } 
  rho.cand.W.list=list()
  for(i in 1:R) { 
    rho.cand.W.list[[i]]=rho.cand[i]*W.list[[i]]
  }
  A.cand=Id-Reduce("+", rho.cand.W.list)
  lndet.cand=determinant(A.cand)$modulus[1]
  Ay.cand=y-colSums((rho.cand*t(Wy)))
  if(log(runif(1))< #accept based on uniform distribution.
     lndet.cand-lndet.curr  
     -(1/(2*s2.curr))*(sum(Ay.cand**2)-2*beta0.cand*sum(Ay.cand)-2*sum(Ay.cand*X.beta.tilde.curr)+g*beta0.cand**2+2*beta0.cand*sum(X.beta.tilde.curr)+sum(X.beta.tilde.curr**2)-rss)   
     +dmvnorm(rho.cand,mean=mu.prior,sigma=Sigma.prior,log=T)-dmvnorm(val.curr[-(R+1)],mean=mu.prior,sigma=Sigma.prior,log=T)
     +dmvnorm(val.curr,mean=mu.cand,sigma=Sigma.cand,log=T)-dmvnorm(val.cand,mean=mu.cand,sigma=Sigma.cand,log=T)
  )
  {val=val.cand;lndet=lndet.cand;Ay=Ay.cand}
  else
  {val=val.curr;lndet=lndet.curr;Ay=Ay.curr}
  parlist=list(val,lndet,Ay)
  names(parlist)=c("value","lndet","Ay")
  return(parlist)
}  

##############################################################################
# Bayesian estimation for len(W) > 1
nam.Bayes=function(y,X,W.list,mu.prior,Sigma.prior,N=100,burnin=0) {
  if(is.matrix(W.list)>0){
    nam.Bayes.1(y,X,W.list,mu.prior,Sigma.prior,N=N,burnin=burnin)
  } else {
    if(is.null(y)>0){stop("The response vector 'y' must not be empty.")}
    if(is.vector(y)<1){stop("The response vector 'y' must be a vector.")}
    if(is.numeric(y)<1){stop("The response vector 'y' must be numeric.")}
    if(is.numeric(X)<1){stop("The covariance matrix 'X' must be numeric.")}
    if(is.matrix(X)<1){stop("The covariance matrix 'X' must be a matrix.")}
    if(min(sapply(W.list,is.numeric))<1){stop("The connectivity matrices 'W' must be numeric.")}
    if(abs(length(y)-nrow(X))>.5){stop("Number of observations in 'X' must match length of 'y'.")}
    g=length(y) # number of nodes 
    R=length(W.list) # number of W matrices 
    if(length(unique(g,sapply(W.list,nrow),sapply(W.list,ncol)))>1){stop("Order(s) of 'W' must match length of 'y'.")}
    if(is.vector(mu.prior)<1){stop("The prior mean 'mu.prior' must be a vector.")}
    if(is.numeric(mu.prior)<1){stop("The prior mean 'mu.prior' must be numeric.")}
    if(is.matrix(Sigma.prior)<1){stop("The prior covariance matrix 'Sigma.prior' must be a matrix.")}
    if(is.numeric(Sigma.prior)<1){stop("The prior covariance matrix 'Sigma.prior' must be numeric.")}
    if(abs(length(mu.prior)-R)>.5){stop("Length of 'mu.prior' must match number of connectivity matrices 'W'.")}
    if(length(unique(R,nrow(Sigma.prior),ncol(Sigma.prior)))>1){stop("Dimensions of 'Sigma.prior' must match number of connectivity matrices 'W'.")}
    if(is.positive.definite(Sigma.prior)<1){stop("The prior covariance matrix 'Sigma.prior' must be positive-definite.")}
    if((N%%1==0)<1|N<1){stop("The number of desired posterior draws 'N' must be a positive integer.")}
    if((burnin%%1==0)<1|burnin<0){stop("The 'burnin' must be a non-negative integer.")}
    
    if(R<3){
      if(max(sapply(W.list,rowSums))==1){
        in.space=in.space.2.rsd
      }else{
        in.space=in.space.2
      }
    }else{
      in.space=in.space.R
    }
    Id=diag(g) # identity matrix
    ones=rep(1,g)
    k=ncol(X) 
    X.tilde=X[,-1]
    XtXi=solve(t(X)%*%X) # transpose of X matrix mult. by X
    XtXiXt=XtXi%*%t(X) # transpose of X matrix mult. by X matrix mult. transpose X
    M=Id-X%*%XtXiXt  # from equation 15 ----- page 185
    yMy=sum(c(M%*%y)**2)
    Wy=matrix(NA,g,R)
    tr.WW=matrix(NA,R,R)
    for(i in 1:R){
      Wy[,i]=c(W.list[[i]]%*%y)
    }
    sumWy=colSums(Wy)
    yWWy=t(Wy)%*%Wy
    for(i in 1:R){
      for(j in i:R){
        tr.WW[i,j]=sum(t(W.list[[i]])*W.list[[j]]); tr.WW[j,i]=tr.WW[i,j] # discussion of this tr() found in A.1. Approximating the Conditional Posterior, page 203. 
      }
    }
    Sigmai.prior=solve(Sigma.prior)
    Sigmaih1=rbind(cbind(yWWy,sumWy),c(sumWy,g))
    Sigmaih2=rbind(cbind(tr.WW+Sigmai.prior,rep(0,R)),rep(0,R+1))
    
    vals.start=lm(y~X.tilde) 
    rho.curr=rep(0,R)
    beta0.curr=vals.start$coefficients[1]
    beta.tilde.curr=vals.start$coefficients[-1]
    s2.curr=sigma(vals.start)**2
    
    rho.curr.W.list=list()
    for(j in 1:R) {
      rho.curr.W.list[[j]]=rho.curr[j]*W.list[[j]]
    }
    A.curr=Id-Reduce("+",rho.curr.W.list)
    lndet.curr=determinant(A.curr)$modulus[1]
    Ay.curr=y-colSums((rho.curr*t(Wy)))
    #print(X.tilde)
    #print(beta.tilde.curr)
    X.beta.tilde.curr=c(X.tilde%*%beta.tilde.curr)
    rss=sum(Ay.curr**2)-2*beta0.curr*sum(Ay.curr)-2*sum(Ay.curr*X.beta.tilde.curr)+g*beta0.curr**2+2*beta0.curr*sum(X.beta.tilde.curr)+sum(X.beta.tilde.curr**2)
    
    rho.samples=matrix(NA,nrow=burnin+N,ncol=R)  
    beta.samples=matrix(NA,nrow=burnin+N,ncol=k)  
    sigma.samples=NULL
    fitted.vals=list()  
    res=list()  
    rho.samples[1,]=rho.curr
    beta.samples[1,]=c(beta0.curr,beta.tilde.curr)
    sigma.samples[1]=sqrt(s2.curr)
    fitted.vals[[1]]=c(solve(A.curr)%*%(X.beta.tilde.curr+rnorm(g,mean=0,sd=sqrt(s2.curr))))
    res[[1]]=y-fitted.vals[[1]]
    
    for (i in 2:(burnin+N)){ #MCMC
      # Draw (rho, beta) pair 
      Sigmai.cand=Sigmaih1/s2.curr+Sigmaih2
      if(is.positive.definite(Sigmai.cand)<1){Sigmai.cand=nearPD(Sigmai.cand)$mat} # if it is NOT positive definite, then compute the nearest positive definite matrix to an approximate one
      Sigma.cand=solve(Sigmai.cand)
      r.tilde=y-c(X.tilde%*%beta.tilde.curr)
      z.rho=c(t(Wy)%*%r.tilde)/s2.curr+c(Sigmai.prior%*%mu.prior)
      z.beta0=sum(r.tilde)/s2.curr
      z=c(z.rho,z.beta0)
      mu.cand=c(Sigma.cand%*%z)
      draw=mh.rho.R.beta0(mu.cand,Sigma.cand,R,W.list,Id,y,Wy,lndet.curr,s2.curr,X.beta.tilde.curr,g,rss,c(rho.curr,beta0.curr),mu.prior,Sigma.prior,Ay.curr,FUN=in.space)  
      rho.curr=draw$value[-(R+1)]
      beta0.curr=draw$value[R+1]
      for(j in 1:R) {
        rho.curr.W.list[[j]]=rho.curr[j]*W.list[[j]]
      }
      A.curr=Id-Reduce("+",rho.curr.W.list)
      lndet.curr=draw$lndet
      Ay.curr=draw$Ay
      rho.samples[i,]=rho.curr
      beta.samples[i,1]=beta0.curr
      
      # Draw beta.tilde
      # #drawn from (k+1) variate normal-- pg 205, step 5 of A.1.2. Sampling Algorithm& section 3.2 Equation 12 
      mu.beta=c(XtXiXt%*%Ay.curr)
      Sigma.beta=s2.curr*XtXi
      Sigma.beta11=Sigma.beta[1:1]
      Sigma.beta12=Sigma.beta[1,2:k]
      Sigma.beta22=Sigma.beta[2:k,2:k]
      mu.beta.tilde=mu.beta[-1]+Sigma.beta12*(beta0.curr-mu.beta[1])/Sigma.beta11
      Sigma.beta.tilde=Sigma.beta22-outer(Sigma.beta12,Sigma.beta12)/Sigma.beta11
      beta.tilde.curr=c(rmvnorm(1,mean=mu.beta.tilde,sigma=Sigma.beta.tilde))
      X.beta.tilde.curr=c(X.tilde%*%beta.tilde.curr)
      # residual sum of squares 
      rss=sum(Ay.curr**2)-2*beta0.curr*sum(Ay.curr)-2*sum(Ay.curr*X.beta.tilde.curr)+g*beta0.curr**2+2*beta0.curr*sum(X.beta.tilde.curr)+sum(X.beta.tilde.curr**2)
      beta.samples[i,-1]=beta.tilde.curr
      
      # Draw Sigma2
      #  #drawn from inverse gamma -- pg 205, step 4 of A.1.2. Sampling Algorithm (Dittrich 2020) & section 3.2 Equation 11 
      s2.curr=1/rgamma(1,shape=g/2,rate=rss/2)
      sigma.samples[i]=sqrt(s2.curr) # Equation 6 (?) 
      fitted.vals[[i]]=c(solve(A.curr)%*%(c(X%*%beta.samples[i,])+rnorm(g,mean=0,sd=sqrt(s2.curr))))
      res[[i]]=y-fitted.vals[[i]]
    }
    
    if(burnin>0){
      rho.samples=rho.samples[-(1:burnin),]
      beta.samples=beta.samples[-(1:burnin),]
      sigma.samples=sigma.samples[-(1:burnin)]
      fitted.vals=fitted.vals[-(1:burnin)]
      res=res[-(1:burnin)]
    }
    if(is.null(colnames(X))>0){colnames(beta.samples)=paste("X",1:k,sep="")}else{colnames(beta.samples)=colnames(X)}
    colnames(rho.samples)=paste("rho",1:R,sep="")
    theta.samples=cbind(rho.samples,beta.samples,sigma.samples)
    colnames(theta.samples)=c(colnames(rho.samples),colnames(beta.samples),"sigma")
    
    o=list()
    o$y=y
    o$X=X
    o$W.list=W.list
    o$mu.prior=mu.prior
    o$Sigma.prior=Sigma.prior
    o$N=N
    o$burnin=burnin
    o$cvm=cov(theta.samples)
    o$rho=rho.samples
    o$beta=beta.samples
    o$sigma=sigma.samples
    o$rho.mean=apply(rho.samples,2,mean)
    o$beta.mean=apply(beta.samples,2,mean)
    o$sigma.mean=mean(sigma.samples)
    names(o$sigma.mean)="sigma"
    o$rho.se=apply(rho.samples,2,sd)
    o$beta.se=apply(beta.samples,2,sd)
    o$sigma.se=sd(sigma.samples)
    o$fitted.values=fitted.vals
    o$residuals=res
    o$call <- match.call()
    class(o) <- c("nam.Bayes")
    o
  }
}
################################################################################

print.nam.Bayes=function(x, digits = max(3, getOption("digits") - 4), ...){
  coefs.mean=c(x$rho.mean,x$beta.mean)
  cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")       
  cat("Coefficients:\n")
  print.default(format(coefs.mean, digits = digits), print.gap = 2, quote = FALSE)
  cat("\n")
}
################################################################################

summary.nam.Bayes=function(object, ...){
  coefs=cbind(object$rho,object$beta)
  coefs.mean=c(object$rho.mean,object$beta.mean)
  coefs.sd=c(object$rho.se,object$beta.se)
  coefs.0=colSums(apply(coefs,2,function(a){a>0}))/object$N
  tab=cbind(coefs.mean,coefs.sd,coefs.0)
  colnames(tab)=c("Estimate","Std. Error","Pr(>0)")
  cat("\nCall:\n", deparse(object$call), "\n\n", sep = "")        
  cat("Coefficients:\n")
  print.default(format(tab, digits = 4), print.gap = 2, quote = FALSE)
  tab=cbind(object$sigma.mean,object$sigma.se)
  rownames(tab)="sigma"
  colnames(tab)=c("Estimate","Std. Error")
  cat("\n---------\n")
  print.default(format(tab, digits = 4), print.gap = 2, quote = FALSE)
}

################################################################################

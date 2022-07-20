# n   : number of regions
# rho : spatial correlation
# nn  : number of replications
# event rate is a linear function of a covariate x: mu=b0+b1*x
# sd_err: error standard deviation
# xMatrix: nn*n by 1 vector of covariates
# X   : nn*n by 2 matrix of covariates and intercept
# nb4rt: neighbors list (Pebesma, Ch 11.5.2, p 335)
simulateSAR_NBCovarGiven=function(n,rho,nn,b0,b1,sd_err,xMatrix,nb4rt){

# represent neighbors as a matrix 
W1=nb2mat(nb4rt)
# replicated neighbors matrix
Wn=W1
i=1
Wnew=NA
while (i<nn){
  Wnew=bdiag(Wn,W1)
  Wn=Wnew
  i=i+1
}
ImBinv=solve(diag(n)-rho*W1)

#### simulate nn event count observations at each of the n regions ####
# initialize data vector and coord matrix
y=NA
corr_err=NA
for (i in 1:nn){
  # independent variables
  x=xMatrix[(n*(i-1)+1):(n*i)]
  mu_x=b0+b1*x
  uncorr_x <- rnorm(n,mean=0,sd=sd_err)
  corr_err[(n*(i-1)+1):(n*i)] <- ImBinv %*% uncorr_x
  y[(n*(i-1)+1):(n*i)]=mu_x+corr_err[(n*(i-1)+1):(n*i)]
}

X=matrix(c(rep(1,n*nn),xMatrix),n*nn,2)
m=length(y)
yvec=matrix(y,n*nn,1)


# return results
result=list(y=yvec,X=X,Wn=Wn,W=W1,xMatrix=xMatrix)
return(result)

}
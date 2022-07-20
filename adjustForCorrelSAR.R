# n: number of regions
# rho: spatial correlation
# nn: number of replications
# event rate is a linear function of a covariate x: mu=b0+b1*x
# sd_err: error standard deviation
adjustForCorrelSAR=function(n,rho,nn,b0,b1,sd_err,y,xMatrix,W1,xL,xU,flag_plot,flag_RiskAdjust){


  # regression and correlation adjusted data
  yAdjMatrix=matrix(NA,nrow=nn,ncol=n)
  xMM=matrix(xMatrix,nrow=nn,byrow=TRUE)
  yMM=matrix(y,nrow=nn,byrow=TRUE)
  ImB=diag(n)-rho*W1
  for (i in 1:nn){
    if(flag_RiskAdjust){
      yAdjMatrix[i,]=(1/as.numeric(sd_err))*ImB%*%(yMM[i,]-b0-b1*xMM[i,])
      }
    else{
      yAdjMatrix[i,]=(1/as.numeric(sd_err))*ImB%*%(yMM[i,]-b0-b1*(xL+xU)*0.5)  # mean of U(5,10) is 7.5
      }
  }
  DistX=as.matrix(dist(X0,method='euclidian',diag=TRUE))
  SigYAdj=as.matrix(cor(yAdjMatrix))
  yAdj=matrix(t(yAdjMatrix),nrow=nn*n,byrow=TRUE)
  
  if (flag_plot){
    #plot(DistX[upper.tri(DistX)],SigYAdj[upper.tri(SigYAdj)],xlab='dist',ylab='correl')
    #title('correl of adjusted data vs distance ')
    # plot y vs risk factor
    par(mfrow=c(1,2))
    plot(y~xMatrix)
    abline(lm(y~xMatrix))
    plot(yAdj~xMatrix)
    abline(lm(yAdj~xMatrix))
    
  }

  
# return results
result=list(yAdj=yAdj)
return(result)

}
getLikelihoodSAR<-function(rho,y,W,X){
  m=length(y)
  B=rho*W
  A=diag(m)-B
  
  beta_hat=solve(t(A%*%X)%*%(A%*%X),t(A%*%X)%*%(A%*%y))
  e_vec=y-X%*%beta_hat
  sig2hat=(1/m)*t(e_vec)%*%(t(A)%*%A)%*%e_vec

  log_likel=as.numeric(log(det(A))-m/2*log(sig2hat))
  return(-log_likel)
}
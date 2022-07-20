# xMatrix: nn*n by 1 vector of covariates
# nb4rt: neighbors list (Pebesma, Ch 11.5.2, p 335)

# average  run length of spatial cusum designed to detect trends or step shifts under trend shifts
get_ARL_SCUSUMCorrel_NeighborCovarGiven<-function(r, X0, b0,b1,rho,sd_err,
                                     l0,delt_c, delt, h_s, time_inj, n_MC,ind_inject,
                                     flag_plot,flag_moran,nb4rt,xMatrix,flagRiskCorrAdjust){

nr=length(r)
ns =dim(X0)[1]
RL=NA
tau_est=NA
center_cluster=matrix(NA,n_MC,2)
rad_cluster=NA
pvalMoranAdj=NA
pvalMoran=NA
n_inj=length(ind_inject)
prmEstVec=matrix(NA,n_MC,3) # cols 1,2,3 b0hat, b1hat,sd_err_hat
j=1
while (j <=n_MC){
  S_prev=matrix(0,nr,ns)
  S_max=0
  Smax_vec=0
  i=1
  infectCount=NA
  tau=matrix(0,nr,ns)
  
  # monitor- online measurement size nnOL=1
  nnOL=1
  while (S_max<h_s| i<=time_inj){
    i=i+1
    res1=simulateSAR_NBCovarGiven(ns,rho,nnOL,b0,b1,sd_err,xMatrix,nb4rt)
    y=res1$y
    Wn=res1$Wn
    W=res1$W
    if (i>time_inj & delt>0){
      b01=b0+delt*sd_err
      # find a way to use the previous xMatrix here
      y2=simulateSAR_NBCovarGiven(ns,rho,nnOL,b01,b1,sd_err,xMatrix,nb4rt)$y
      y[ind_inject,]=y2[ind_inject,]  
    }
    k=i-1

    if(flagRiskCorrAdjust){
      flag_RiskAdjust=TRUE;xL=NA;xU=NA
      yAdj=adjustForCorrelSAR(ns,rho,nnOL,b0,b1,sd_err,y,xMatrix,W,xL,xU,flag_plot,flag_RiskAdjust)$yAdj
    }else{
      yAdj=(y-mean(y))/sd(y)
    }
    
    
    # sustained shift cusum for adjusted data
    res_S=get_spaceTimeCUSUMNormalStandard(yAdj,X0,S_prev,delt_c/2,r)

    if (flag_moran){
      # moran test on residuals
      pvalMoranAdj=cbind(pvalMoranAdj,moran.test(yAdj, listw=mat2listw(Wn))$p.value)
      # moran test data
      pvalMoran=cbind(pvalMoran,moran.test(y, listw=mat2listw(Wn))$p.value)
    }
    
    # reset tau (change point estimator for each subregion) to k if cusum = 0
    S_prev=res_S$S_new
    tau[which(S_prev==0)]=k    
    # read cusum stat
    S_max=res_S$S_max
    #S_prev=res_S$S_new
    Smax_vec=c(Smax_vec,S_max)
    # for plotting, count in infected regions over time
    infectCount[i]=sum(y[ind_inject])
    #print(i)
    #print(c(S_max,h_s))
  }
  if (i>time_inj){
    print(paste("MC run ",j,sep=""))
    RL[j]=i-time_inj
    # cluster center 
    center_cluster[j,]=as.numeric(X0[res_S$region_max,])
    # cluster radius
    rad_cluster[j]=as.numeric(r[res_S$radius_max])
    # change point estimator
    ind_max=which(res_S$S_new==S_max)
    tau_est[j]=tau[ind_max[1]]
    j=j+1
  }
}

ARL=mean(RL)
seRL=sqrt(var(RL)/n_MC)

# bias and mse of change point estimates
bias=sum(tau_est-time_inj)/n_MC
mse=1/n_MC*(t(tau_est-time_inj)%*%(tau_est-time_inj))


# return results
result=list(ARL=ARL,seRL=seRL,tau_est=tau_est,bias=bias,
            mse=mse,Smax_vec=Smax_vec,infectCount=infectCount,
            pvalMoranAdj=pvalMoranAdj,pvalMoran=pvalMoran,prmEstVec=prmEstVec)
return(result)

}
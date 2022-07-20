# Space time circular CUSUM for sustained shift in iid standard normal data
# Author:   Arda Vanli
# Version:  3/8/2021

# y~N(0,1) when in control, y~N(2k,1) when out of control
# find space time CUSUM for time t
# y is (n by 1) vector of counts where n is number of regions
# X0 is (n by 2) vector of region coordinates
get_spaceTimeCUSUMNormalStandard<-function(y,X0,S_prev,k,r){
  # number of space points
  ns =dim(X0)[1]
  
  Xs=X0
  colnames(Xs)=c("x1","x2")

  # 3 radius values, 0, 0.167 and 0.333 (for 6 by 6 grid)
  #r=seq(0,2*step_size,by=step_size)
  nr=length(r)

  S=array(0,dim=c(nr,ns))
  
  # number of subregions scanned for each S[i,j]
  nc_mat=matrix(NA,nr,ns)
  
  # total number of events
  Nts=sum(y)
  for (i in 1:nr) {
    #print(paste('radius no = ',i,sep=""))
    for (j in 1:ns){
      fl_withinCirc=logical(length=ns)
      # find rows within radius of r
      for (j1 in 1:ns){
        distR=sqrt((Xs[j,1]-X0[j1,1])^2+(Xs[j,2]-X0[j1,2])^2)
        #fl_withinCirc[j1]=distR<r[i]&distR!=0
        fl_withinCirc[j1]=distR<=(r[i]+0.001)
      }          
      # number of events within radius r
      nZts=sum(y[fl_withinCirc])
      # number of regions within circle
      n_c=sum(fl_withinCirc)
      nc_mat[i,j]=n_c
      # new cusum for this circle
      S[i,j]=max(0,S_prev[i,j]+(nZts-n_c*k))
    }
  }
  
  # max over all regions and radii 
  S_max=max(S)
  ind1=which(S==S_max,arr.ind = TRUE)
  
  # second and third largest elements
  S_max2=S[which(rank(-S,ties.method = 'first')==2)]  # may have to use larger numbers to avoid overlapping clusters
  S_max3=S[which(rank(-S,ties.method = 'first')==3)]
  ind2=which(S==S_max2,arr.ind = TRUE)
  ind3=which(S==S_max3,arr.ind = TRUE)
  # number of subregions scanned for S_max
  nc_max=nc_mat[ind1[1]]
  #find radius and region number of the max 
  radius_max=ind1[1,1]
  region_max=ind1[1,2] 
  
  # radius and region number of the second and third largest  
  RadiusRegionTop3=matrix(c(ind1[1,1],ind2[1,1],ind3[1,1],ind1[1,2],ind2[1,2],ind3[1,2]),byrow=TRUE,nrow=2)
  
  # return results
  result=list(S_max_t=S_max, S_new=S, 
              radius_max=radius_max,region_max=region_max,nc_max=nc_max,RadiusRegionTop3=RadiusRegionTop3)
  return(result)
  
}

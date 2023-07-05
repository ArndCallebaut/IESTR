

print_tensor = function(tensor,height_map){
  rvb_tensor[is.na(height_map)]=1
  rvb_tensor2 = round(rvb_tensor*255)
  raster_RGB = stack(raster(rvb_tensor2[,,1]),raster(rvb_tensor2[,,2]),raster(rvb_tensor2[,,3]))
  plt1 <- plotRGB(flip(raster_RGB,1))
}


# 0-preprocessing
cost = cost1

# 0.1 : evaluate everything
evaluation_sites = matrix(0,N_cycles,length(gsc[,1]))
for (i in 1:N_cycles){
  cm_loc = cm[[i]]
  for (j in 1:length(gsc[,1])){
    evaluation_sites[i,j] = sum(cm_loc[,j])
  }
}

# 0.2 : evaluate everything
library(ramify)
max_each_site = apply(evaluation_sites,FUN=max,2)
ind_each_site = (argmax(t(evaluation_sites)))

# diag_cms = matrix(0,N_cycles,length(gsc[,1]))
# for (i in 1:N_cycles){
#   cm_loc = cm[[i]]
#   diag_cms[i,] = apply(cm_loc,2,sum)
# }
# 
# maxi_diags = apply(diag_cms,FUN=max,2)


efficienty3 = max_each_site / cost[gsc[,c(1,2)]+1]
eff3 = efficienty3 / max(efficienty3)


val_matrix = array(0,c(nrow,ncol))
ind_matrix = array(0,c(nrow,ncol))
for (i in 1:length(efficienty3)){
  val_tensor[gsc[i,1]+1,gsc[i,2]+1]= eff3[i]
  ind_tensor[gsc[i,1]+1,gsc[i,2]+1]= ind_each_site[i]
}


rvb_tensor = array(0.9,c(nrow,ncol,3))
rvb_tensor[is.na(height_map)]=1

rvb_tensor[,,1] = val_tensor


external_coord = which(is.na(height_map[2:(nrow-1),2:(ncol-1)]),arr.ind=TRUE)
is_on_border = apply(external_coord,1,FUN = function(x) 0!=sum(!is.na(height_map[(x[1]):(x[1]+2),(x[2]):(x[2]+2)])) )
coord_border = external_coord[is_on_border,]
rvb_tensor[,,1][cbind(coord_border[,1]+1,coord_border[,2]+1)] = 0.2
rvb_tensor[,,2][cbind(coord_border[,1]+1,coord_border[,2]+1)]= 0.2
rvb_tensor[,,3][cbind(coord_border[,1]+1,coord_border[,2]+1)] = 0.2 

# 1 - Structure of the RVB tensor (use : height maps, nrow, ncol)

print_tensor(tensor,height_map)


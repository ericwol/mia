rasterize.voi <- function(voi,def=NA){
# In: impro.R
# Input:
#   voi is in list-format (u,w,x,y,z) as in AMIDE
# Value:
#   returns a raster form of the voi
# The output raster is filled with 0's where values 
# are not found in voi.
	xg = sort(unique(voi[,3]))
	yg = sort(unique(voi[,4]))
	zg = sort(unique(voi[,5]))
	NR = length(yg)
	NC = length(xg)
	NS = length(zg)
	ras = array(def,dim=c(NR,NC,NS))
	for(v in 1:nrow(voi)){
		i = which(yg==voi[v,4])
		j = which(xg==voi[v,3])
		k = which(zg==voi[v,5])
		ras[i,j,k] = voi[v,1]
	}
	return(ras)
}

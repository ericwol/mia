#' @export
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

#' @export
scan.views <- function(tsv,mlab="",LWD=.8){
# Displays transverse, sagittal and coronal views of input 3D scan tsv.
# mlab is a patch string for more detailed plot main title.
# LWD indicates width of crosshairs (default .8).
	rtsv = rasterize.voi(tsv)
	dd = dim(rtsv)
	s = floor(dd/2)	
	graphics::image(t(rtsv[,,s[3]]),axes=FALSE,col=grDevices::grey(c(0:255)/255),
			main=paste(mlab,"(Transverse)")) # mid-volume transverse slice
	graphics::abline(h=s[1]/dd[1],col=4,lwd=LWD,lty=1)
	graphics::abline(v=s[2]/dd[2],col=4,lwd=LWD,lty=1)
	graphics::image(t(rtsv[,s[2],]),axes=FALSE,col=grDevices::grey(c(0:255)/255),
			main=paste(mlab,"(Sagittal)")) # mid-volume sagittal slice
	graphics::abline(h=s[1]/dd[1],col=4,lwd=LWD,lty=1)
	graphics::abline(v=s[3]/dd[3],col=4,lwd=LWD,lty=1)
	graphics::image(t(rtsv[s[1],,]),axes=FALSE,col=grDevices::grey(c(0:255)/255),
			main=paste(mlab,"(Coronal)")) # mid-volume coronal slice
	graphics::abline(h=s[2]/dd[2],col=4,lwd=LWD,lty=1)
	graphics::abline(v=s[3]/dd[3],col=4,lwd=LWD,lty=1)
}


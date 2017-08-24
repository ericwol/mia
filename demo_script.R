#! /usr/bin/Rscript --vanilla
# while( dev.next()>1 ){ dev.off() } # close all
# rm(list=ls(pos=.GlobalEnv), pos=.GlobalEnv) # clear all

######################################################################################
# Note: uptake values are expressed in SUV scale for these test datasets
# CaseA.ClinicalSUVmax = 2.1  # Female
# CaseB.ClinicalSUVmax = 12.4 # Male
######################################################################################

source("impro/impro.R") # use our own image/raster-handling functions for display
source("ellmod/hetE.R")  # implementation of (ellipsoidal) structural analysis tools

# ------------------------------------------------------------------------ 
# initialisations

ptids = paste("Case_",c("A","B"),sep="")
cqs = c(.25,.30,.40,.50,.60,.70,.75,.80,.85,.90,.95,1) # quantiles for gradients
Q = length(cqs)

het0 = het1 = het0.b = het1.b = numeric(2)
grad.b = rgrad.b = wgrad.b = matrix(nr=2,nc=Q)
grad = rgrad = wgrad = matrix(nr=2,nc=Q)

# ------------------------------------------------------------------------ 
# display the two demo cases

par(mfrow=c(2,3),font=2,font.lab=2,font.axis=2)
for(i in 1:2){ 
	ptid = ptids[i]
	
	# bounding box input VOI
	tsv = matrix(scan(file=paste("ellmod/data/",ptid,"_tumor.tsv",sep=""),skip=8),nc=5,byrow=TRUE)
	tsv[which(tsv[,1]<0),1] = abs(tsv[which(tsv[,1]<0),1])	
	rtsv = rasterize.voi(tsv)
	dd = dim(rtsv)
	s = floor(dd/2)	
	LWD = .8 # width of crosshairs
	image(t(rtsv[,,s[3]]),axes=FALSE,col=gray(c(0:255)/255),
			main=paste(ptid,"(Transverse)")) # mid-volume transverse slice
	abline(h=s[1]/dd[1],col=4,lwd=LWD,lty=1)
	abline(v=s[2]/dd[2],col=4,lwd=LWD,lty=1)
	image(t(rtsv[,s[2],]),axes=FALSE,col=gray(c(0:255)/255),
			main=paste(ptid,"(Sagittal)")) # mid-volume sagittal slice
	abline(h=s[1]/dd[1],col=4,lwd=LWD,lty=1)
	abline(v=s[3]/dd[3],col=4,lwd=LWD,lty=1)
	image(t(rtsv[s[1],,]),axes=FALSE,col=gray(c(0:255)/255),
			main=paste(ptid,"(Coronal)")) # mid-volume coronal slice
	abline(h=s[2]/dd[2],col=4,lwd=LWD,lty=1)
	abline(v=s[3]/dd[3],col=4,lwd=LWD,lty=1)
}

# ------------------------------------------------------------------------ 
# analyse the two demo cases
for(i in 1:2){ 
	ptid = ptids[i]
	
	# bounding box input VOI
	tsv = matrix(scan(file=paste("data/",ptid,"_tumor.tsv",sep=""),skip=8),nc=5,byrow=TRUE)
	tsv[which(tsv[,1]<0),1] = abs(tsv[which(tsv[,1]<0),1])	
	rtsv = rasterize.voi(tsv)
	
	# retrieve voxel dimensions
	tsvz = sort(unique(tsv[,5]))
	tsvx = sort(unique(tsv[,3]))
	tsvy = sort(unique(tsv[,4]))
	dimz = round(sort(unique(diff(tsvz)))[1],2)
	dimx = round(sort(unique(diff(tsvx)))[1],2)
	dimy = round(sort(unique(diff(tsvy)))[1],2)

	# re-arrange VOI columns
	troi = tsv[,c(2:5,1)]
	
	# run ellipsoidal analyses
	outH = hetE(troi)
	# final het + gradient quants:
	outSQ = struct.quant(outH,zdim=dimz,xdim=dimx,ddf=10)
	het0[i] = outH$het0
	het1[i] = outH$het1
	dg.r = outSQ$gradients
	dg = outSQ$normalized.gradients
	dg.w = outSQ$w.gradients
	rgrad[i,] = quantile(dg.r,cqs,na.rm=TRUE)
	grad[i,] = quantile(dg,cqs,na.rm=TRUE)
	wgrad[i,] = quantile(dg.w,cqs,na.rm=TRUE) 
	# plot all this:
	if(.Platform$OS.type=="unix"){ quartz() } else { x11() }
	plot.profile(outH,outSQ,sdf=10,mpatch=paste(", ",ptid))
}

het.df = data.frame(het0,het1)
row.names(het.df) = ptids

gradients.df = data.frame(quantile=cqs,t(rgrad),t(grad),t(wgrad))
names(gradients.df) = c("quantile",paste("raw",ptids),paste("normalized",ptids),paste("weighted",ptids))

"Heterogeneity summaries:"
round(het.df,2)

"Gradient summaries:"
round(gradients.df,2)
#' @import utils
#' @import stats
#' @import graphics
#' @import grDevices
#' @export
demo.ellmod <- function(){
# Demonstrates the use of hetE output. 
# Runs heterogeneity and gradient analysis for two sample datasets
# Uses "impro.R", our own image/raster-handling functions for display
# Uses "hetE.R", an implementation of (ellipsoidal) structural analysis tools
# Note: uptake values are expressed in SUV scale for these test datasets
# CaseA.ClinicalSUVmax = 2.7
# CaseB.ClinicalSUVmax = 7.8
#
	# ------------------------------------------------------------------------ 
	# internal functions
	analyse.and.plot <- function(tsv,ptid){
	# Internal function; analyses input VOI tsv and plots uptake profiles.
	# Value:
	# 	a list containing outputs of hetE() and struct.quant()
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
		outSQ = struct.quant(outH,zdim=dimz,xdim=dimx)
		het0 = outH$het0
		het1 = outH$het1
		dg.r = outSQ$gradients
		dg = outSQ$normalized.gradients
		dg.w = outSQ$w.gradients
		rgrad = stats::quantile(dg.r,cqs,na.rm=TRUE)
		grad = stats::quantile(dg,cqs,na.rm=TRUE)
		wgrad = stats::quantile(dg.w,cqs,na.rm=TRUE) 
		# plot all this:
		if(.Platform$OS.type=="unix"){ grDevices::quartz() } else { grDevices::x11() }
		graphics::par(font=2,font.axis=2,font.lab=2)
		plot.profile(outH,outSQ,sdf=10,mpatch=paste(", ",ptid))

		# return
		return(list(outH=outH,outSQ=outSQ,het0=het0,het1=het1,
						grad=grad,rgrad=rgrad,wgrad=wgrad))
	}
	
	# ------------------------------------------------------------------------ 
	# initialisations
	ptids = paste("Case_",c("A","B"),sep="")
	cqs = c(.25,.30,.40,.50,.60,.70,.75,.80,.85,.90,.95,1) # quantiles for gradients
	Q = length(cqs)
	cat("This script generates three figures:
			1) a visualisation of the two demo sarcoma datasets
			2) the uptake profile fit and corresponding gradients for Case A
			3) the uptake profile fit and corresponding gradients for Case B\n")
		
	# ------------------------------------------------------------------------ 
	# load + display the two demo cases
	grDevices::dev.new()
	op <- graphics::par(mfrow=c(2,3),font=2,font.lab=2,font.axis=2)
	Case_A_ROI <- Case_B_ROI <- NULL
	load("data/Case_A_ROI.rda")
	load("data/Case_B_ROI.rda")
	scan.views(Case_A_ROI,mlab=ptids[1])
	scan.views(Case_B_ROI,mlab=ptids[2])
	
	# ------------------------------------------------------------------------ 
	# analyse the two demo cases
	outA = analyse.and.plot(Case_A_ROI,ptids[1])
	outB = analyse.and.plot(Case_B_ROI,ptids[2])

	# ------------------------------------------------------------------------ 
	# present results
	het.df = data.frame(H0=c(outA$het0,outB$het0),H1=c(outA$het1,outB$het1))
	row.names(het.df) = ptids
	gradients.df = data.frame(quantile=cqs,
						(outA$rgrad),(outB$rgrad),
						(outA$grad),(outB$grad),
						(outA$wgrad),(outB$wgrad))
	names(gradients.df) = c("quantile",paste("raw",ptids),
										paste("norm'd",ptids),
										paste("w'ed",ptids))
	
	cat("\nHeterogeneity summaries:\n")
	print(round(het.df,2))
	
	cat("\nGradient summaries:\n")
	print(round(gradients.df,2))
	
	rm(Case_A_ROI)	
	rm(Case_B_ROI)
	graphics::par(op) # reset par()'s initial settings
	return(1)
}

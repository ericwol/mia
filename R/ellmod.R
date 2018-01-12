# ----------------------------------------------------------------------------------- 
# utils

#' Uptake-weighted PCA of 3D coordinates x... y is the uptake.
#'
#' @import stats
#' @export
project <- function(x,y){
	# u0 = mean of the xi weighted by yi
	wy = ifelse(y>0,y,mean(y)/100)	
	# weighted means of coordinates
	u0 = c(NA,NA,NA)
	for (i in 1:3){ 
		u0[i] = sum(wy*x[,i])/sum(wy) 
	}
	# get weighted mean and covariance, Sw
	S = cov.wt(x,wt=wy,cor=F,center=TRUE,method="unbiased")
	# u0=S$center - this is not right so calculate as above
	Sw = S$cov
	# Transform xi to simplify parameters - new matrix is xstar
	o = eigen(Sw)
	G = o$vectors; d=o$values; D=diag(d**(-.5))
	# S = G %*% D %*% t(G)
	cent =  t(apply(x, 1, function(x) x-u0))
	return(t(D %*% t(G) %*% t(cent)))
}

#' Local function, plots output of hetE and struct.quant analyses...
#' 	hetEo: 	output of hetE()
#'	sqo: 	output from struct.quant()
#'
#' @import stats
#' @import graphics
#' @export
view.profile <- function(hetEo,sqo,mpatch=""){
	v = hetEo$iso.reg.out
	yh0 = v[,2] 
	is = order(v[,1])
	u = v[is,1]
	dg = sqo$normalized.gradients
	dg.w = sqo$w.gradients
	graphics::par(mfrow=c(2,1))
	# uptake profiles
	ylimi=range(c(v[is,c(3,2)]))
	graphics::matplot(u,v[is,c(3,2)],
		xlab="u (radial ellipsoidal coordinates)",ylab="uptake",
		type="pl",lty=c(1),pch=20,ylim=ylimi,lwd=3,cex=.5,
		main=paste("Initial data and smooth regressions",mpatch,sep=""))
	graphics::points(u,sqo$yh.i,ylim=ylimi,cex=1.2,col=2)
	graphics::points(u,sqo$yh.b,col=4,pch=20,cex=1)
	graphics::abline(h=0,col=8,lwd=.8)	
	graphics::legend("topright",legend=c("isotonic","bitonic"),lty=1,
		lwd=c(5,4),col=c(2,4))
	# gradients
	graphics::plot(u,dg.w,pch=20,main=paste("gradients",mpatch,sep=""),
		xlab="u (radial ellipsoidal coordinates)",
		ylab="gradient summaries")
	graphics::points(u,dg,pch=20,col=8,cex=1.2)	
	graphics::abline(h=0,col=8,lwd=.8)	
	graphics::legend("topright",legend=c("weighted","normalized"),pch=20,col=c(1,2))
}

# ----------------------------------------------------------------------------------- 
# monotonic regression

#' Computes monotonically decreasing nonparametric least squares regression.
#' Input:
#'	 xstar:		projected coordinates
#'	 y:			uptake
#'	 uu1--uu3:	positioning of ellipse (center point)
#'	 a1--a5:		makings of covar matrix (ellipse)
#' Returns:
#' 	3-col array [tx,ytx,y] where
#' 	 tx:	 	ellipsoidal coordinates,
#' 	 ytx:		evaluation of monotnoic regression on tx grid,
#' 	 y: 		input uptake y.
#'
#' @import stats
#' @export
iso.reg <- function(xstar,y,uu1,uu2,uu3,a1,a2,a3,a4,a5) {
	# conversion matrix
	y = y*(-1) # because code fits a monotonic-increasing curve
	# compute ellipsoidal radii tx form xstar
	nr = length(y)
	R = cbind( c(a1,0,0) , c(a3,a2,0), c(a4,a5,1) ) 
	sigma = R %*% t(R)
	uu0 = cbind(uu1,uu2,uu3)
	beta = c(uu1,uu2,uu3,a1,a2,a3,a4,a5)
	tx = rep(NA,nr)
	for (i in 1:nr) {
		tx[i] = (xstar[i,]-uu0) %*% sigma %*% t(xstar[i,]-uu0)
	}
	# run isotonic regression
	iso = isoreg(y ~ tx)
	iso$yf = iso$yf*(-1); iso$y=iso$y*(-1); iso$yc=iso$yc*(-1)
	oo = order(iso$ord)
	keep = cbind(tx,stats::approx(iso$x,iso$yf[oo],xout=tx,rule=2)$y,-y)
	return(keep)
}

#' Makes a gradient to use in nls
#'
#' @export
grad <- function(xstar,y,uu1,uu2,uu3,a1,a2,a3,a4,a5){
	yh = iso.reg(xstar,y,uu1,uu2,uu3,a1,a2,a3,a4,a5)[,2]
	dev = (y-yh)
	chg = .01  
	nr = length(y)
	keeppar = matrix(nrow=nr,ncol=8)
	for (j in 1:8) {
		add = rep(0,8)
		add[j] = chg 
		yh1 = iso.reg(xstar,y,uu1+add[1],uu2+add[2],uu3+add[3],a1+add[4],	
					a2+add[5],a3+add[6],a4+add[7],a5+add[8])[,2]
		keeppar[,j] = (yh-yh1)/chg
	}
	nr = length(y)
	sds = sqrt(apply(keeppar,2,var))
	sdl = max(sds)/10^5 
	for(j in 1:8) { 
		if(sds[j]<sdl){ 
			keeppar[,j] = rnorm(nr)*sdl 
		}
	}	
	attr(dev, "gradient") = -1*(y-yh)*keeppar
	return(dev)
}

# ----------------------------------------------------------------------------------- 
# monotonic analyzer

#' Runs heterogeneity analysis based on isotonic regression, as in 
#'   Eary JF, O’Sullivan F, O’Sullivan J, Conrad EU. 
#'   Spatial heterogeneity in sarcoma 18F-FDG uptake as a predictor of patient outcome. 
#'   J Nucl Med. 2008; 49:1973–1979
#' Arguments
#' z:	input ROI
#'		y=z[,5]; x=z[,2:4]   #x is location and y is uptake
#' 		eg: read ROID; output is column 5; location is col 2:4
#' 			nr=scan('ROID',n=1);
#'			z=matrix(scan('ROID',skip=1),nrow=nr,byrow=T)
#' par0: initial parameters for nls
#'		uu1=par0[1]; uu2=par0[2]; uu3=par0[3]; 
#'		a1=par0[4] ; a2=par0[5] ; a3=par0[6] ; a4=par0[7] ; a5=par0[8]
#' Values
#' 	v:		 output of iso.reg:
#' 			   - fitted monotonic regression line: yh=v[,2] 
#'			   - residuals: res = (v[,3]-v[,2])
#' 	het0: 	 100*mean(res^2)/mean(y^2)	
#' 	het1: 	 100*mean(res^2)/var(y)
#'	vary:	 var(y)
#'	hetvals: round(c(het0,het1,var(y)),4) [for retro-compatibility purposes]
#' Updates:
#'	12 Dec 2016: - restricted call to nls to case where sd>mv, where 
#'					'sd' = sqrt(var(yh0)) (yh0 = fitted isotonic values) 
#'					'mv' = .0001*abs(mean(yh0))
#'				 - more robust call to nls
#'
#'
#' @import graphics
#' @export
hetE <- function(z,par0=NULL,doplot=FALSE){
	# Standardise	
	y=z[,5]; x=z[,2:4]   #x is location and y is uptake
	xstar = project(x,y)

	# Optim..................... 
	if(is.null(par0)){
		uu1=0; uu2=0; uu3=0; 
		a1=1 ; a2=1 ; a3=0 ; a4=0 ; a5=0
	} else {
		uu1=par0[1]; uu2=par0[2]; uu3=par0[3]; 
		a1=par0[4] ; a2=par0[5] ; a3=par0[6] ; a4=par0[7] ; a5=par0[8]
	}
	# run isotonic regression
	v = iso.reg(xstar,y,uu1,uu2,uu3,a1,a2,a3,a4,a5)
	if(doplot){
		graphics::matplot(sort(v[,1]),v[order(v[,1]),c(3,2)],type="pl",pch=".")
	}
	yh0 = v[,2] 
	nls.cv='FALSE' 
	aa=NULL 
	sd=sqrt(var(yh0)) 
	mv=.0001*abs(mean(yh0))
	if(sd>mv){
		# datalist=list(y=y,xstar=xstar) ; 
		iters=25 
		coef=c(uu1,uu2,uu3,a1,a2,a3,a4,a5) 
		nls.cv='FALSE'
		result=try(aa<-nls(~ grad(xstar,y,uu1,uu2,uu3,a1,a2,a3,a4,a5),
			control=nls.control(maxiter=iters,warnOnly=TRUE,minFactor=1/(4096*256)),
			trace=FALSE,
			#data=datalist,
			start=list(uu1=uu1,uu2=uu2,uu3=uu3,a1=a1,a2=a2,a3=a3,a4=a4,a5=a5)),silent=TRUE)
		if(!inherits(result,"try-error")){ 
			nls.cv = summary(aa)$convInfo$isConv 
			coef = summary(aa)$coef[,1] 
		} # mute if warnonly is TRUE 
		### need to make provision for non-convergent	
	}
	if (nls.cv=='TRUE'){
		uu1=coef[1] ; uu2=coef[2] ; uu3=coef[3]
		a1=coef[4]; a2=coef[5]; a3=coef[6]; a4=coef[7]; a5=coef[8]
	}
	v=iso.reg(xstar,y,uu1,uu2,uu3,a1,a2,a3,a4,a5)
	yh=v[,2] 
	res = (v[,3]-v[,2])
	het1=100*mean(res^2)/var(y)
	het0=100*mean(res^2)/mean(y^2)	
	hetvals=round(c(het0,het1,var(y)),4) 
	return(list(iso.reg.out=v,
				het0=het0,
				het1=het1,
				hetvals=hetvals,
				vary=var(y),
				nls.out=aa,
				nls.cv=nls.cv,
				nls.coef=coef))
}

# ----------------------------------------------------------------------------------- 
# Unimodal (bitonic) regression

#' Unimodal data fit when mode is selected via a simple rule.
#' Returns cbind(x,yhat), i.e. x values and regression fit y-values.
#' Version written 12/12/2016
#'
#' @import stats
#' @import quadprog
#' @export
unismooth <- function(x,y,w=rep(1,length(x))){ 
	# Mode point -> mx
	out = stats::smooth.spline(x,y,w=w,spar=.75) 
	sx = out$x 
	sy = out$y
	low = stats::quantile(x,.025) 
	high = stats::quantile(x,.975) 
	ssx = sx[sx>=low&sx<=high]
	ssy = sy[sx>=low&sx<=high]
	mx = ssx[order(-ssy)[1]]
	# Set breakpoints
 	ux = unique(sort(x)) 
 	ns = length(ux)
 	if(ns<100){ 
 		xs = unique(sort(c(ux,mx))) 
 		ns = length(xs) 
 	}
 	if(ns>100){
 		xs = unique(sort(c(ux[round(1+c(0:99)*(ns-1)/99)],mx)))
 		ns = length(xs) 
 	}
 	# Construct basis elements for unimodal smooth 
 	# (piecewise linear with slopes positive<mode and negative>mode)
	n = length(x)
	xx0 = rep(1,n)
	for(j in 1:(ns-1)){ 
		ds = xs[j+1]-xs[j]
		xj = ifelse(x>=xs[j+1],1,0)
		xj[(x>=xs[j])&(x<xs[j+1])] = (x[(x>=xs[j])&(x<xs[j+1])]-xs[j])/ds 
		xx0=cbind(xx0,xj) 
	}
 	# Design matrix and regularization
 	xmat = sqrt(w)*xx0 
 	xmx = t(xmat)%*%xmat 
  	o = eigen(xmx)
  	lvals = o$values
  	lvals = unique(c(max(lvals)/10^6,lvals[lvals>max(lvals)/10^20]))
 	lamn = NULL
 	dfns = NULL
 	M = diag(rep(1,ncol(xmat)))
 	for(lam in lvals){
 		o = eigen(xmx+lam*M)
 		v = o$vectors
 		l = o$values 
	 	if(min(l)>(.1e-9*max(l))){
	 		lamn = c(lamn,lam) 
	 		dfns = c(dfns, sum(diag(solve(xmx+lam*M,xmx)))) 
	 	}
  	}
 	nx = ncol(xmat)
 	lam = median(lamn,na.rm=T) 
 	if(length(dfns)>1){
 		lam = stats::approx(dfns,lamn,xout=min(n,nx)-.001,rule=2)$y
 	}
	AM = diag(rep(1,ncol(xmx)))
	kdel = order(abs(xs-mx))[1]
	xx = xx0
	if(kdel<ns){ 
		xx[,c((kdel+1):ns)] = 0-xx[,c((kdel+1):ns)] 
	}
	xmat = sqrt(w)*xx
	xmx = t(xmat)%*%xmat 
	DM = xmx+lam*M 
	dvec = t(xmat)%*%c(sqrt(w)*y) 
	o = quadprog::solve.QP(DM,dvec,AM) 
	bhatu = o$unconstrained.solution
	bhat = o$solution
	beta = ifelse(bhat>0,bhat,0)
	yhat = xmat%*%beta/sqrt(w) 
	yhat0 = xmat%*%bhatu/sqrt(w)
	wrssv = sum(w*(y-yhat)^2)  
 	return(cbind(x,yhat))
}

#' Returns a smoothed version of unismooth(u,y)'s bitonic curve output,
#' using a cubic smoothing spline.
#' ddf controls the number of degrees of freedom for smooth.spline.
#' Value
#' 	u: same as input u (radial coordinate of voxels)
#'	y: same as input y (uptake values at those voxels)
#'	breg:typical approx() output y, i.e. smoothed y-values, same length as u.
#'
#' @import stats
#' @export
unismooth.sp <- function(u,y,ddf=50){
	ur = unismooth(u,y)
	bs = stats::smooth.spline(ur,df=ddf)
	gb = stats::approx(bs$x,bs$y,xout=u,rule=2)$y 
	return(list(u=u,y=y,breg=gb))
}


#----------------------------------------------------------------------- final quantitation

#' Structural quantitation summaries of tumor uptake.
#' Arguments:
#'	out = output of hetE() 
#' 	zdim = slice thickness (defaults to 1)
#' 	xdim = voxel transverse dimensions (defaults to 1)
#' 			Ex: zdim=xx$Thickness[kk]; xdim=xx$InPlaneDim[kk]
#' 	bndT = cutoff quantile threshold defining the start of the boundary rim 
#'			(defaults to 75th quantile of VOI radial coordinates) 
#'	re.sign = if TRUE, flips gradient curve sign for better interpretation
#'			(a negative output gradient value would then indicate a decreasing rate)
#'  sm.reg = if TRUE, gradients are computed from a smoothed unimodal regression fit
#'  ddf = corresponding number of smooth.spline df's (between 5-10 usually works well)
#' Values:
#'	u:		radial coordinates (sorted in increasing order)
#'	yh.i: 	isotonic fit for u values
#'	yh.b: 	bitonic fit for u values
#' 	fractional.void = F1
#' 	volume = F2
#' 	contrast.ratio = F3
#' 	max = F4
#' 	boundary.gradient = F5
#' 	cv = F6
#' 	y = model uptake (where y = out$v[order(v[,1],3])
#' 	gradients = all gradient values evaluated at locations u 
#' 			(expressed in the units of the input uptake values)
#' 			(where u is a normalized version of sort(out$v[,1]))
#' 	normalized.gradients = gradients/F4 (unitless, normalized to (-1,1))
#' 	w.gradients = uptake-weighted normalized gradients,
#' 	bnd.gradients = gradients at boundary
#' 	normalized.bnd.gradients = normalized gradients at boundary
#' 	w.bnd.gradients = uptake-weighted (normalized) gradients at boundary
#' 
#' @import stats
#' @export
struct.quant <- function(out,zdim=1,xdim=1,bndT=.75,re.sign=TRUE,sm.reg=TRUE,ddf=5){
	vb = out$iso.reg.out
	o = order(vb[,1])
	y = vb[o,3]
	u = vb[o,1] 
	# u = u/median(ifelse(u>.001,u,.001))
	if(!sm.reg){ # just a switch:
		# ..... one could switch to this if unsmoothed regression cuve is preferred
		uu = unismooth(u,y)
	} else {
		uus = unismooth.sp(u,y,ddf)
		uu = cbind(uus$u,uus$breg)
	}
	yh = uu[,2]
	um = u[order(-yh)[1]]   # mode of unimodal uptake profile
	ym = min(yh) 			# uptake value at mode
	nu = length(u)
	maxyh = max(yh)
	du = (ifelse(u>.01,u,.01)*.2)
	grad = (stats::approx(u,yh,xout=u*1.1,rule=2)$y-stats::approx(u,yh,xout=u*.9,rule=2)$y)/du
	if(re.sign){grad = -grad} # flip curve for better interpretation
	wgrad = (grad/maxyh)*y
	A = mean(yh[u<=um])-ym 
	B = mean(yh[u>um])-ym 
	F1 = A/(A+B)										# Fractional void (voxel weighted)
	F2 = nu*zdim*xdim*xdim								# Volume of tumor
	F3 = 1-min(yh)/maxyh	    		   				# Contrast ratio
	F4 = maxyh									 		# Max
	F5 = 0-mean(grad[u>=stats::quantile(u,bndT)])/maxyh # (Mean) boundary gradient
	F6 = mean( abs(y-yh) )/maxyh						# CoV about the tumor profile
	res = list( fractional.void=F1,
				volume=F2,
				contrast.ratio=F3,
				max=F4,
				boundary.gradient=F5,
				cv=F6,
				y=y,
				u=u,
				yh.i=vb[o,2], # isotonic
				yh.b=uu[,2], # bitonic
				gradients=grad,
				normalized.gradients=grad/maxyh,
				w.gradients=wgrad,
				bnd.gradients=grad[u>=stats::quantile(u,bndT)],
				normalized.bnd.gradients=grad[u>=stats::quantile(u,bndT)]/maxyh,
				w.bnd.gradients=wgrad[u>=stats::quantile(u,bndT)])
	return(res)
}

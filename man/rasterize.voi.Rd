% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/impro.R
\name{rasterize.voi}
\alias{rasterize.voi}
\title{Turns a list-format VOI into a raster form (i.e. a 3D array).
Input:
  voi is in list-format (u,w,x,y,z) as in AMIDE
Value:
  returns a raster form of the voi
The output raster is filled with NA's where values 
are not found in voi.}
\usage{
rasterize.voi(voi, def = NA)
}
\description{
Turns a list-format VOI into a raster form (i.e. a 3D array).
Input:
  voi is in list-format (u,w,x,y,z) as in AMIDE
Value:
  returns a raster form of the voi
The output raster is filled with NA's where values 
are not found in voi.
}

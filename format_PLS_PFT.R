#Format the PLS composition MCMC output provided by Chris Paciorek.
#Input is 8km Albers projection of PLS statistical output at taxon level. 
#Output is netCDF file of deciduous/evergreen PFTs in lon,lat 0.5-degree grid and 
#an up-scaled to 2.5-degree grid, downscaled to 0.5-degree grid version for resolution analysis. 
#Jaclyn Hatala Matthes, 5/1/15

library(ncdf)
library(fields)
library(maps)
library(sp)
library(rgdal)
library(raster)
library(abind)

#local functions for matrix manipulation
mirror.matrix <- function(x) {
  xx <- as.data.frame(x);
  xx <- rev(xx);
  xx <- as.matrix(xx);
  xx;
}
rotate180.matrix <- function(x) {
  xx <- rev(x);
  dim(xx) <- dim(x);
  xx;
}
flip.matrix <- function(x) {
  mirror.matrix(rotate180.matrix(x))
}

# #savanna raster - REVISIT THIS LATER
# str_name <- "/Users/jmatthes/Dropbox/BU/anderson.tif" 
# savanna_raster=raster(str_name)
# savanna_raster[savanna_raster>1] <- 0 #mask prairie and savanna types
# savanna.mask <- flip.matrix(rotate180.matrix(t(as.matrix(savanna_raster))))
# savanna.mask[savanna.mask==1] <- 10
# savanna.mask[is.na(savanna.mask)] <- 1

#load data
pls.path <- "/Users/jmatthes/Dropbox/BU/"
pls.nc   <- open.ncdf(paste(pls.path,"composition_v0.3.nc",sep=""))
pls.x <- get.var.ncdf(pls.nc,"x")
pls.y <- get.var.ncdf(pls.nc,"y")
pls.mask <- open.ncdf(paste(pls.path,"paleonmask.nc",sep=""))

#load PLS domain raster for Albers projection reference
pls.ref.raster <- raster(paste(pls.path,"paleon_full_alb_v0.1.tif",sep=""))

#define PFTs
decid.pft <- c("Ash","Basswood","Beech","Birch","Chestnut","Hickory","Ironwood",
               "Maple","Oak","Tamarack","Walnut","Other hardwood")
conif.pft <- c("Atlantic white cedar","Fir","Hemlock","Pine","Spruce")

#aggregate deciduous taxa into single variable
for(p in 1:length(decid.pft)){
  pft.tmp <- apply(get.var.ncdf(pls.nc,decid.pft[p]),c(1,2),median)
  pft.quan  <- apply(get.var.ncdf(pls.nc,decid.pft[p]),c(1,2),quantile,c(0.05,0.95),na.rm=TRUE)
  
  pft.tmp[pft.tmp==1e+30] <- 0
  if(p==1){
    dec.pft <- pft.tmp
  } else {
    #add PFTs while ignoring NA values in matrix since some taxa are only in East or West domains
    dec.pft <- ifelse(is.na(dec.pft), ifelse(is.na(pft.tmp), NA, pft.tmp), ifelse(is.na(pft.tmp), dec.pft, dec.pft + pft.tmp))
  }
}

#aggregate evergreen taxa into single variable
for(p in 1:length(conif.pft)){
  pft.tmp <- apply(get.var.ncdf(pls.nc,conif.pft[p]),c(1,2),median)
  pft.tmp[pft.tmp==1e+30] <- 0
  if(p==1){
    evg.pft <- pft.tmp
  } else {
    #add PFTs while ignoring NA values in matrix since some taxa are only in East or West domains    
    evg.pft <- ifelse(is.na(evg.pft), ifelse(is.na(pft.tmp), NA, pft.tmp), ifelse(is.na(pft.tmp), evg.pft, evg.pft + pft.tmp))
  }
}

#convert matrix to geo-referenced raster with defined Albers projection
rast.dec <- raster(mirror.matrix(rotate180.matrix(t(dec.pft))),
                   xmn=range(nc.x)[1], xmx=range(nc.x)[2],
                   ymn=range(nc.y)[1], ymx=range(nc.y)[2], crs=crs(pls.ref.raster))
rast.evg <- raster(mirror.matrix(rotate180.matrix(t(evg.pft))),
                   xmn=range(nc.x)[1], xmx=range(nc.x)[2],
                   ymn=range(nc.y)[1], ymx=range(nc.y)[2], crs=crs(pls.ref.raster))

#convert projection to 0.5-degree lat-lon base grid for the whole domain:
base.rast <- raster(xmn = -99.75, xmx = -60.25, ncols=80,
                    ymn = 35.25,  ymx = 49.75, nrows = 30,
                    crs = "+proj=longlat +datum=WGS84")
dec.pft.lonlat <- projectRaster(from=rast.dec, to=base.rast, method="bilinear")
evg.pft.lonlat <- projectRaster(from=rast.evg, to=base.rast, method="bilinear")

#convert projection to 2.5-degree lat-lon base grid for the whole domain:
base.rast.lo <- raster(xmn = -99.75, xmx = -60.25, ncols=16,
                    ymn = 35.25,  ymx = 49.75, nrows = 6,
                    crs = "+proj=longlat +datum=WGS84")
dec.pft.lonlat.tmp <- projectRaster(from=rast.dec, to=base.rast.lo, method="bilinear")
evg.pft.lonlat.tmp <- projectRaster(from=rast.evg, to=base.rast.lo, method="bilinear")
dec.pft.lonlat.lo <- projectRaster(from=dec.pft.lonlat.tmp, to=base.rast, method="bilinear")
evg.pft.lonlat.lo <- projectRaster(from=evg.pft.lonlat.tmp, to=base.rast, method="bilinear")

#export the lon,lat 0.5-degree deciduous/evergreen netcdf files 
#set up files/variables - had trouble getting rasters to netCDF, some dumb formatting in here
mv <- 1e30    #missing value for output .nc file
dec.rnc <- writeRaster(dec.pft.lonlat,filename=paste(pls.path,"dec_out.nc",sep=""),format="CDF",overwrite=TRUE)
evg.rnc <- writeRaster(evg.pft.lonlat,filename=paste(pls.path,"evg_out.nc",sep=""),format="CDF",overwrite=TRUE)
dec.lo.rnc <- writeRaster(dec.pft.lonlat.lo,filename=paste(pls.path,"dec_lo_out.nc",sep=""),format="CDF",overwrite=TRUE)
evg.lo.rnc <- writeRaster(evg.pft.lonlat.lo,filename=paste(pls.path,"evg_lo_out.nc",sep=""),format="CDF",overwrite=TRUE)

dec.rnc <- open.ncdf(paste(pls.path,"dec_out.nc",sep=""))
dec.dat <- mirror.matrix(get.var.ncdf(dec.rnc,"layer"))
evg.rnc <- open.ncdf(paste(pls.path,"evg_out.nc",sep=""))
evg.dat <- mirror.matrix(get.var.ncdf(evg.rnc,"layer"))

dec.lo.rnc <- open.ncdf(paste(pls.path,"dec_lo_out.nc",sep=""))
dec.lo.dat <- mirror.matrix(get.var.ncdf(dec.lo.rnc,"layer"))
evg.lo.rnc <- open.ncdf(paste(pls.path,"evg_lo_out.nc",sep=""))
evg.lo.dat <- mirror.matrix(get.var.ncdf(evg.lo.rnc,"layer"))

lon <- seq(-99.75,-60.25,by=0.5) 
lat <- seq(35.25, 49.75, by=0.5)

ncname <- paste(pls.path,"PLS_pft_gridded.nc",sep="")
nc_variable_units <- "fraction cover"
dimX <- dim.def.ncdf( "lon", "longitude: degrees", lon )
dimY <- dim.def.ncdf( "lat", "latitude: degrees", lat )
dimZ <- dim.def.ncdf( "functional type", "plant (tree) functional type: 1 = evergreen, 2 = deciduous",1:2)
var.nc <- var.def.ncdf("pft",nc_variable_units, list(dimX,dimY,dimZ),mv) #set up variables: lon, lat, pft
nc <- create.ncdf(ncname, list(var.nc) ) #create the file
put.var.ncdf(nc, var.nc, abind(evg.dat,dec.dat,along=3)) #write data to the file
close.ncdf(nc)

ncname.lo <- paste(pls.path,"PLS_pft_gridded_lo.nc",sep="")
nc.lo <- create.ncdf(ncname.lo, list(var.nc) ) #create the file
put.var.ncdf(nc.lo, var.nc, abind(evg.lo.dat,dec.lo.dat,along=3)) #write data to the file
close.ncdf(nc.lo)

# #look at 95% CI for deciduous taxa 
# for(p in 1:length(decid.pft)){
#   pft.tmp <- apply(get.var.ncdf(pls.nc,decid.pft[p]),c(1,2),quantile,c(0.05,0.95),na.rm=TRUE)
# #  pft.tmp <- get.var.ncdf(pls.nc,decid.pft[p])
#   pft.tmp[pft.tmp==1e+30] <- 0
#   if(p==1){
#     dec.pft <- pft.tmp
#   } else {
#     #add PFTs while ignoring NA values in matrix since some taxa are only in East or West domains
#     dec.pft <- ifelse(is.na(dec.pft), ifelse(is.na(pft.tmp), NA, pft.tmp), ifelse(is.na(pft.tmp), dec.pft, dec.pft + pft.tmp))
#   }
# }
# 


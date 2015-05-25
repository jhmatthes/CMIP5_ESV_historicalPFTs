#Aggregate CRUNCEP monthly files at 6-hourly resolution into a 30-year climatology
#for mean annual temperature (MAT) and mean total annual precipitation (MAP).
#Input:
#cru.path = directory of CRUNCEP data, where tas/ and pr/ correspond to monthly 6-hour mean air temperature & precipitation
#pls.mask = mask for the extent of the PLS (ESV) dataset.
#Jaclyn Hatala Matthes, 5/6/15

cruncep_climatology <- function(cru.path,pls.mask){

  #load & format 30-year CRUNCEP temperature/precip climatology 1980-2010
  tas.files <- list.files(paste(cru.path,"tas/",sep=""))[c(grep("198",list.files(paste(cru.path,"tas/",sep=""))),grep("199",list.files(paste(cru.path,"tas/",sep=""))),
                                                           grep("200",list.files(paste(cru.path,"tas/",sep=""))),grep("201",list.files(paste(cru.path,"tas/",sep=""))))]
  pr.files <- list.files(paste(cru.path,"pr/",sep=""))[c(grep("198",list.files(paste(cru.path,"pr/",sep=""))),grep("199",list.files(paste(cru.path,"pr/",sep=""))),
                                                         grep("200",list.files(paste(cru.path,"pr/",sep=""))),grep("201",list.files(paste(cru.path,"pr/",sep=""))))]
  
  #CRUNCEP mean annual temperature, 1981-2010
  tas.cru <- list()
  for(f in 1:length(tas.files)){
    ncf <- open.ncdf(paste(cru.path,"tas/",tas.files[f],sep=""))
    cm.dat    <- get.var.ncdf(ncf,'tair')
    fillvalue <- att.get.ncdf(ncf,'tair','_FillValue')
    cm.dat[cm.dat==fillvalue$value] <- NA
    close.ncdf(ncf)
    
    #average monthly temp 
    tas.cru[[f]] <- (apply(cm.dat,c(1,2),mean)-273.15)*pls.mask #mask from PLS data
    if(f==1){
      tas.cru.clim <- tas.cru[[f]]
    }else{
      tas.cru.clim <- abind(tas.cru.clim,tas.cru[[f]],along=3)
    }
  }
  tas.cru.mean <- apply(tas.cru.clim,c(1,2),mean)     
  
  pr.cru <- list()
  yr.ind <- seq(1,length(pr.files)+1,by=12)
  for(f in 1:length(pr.files)){
    ncf <- open.ncdf(paste(cru.path,"pr/",pr.files[f],sep=""))
    cm.dat    <- get.var.ncdf(ncf,"precipf")
    fillvalue <- 1.00000001504747e+30
    cm.dat[cm.dat==fillvalue] <- NA
    close.ncdf(ncf)
    
    #monthly precip sum into yearly sum
    pr.cru[[f]] <- (apply(cm.dat,c(1,2),sum,na.rm=TRUE)*(60*60*6))*pls.mask #mask from PLS data
    if(f %in% yr.ind){ #sum along months of year
      pr.cru.ann <- pr.cru[[f]]
    }else{
      pr.cru.ann <- pr.cru.ann + pr.cru[[f]]
    }
    
    if(f %in% (yr.ind-1)){ #if it's the end of the year
      if(f==12){
        pr.cru.clim <- pr.cru.ann
      } else {
        pr.cru.clim <- abind(pr.cru.clim,pr.cru.ann,along=3)
      }
    }
  }
  pr.cru.mean <- apply(pr.cru.clim[,,2:30],c(1,2),mean) 
  pr.cru.mean[pr.cru.mean==0] <- mean(pr.cru.mean, na.rm=TRUE)
  
  cru.clim     <- list()
  cru.clim$tas <- tas.cru.mean 
  cru.clim$pr  <- pr.cru.mean
  cru.clim
}
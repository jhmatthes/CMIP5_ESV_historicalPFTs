#Aggregate CMIP5 output at monthly resolution into a 30-year climatology
#for mean annual temperature (MAT) and mean total annual precipitation (MAP).
#Input: 
#cmip5.path = CMIP5 output directory, where tas/ and pr/ each hold the files for all models 
#             for monthly mean air temperature & precipitation
#pls.mask = mask for the extent of the PLS (ESV) dataset
#model.full = list of full model names from PFT list to match model names to tas and pr files
#Jaclyn Hatala Matthes, 5/6/15

cmip5_climatology <- function(cmip5.path,pls.mask,model.full){

  #Aggregate CMIP5 30-year climatology
  tas.gcm.files <- list.files(paste(cmip5.path,"tas/",sep=""),pattern = "\\.nc$")
  pr.gcm.files  <- list.files(paste(cmip5.path,"pr/",sep=""),pattern = "\\.nc$")
  
  #match models to tas & pr files
  tas.fmatch <- pr.fmatch <- vector(length=length(model.full))
  for(m in 1:length(model.full)){
    match <- grep(model.full[m],tas.gcm.files)
    tas.fmatch[m] <- tas.gcm.files[match]
    match <- grep(model.full[m],pr.gcm.files)
    pr.fmatch[m]  <- pr.gcm.files[match]
  }
  
  #calculate first 30-year annual temp mean 
  tas.mean <- list()
  for(f in 1:length(tas.fmatch)){
    nc  <- open.ncdf(paste(cmip5.path,"tas/",tas.fmatch[f],sep=""))
    dat <- get.var.ncdf(nc,"tas")
    close.ncdf(nc)
    
    #30-year mean annual temperature 
    tas.mean[[f]] <- (apply(dat[,,1:(12*30)],c(1,2),mean)-273.15)*pls.mask
  }
  
  #calculate 30-year annual precip mean 
  pr.mean <- list()
  for(f in 1:length(pr.fmatch)){
    nc  <- open.ncdf(paste(cmip5.path,"pr/",pr.fmatch[f],sep=""))
    dat <- get.var.ncdf(nc,"pr")
    close.ncdf(nc)
    
    for(y in 1:30){
      yr.ind <- ((y-1)*12+1):(y*12)
      #convert mean monthly precip to total annual precip
      pr.tmp <- (apply(dat[,,yr.ind],c(1,2),sum)*(60*60*24*365/12))*pls.mask
      if(y==1){
        pr.mean[[f]] <- pr.tmp
      } else {
        pr.mean[[f]] <- abind(pr.mean[[f]],pr.tmp,along=3)
      }
    }
    #30-year mean annual precip
    pr.mean[[f]] <- apply(pr.mean[[f]],c(1,2),mean) 
  }

  cmip5.clim     <- list()
  cmip5.clim$tas <- tas.mean
  cmip5.clim$pr  <- pr.mean
  cmip5.clim$mod <- model.full
  cmip5.clim
}
#Aggregate CMIP5 output within the regional domain from monthly mean resolution to the 30-year mean total annual
#net primary productivity (npp), mean total annual evapotranspiration (tran), and mean annual albedo (rsds/rsus).
#Input: 
#cmip5.path = CMIP5 output directory, where npp/, tran/, rsds/, and rsus/ each hold the files for all models 
#             for monthly mean npp, evapotranspiration, down-welling radiation (rsds), and up-welling radiation (rsus).
#pls.mask = mask for the extent of the PLS (ESV) dataset
#model.full = list of full model names from PFT list to match model names to tas and pr files
#Jaclyn Hatala Matthes, 5/6/15

cmip5_land_fluxes <- function(cmip5.path,pls.mask,model.full){
  
  #match models to ecosystem flux files
  tran.files <- list.files(paste(cmip5.path,"tran/",sep=""),pattern = "\\.nc$")
  npp.files <- list.files(paste(cmip5.path,"npp/",sep=""),pattern = "\\.nc$")
  rsds.files  <- list.files(paste(cmip5.path,"rsds/",sep=""),pattern = "\\.nc$")
  rsus.files  <- list.files(paste(cmip5.path,"rsus/",sep=""),pattern = "\\.nc$")
  
  tran.fmatch <- npp.fmatch <- rsds.fmatch <- rsus.fmatch <- vector()
  for(m in 1:length(model.full)){
    match <- grep(model.full[m],tran.files)
    if(length(match)>0){tran.fmatch[m] <- tran.files[match]}
    match <- grep(model.full[m],npp.files)
    if(length(match)>0){npp.fmatch[m]  <- npp.files[match]}
    match <- grep(model.full[m],rsds.files)
    if(length(match)>0){rsds.fmatch[m]  <- rsds.files[match]}
    match <- grep(model.full[m],rsus.files)
    if(length(match)>0){rsus.fmatch[m]  <- rsus.files[match]}
  }
  
  #calculate 30-year mean of annual sum npp and tran and mean annual albedo
  #have to do each of these separately for each variable since some models are missing flux output variables
  npp.mean <- tran.mean <- alb.mean <- list()
  for(f in 1:length(npp.fmatch)){
    if(!is.na(npp.fmatch[f])){
      nc  <- open.ncdf(paste(cmip5.path,"npp/",npp.fmatch[f],sep=""))
      dat <- get.var.ncdf(nc,"npp")
      close.ncdf(nc)
      
      for(y in 1:30){
        yr.ind <- ((y-1)*12+1):(y*12)
        npp.tmp <- (apply(dat[,,yr.ind],c(1,2),sum)*(60*60*24*365/12))*pls.mask
        if(y==1){
          npp.mean[[f]] <- npp.tmp
        } else {
          npp.mean[[f]] <- abind(npp.mean[[f]],npp.tmp,along=3)
        }
      }
      npp.mean[[f]] <- apply(npp.mean[[f]],c(1,2),mean) 
    } 
    else { npp.mean[[f]] <- NA }
  }
  
  #regional domain TRANSPIRATION
  for(f in 1:length(tran.fmatch)){
    if(!is.na(tran.fmatch[f])){
      nc  <- open.ncdf(paste(cmip5.path,"tran/",tran.fmatch[f],sep=""))
      dat <- get.var.ncdf(nc,"tran")
      close.ncdf(nc)
      
      for(y in 1:30){
        yr.ind <- ((y-1)*12+1):(y*12)
        tran.tmp <- (apply(dat[,,yr.ind],c(1,2),sum)*(60*60*24*365)/12)*pls.mask
        if(y==1){
          tran.mean[[f]] <- tran.tmp
        } else {
          tran.mean[[f]] <- abind(tran.mean[[f]],tran.tmp,along=3)
        }
      }
      tran.mean[[f]] <- apply(tran.mean[[f]],c(1,2),mean) 
    } 
    else { tran.mean[[f]] <- NA }
  }
  
  #regional ALBEDO
  for(f in 1:length(rsus.fmatch)){
    if(!is.na(rsds.fmatch[f])){
      nc  <- open.ncdf(paste(cmip5.path,"rsds/",rsds.fmatch[f],sep=""))
      dat.d <- get.var.ncdf(nc,"rsds")
      close.ncdf(nc)
      
      nc  <- open.ncdf(paste(cmip5.path,"rsus/",rsus.fmatch[f],sep=""))
      dat.u <- get.var.ncdf(nc,"rsus")
      close.ncdf(nc)
      
      dat <- dat.u/dat.d
      
      for(y in 1:30){
        yr.ind <- ((y-1)*12+1):(y*12)
        alb.tmp <- apply(dat[,,yr.ind],c(1,2),mean)*pls.mask
        if(y==1){
          alb.mean[[f]] <- alb.tmp
        } else {
          alb.mean[[f]] <- abind(alb.mean[[f]],alb.tmp,along=3)
        }
      }
      alb.mean[[f]] <- apply(alb.mean[[f]],c(1,2),mean) 
    } 
    else { alb.mean[[f]] <- NA }
  }
  
  cmip5.reg.fluxes      <- list()
  cmip5.reg.fluxes$npp  <- npp.mean
  cmip5.reg.fluxes$tran <- tran.mean
  cmip5.reg.fluxes$alb  <- alb.mean
  cmip5.reg.fluxes$mod  <- model.full
  cmip5.reg.fluxes
}
#Aggregate CMIP5 output within the regional domain from monthly mean resolution to the 30-year mean total annual
#net primary productivity (npp), mean total annual evapotranspiration (tran), and mean annual albedo (rsds/rsus).
#Input: 
#cmip5.path = CMIP5 output directory, where npp/, tran/, rsds/, and rsus/ each hold the files for all models 
#             for monthly mean npp, evapotranspiration, down-welling radiation (rsds), and up-welling radiation (rsus).
#pls.mask = mask for the extent of the PLS (ESV) dataset
#model.full = list of full model names from PFT list to match model names to tas and pr files
#model.code = binary (1/0) indicator of PFT file dimension
#Jaclyn Hatala Matthes, 5/6/15

format_cmip5_pft <- function(cmip5.path,model.full,mod.pft,pls.mask){
  
  lcf.files <- list.files(paste(cmip5.path,"landCoverFrac/",sep=""),pattern = "\\.nc$")
  
  lcf.evg <- lcf.dec <- lcf.evg.map <- lcf.dec.map <- list()
  for(f in 1:length(lcf.files)){
    
    #get model-specific PFT codes
    dec.ind <- unique(as.numeric(mod.pft[grep(model.full[f],paste(as.character(mod.pft[,1]),"_",sep="")),3]))
    evg.ind <- unique(as.numeric(mod.pft[grep(model.full[f],paste(as.character(mod.pft[,1]),"_",sep="")),2]))
    
    #open .nc landCoverFrac data for each model
    nc  <- open.ncdf(paste(cmip5.path,"landCoverFrac/",lcf.files[f],sep=""))
    dat <- get.var.ncdf(nc,"landCoverFrac")
    close.ncdf(nc)
    
    print(f)
    
    #aggregate model landCoverFrac into deciduous and evergreen PFTs
    if(length(dim(dat))==3){ #if file only has 3 dims - no time
      lcf.dec.map[[f]] <- apply(dat[,,dec.ind]/100,c(1,2),sum)/
        apply(dat[,,c(dec.ind,evg.ind)]/100,c(1,2),sum)*pls.mask
      lcf.dec[[f]] <- as.vector(lcf.dec.map[[f]])
      lcf.dec[[f]] <- lcf.dec[[f]][complete.cases(lcf.dec[[f]])]
      
      lcf.evg.map[[f]] <- apply(dat[,,evg.ind]/100,c(1,2),sum)/
        apply(dat[,,c(dec.ind,evg.ind)]/100,c(1,2),sum)*pls.mask
      lcf.evg[[f]] <- as.vector(lcf.evg.map[[f]])
      lcf.evg[[f]] <- lcf.evg[[f]][complete.cases(lcf.evg[[f]])]
      
    } else if(length(dim(dat))>3) { 
      #grab first model timepoint 
      lcf.dec.map[[f]] <- apply(dat[,,dec.ind,1]/100,c(1,2),sum)/
        (apply(dat[,,c(dec.ind,evg.ind),1]/100,c(1,2),sum))*pls.mask
      lcf.dec[[f]] <- as.vector(lcf.dec.map[[f]])
      lcf.dec[[f]] <- lcf.dec[[f]][complete.cases(lcf.dec[[f]])]
      
      lcf.evg.map[[f]] <- apply(dat[,,evg.ind,1]/100,c(1,2),sum)/
        (apply(dat[,,c(evg.ind,dec.ind),1]/100,c(1,2),sum))*pls.mask
      lcf.evg[[f]] <- as.vector(lcf.evg.map[[f]])
      lcf.evg[[f]] <- lcf.evg[[f]][complete.cases(lcf.evg[[f]])]
    }
  }
  
  cmip5.pft         <- list()
  cmip5.pft$evg     <- lcf.evg
  cmip5.pft$dec     <- lcf.dec
  cmip5.pft$evg.map <- lcf.evg.map
  cmip5.pft$dec.map <- lcf.dec.map
  cmip5.pft$mod     <- model.full
  cmip5.pft
}
#Calculate the climate niche correspondence for deciduous/evergreen PFTs in the ESV-CRUNCEP data products and CMIP5 models. 
#Requires 'phyloclim' R library to calculate niche overlap in climate space. 
#Input:
#esv.evg/.dec = ESV PFT fraction evergreen and deciduous
#cmip5.pft    = CMIP5 PFT fraction for evergreen and deciduous 
#cru.clim     = CRUNCEP 30-year climatology for tas (mean annual temp) and pr (mean annual precipitation)
#cmip5.clim   = CMIP5 30-year climatology for tas (mean annual temp) and pr (mean annual precipitation)
#samp         = number of samples for random niche bootstrapping
#
#Output:
#stats = table with 
#Jaclyn Hatala Matthes, 5/18/15
esv_cmip5_niche_overlap <- function(esv.evg, esv.dec, cmip5.pft, cru.clim, cmip5.clim, model.full, samp, out.path){

  #format ESV dataset
  z1.esv <- as.vector(esv.evg)
  z2.esv <- as.vector(esv.dec)
  x.esv  <- as.vector(cru.clim$tas)
  y.esv  <- as.vector(cru.clim$pr)
  esv.comp.set1 <- c(!is.na(z1.esv) & !is.na(x.esv) & !is.na(y.esv))
  esv.comp.set2 <- c(!is.na(z2.esv) & !is.na(x.esv) & !is.na(y.esv))
  
  #interpolate on same MAT/MAP grid to evaluate CMIP5 vs ESV dataset niche overlaps 
  tas.grid <- seq(0,15,by=0.5)
  pr.grid  <- seq(475,1550,by=50)
  
  evgdat.esv <- cbind(x.esv[esv.comp.set1],y.esv[esv.comp.set1],z1.esv[esv.comp.set1])
  colnames(evgdat.esv) <- c("x.esv","y.esv","z1.esv")
  evgdat.esv <- data.frame(evgdat.esv)
  evg.esv <- with(evgdat.esv,interp(x.esv,y.esv,z1.esv,xo=tas.grid,yo=pr.grid),length=20)
  evg.esv <- evg.esv$z
  evg.esv[is.na(evg.esv)] <- 0
  
  decdat.esv <- cbind(x.esv[esv.comp.set2],y.esv[esv.comp.set2],z2.esv[esv.comp.set2])
  colnames(decdat.esv) <- c("x.esv","y.esv","z2.esv")
  decdat.esv <- data.frame(decdat.esv)
  dec.esv <- with(decdat.esv,interp(x.esv,y.esv,z2.esv,xo=tas.grid,yo=pr.grid),length=20)
  dec.esv <- dec.esv$z
  dec.esv[is.na(dec.esv)] <- 0
  
  #loop through CMIP5 models and calculate climate space niches
  for(f in 1:length(cmip5.clim$tas)){
    
    #set up each CMIP5 model-specific comparison set
    z1 <- as.vector(cmip5.pft$evg.map[[f]])
    z2 <- as.vector(cmip5.pft$dec.map[[f]])
    x <- as.vector(cmip5.clim$tas[[f]])
    y <- as.vector(cmip5.clim$pr[[f]])
    comp.set1 <- c(!is.na(z1) & !is.na(x) & !is.na(y))
    comp.set2 <- c(!is.na(z2) & !is.na(x) & !is.na(y))
    gcm.evg <- data.frame(cbind(x[comp.set1],y[comp.set1],z1[comp.set1]))
    gcm.dec <- data.frame(cbind(x[comp.set1],y[comp.set1],z2[comp.set2]))
    colnames(gcm.evg) <- c("x.gcm","y.gcm","z.gcm")
    colnames(gcm.dec) <- c("x.gcm","y.gcm","z.gcm")
    
    #interpolate response in temp/precip space
    evg.gcm <- with(gcm.evg,interp(x.gcm,y.gcm,z.gcm,xo=tas.grid,yo=pr.grid),length=20)
    evg.den <- evg.gcm$z
    evg.den[is.na(evg.den)] <- 0
    
    dec.gcm <- with(gcm.dec,interp(x.gcm,y.gcm,z.gcm,xo=tas.grid,yo=pr.grid),length=20)
    dec.den <- dec.gcm$z
    dec.den[is.na(dec.den)] <- 0
    
    if(f==1){
      evg.all <- cbind(rep(1,length(as.vector(evg.esv))),as.vector(evg.esv),as.vector(evg.den))
      dec.all <- cbind(rep(1,length(as.vector(dec.esv))),as.vector(dec.esv),as.vector(dec.den)) 
    } else{
      evg.all <- cbind(evg.all,as.vector(evg.den))
      dec.all <- cbind(dec.all,as.vector(dec.den))
    }
  }
  colnames(dec.all) <- c("dummy","PLS",model.full)
  colnames(evg.all) <- c("dummy","PLS",model.full)
  dec.all <- data.frame(dec.all)
  evg.all <- data.frame(evg.all)
  
  no.evg <- niche.overlap(evg.all)
  no.dec <- niche.overlap(dec.all)
  
  #bootstrap data to find distributions of random niche correlation
  hel.dist.evg <- hel.dist.dec <- matrix(ncol=ncol(evg.all)-1,nrow=samp)
  for(s in 1:samp){
    samp.mat <- matrix(ncol=(ncol(evg.all)-1),nrow=nrow(evg.all))
    for(m in 1:(ncol(evg.all)-1)){
      samp.mat[,m] <- sample(evg.all[,(m+1)],replace=TRUE)
    }
    samp.mat <- data.frame(cbind(evg.all[,1],samp.mat))
    colnames(samp.mat) <- colnames(evg.all)
    
    no.samp <- niche.overlap(samp.mat)
    hel.dist.evg[s,] <- no.samp[,1]
    
    samp.mat <- matrix(ncol=(ncol(dec.all)-1),nrow=nrow(dec.all))
    for(m in 1:(ncol(dec.all)-1)){
      samp.mat[,m] <- sample(dec.all[,(m+1)],replace=TRUE)
    }
    samp.mat <- data.frame(cbind(dec.all[,1],samp.mat))
    colnames(samp.mat) <- colnames(dec.all)
    
    no.samp <- niche.overlap(samp.mat)
    hel.dist.dec[s,] <- no.samp[,1]
  }
  evg.sd <- apply(hel.dist.con[,2:13],2,sd)
  evg.mean <- apply(hel.dist.con[,2:13],2,mean)
  evg.95 <- apply(hel.dist.con[,2:13],2,quantile, probs = c(0.05, 0.95),  na.rm = TRUE) 
  
  dec.sd <- apply(hel.dist.dec[,2:13],2,sd)
  dec.mean <- apply(hel.dist.dec[,2:13],2,mean)
  dec.95 <- apply(hel.dist.dec[,2:13],2,quantile, probs = c(0.05, 0.95),  na.rm = TRUE) 
  
  #organize & export summary statistics to table
  stats.out <- cbind(no.evg[2:13,1],evg.mean,evg.sd,t(evg.95),no.dec[2:13,1],dec.mean,dec.sd,t(dec.95))
  colnames(stats.out) <- c("evg_niche_ovlp","evg_rand_mean","evg_rand_sd","evg_rand_05","evg_rand_95",
                           "dec_niche_ovlp","dec_rand_mean","dec_rand_sd","dec_rand_05","dec_rand_95")
  rownames(stats.out) <- model.full
  
  write.csv(stats.out,file=paste(out.path,"niche_comp_stats.csv",sep=""),row.names=TRUE)
  
  #output variables
  niche.comp <- list()
  niche.comp$stats  <- stats.out
  niche.comp$no.evg <- no.evg
  niche.comp$no.dec <- no.dec
  niche.comp
}

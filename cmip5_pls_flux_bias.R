#Calculate bias between CMIP5 models and PLS by mapping global NPP/transpiration/albedo for high-PFT
#grid cells, then calculating what variable should have been given PLS mix of dec/evg.
#Input: 
#.gcm vars = CMIP5 variables for PFT fraction, MAT, MAP, and the flux variable 
#.reg. = CMIP5 output re-gridded to 0.5-deg resolution and clipped to ESV extent 
#pls.evg/.dec = PLS fraction evergreen and deciduous
#out.path = path to save figure
#var = string for "npp"/"trn"/"alb" to save file output
#
#Output:
#bias = absolute magnitude of the difference between CMIP5 and the PLS-reconstructed variable
#val  = PLS-reconstructed variable
#perc = percent bias (bias/original value)
#stats = table with mean, SD for absolute and percent biases
#Jaclyn Hatala Matthes, 5/14/15

cmip5_pls_flux_bias <- function(model.names, cmip5.global.pft, cmip5.global.clim, 
                                var.global.gcm, var.reg.gcm, cmip5.reg.clim, pls.evg, pls.dec, out.path, var){
  
  var.bias <- var.calc <- var.bias.frac <- list()
  var.val.mean <- var.val.sd <- var.perc.mean <- var.perc.sd <- vector()
  
  for(f in 1:length(model.names)){
  
    if(length(var.global.gcm[[f]])==1){ #not all models have flux variables, so there is a single NA value
      var.bias[[f]] <- NA
      var.calc[[f]] <- NA
      var.val.mean[f] <- var.val.sd[f] <- var.perc.mean[f] <- var.perc.sd[f] <- NA
    } else { #if there are values
      
      #convert global matrices into vectors for easier calculations
      tas.gcm <- as.vector(cmip5.global.clim$tas[[f]])
      pr.gcm  <- as.vector(cmip5.global.clim$pr[[f]])
      var.gcm <- as.vector(var.global.gcm[[f]]) 
      dec.gcm <- as.vector(cmip5.global.pft$dec.map[[f]])
      evg.gcm <- as.vector(cmip5.global.pft$evg.map[[f]])
      tas.reg <- cmip5.reg.clim$tas[[f]]
      pr.reg  <- cmip5.reg.clim$pr[[f]]
      var.reg <- var.reg.gcm[[f]]
      
      #find global high evergreen/deciduous grid cells
      hi.evg  <- which((evg.gcm/(evg.gcm+dec.gcm))>0.8 & !is.na(var.gcm))
      hi.dec  <- which((dec.gcm/(evg.gcm+dec.gcm))>0.8 & !is.na(var.gcm))
      
      #interpolate relationship between MAT, MAP, and flux variable for high evg/dec cells
      gcm.evg <- cbind(tas.gcm[hi.evg],pr.gcm[hi.evg],var.gcm[hi.evg])
      gcm.evg <- data.frame(gcm.evg)
      im.evg <- with(gcm.evg,interp(tas.gcm[hi.evg],pr.gcm[hi.evg],var.gcm[hi.evg],
                                    xo=seq(min(tas.gcm[hi.evg]),max(tas.gcm[hi.evg]),length=20),
                                    yo=seq(min(pr.gcm[hi.evg]),max(pr.gcm[hi.evg]),length=20)))
      image.plot(im.evg,xlab="MAT",ylab="MAP",main=paste("Evergreen: ",model.names[f],sep=""))
      
      gcm.dec <- cbind(tas.gcm[hi.dec],pr.gcm[hi.dec],var.gcm[hi.dec])
      gcm.dec <- data.frame(gcm.dec)
      im.dec <- with(gcm.dec,interp(tas.gcm[hi.dec],pr.gcm[hi.dec],var.gcm[hi.dec],
                                    xo=seq(min(tas.gcm[hi.dec]),max(tas.gcm[hi.dec]),length=20),
                                    yo=seq(min(pr.gcm[hi.dec]),max(pr.gcm[hi.dec]),length=20)))
      image.plot(im.dec,xlab="MAT",ylab="MAP",main=paste("Deciduous: ",model.names[f],sep=""))
      print(f)
      var.bias[[f]] <- matrix(nrow=80,ncol=30)
      var.calc[[f]] <- matrix(nrow=80,ncol=30)
      
      #loop over PLS points and calculate model flux variable with sett veg PFTs
      for(p in 1:length(pls.evg)){ 
        if(!is.na(pls.evg[p])){
          #find nearest value of the global interpolated MAT/MAP/variable relationship to the GCM MAT & MAP
          t.ind.d <- which.min(abs(im.dec$x - tas.reg[p]))
          p.ind.d <- which.min(abs(im.dec$y - pr.reg[p]))
          var.d   <- pls.dec[p]*im.dec$z[t.ind.d,p.ind.d]
          
          t.ind.e <- which.min(abs(im.evg$x - tas.reg[p]))
          p.ind.e <- which.min(abs(im.evg$y - pr.reg[p]))
          var.e   <- pls.evg[p]*im.evg$z[t.ind.e,p.ind.e]
          var.bias[[f]][p] <- var.reg[p] - (var.e+var.d)          
          var.calc[[f]][p] <- var.e + var.d
        } else {
          var.bias[[f]][p] <- NA
          var.calc[[f]][p] <- NA
        }
      }
      var.val.mean[f]  <- median(as.vector(var.bias[[f]]),na.rm=TRUE)
      var.val.sd[f]    <- sd(as.vector(var.bias[[f]]),na.rm=TRUE)
      var.bias.frac[[f]]    <- var.bias[[f]]/var.reg #map of fraction (%) bias
      var.bias.frac[[f]][!is.finite(var.bias.frac[[f]])] <- NA #remove infinite values
      var.bias.frac[[f]][var.bias.frac[[f]] < -1] <- NA
      var.perc.mean[f] <- median(as.vector(var.bias.frac[[f]]),na.rm=TRUE)
      var.perc.sd[f]   <- sd(as.vector(var.bias.frac[[f]]),na.rm=TRUE)
    } 
  }
  
  #data table - mean stats
  dat <- cbind(model.names,var.val.mean,var.val.sd,var.perc.mean*100,var.perc.sd*100)
  write.csv(dat,paste(out.path,"bias_stats_",var,".csv",sep=""))
  
  cmip5.bias <- list()
  cmip5.bias$bias <- var.bias
  cmip5.bias$val  <- var.calc
  cmip5.bias$frac <- var.bias.frac
  cmip5.bias$stats <- dat
  cmip5.bias
}
#Calculate the climate similarity between the CRUNCEP data product and CMIP5 models. 
#Requires 'phyloclim' R library to calculate niche overlap in climate space. 
#Input:
#esv.evg/.dec = ESV PFT fraction evergreen and deciduous
#cmip5.pft    = CMIP5 PFT fraction for evergreen and deciduous 
#cru.clim     = CRUNCEP 30-year climatology for tas (mean annual temp) and pr (mean annual precipitation)
#cmip5.clim   = CMIP5 30-year climatology for tas (mean annual temp) and pr (mean annual precipitation)
#
#Output:
#bias = absolute magnitude of the difference between CMIP5 and the esv-reconstructed variable
#val = esv-reconstructed variable
#Jaclyn Hatala Matthes, 5/18/15
cruncep_cmip5_clim_overlap <- function(esv.evg, esv.dec, cmip5.pft, cru.clim, cmip5.clim, model.full, out.path){
  
  #format ESV dataset as vectors and to remove NAs
  z1.esv <- as.vector(esv.evg)
  z2.esv <- as.vector(esv.dec)
  x.esv  <- as.vector(cru.clim$tas)
  y.esv  <- as.vector(cru.clim$pr)
  esv.comp.set1 <- c(!is.na(z1.esv) & !is.na(x.esv) & !is.na(y.esv))
  esv.comp.set2 <- c(!is.na(z2.esv) & !is.na(x.esv) & !is.na(y.esv))
  
  #calculate climate correlation between CRUNCEP and CMIP5 models
  r.tas <- r.pr <- t.bias <- t.sd <- p.bias <- p.sd <- vector()
  coef.tas <- coef.pr <-  matrix(nrow=length(cmip5.clim$tas),ncol=4)
  for(f in 1:length(cmip5.clim$tas)){
    z1 <- as.vector(cmip5.pft$evg.map[[f]])
    z2 <- as.vector(cmip5.pft$dec.map[[f]])
    x <- as.vector(cmip5.clim$tas[[f]])
    y <- as.vector(cmip5.clim$pr[[f]])
    comp.set1 <- c(!is.na(z1) & !is.na(x) & !is.na(y))
    comp.set2 <- c(!is.na(z2) & !is.na(x) & !is.na(y))
    
    #calculate correlation between tas and pr output
    r.tas[f] <- cor(x.esv[esv.comp.set1],as.vector(cmip5.clim$tas[[f]])[esv.comp.set1]) 
    r.pr[f]  <- cor(y.esv[esv.comp.set1],as.vector(cmip5.clim$pr[[f]])[esv.comp.set1]) 
    
    #calculate bias as mean difference in tas and pr
    t.bias[f] <- mean(as.vector(cmip5.clim$tas[[f]])[esv.comp.set1]-x.esv[esv.comp.set1],na.rm=TRUE)
    t.sd[f]   <- sd(as.vector(cmip5.clim$tas[[f]])[esv.comp.set1]-x.esv[esv.comp.set1],na.rm=TRUE)
    
    p.bias[f] <- mean(as.vector(cmip5.clim$pr[[f]])[esv.comp.set1]-y.esv[esv.comp.set1],na.rm=TRUE)
    p.sd[f]   <- sd(as.vector(cmip5.clim$pr[[f]])[esv.comp.set1]-y.esv[esv.comp.set1],na.rm=TRUE)
    
    #also save linear regression coefs for bias comparison
    f.tas <- lm(as.vector(cmip5.clim$tas[[f]])[esv.comp.set1] ~ x.esv[esv.comp.set1])
    f.pr  <- lm(as.vector(cmip5.clim$pr[[f]])[esv.comp.set1] ~ y.esv[esv.comp.set1])
    coef.tas[f,] <- c(f.tas$coef[1],summary(f.tas)$coefficients[[3]],f.tas$coef[2],summary(f.tas)$coefficients[[4]])
    coef.pr[f,]  <- c(f.pr$coef[1],summary(f.pr)$coefficients[[3]],f.pr$coef[2],summary(f.pr)$coefficients[[4]])
  }
  
  #organize & export summary statistics to table
  stats.out <- cbind(coef.tas,coef.pr,r.tas,r.pr,t.bias,t.sd,p.bias,p.sd)
  colnames(stats.out) <- c("temp_int","temp_int_se","temp_sl","temp_sl_se","prep_int",
                           "prep_int_se","prep_sl","prep_sl_se","r_temp","r_prep",
                           "t_bias","t_bias_sd","p_bias","p_bias_sd")
  rownames(stats.out) <- cmip5.clim$mod
  
  write.csv(stats.out,file=paste(out.path,"clim_comp_stats.csv",sep=""),col.names=TRUE,row.names=TRUE)
  
  #output variables
  clim.comp <- list()
  clim.comp$stats  <- stats.out
  clim.comp
}
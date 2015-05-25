#This workflow should be executed after running format_PLS_PFT.R to translate the resolution of the PLS data. 
#This code conducts comparisons between the ESV dataset and CMIP5 models by examining:
#1. PFT spatial extent (basic mapping & Taylor plots, taylor.diagram.2.R)
#2. PFT abundance in climate space (cruncep_climatology.R and cmip5_climatology.R) using Helliker distance (niche_test.R)
#3. Calculation of bias in land-atmosphere exchange as a result of mis-mapping PFTs (cmip5_pls_flux_bias.R).
#4. Calculation of sensitivity between modeled variables (elasticity.R).
#Jaclyn Hatala Matthes, 5/11/15

library(ncdf)
library(spam)
library(maps)
library(fields)
library(plotrix)
library(abind)
library(akima)
library(plot3D)
library(phyloclim)

base.path   <- ""
pls.path    <- "path to ESV dataset"
cmip5.path  <- "path to CMIP5 regional re-gridded output"
cmip5.global.path <- "path to CMIP5 global original output"
cru.path    <- "path to CRUNCEP data"
out.path    <- "path to output tables/figures"

#set up CMIP5 model PFT codes
mod.pft <- data.frame(read.csv(paste(pls.path,"cmip5_model_pft_key.csv",sep=""),header=TRUE))

#load the ESV PFT data & normalize fractions to sum to 1
pls.nc    <- open.ncdf(paste(pls.path,"PLS_pft_gridded.nc",sep=""))
pls.pft   <- get.var.ncdf(pls.nc,"pft")
pls.lon   <- get.var.ncdf(pls.nc,"lon")
pls.lat   <- get.var.ncdf(pls.nc,"lat")
close.ncdf(pls.nc)
pls.mv    <- 1e30
pls.pft[pls.pft>=pls.mv] <- NA

pls.dec <- pls.pft[,,2]/(pls.pft[,,1]+pls.pft[,,2]) #ESV deciduous fraction
pls.evg <- pls.pft[,,1]/(pls.pft[,,1]+pls.pft[,,2]) #ESV evergreen fraction

#make mask for ESV extent of lo-res versions 
pls.mask <- pls.dec
pls.mask[pls.mask>0] <- 1

#load low-res ESV PFT data & normalize fractions to sum to 1
pls.lo.nc    <- open.ncdf(paste(pls.path,"PLS_pft_gridded_lo.nc",sep=""))
pls.lo.pft   <- get.var.ncdf(pls.lo.nc,"pft")
close.ncdf(pls.lo.nc)
pls.mv    <- 1e30
pls.lo.pft[pls.pft>=pls.mv] <- NA

pls.lo.dec <- (pls.lo.pft[,,2]/(pls.lo.pft[,,1]+pls.lo.pft[,,2])) #ESV lo-res deciduous fraction
pls.lo.dec[is.na(pls.lo.dec)] <- 0 
pls.lo.dec <- pls.lo.dec*pls.mask
pls.lo.evg <- (pls.lo.pft[,,1]/(pls.lo.pft[,,1]+pls.lo.pft[,,2])) #ESV lo-res evergreen fraction
pls.lo.evg[is.na(pls.lo.evg)] <- 0
pls.lo.evg <- pls.lo.evg*pls.mask

#get set of CMIP5 models with PFT (landCoverFrac variable)
lcf.files <- list.files(paste(cmip5.path,"landCoverFrac/",sep=""),pattern = "\\.nc$")
lcf.split <- strsplit(lcf.files,"_")
models <- model.names <- model.full <- vector() #get list of CMIP5 models
for(m in 1:length(lcf.split)){
  tmp       <- strsplit(lcf.split[[m]][3],"-")
  models[m] <- tmp[[1]][1]
  model.names[m] <- paste(lcf.split[[m]][3],lcf.split[[m]][4],sep="_")
  model.full[m] <- paste(strsplit(lcf.split[[m]],"_")[[3]][1],"_",sep="")
}

#aggregate CMIP5 PFTs into deciduous/evergreen
cmip5.pft <- format_cmip5_pft(cmip5.path,model.full,mod.pft)

#2. Analyze PFT fractions in climate space
#aggregate monthly CRUNCEP product (6-hour resolution) into 30-year mean annual temperature and mean annual total precipitaiton 
cru.clim   <- cruncep_climatology(cru.path,pls.mask)

#aggregate CMIP5 files (mean monthly resolution) into 30-year mean annual temperature and mean annual total precipitaiton 
cmip5.clim <- cmip5_climatology(cmip5.path,pls.mask,model.full)

#3. Compare climate niches for ESV-CRUNCEP and CMIP5 modeled PFTs
#niche comparison between ESV PFTs in CRUNCEP climate space and CMIP5 PFTs in climate space
niche.comp <- esv_cmip5_niche_overlap(pls.evg, pls.dec, cmip5.pft, cru.clim, cmip5.clim, model.full, samp=5000, out.path)

#compare CRUNCEP and CMIP5 MAT/MAP 
clim.comp <- cruncep_cmip5_clim_overlap(pls.evg, pls.dec, cmip5.pft, cru.clim, cmip5.clim, model.full, out.path)

#4. Compare CMIP5 modeled NPP, transpiration, albedo by global MAT/MAP deciduous/evergreen PFT relationships
#aggregate regional mean annual net primary productivity (npp), transpiration (tran), and albedo (alb) to match climatology
cmip5.reg.fluxes <- cmip5_land_fluxes(cmip5.path,pls.mask,model.full)

#aggregate CMIP5 global PFTs
cmip5.global.pft <- format_cmip5_pft(cmip5.global.path,model.full,mod.pft,1)

#aggregate CMIP5 global climatology (original: mean monthly resolution) into 30-year MAT, MAP
cmip5.global.clim <- cmip5_climatology(cmip5.global.path,1,model.full)

#aggregate CMIP5 global fluxes (original: mean monthly resolution) into annual NPP, transpiration, and mean albedo
cmip5.global.fluxes <- cmip5_land_fluxes(cmip5.global.path,1,model.full)

#calculate bias due to mis-mapping PFTs in CMIP5 models
cmip5.npp.bias <- cmip5_pls_flux_bias(model.names, cmip5.global.pft, cmip5.global.clim, var.global.gcm=cmip5.global.fluxes$npp, 
                                      var.reg.gcm=cmip5.reg.fluxes$npp, cmip5.reg.clim=cmip5.clim, pls.evg, pls.dec, out.path, "npp")

cmip5.trn.bias <- cmip5_pls_flux_bias(model.names, cmip5.global.pft, cmip5.global.clim, var.global.gcm=cmip5.global.fluxes$tran, 
                                      var.reg.gcm=cmip5.reg.fluxes$tran, cmip5.reg.clim=cmip5.clim, pls.evg, pls.dec, out.path, "tran")

cmip5.alb.bias <- cmip5_pls_flux_bias(model.names, cmip5.global.pft, cmip5.global.clim, var.global.gcm=cmip5.global.fluxes$alb, 
                                      var.reg.gcm=cmip5.reg.fluxes$alb, cmip5.reg.clim=cmip5.clim, pls.evg, pls.dec, out.path, "alb")

#calculate post-hoc sensitivities between CMIP5 modeled variables
dec.tas.sens <- elasticity(cmip5.clim$tas,cmip5.pft$dec.map)
dec.pr.sens  <- elasticity(cmip5.clim$pr,cmip5.pft$dec.map)

tas.alb.sens <- elasticity(cmip5.clim$tas,cmip5.reg.fluxes$alb)
tas.npp.sens <- elasticity(cmip5.clim$tas,cmip5.reg.fluxes$npp)
tas.trn.sens <- elasticity(cmip5.clim$tas,cmip5.reg.fluxes$tran)

pr.alb.sens <- elasticity(cmip5.clim$pr,cmip5.reg.fluxes$alb)
pr.npp.sens <- elasticity(cmip5.clim$pr,cmip5.reg.fluxes$npp)
pr.trn.sens <- elasticity(cmip5.clim$pr,cmip5.reg.fluxes$tran)

dec.alb.sens <- elasticity(cmip5.pft$dec.map,cmip5.reg.fluxes$alb)
dec.npp.sens <- elasticity(cmip5.pft$dec.map,cmip5.reg.fluxes$npp)
dec.trn.sens <- elasticity(cmip5.pft$dec.map,cmip5.reg.fluxes$tran)

#make elasticity matrix across models and variables
dat <- cbind(dec.tas.sens$elas,dec.pr.sens$elas,tas.alb.sens$elas,tas.npp.sens$elas,tas.trn.sens$elas,
             pr.alb.sens$elas,pr.npp.sens$elas,pr.trn.sens$elas,dec.alb.sens$elas,dec.npp.sens$elas,dec.trn.sens$elas)
colnames(dat) <- c("tas-dec","pr-dec","tas-alb","tas-npp","tas-trn","pr-alb","pr-npp","pr-trn","dec-alb","dec-npp","dec-trn")
rownames(dat) <- mod.lab
write.csv(dat,paste(out.path,"cmip5_elasticity.csv",sep=""))

#aggregate bias stats into single matrix
n <- 0
for(f in 1:length(cmip5.npp.bias$bias)){
  if(length(cmip5.npp.bias$bias[[f]])>1){
    if(n==0){
      npp.all <- as.vector(cmip5.npp.bias$bias[[f]])
      npp.all.per <- as.vector(cmip5.npp.bias$frac[[f]])
      n <- n+1
    } else{
      npp.all <- cbind(npp.all,as.vector(cmip5.npp.bias$bias[[f]]))
      npp.all.per <- cbind(npp.all.per,as.vector(cmip5.npp.bias$frac[[f]]))
      n <- n+1
    }
  }
}

n <- 0
for(f in 1:length(cmip5.trn.bias$bias)){
  if(length(cmip5.trn.bias$bias[[f]])>1){
    if(n==0){
      trn.all <- as.vector(cmip5.trn.bias$bias[[f]])
      trn.all.per <- as.vector(cmip5.trn.bias$frac[[f]])
      n <- n+1
    } else{
      trn.all <- cbind(trn.all,as.vector(cmip5.trn.bias$bias[[f]]))
      trn.all.per <- cbind(trn.all.per,as.vector(cmip5.trn.bias$frac[[f]]))
      n <- n+1
    }
  }
}

n <- 0
for(f in 1:length(cmip5.alb.bias$bias)){
  if(length(cmip5.alb.bias$bias[[f]])>1){
    if(n==0){
      alb.all <- as.vector(cmip5.alb.bias$bias[[f]])
      alb.all.per <- as.vector(cmip5.alb.bias$frac[[f]])
      n <- n+1
    } else{
      alb.all <- cbind(alb.all,as.vector(cmip5.alb.bias$bias[[f]]))
      alb.all.per <- cbind(alb.all.per,as.vector(cmip5.alb.bias$frac[[f]]))
      n <- n+1
    }
  }
}

##### PLOT FIGURES
mod.lab <- c("ACCESS1-0","ACCESS1-3","HadGEM2-CC","HadGEM2-ES",
             "IPSL-CM5A-LR","IPSL-CM5A-MR","IPSL-CM5B-LR","MIROC-ESM-CHEM","MIROC-ESM",
             "MPI-ESM-LR","MPI-ESM-MR","MPI-ESM-P")

#FIGURE 1: plot PLS PFT maps
pdf(paste(out.path,"ESV_PFT_map.pdf",sep=""),height=10,width=10)
par(mfrow=c(2,1))
par(mar=c(5,5,1,2))
image.plot(pls.lon,pls.lat,pls.evg,xlab="",ylab="",main="ESV Fraction Evergreen",col=rev(terrain.colors(20)))
map("state",xlim=range(pls.lon),ylim=range(pls.lat),add=TRUE)
map("world",xlim=range(pls.lon),ylim=range(pls.lat),add=TRUE)
par(mar=c(5,5,1,2))
image.plot(pls.lon,pls.lat,pls.dec,xlab="",ylab="",main="ESV Fraction Deciduous",col=rev(terrain.colors(20)))
map("state",xlim=range(pls.lon),ylim=range(pls.lat),add=TRUE)
map("world",xlim=range(pls.lon),ylim=range(pls.lat),add=TRUE)
dev.off()

#FIGURE S1: plot CMIP5 PFT maps 
pdf(paste(out.path,"CMIP5_map.pdf",sep=""),height=16,width=12)
for(p in 1:length(lcf.evg.map)){
  par(mfrow=c(2,1))
  image.plot(pls.lon,pls.lat,cmip5.pft$evg.map[[p]],xlab="",ylab="",main=paste(model.full[p],": Fraction Evergreen",sep=""),
             col=rev(terrain.colors(20)),zlim=c(0,1))
  map("state",xlim=range(pls.lon),ylim=range(pls.lat),add=TRUE)
  map("world",xlim=range(pls.lon),ylim=range(pls.lat),add=TRUE)
  
  image.plot(pls.lon,pls.lat,cmip5.pft$dec.map[[p]],xlab="",ylab="",main=paste(model.full[p],": Fraction Deciduous",sep=""),
             col=rev(terrain.colors(20)),zlim=c(0,1))
  map("state",xlim=range(pls.lon),ylim=range(pls.lat),add=TRUE)
  map("world",xlim=range(pls.lon),ylim=range(pls.lat),add=TRUE)
}
dev.off()

#FIGURE 2: Make Taylor plots to compare the PLS and model PFTs
#format PLS data for Taylor plots
ref.dec <- as.vector(pls.dec*pls.mask)
ref.dec <- ref.dec[complete.cases(ref.dec)]
ref.evg <- as.vector(pls.evg*pls.mask)
ref.evg <- ref.evg[complete.cases(ref.evg)]
ref.lo.dec <- as.vector(pls.lo.dec)
ref.lo.dec <- ref.lo.dec[complete.cases(ref.lo.dec)]
ref.lo.evg <- as.vector(pls.lo.evg)
ref.lo.evg <- ref.lo.evg[complete.cases(ref.lo.evg)]

#CMIP5-PLS Taylor plots of PFT comparison
pdf(paste(out.path,"taylor_legend_v2.pdf",sep=""),height=5,width=20)
plot(1:10,1:10,pch="")
leg.text <- c("ACCESS1-0","ACCESS1-3","HadGEM2-CC","HadGEM2-ES",
              "IPSL-CM5A-LR","IPSL-CM5A-MR","MIROC-ESM-CHEM","MIROC-ESM",
              "MPI-ESM-LR","MPI-ESM-MR","MPI-ESM-P","Re-scaled ESV")
legend(3,10,leg.text,pch=c(1:length(diff.use),16),cex=1.3,ncol=3)
dev.off()

pdf(paste(out.path,"taylor_dec_v2.pdf",sep=""),height=6,width=8)
diff.use <- c(1,2,3,4,5,6,8,9,10,11,12)
taylor.diagram.2(ref.dec[which(!is.na(cmip5.pft$dec[[diff.use[1]]]))],
                 cmip5.pft$dec[[1]][which(!is.na(cmip5.pft$dec[[diff.use[1]]]))],
                 col=1,pch=1,pos.cor=FALSE,show.gamma=4,
                 main="Pre-industrial fraction deciduous",pcex=1.3)
for(p in 2:length(diff.use)){
  taylor.diagram.2(ref.dec[which(!is.na(cmip5.pft$dec[[diff.use[p]]]))],
                   cmip5.pft$dec[[diff.use[p]]][which(!is.na(cmip5.pft$dec[[diff.use[p]]]))],
                   add=TRUE,col=1,pch=p,pos.cor=FALSE,pcex=1.3)
}
taylor.diagram.2(ref.dec[which(!is.na(ref.lo.dec))],ref.lo.dec[which(!is.na(ref.lo.dec))],col=1,pos.cor=FALSE,pcex=1.3,add=TRUE)
dev.off()

#FIGURE S2: plot ESV/CMIP5 PFT distributions in MAT/MAP space
z1.pls <- as.vector(pls.evg)
z2.pls <- as.vector(pls.dec)
x.pls  <- as.vector(cru.clim$tas)
y.pls  <- as.vector(cru.clim$pr)
comp.set <- c(!is.na(z1.pls) & !is.na(z2.pls) & 
                !is.na(x.pls) & !is.na(y.pls))
x.pls <- x.pls[comp.set]
y.pls <- y.pls[comp.set]
z1.pls <- z1.pls[comp.set]
z2.pls <- z2.pls[comp.set]

#plot ESV & CMIP5 PFT-climate relationships
pls.scale <- 10
pdf(paste(out.path,"CMIP5_climspace.pdf",sep=""),width=22,height=12)
par(mfrow=c(1,2),mar=rep(7,4))

scatter2D(x.pls,y.pls,colvar=z1.pls,main="ESV evergreen fraction",
          xlab="Mean annual temperature [C]", cex.axis=2,cex.lab=2,
          ylab="Annual precipitation [mm/yr]",pch=16,cex.main=2,
          ylim=c(400,1750),xlim=c(-2,18),clim=c(0,1.0))
scatter2D(x.pls,y.pls,colvar=z2.pls,main="ESV deciduous fraction",
          xlab="Mean annual temperature [C]",
          ylab="",pch=16, cex.axis=2,cex.lab=2,cex.main=2,
          clim=c(0,1.0),ylim=c(400,1750),xlim=c(-2,18))

for(f in 1:12){
  z1 <- as.vector(cmip5.pft$evg.map[[f]])
  z2 <- as.vector(cmip5.pft$dec.map[[f]])
  x <- as.vector(cmip5.clim$tas[[f]])
  y <- as.vector(cmip5.clim$pr[[f]])
  comp.set <- c(!is.na(z1) & !is.na(z2) & !is.na(x) & !is.na(y))
  x  <- x[comp.set]
  y  <- y[comp.set]
  z1 <- z1[comp.set]
  z2 <- z2[comp.set]
  
  #evergreen
  image.plot(im.pls$x,im.pls$y,im.pls$z,col = rev(grey(seq(0, 1, length = pls.scale))),
             xlab="Mean annual temperature [C]",cex.axis=2,cex.lab=2, cex.main=2,
             #             legend.width=0,
             main=paste(model.full[f],": evergreen fraction",sep=""),
             ylab="Annual precipitation [mm/yr]",ylim=c(400,1750),
             xlim=c(-2,18),zlim=c(0,1.0))
  scatter2D(x,y,colvar=z1,main=paste(model.names[f],": % evergreen",sep=""),
            xlab="Mean annual temperature [C]",cex.axis=2,cex.lab=2,
            ylab="Annual precipitation [mm/yr]",pch=16,add=TRUE,
            colkey=FALSE,clim=c(0,1.0))
  
  #deciduous          
  image.plot(im.plsd$x,im.plsd$y,im.plsd$z,col = rev(grey(seq(0, 1, length = pls.scale))),
             xlab="Mean annual temperature [C]",cex.axis=2,
             #             legend.width=0,
             cex.lab=2,cex.main=2,
             main=paste(model.full[f],": deciduous fraction",sep=""),
             ylab="",ylim=c(400,1750),xlim=c(-2,18),zlim=c(0,1.0))
  scatter2D(x,y,colvar=z2,main=paste(model.names[f],": % deciduous"),
            cex.lab=2,
            xlab="Mean annual temperature [C]",cex.axis=2,add=TRUE,
            ylab="",pch=16,clim=c(0,1.0),colkey=FALSE) 
}
dev.off()

#FIGURE 3: plot climate overlap
pdf(paste(out.path,"climate_overlap.pdf",sep=""),width=8,height=12)
par(mfrow=c(2,1))
par(mar=c(10,5,1,2))
#tas bias
plot(clim.comp$stats[,11],ylim=c(-4,5),ylab="MAT difference [C]",xaxt="n",xlab="",cex.lab=1.5,cex.axis=1.25)
for(p in 1:nrow(clim.comp$stats)){
  lines(rep(p,2),c(clim.comp$stats[p,11]-clim.comp$stats[p,12],clim.comp$stats[p,11]+clim.comp$stats[p,12]),col="gray")
}
lines(0:13,rep(0,14),lty=2)
axis(1, at=seq(1, 12, by=1), labels = FALSE)
text(0.5,4.75,"a)",cex=1.5,pos=4)
text(seq(1, 12, by=1) + par("usr")[1]+0.25, par("usr")[3] - 0.3, labels = mod.lab, 
     srt = 45, pos = 2, xpd = TRUE, cex=1.5)
#pr bias
par(mar=c(10,5,1,2))
plot(clim.comp$stats[,13],ylim=c(-50,425),ylab=expression(paste("MAP difference [mm yr"^"-1"~"]")),
     xaxt="n",xlab="",cex.lab=1.5,cex.axis=1.25)
for(p in 1:nrow(clim.comp$stats)){
  lines(rep(p,2),c(clim.comp$stats[p,13]-clim.comp$stats[p,14],clim.comp$stats[p,13]+clim.comp$stats[p,14]),col="gray")
}
lines(0:13,rep(0,14),lty=2)
axis(1, at=seq(1, 12, by=1), labels = FALSE)
text(0.5,410,"b)",cex=1.5,pos=4)
text(seq(1, 12, by=1) + par("usr")[1]+0.25, par("usr")[3] - 20.0, labels = mod.lab, 
     srt = 45, pos = 2, xpd = TRUE, cex=1.5)
dev.off()

#FIGURE 4: plot niche overlap
pdf(paste(out.path,"niche_overlap.pdf",sep=""),width=8,height=12)
par(mfrow=c(2,1))
par(mar=c(10,5,1,2))
#evergreen overlap
plot(niche.comp$stats[,2],ylim=c(0,1),ylab="Niche overlap",xaxt="n",xlab="",cex.lab=1.5,cex.axis=1.25)
for(p in 1:nrow(niche.comp$stats)){
  lines(rep(p,2),c(niche.comp$stats[p,4],niche.comp$stats[p,5]),col="gray")
}
points(1:12,niche.comp$stats[,1],pch="+",col="black",cex=2)
axis(1, at=seq(1, 12, by=1), labels = FALSE)
text(0.5,0.95, "a) Evergreen",cex=1.5,pos=4)
text(seq(1, 12, by=1) + par("usr")[1]+0.25, par("usr")[3] - 0.06, labels = mod.lab, 
     srt = 45, pos = 2, xpd = TRUE, cex=1.5)
legend(0.75,0.3,c("Random overlap","Actual overlap"),pch=c(1,3),col=1)

#deciduous overlap
par(mar=c(10,5,1,2))
plot(niche.comp$stats[,7],ylim=c(0,1),ylab="Niche overlap",xaxt="n",xlab="",cex.lab=1.5,cex.axis=1.25)
for(p in 1:nrow(niche.comp$stats)){
  lines(rep(p,2),c(niche.comp$stats[p,9],niche.comp$stats[p,10]),col="gray")
}
points(1:12,niche.comp$stats[,6],pch="+",col="black",cex=2)
axis(1, at=seq(1, 12, by=1), labels = FALSE)
text(0.5,0.95,"b) Deciduous",cex=1.5,pos=4)
text(seq(1, 12, by=1) + par("usr")[1]+0.25, par("usr")[3] - 0.06, labels = mod.lab, 
     srt = 45, pos = 2, xpd = TRUE, cex=1.5)
dev.off()

#FIGURE 6/S3: plots of albedo/NPP/transpiration bias across models
pdf(paste(out.path,"bias_plots_val.pdf",sep=""),width=8,height=16)
par(mfrow=c(3,1))
par(mar=c(10,5,3,2))
boxplot(alb.all,cex.lab=1.5,cex.axis=1.5,
        ylab="Absolute bias [fraction]",xaxt="n")
axis(1, at=seq(1, 12, by=1), labels = FALSE)
text(seq(1, 12, by=1) + par("usr")[1]+0.5, par("usr")[3] - 0.02, labels = mod.lab, 
     srt = 45, pos = 2, xpd = TRUE, cex=1.5)
text(1,0.14,"a) Albedo",cex=2)
lines(0:13,rep(0,14),lty=5)

par(mar=c(10,5,3,2))
boxplot(cbind(rep(NA,nrow(npp.all)),rep(NA,nrow(npp.all)),npp.all),
        cex.lab=1.5,cex.axis=1.5,ylab=expression(paste("Absolute bias [kg-C m"^"-2"~"yr"^"-1"~"]")),xaxt="n")
axis(1, at=seq(1, 12, by=1), labels = FALSE)
text(seq(1, 12, by=1) + par("usr")[1]+0.5, par("usr")[3] - 0.075, labels = mod.lab, 
     srt = 45, pos = 2, xpd = TRUE, cex=1.5)
text(2.5,0.55,"b) Net primary productivity",cex=2)
lines(0:13,rep(0,14),lty=5)

par(mar=c(10,5,3,2))
boxplot(cbind(rep(NA,nrow(trn.all)),rep(NA,nrow(trn.all)),
              rep(NA,nrow(trn.all)),rep(NA,nrow(trn.all)),trn.all),
        cex.lab=1.5,cex.axis=1.5,ylab=expression(paste("Absolute bias [mm yr"^"-1"~"]")),xaxt="n")
axis(1, at=seq(1, 12, by=1), labels = FALSE)
text(seq(1, 12, by=1) + par("usr")[1]+0.5, par("usr")[3] - 20.0, labels = mod.lab, 
     srt = 45, pos = 2, xpd = TRUE, cex=1.5)
text(1.5,250,"c) Transpiration",cex=2)
lines(0:13,rep(0,14),lty=5)
dev.off()

pdf(paste(out.path,"bias_plots_percent.pdf",sep=""),width=8,height=16)
par(mfrow=c(3,1))
par(mar=c(10,5,3,2))
boxplot(alb.all.per*100,ylim=c(-100,100),cex.lab=1.5,cex.axis=1.5,
        ylab="Percent bias",xaxt="n")
axis(1, at=seq(1, 12, by=1), labels = FALSE)
text(seq(1, 12, by=1) + par("usr")[1]+0.5, par("usr")[3] - 8.0, labels = mod.lab, 
     srt = 45, pos = 2, xpd = TRUE, cex=1.5)
text(1,100,"a) Albedo",cex=2)
lines(0:13,rep(0,14),lty=5)
par(mar=c(10,5,3,2))
boxplot(cbind(rep(NA,nrow(npp.all.per)),rep(NA,nrow(npp.all.per)),npp.all.per)*100,
        cex.lab=1.5,cex.axis=1.5,ylim=c(-100,100),ylab="Percent bias",xaxt="n")
axis(1, at=seq(1, 12, by=1), labels = FALSE)
text(seq(1, 12, by=1) + par("usr")[1]+0.5, par("usr")[3] - 8.0, labels = mod.lab, 
     srt = 45, pos = 2, xpd = TRUE, cex=1.5)
text(2.5,100,"b) Net primary productivity",cex=2)
lines(0:13,rep(0,14),lty=5)
par(mar=c(10,5,3,2))
boxplot(cbind(rep(NA,nrow(trn.all.per)),rep(NA,nrow(trn.all.per)),
              rep(NA,nrow(trn.all.per)),rep(NA,nrow(trn.all.per)),trn.all.per)*100,
        cex.lab=1.5,cex.axis=1.5,ylim=c(-100,100),ylab="Percent bias",xaxt="n")
axis(1, at=seq(1, 12, by=1), labels = FALSE)
text(seq(1, 12, by=1) + par("usr")[1]+0.5, par("usr")[3] - 8.0, labels = mod.lab, 
     srt = 45, pos = 2, xpd = TRUE, cex=1.5)
text(1.5,100,"c) Transpiration",cex=2)
lines(0:13,rep(0,14),lty=5)
dev.off()

#SUPP FIGURE S2: map the modeled, ESV-reconstructed, and bias in the CMIP5 fluxes
pdf(paste(out.path,"map_bias_plots_albedo.pdf",sep=""),width=10,height=16)
par(mfrow=c(3,1))
for(f in 1:length(model.full)){
  image.plot(pls.lon,pls.lat,cmip5.reg.fluxes$alb[[f]],main=paste("Modeled albedo: ",substr(model.full[f],1,nchar(model.full[f])-1),sep=""),
             zlim=c(0,0.35),xlab="Longitude",ylab="Latitude")
  map("state",xlim=range(pls.lon),ylim=range(pls.lat),add=TRUE)
  map("world",xlim=range(pls.lon),ylim=range(pls.lat),add=TRUE) 
  image.plot(pls.lon,pls.lat,cmip5.alb.bias$val[[f]],main=paste("ESV reconstructed albedo: ",substr(model.full[f],1,nchar(model.full[f])-1),sep=""),
             zlim=c(0,0.35),xlab="Longitude",ylab="Latitude")
  map("state",xlim=range(pls.lon),ylim=range(pls.lat),add=TRUE)
  map("world",xlim=range(pls.lon),ylim=range(pls.lat),add=TRUE)
  image.plot(pls.lon,pls.lat,cmip5.alb.bias$bias[[f]],main=paste("Albedo bias: ",substr(model.full[f],1,nchar(model.full[f])-1),sep=""),
             zlim=c(-0.2,0.2),xlab="Longitude",ylab="Latitude")
  map("state",xlim=range(pls.lon),ylim=range(pls.lat),add=TRUE)
  map("world",xlim=range(pls.lon),ylim=range(pls.lat),add=TRUE)
}
dev.off()

pdf(paste(out.path,"map_bias_plots_npp.pdf",sep=""),width=10,height=16)
par(mfrow=c(3,1))
for(f in 1:length(model.full)){
  if(length(cmip5.reg.fluxes$npp[[f]])>1){
    image.plot(pls.lon,pls.lat,cmip5.reg.fluxes$npp[[f]],main=paste("Modeled NPP: ",substr(model.full[f],1,nchar(model.full[f])-1),sep=""),
               #             zlim=c(0,0.35),
               xlab="Longitude",ylab="Latitude")
    map("state",xlim=range(pls.lon),ylim=range(pls.lat),add=TRUE)
    map("world",xlim=range(pls.lon),ylim=range(pls.lat),add=TRUE) 
    image.plot(pls.lon,pls.lat,cmip5.npp.bias$val[[f]],main=paste("ESV reconstructed NPP: ",substr(model.full[f],1,nchar(model.full[f])-1),sep=""),
               #             zlim=c(0,0.35),
               xlab="Longitude",ylab="Latitude")
    map("state",xlim=range(pls.lon),ylim=range(pls.lat),add=TRUE)
    map("world",xlim=range(pls.lon),ylim=range(pls.lat),add=TRUE)
    image.plot(pls.lon,pls.lat,cmip5.npp.bias$bias[[f]],main=paste("NPP bias: ",substr(model.full[f],1,nchar(model.full[f])-1),sep=""),
               #             zlim=c(-0.2,0.2),
               xlab="Longitude",ylab="Latitude")
    map("state",xlim=range(pls.lon),ylim=range(pls.lat),add=TRUE)
    map("world",xlim=range(pls.lon),ylim=range(pls.lat),add=TRUE)
  }
}
dev.off()

pdf(paste(out.path,"map_bias_plots_tran.pdf",sep=""),width=10,height=16)
par(mfrow=c(3,1))
for(f in 1:length(model.full)){
  if(length(cmip5.reg.fluxes$tran[[f]])>1){
    image.plot(pls.lon,pls.lat,cmip5.reg.fluxes$tran[[f]],main=paste("Modeled Transpiration: ",substr(model.full[f],1,nchar(model.full[f])-1),sep=""),
               #             zlim=c(0,0.35),
               xlab="Longitude",ylab="Latitude")
    map("state",xlim=range(pls.lon),ylim=range(pls.lat),add=TRUE)
    map("world",xlim=range(pls.lon),ylim=range(pls.lat),add=TRUE) 
    image.plot(pls.lon,pls.lat,cmip5.trn.bias$val[[f]],main=paste("ESV reconstructed Transpiration: ",substr(model.full[f],1,nchar(model.full[f])-1),sep=""),
               #             zlim=c(0,0.35),
               xlab="Longitude",ylab="Latitude")
    map("state",xlim=range(pls.lon),ylim=range(pls.lat),add=TRUE)
    map("world",xlim=range(pls.lon),ylim=range(pls.lat),add=TRUE)
    image.plot(pls.lon,pls.lat,cmip5.trn.bias$bias[[f]],main=paste("Transpiration bias: ",substr(model.full[f],1,nchar(model.full[f])-1),sep=""),
               #             zlim=c(-0.2,0.2),
               xlab="Longitude",ylab="Latitude")
    map("state",xlim=range(pls.lon),ylim=range(pls.lat),add=TRUE)
    map("world",xlim=range(pls.lon),ylim=range(pls.lat),add=TRUE)
  }
}
dev.off()


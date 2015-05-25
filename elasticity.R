#Calculate the sensitivity between two variables, var1 & var2 from a list of models.
#Output a vector with sensitivity and elasticity between the two variables for each model. 
#Jaclyn Hatala Matthes, 5/18/15

elasticity <- function(var1,var2){
  for(f in 1:length(var1)){
    if(length(var2[[f]])==1){
      sens[f] <- NA
      elas[f] <- NA
    } else{
      use <- which(!is.na(as.vector(var1[[f]]))&!is.na(as.vector(var2[[f]])))
      v1  <- as.vector(var1[[f]])[use]
      q1    <- quantile(v1)
      v2  <- as.vector(var2[[f]])[use]
      q2    <- quantile(v2)
      sens.spline <- splinefun(q1, q2, method = "monoH.FC")
      sens[f]     <- sens.spline(median(v1), 1)
      elas[f]     <- sens[f] / (median(v2) / median(v1))
    }
  }
  output <- list()
  output$sens <- sens
  output$elas <- elas
  output
}
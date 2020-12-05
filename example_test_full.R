




#full_D
pall <- c(0.6,0.45,0.12,0.12,0.02,0.02,0.02,0.01)
source("util_FC.R")
geno <- sim_geno_auto(allpar=pall,n=1000)

res <- geno_est_auto(geno=geno)
res

allres <- c()
for(i in 1:1000){
  geno <- sim_geno_auto(allpar=pall,n=1000)
  
  res <- geno_est_auto(geno=geno)
  allres <- rbind(allres,res)
  
}
colMeans(allres)

#full_H
pall <- c(0.6,0.45,0.12)
source("util_FC_m.R")

geno <- sim_geno_auto_m(allpar=pall,n=1000)

res <- geno_est_auto_m(geno=geno)
res

allres <- c()
for(i in 1:1000){
  geno <- sim_geno_auto_m(allpar=pall,n=1000)
  
  res <- geno_est_auto_m(geno=geno)
  allres <- rbind(allres,res)
  
}
colMeans(allres)

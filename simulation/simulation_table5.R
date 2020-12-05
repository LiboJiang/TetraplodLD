

#Table4_D_P
#data simualtion by model1
pall <- c(0.6,0.45,0,0,0,0,0,0)
source("util2.R")
source("util2_1.R")


#H_100
pv_H_100 <- c()
for(i in 1:1000){
  geno <- sim_geno_auto(allpar=pall,n=100)
  res <- geno_est_auto(geno=geno)
  pv_H_100 <- rbind(pv_H_100,res)
}
colMeans(pv_H_100)


#H_200
pv_H_200 <- c()
for(i in 1:1000){
  geno <- sim_geno_auto(allpar=pall,n=200)
  res <- geno_est_auto(geno=geno)
  pv_H_200 <- rbind(pv_H_200,res)
}
colMeans(pv_H_200)

#H_400
pv_H_400 <- c()
for(i in 1:1000){
  geno <- sim_geno_auto(allpar=pall,n=400)
  res <- geno_est_auto(geno=geno)
  pv_H_400 <- rbind(pv_H_400,res)
}
colMeans(pv_H_400)



#HP_100
pv_H_100 <- c()
for(i in 1:1000){
  geno <- sim_geno_autoP(allpar=pall,n=100)
  res <- geno_est_autoP(geno=geno)
  pv_H_100 <- rbind(pv_H_100,res)
}
colMeans(pv_H_100)


#HP_200
pv_H_200 <- c()
for(i in 1:1000){
  geno <- sim_geno_autoP(allpar=pall,n=200)
  res <- geno_est_autoP(geno=geno)
  pv_H_200 <- rbind(pv_H_200,res)
}
colMeans(pv_H_200)

#HP_400
pv_H_400 <- c()
for(i in 1:1000){
  geno <- sim_geno_auto(allpar=pall,n=400)
  res <- geno_est_auto(geno=geno)
  pv_H_400 <- rbind(pv_H_400,res)
}
colMeans(pv_H_400)



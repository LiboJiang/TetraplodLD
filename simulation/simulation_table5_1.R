

#Table4_H_F
#data simualtion by model1
pall <- c(0.6,0.45,0)
source("util2.R")
source("util2_1.R")


#H_100
pv_H_100 <- c()
for(i in 1:1000){
  geno <- sim_geno_auto_m(allpar=pall,n=100)
  res <- geno_est_auto_m(geno=geno)
  pv_H_100 <- rbind(pv_H_100,res)
}

colMeans(pv_H_100)

#H_200
pv_H_200 <- c()
for(i in 1:1000){
  geno <- sim_geno_auto_m(allpar=pall,n=200)
  res <- geno_est_auto_m(geno=geno)
  pv_H_200 <- rbind(pv_H_200,res)
}

colMeans(pv_H_200)


#H_400
pv_H_400 <- c()
for(i in 1:1000){
  geno <- sim_geno_auto_m(allpar=pall,n=400)
  res <- geno_est_auto_m(geno=geno)
  pv_H_400 <- rbind(pv_H_400,res)
}

colMeans(pv_H_400)





#Table4_H_P
#data simualtion by model1
pall <- c(0.6,0.45,0)
source("util3.R")
source("util3_1.R")


#H_100
pv_H_100 <- c()
for(i in 1:1000){
  geno <- sim_geno_autoP_m(allpar=pall,n=100)
  res <- geno_est_autoP_m(geno=geno)
  pv_H_100 <- rbind(pv_H_100,res)
}
colMeans(pv_H_100)


#H_200
pv_H_200 <- c()
for(i in 1:1000){
  geno <- sim_geno_autoP_m(allpar=pall,n=200)
  res <- geno_est_autoP_m(geno=geno)
  pv_H_200 <- rbind(pv_H_200,res)
}
colMeans(pv_H_200)




#H_400
pv_H_400 <- c()
for(i in 1:1000){
  geno <- sim_geno_autoP_m(allpar=pall,n=400)
  res <- geno_est_autoP_m(geno=geno)
  pv_H_400 <- rbind(pv_H_400,res)
}
colMeans(pv_H_400)


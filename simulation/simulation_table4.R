

#Table4 

#data simualtion by model1
pall <- c(0.6,0.45,0.12)
source("../util_PC.R")
source("../util_PC_m.R")
#H_100
cret_H_100 <- c()
for(i in 1:1000){
  geno <- sim_geno_autoP_m(allpar=pall,n=100)
  res <- geno_est_autoP_m(geno=geno)
  pvr <- Power_cal_mp(res=res,geno=geno)
  cret_H_100 <- rbind(cret_H_100,c(res,pvr))
}
colMeans(cret_H_100)
apply(cret_H_100,2,sd)



#H_200
cret_H_200 <- c()
for(i in 1:1000){
  geno <- sim_geno_autoP_m(allpar=pall,n=200)
  res <- geno_est_autoP_m(geno=geno)
  pvr <- Power_cal_mp(res=res,geno=geno)
  cret_H_200 <- rbind(cret_H_200,c(res,pvr))
}
colMeans(cret_H_200)
apply(cret_H_200,2,sd)



#H_400
cret_H_400 <- c()
for(i in 1:1000){
  geno <- sim_geno_autoP_m(allpar=pall,n=400)
  res <- geno_est_autoP_m(geno=geno)
  pvr <- Power_cal_mp(res=res,geno=geno)
  cret_H_400 <- rbind(cret_H_400,c(res,pvr))
}
colMeans(cret_H_400)
apply(cret_H_400,2,sd)








#data simualtion by model1
pall <- c(0.6,0.45,0.12,0.12,0.02,0.02,0.02,0.01)
source("../util_PC.R")
source("../util_PC_m.R")
#D_100
cret_D_100 <- c();allp_100 <- c()
for(i in 1:1000){
  geno <- sim_geno_autoP(allpar=pall,n=100)
  res <- geno_est_autoP(geno=geno)
  allp_100 <- rbind(allp_100,Power_cal_p(res=res,geno=geno))
  cret_D_100 <- rbind(cret_D_100,res)
}
colMeans(cret_D_100)
apply(cret_D_100,2,sd)



#D_200
cret_D_200 <- c();allp_200 <- c()
for(i in 1:1000){
  geno <- sim_geno_autoP(allpar=pall,n=200)
  res <- geno_est_autoP(geno=geno)
  allp_200 <- rbind(allp_200,Power_cal_p(res=res,geno=geno))
  cret_D_200 <- rbind(cret_D_200,res)
}
colMeans(cret_D_200)
apply(cret_D_200,2,sd)


#D_400
cret_D_400 <- c();allp_400 <- c()
for(i in 1:1000){
  geno <- sim_geno_autoP(allpar=pall,n=400)
  res <- geno_est_autoP(geno=geno)
  allp_400 <- rbind(allp_400,Power_cal_p(res=res,geno=geno))
  cret_D_400 <- rbind(cret_D_400,res)
}
colMeans(cret_D_400)
apply(cret_D_400,2,sd)




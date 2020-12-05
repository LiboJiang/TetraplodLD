
#full_mix
pall <- c(0.6,0.45,0.12)
source("../util_FC.R")
source("../util_FC_m.R")
#H_simulation
dih_100 <- c()
for(ii in 1:20){
  cret_H_100 <- c()
  for(i in 1:1000){
    geno <- sim_geno_auto_m(allpar=pall,n=100)
    cret_H_100 <- rbind(cret_H_100,mix_est_f(geno=geno))
  }
  dih_100 <- c(dih_100,length(which(cret_H_100[,2]>0.05))/1000)
  cat("ii=",ii,"\n")
}
mean(dih_100)


dih_200 <- c()
for(ii in 1:20){
  cret_H_200 <- c()
  for(i in 1:1000){
    geno <- sim_geno_auto_m(allpar=pall,n=200)
    cret_H_200 <- rbind(cret_H_200,mix_est_f(geno=geno))
  }
  dih_200 <- c(dih_200,length(which(cret_H_200[,2]>0.05))/1000)
  cat("ii=",ii,"\n")
}
mean(dih_200)


dih_400 <- c()
for(ii in 1:20){
  cret_H_400 <- c()
  for(i in 1:1000){
    geno <- sim_geno_auto_m(allpar=pall,n=400)
    cret_H_400 <- rbind(cret_H_400,mix_est_f(geno=geno))
  }
  dih_400 <- c(dih_400,length(which(cret_H_400[,2]>0.05))/1000)
  cat("ii=",ii,"\n")
}
mean(dih_400)

#D_simulation
pall <- c(0.6,0.45,0.12,0.12,0.02,0.02,0.02,0.01)
source("../util_FC.R")
source("../util_FC_m.R")
did_100 <- c()
for(ii in 1:20){
  cret_D_100 <- c()
  for(i in 1:1000){
    geno <- sim_geno_auto(allpar=pall,n=100)
    cret_D_100 <- rbind(cret_D_100,mix_est_f(geno=geno))
  }
  did_100 <- c(did_100,length(which(cret_D_100[,2]>0.05))/1000)
  cat("ii=",ii,"\n")
}
mean(did_100)


did_200 <- c()
for(ii in 1:20){
  cret_D_200 <- c()
  for(i in 1:1000){
    geno <- sim_geno_auto(allpar=pall,n=200)
    cret_D_200 <- rbind(cret_D_200,mix_est_f(geno=geno))
  }
  did_200 <- c(did_200,length(which(cret_D_200[,2]>0.05))/1000)
  cat("ii=",ii,"\n")
}
mean(did_200)


did_400 <- c()
for(ii in 1:20){
  cret_D_400 <- c()
  for(i in 1:1000){
    geno <- sim_geno_auto(allpar=pall,n=400)
    cret_D_400 <- rbind(cret_D_400,mix_est_f(geno=geno))
  }
  did_400 <- c(did_400,length(which(cret_D_400[,2]>0.05))/1000)
  cat("ii=",ii,"\n")
}
mean(did_400)






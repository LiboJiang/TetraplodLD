

#full_D
pall <- c(0.6,0.45,0.12,0.12,0.02,0.02,0.02,0.01)
source("util_PC.R")
geno <- sim_geno_autoP(allpar=pall,n=1000)

res <- geno_est_autoP(geno=geno)
res

allres <- c()
for(i in 1:1000){
  geno <- sim_geno_autoP(allpar=pall,n=400)
  
  res <- geno_est_autoP(geno=geno)
  allres <- rbind(allres,res)
  
}
colMeans(allres)

#full_H
pall <- c(0.6,0.45,0.12)
source("util_PC_m.R")

geno <- sim_geno_autoP_m(allpar=pall,n=1000)

res <- geno_est_autoP_m(geno=geno)
res

allres <- c()
for(i in 1:1000){
  geno <- sim_geno_autoP_m(allpar=pall,n=1000)
  
  res <- geno_est_autoP_m(geno=geno)
  allres <- rbind(allres,res)
  
}
colMeans(allres)



#work_autoP

source("util_PC_m.R")
source("util_PC.R")

snpinfo <-readRDS("./data/annot.rds")
del <- unique(c(which(is.na(snpinfo$snp)),which(is.na(snpinfo$pos))))
minfo <- snpinfo[-del,]

load("./data/geno_up2_version2.RData")
load("./data/geno_low2_version2.RData")

ii <- match(colnames(geno_up2),colnames(geno_low2))
geno_up21 <- geno_up2[,which(!is.na(ii))]
geno_low21 <- geno_low2[,ii[which(!is.na(ii))]]

INFO1 <- match(colnames(geno_low21),minfo[,1])
geno_up21_1 <- geno_up21[,which(!is.na(INFO1))]
geno_low21_1 <- geno_low21[,which(!is.na(INFO1))]

minfo1 <- minfo[INFO1[which(!is.na(INFO1))],]


chr03a_i <- which(minfo1$chr=="Chr03a")

chr03a_m <- geno_low21_1[,chr03a_i]

LD <- work_test1_mix(M=chr03a_m,mn=minfo1[chr03a_i,1])
LDr2 <- r2c(pA=LD[,6],pB=LD[,7],D=LD[,10:13])
#write.csv(LD,file="LD_chr3a_400.csv")



boxplot(LD[,18])

#ti <- seq(1,dim(geno_up21_1)[2],2)
#geno_up22 <- geno_up21_1[,ti]
#geno_low22 <- geno_low21_1[,ti]

#minfo2 <- minfo1[ti,]
#dat <-list(up=geno_up22,low=geno_low22,info=minfo2)#5424
#save(dat,file="dat.RData",version = 2)







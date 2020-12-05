

source("util_PC_m.R")
source("util_PC.R")
source("LD.R")
load("dat.RData")



chr03a_i <- which(dat$info$chr=="Chr03a")
chr04a_i <- which(dat$info$chr=="Chr04a")

chr03a_mf <- dat$info[chr03a_i,]
chr04a_mf <- dat$info[chr04a_i,]
chr <- rbind(chr03a_mf,chr04a_mf)

chr03a_ml <- dat$low[,chr03a_i]
chr04a_ml <- dat$low[,chr04a_i]
chr_l <- cbind(chr03a_ml,chr04a_ml)

chr03a_mu <- dat$up[,chr03a_i]
chr04a_mu <- dat$up[,chr04a_i]
chr_u <- cbind(chr03a_mu,chr04a_mu)


LD_03_l <- work_test1_mix(M=chr03a_ml,mn=chr03a_mf[,1])
LD_04_l <- work_test1_mix(M=chr04a_ml,mn=chr04a_mf[,1])
LD_c_l <- work_test1_mix(M=chr_l,mn=chr[,1])


LD_03_u <- work_test1_mix(M=chr03a_mu,mn=chr03a_mf[,1])
LD_04_u <- work_test1_mix(M=chr04a_mu,mn=chr04a_mf[,1])
LD_c_u <- work_test1_mix(M=chr_u,mn=chr[,1])




#normalize
LDn2_03_ls <- LD_chr(dn=LD_03_l[,16:19],mi=chr03a_mf)
LDn2_04_ls <- LD_chr(dn=LD_04_l[,16:19],mi=chr04a_mf)

LDn2_c_ls <- LD_chr_inv(dn=LD_c_l[,16:19],mi=chr)
LDn2_c_ls1 <- Lcs(l1=LDn2_03_ls,l2=LDn2_04_ls,l3=LDn2_c_ls)


LDn2_03_us <- LD_chr(dn=LD_03_u[,16:19],mi=chr03a_mf)
LDn2_04_us <- LD_chr(dn=LD_04_u[,16:19],mi=chr04a_mf)

LDn2_c_us <- LD_chr_inv(dn=LD_c_u[,16:19],mi=chr)
LDn2_c_us1 <- Lcs(l1=LDn2_03_us,l2=LDn2_04_us,l3=LDn2_c_us)






y1=LDn2_03_ls
y2=LDn2_04_ls
y3=LDn2_c_ls1

x1=LDn2_03_us
x2=LDn2_04_us
x3=LDn2_c_us1



x1_1 <- x1[which(x1[,1]<1),];x1_2 <- x1[which(x1[,1]>1),];
x2_1 <- x2[which(x2[,1]<1),];x2_2 <- x2[which(x2[,1]>1),];

y1_1 <- y1[which(y1[,1]<1),];y1_2 <- y1[which(y1[,1]>1),];
y2_1 <- y2[which(y2[,1]<1),];y2_2 <- y2[which(y2[,1]>1),];

Deab_u <- c(mean(x1_1[,2]),mean(x1_2[,2]),mean(x3[,2]),mean(x2_1[,2]),mean(x2_2[,2]))
DAb_u <- c(mean(x1_1[,3]),mean(x1_2[,3]),mean(x3[,3]),mean(x2_1[,3]),mean(x2_2[,3]))
DaB_u <- c(mean(x1_1[,4]),mean(x1_2[,4]),mean(x3[,4]),mean(x2_1[,4]),mean(x2_2[,4]))
DAB_u <- c(mean(x1_1[,5]),mean(x1_2[,5]),mean(x3[,5]),mean(x2_1[,5]),mean(x2_2[,5]))

Deab_l <- c(mean(y1_1[,2]),mean(y1_2[,2]),mean(y3[,2]),mean(y2_1[,2]),mean(y2_2[,2]))
DAb_l <- c(mean(y1_1[,3]),mean(y1_2[,3]),mean(y3[,3]),mean(y2_1[,3]),mean(y2_2[,3]))
DaB_l <- c(mean(y1_1[,4]),mean(y1_2[,4]),mean(y3[,4]),mean(y2_1[,4]),mean(y2_2[,4]))
DAB_l <- c(mean(y1_1[,5]),mean(y1_2[,5]),mean(y3[,5]),mean(y2_1[,5]),mean(y2_2[,5]))


pdf("Figure_res.pdf",height = 6,width=8)
par(mar=c(0.5,3,0.8,0.5),oma=c(2,1,0,0),mfrow=c(2,2))
plot(NA,NA,pch=16,type="n",xlab=" ",ylab=" ",xlim=c(-0.1,4.6),ylim=c(-0.65,0.65),xaxt="n",yaxt="n",xaxs="i", yaxs="i")
#A
rect(0.2,0,0.8,Deab_u[1],border="grey",col="red")
rect(0.9,0,1.5,Deab_u[2],border="grey",col="blue")
rect(2.0,0,2.6,Deab_u[3],border="grey",col="yellow")
rect(3.1,0,3.7,Deab_u[4],border="grey",col="red")
rect(3.8,0,4.4,Deab_u[5],border="grey",col="blue")
segments(-100,0,100,0,col="grey",lwd=2)
rect(0.2,0,0.8,-Deab_l[1],border="grey",col="red")
rect(0.9,0,1.5,-Deab_l[2],border="grey",col="blue")
rect(2.0,0,2.6,-Deab_l[3],border="grey",col="yellow")
rect(3.1,0,3.7,-Deab_l[4],border="grey",col="red")
rect(3.8,0,4.4,-Deab_l[5],border="grey",col="blue")

text(2.25,0.65*0.9,"Upland",cex=1.5)
text(2.25,-0.65*0.9,"Lowland",cex=1.5)

axis(2,seq(0,0.6,0.2),seq(0,0.6,0.2),cex.axis=1.2,lwd=1,las=1)
axis(2,seq(-0.2,-0.6,-0.2),seq(0.2,0.6,0.2),cex.axis=1.2,lwd=1,las=1)
axis(1,c(0.85,2.3,3.75),c(" "," "," "),cex.axis=1.2,lwd=1,las=1)
text(0.2,0.65*0.9,"A",cex=1.5)
mtext(expression(italic(D)^"'"),2,cex=1.5,line=2)

par(mar=c(0.5,3,0.8,0.5))
plot(NA,NA,pch=16,type="n",xlab=" ",ylab=" ",xlim=c(-0.1,4.6),ylim=c(-0.41,0.41),xaxt="n",yaxt="n",xaxs="i", yaxs="i")
#B
rect(0.2,0,0.8,DAb_u[1],border="grey",col="red")
rect(0.9,0,1.5,DAb_u[2],border="grey",col="blue")
rect(2.0,0,2.6,DAb_u[3],border="grey",col="yellow")
rect(3.1,0,3.7,DAb_u[4],border="grey",col="red")
rect(3.8,0,4.4,DAb_u[5],border="grey",col="blue")
segments(-100,0,100,0,col="grey",lwd=2)
rect(0.2,0,0.8,-DAb_l[1],border="grey",col="red")
rect(0.9,0,1.5,-DAb_l[2],border="grey",col="blue")
rect(2.0,0,2.6,-DAb_l[3],border="grey",col="yellow")
rect(3.1,0,3.7,-DAb_l[4],border="grey",col="red")
rect(3.8,0,4.4,-DAb_l[5],border="grey",col="blue")

text(2.25,0.41*0.9,"Upland",cex=1.5)
text(2.25,-0.41*0.9,"Lowland",cex=1.5)

axis(2,seq(0,0.4,0.1),seq(0,0.4,0.1),cex.axis=1.2,lwd=1,las=1)
axis(2,seq(-0.1,-0.4,-0.1),seq(0.1,0.4,0.1),cex.axis=1.2,lwd=1,las=1)
axis(1,c(0.85,2.3,3.75),c(" "," "," "),cex.axis=1.2,lwd=1,las=1)
text(0.2,0.41*0.9,"B",cex=1.5)


par(mar=c(0.5,3,0.8,0.5))
plot(NA,NA,pch=16,type="n",xlab=" ",ylab=" ",xlim=c(-0.1,4.6),ylim=c(-0.41,0.41),xaxt="n",yaxt="n",xaxs="i", yaxs="i")
#C
rect(0.2,0,0.8,DaB_u[1],border="grey",col="red")
rect(0.9,0,1.5,DaB_u[2],border="grey",col="blue")
rect(2.0,0,2.6,DaB_u[3],border="grey",col="yellow")
rect(3.1,0,3.7,DaB_u[4],border="grey",col="red")
rect(3.8,0,4.4,DaB_u[5],border="grey",col="blue")
segments(-100,0,100,0,col="grey",lwd=2)
rect(0.2,0,0.8,-DaB_l[1],border="grey",col="red")
rect(0.9,0,1.5,-DaB_l[2],border="grey",col="blue")
rect(2.0,0,2.6,-DaB_l[3],border="grey",col="yellow")
rect(3.1,0,3.7,-DaB_l[4],border="grey",col="red")
rect(3.8,0,4.4,-DaB_l[5],border="grey",col="blue")

text(2.25,0.41*0.9,"Upland",cex=1.5)
text(2.25,-0.41*0.9,"Lowland",cex=1.5)

axis(2,seq(0,0.4,0.1),seq(0,0.4,0.1),cex.axis=1.2,lwd=1,las=1)
axis(2,seq(-0.1,-0.4,-0.1),seq(0.1,0.4,0.1),cex.axis=1.2,lwd=1,las=1)
axis(1,c(0.85,2.3,3.75),c("chr03a","chr03a¡Áchr04a","chr04a"),cex.axis=1.2,lwd=1,las=1)
text(0.2,0.41*0.9,"C",cex=1.5)
mtext(expression(italic(D)^"'"),2,cex=1.5,line=2)

par(mar=c(0.5,3,0.8,0.5))
plot(NA,NA,pch=16,type="n",xlab=" ",ylab=" ",xlim=c(-0.1,4.6),ylim=c(-0.35,0.35),xaxt="n",yaxt="n",xaxs="i", yaxs="i")
#D
rect(0.2,0,0.8,DAB_u[1],border="grey",col="red")
rect(0.9,0,1.5,DAB_u[2],border="grey",col="blue")
rect(2.0,0,2.6,DAB_u[3],border="grey",col="yellow")
rect(3.1,0,3.7,DAB_u[4],border="grey",col="red")
rect(3.8,0,4.4,DAB_u[5],border="grey",col="blue")
segments(-100,0,100,0,col="grey",lwd=2)
rect(0.2,0,0.8,-DAB_l[1],border="grey",col="red")
rect(0.9,0,1.5,-DAB_l[2],border="grey",col="blue")
rect(2.0,0,2.6,-DAB_l[3],border="grey",col="yellow")
rect(3.1,0,3.7,-DAB_l[4],border="grey",col="red")
rect(3.8,0,4.4,-DAB_l[5],border="grey",col="blue")

text(2.25,0.35*0.9,"Upland",cex=1.5)
text(2.25,-0.35*0.9,"Lowland",cex=1.5)

axis(2,seq(0,0.3,0.1),seq(0,0.3,0.1),cex.axis=1.2,lwd=1,las=1)
axis(2,seq(-0.1,-0.3,-0.1),seq(0.1,0.3,0.1),cex.axis=1.2,lwd=1,las=1)
axis(1,c(0.85,2.3,3.75),c("chr03a","chr03a¡Áchr04a","chr04a"),cex.axis=1.2,lwd=1,las=1)
text(0.2,0.35*0.9,"D",cex=1.5)
dev.off()







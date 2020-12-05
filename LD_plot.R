

#chr03a


load("tmp_res/tmp_chr03a.RData")
source("LD.R")


LD_dr <- LD_chr(dn=LDr2,mi=minfo1[chr03a_i,])
LD_dn <- LD_chr(dn=LD[,10:13],mi=minfo1[chr03a_i,])
#points

LD_points_chr(res=LD_dr[,c(1,2)],xll=seq(0,50,10),yll=seq(0,1,0.2),
                          h=5,w=7,filen="LD1_point_chr03a_Deab_r2_low.pdf")
LD_points_chr(res=LD_dr[,c(1,3)],xll=seq(0,50,10),yll=seq(0,0.08,0.02),
              h=5,w=7,filen="LD2_point_chr03a_DAb_r2_low.pdf")
LD_points_chr(res=LD_dr[,c(1,4)],xll=seq(0,50,10),yll=seq(0,0.08,0.02),
              h=5,w=7,filen="LD3_point_chr03a_DaB_r2_low.pdf")
LD_points_chr(res=LD_dr[,c(1,5)],xll=seq(0,50,10),yll=seq(0,0.3,0.05),
              h=5,w=7,filen="LD4_point_chr03a_DAB_r2_low.pdf")

#normalize

LD_pointsN_chr(res=LD_dn[,c(1,2)],xll=seq(0,50,10),yll=seq(0,1,0.2),
              h=5,w=7,filen="LD1_point_chr03a_Deab_N_low.pdf")
LD_pointsN_chr(res=LD_dn[,c(1,3)],xll=seq(0,50,10),yll=seq(0,0.1,0.02),
              h=5,w=7,filen="LD2_point_chr03a_DAb_N_low.pdf")
LD_pointsN_chr(res=LD_dn[,c(1,4)],xll=seq(0,50,10),yll=seq(0,0.1,0.02),
              h=5,w=7,filen="LD3_point_chr03a_DaB_N_low.pdf")
LD_pointsN_chr(res=LD_dn[,c(1,5)],xll=seq(0,50,10),yll=seq(0,0.1,0.02),
              h=5,w=7,filen="LD4_point_chr03a_DAB_N_low.pdf")



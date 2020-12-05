
DF10 <- function(pall){
  
  
  pA <- pall[1]
  pa <- 1-pA
  pB <- pall[2]
  pb <- 1-pB
  DA <- pall[3]
  DB <- pall[4]
  Deab <- pall[5]
  DAb <- pall[6]
  DaB <- pall[7]
  DAB <- pall[8]
  
  PAABB <- pA^2*pB^2+pA^2*DB+pB^2*DA+2*pA*pB*Deab+2*pB*DAb+2*pA*DaB+DA*DB+Deab^2+DAB
  
  PAABb <- 2*(pA^2*pB*pb-pA^2*DB+pB*pb*DA+(pA*pb-pA*pB)*Deab+(pb-pB)*DAb-2*pA*DaB)-2*(DA*DB+Deab^2+DAB)
  
  PAAbb <- pA^2*pb^2+pA^2*DB+pb^2*DA-2*pA*pb*Deab-2*pb*DAb+2*pA*DaB+DA*DB+Deab^2+DAB
  
  PAaBB <- 2*(pA*pa*pB^2+pA*pa*DB-pB^2*DA+(pa*pB-pA*pB)*Deab-2*pB*DAb+(pa-pA)*DaB)-2*(DA*DB+Deab^2+DAB)
  
  PAaBb <- 2*(pA*pa*pB*pb-pA*pa*DB-pB*pb*DA+(pA*pB+pa*pb)*Deab-(pA*pb+pa*pB)*Deab+(pB-pb)*DAb+(pA-pa)*DaB)+
    2*(pA*pa*pB*pb-pA*pa*DB-pB*pb*DA+(pB-pb)*DAb+(pA-pa)*DaB)+4*(DA*DB+Deab^2+DAB)
  
  PAabb <- 2*(pA*pa*pb^2+pA*pa*DB-pb^2*DA+(pA*pb-pa*pb)*Deab+2*pb*DAb+(pa-pA)*DaB)-2*(DA*DB+Deab^2+DAB)
  
  PaaBB <- pa^2*pB^2+pa^2*DB+pB^2*DA-2*pa*pB*Deab+2*pB*DAb-2*pa*DaB +DA*DB+Deab^2+DAB
  PaaBb <- 2*(pa^2*pB*pb-pa^2*DB+pB*pb*DA+(pa*pB-pa*pb)*Deab+(pb-pB)*DAb+2*pa*DaB)-2*(DA*DB+Deab^2+DAB)
  Paabb <- pa^2*pb^2+pb^2*DA+pa^2*DB+2*pa*pb*Deab-2*pb*DAb-2*pa*DaB+DA*DB+Deab^2+DAB
  
  tmp <- c(PAABB,PAABb,PAAbb,PAaBB,PAaBb,PAabb,PaaBB,PaaBb,Paabb)
  #sum(tmp)
  return(tmp)
}



autoP_DF10 <- function(pall){
  
  
  gp <- DF10(pall=pall)
  
  PAABB <- gp[1]
  PAABb <- gp[2]
  PAAbb <- gp[3] 
  PAaBB <- gp[4] 
  PAaBb <- gp[5] 
  PAabb <- gp[6] 
  PaaBB <- gp[7] 
  PaaBb <- gp[8] 
  Paabb <- gp[9] 
  
  zy1 <- (PAABB)^2 
  zy2 <- 2*PAABB*PAABb + 2*PAABB*PAAbb+(PAABb)^2 + 2*PAABb*PAAbb
  zy3 <- (PAAbb)^2
  
  zy4 <- 2*PAABB*PAaBB + (PAaBB)^2+2*PAABB*PaaBB + 2*PAaBB*PaaBB
  
  
  zy5 <- 2*PAABB*PAaBb+2*PAaBB*PAABb + 2*PAABb*PAaBb+2*PAaBB*PAAbb+2*PAabb*PAABB + 2*PAAbb*PAaBb+2*PAabb*PAABb+
    2*PAaBB*PAaBb+2*PAABB*PaaBb+2*PAABb*PaaBB+ (PAaBb)^2+2*PAABb*PaaBb+2*PAABB*Paabb+2*PAAbb*PaaBB+2*PAabb*PAaBB+
    2*PAabb*PAaBb+2*PAABb*Paabb+2*PAAbb*PaaBb+2*PaaBB*PAaBb+2*PaaBb*PAaBB+2*PaaBb*PAaBb+2*Paabb*PAaBB+2*PaaBB*PAabb+
    2*Paabb*PAaBb+2*PaaBb*PAabb
  
  
  zy6 <- 2*PAabb*PAAbb + (PAabb)^2+2*PAAbb*Paabb + 2*Paabb*PAabb
  
  zy7 <- (PaaBB)^2
  zy8 <- 2*PaaBB*PaaBb + (PaaBb)^2+2*PaaBB*Paabb + 2*PaaBb*Paabb
  zy9 <- (Paabb)^2
  
  tmp <- c(zy1,zy2,zy3,zy4,zy5,zy6,zy7,zy8,zy9)
  #sum(tmp)
  return(tmp)
}




sim_geno_autoP  <- function(allpar,n=1000){
  
  
  conp <- autoP_DF10(pall=allpar)
  
  id <- sample(1:9,n,replace = T,prob=conp)
  
  return(id)
}




geno_est_autoP <- function(geno){
  
  nts <- rep(0,9)
  nt <- table(geno)
  nts[as.numeric(names(nt))] <- as.numeric(nt)
  nn <- sum(nts)
  
  res4 <- EM_CP(nts=nts)
  
  #DA;DB
  pAA <- sum(res4[1:3])
  pAa <- sum(res4[4:6])
  paa <- sum(res4[7:9])
  
  pBB <- sum(res4[c(1,4,7)])
  pBb <- sum(res4[c(2,5,8)])
  pbb <- sum(res4[c(3,6,9)])
  
  pA <- pAA+pAa/2;pa <- paa+pAa/2;
  pB <- pBB+pBb/2;pb <- pbb+pBb/2;
  
  DA <- pAA-pA^2;DB <- pBB-pB^2;
  
  #Deab
  pAB <- (2*res4[1] + sum(res4[c(2,4)])+res4[5]/2)/2
  pAb <- (2*res4[3] + sum(res4[c(2,6)])+res4[5]/2)/2
  paB <- (2*res4[7] + sum(res4[c(4,8)])+res4[5]/2)/2
  pab <- (2*res4[9] + sum(res4[c(6,8)])+res4[5]/2)/2
  Deab <- 2*(pAB-pA*pB)
  
  #DAb
  DAb <- res4[1]+res4[2]/2-pA*Deab-pB*DA-pA^2*pB
  #-(res4[3]+res4[2]/2)-pA*Deab+pb*DA+pA^2*pb
  #-(res4[4]+res4[5]/2)/2-(pA-pa)/2*Deab-pB*DA+pA*pa*pB
  #(res4[6]+res4[5]/2)/2-(pA-pa)/2*Deab+pb*DA-pA*pa*pb
  #res4[7]+res4[8]/2+pa*Deab-pB*DA-pa^2*pB
  #-(res4[9]+res4[8]/2)+pa*Deab+pb*DA+pa^2*pb
  #DaB
  DaB <- res4[1]+res4[4]/2-pB*Deab-pA*DB-pA*pB^2
  #-(res4[7]+res4[4]/2)-pB*Deab+pa*DB+pa*pB^2
  #-(res4[2]+res4[5]/2)/2-(pB-pb)/2*Deab-pA*DB+pA*pB*pb
  #(res4[8]+res4[5]/2)/2-(pB-pb)/2*Deab+pa*DB-pa*pB*pb
  #(res4[3]+res4[6]/2)+pb*Deab-pA*DB-pA*pb^2
  #-(res4[9]+res4[6]/2)+pb*Deab+pa*DB+pa*pb^2
  
  #DAB
  DAB <- res4[1]-(pA^2*pB^2+pB^2*DA+pA^2*DB+2*pB*DAb+2*pA*DaB+2*pA*pB*Deab)-DA*DB-Deab^2
  #pA^2*pB*pb-pA^2*DB+pB*pb*DA+(pA*pb-pA*pB)*Deab+(pb-pB)*DAb-2*pA*DaB-res4[2]/2-DA*DB-Deab^2
  estp <- as.numeric(c(pA,pB,DA,DB,Deab,DAb,DaB,DAB))
  return(estp)
}







EM_DP <- function(ntsAA){
  
  p4 <- rep(1/3,3)
  nn <- sum(ntsAA)
  iter <- 0
  
  while(1){
    
    p41 <- p4
    PAM <- c(2*p4[1]*p4[2],2*p4[1]*p4[3],p4[2]^2,2*p4[2]*p4[3])
    phi <- PAM/sum(PAM)
    p1 <- (2*ntsAA[1]+ntsAA[2]*(phi[1]+phi[2]))/(2*nn)
    p2 <- (ntsAA[2]*(1-phi[2]+phi[3]))/(2*nn)
    p3 <- ((1-phi[1]-phi[3])*ntsAA[2]+2*ntsAA[3])/(2*nn)
    
    
    p4 <- c(p1,p2,p3)
    iter <- iter + 1
    if(max(abs(p4-p41))<1e-5)
      break
  }
  return(p4)
  
}


EM_abP <- function(nts){
  
  
  
  p4 <- rep(0.25,4)
  nn <- sum(nts)
  iter <- 0
  
  while(1){
    
    p41 <- p4
    emm <- EM_EP(p=p4)
    p1 <- sum(emm[1,]*nts)/(4*nn)
    p2 <- sum(emm[2,]*nts)/(4*nn)
    p3 <- sum(emm[3,]*nts)/(4*nn)
    p4 <- sum(emm[4,]*nts)/(4*nn)
    
    p4 <- c(p1,p2,p3,p4)
    iter <- iter + 1
    if(max(abs(p4-p41))<1e-5)
      break
  }
  return(p4)
}

EM_EP <- function(p){
  
  PAB <- p[1]
  PAb <- p[2]
  PaB <- p[3]
  Pab <- p[4]
  
  T1 <- c(4*PAB^3*PAb,6*PAB^2*PAb^2,4*PAB*PAb^3)
  phi_11 <- (4*PAB^3*PAb)/sum(T1)
  phi_12 <- (6*PAB^2*PAb^2)/sum(T1)
  
  T2 <- c(4*PAB^3*PaB,6*PAB^2*PaB^2,4*PAB*PaB^3)
  phi_21 <- (4*PAB^3*PaB)/sum(T2)
  phi_22 <- (6*PAB^2*PaB^2)/sum(T2)
  
  T3 <- c(4*PAB^3*Pab,12*PAB^2*PAb*PaB,12*PAB^ 2*PAb*Pab,12*PAB*PAb^2*PaB,12*PAB*PAb^2*Pab,4*PAb^3*PaB,12*PAB^2*PaB*Pab,
          12*PAB*PAb*PaB^2,6*PAB^2*Pab^2,24*PAB*PAb*PaB*Pab,6*PAb^2*PaB^2,12*PAB*PAb*Pab^2,12*PAb^2*PaB*Pab,
          12*PAB*PaB^2*Pab,4*PAb*PaB^3,12*PAB*PaB*Pab^2,12*PAb*PaB^2*Pab,4*PAB*Pab^3,12*PAb*PaB*Pab^2)
  
  phi3 <- T3/sum(T3)
  i11 <- c(3,2,2,1,1,2,1,2,1,1,1,1,1)
  i12 <- c(1,2,3,4,5,7,8,9,10,12,14,16,18)
  phi3_AB <- sum(i11*phi3[i12])
  i21 <- c(1,1,2,2,3,1,1,2,1,2,1,1,1)
  i22 <- c(2,3,4,5,6,8,10,11,12,13,15,17,19)
  phi3_Ab <- sum(i21*phi3[i22])
  
  i31 <- c(1,1,1,1,2,1,2,1,2,3,1,2,1)
  i32 <- c(2,4,6,7,8,10,11,13,14,15,16,17,19)
  phi3_aB <- sum(i31*phi3[i32])
  
  i41 <- c(1,1,1,1,2,1,2,1,1,2,1,3,2)
  i42 <- c(1,3,5,7,9,10,12,13,14,16,17,18,19)
  phi3_ab <- sum(i41*phi3[i42])
  
  T4 <- c(4*PAb^3*Pab,6*PAb^2*Pab^2,4*PAb*Pab^3)
  phi_41 <- (4*PAb^3*Pab)/sum(T4)
  phi_42 <- (6*PAb^2*Pab^2)/sum(T4)
  
  T5 <- c(4*PaB^3*Pab,6*PaB^2*Pab^2,4*PaB*Pab^3)
  phi_51 <- (4*PaB^3*Pab)/sum(T5)
  phi_52 <- (6*PaB^2*Pab^2)/sum(T5)
  
  epAB <- c(4,1+2*phi_11+phi_12,0,1+2*phi_21+phi_22,phi3_AB,0,0,0,0)
  epAb <- c(0,3-2*phi_11-phi_12,4,0,phi3_Ab,1+2*phi_41+phi_42,0,0,0)
  epaB <- c(0,0,0,3-2*phi_21-phi_22,phi3_aB,0,4,1+2*phi_51+phi_52,0)
  epab <- c(0,0,0,0,phi3_ab,3-2*phi_41-phi_42,0,3-2*phi_51-phi_52,4)
  
  ep <- rbind(epAB,epAb,epaB,epab)
  ep[which(is.nan(ep))] <- 0
  return(ep)
}

EM_CP <- function(nts){
  
  
  
  p4 <- rep(1/9,9)
  nn <- sum(nts)
  iter <- 0
  nts11 <- matrix(rep(nts,9),nrow=9,byrow=T)
  while(1){
    
    p41 <- p4
    emm <- CEMP(p=p4)
    p4 <- rowSums(emm*nts11)/(2*nn)
    
    iter <- iter + 1
    if(max(abs(p4-p41))<1e-5)
      break
  }
  return(p4)
}



CEMP <- function(p){
  
  PAABB <- p[1]
  PAABb <- p[2]
  PAAbb <- p[3]
  PAaBB <- p[4]
  PAaBb <- p[5]
  PAabb <- p[6]
  PaaBB <- p[7]
  PaaBb <- p[8]
  Paabb <- p[9]
  
  
  
  
  T1 <- c(2*PAABB*PAABb,2*PAABB*PAAbb,(PAABb)^2,2*PAABb*PAAbb)
  phi_11 <- (2*PAABB*PAABb)/sum(T1)
  phi_12 <- (2*PAABB*PAAbb)/sum(T1)
  phi_13 <- ((PAABb)^2)/sum(T1)
  
  T2 <- c(2*PAABB*PAaBB , (PAaBB)^2,2*PAABB*PaaBB , 2*PAaBB*PaaBB)
  phi_21 <- (2*PAABB*PAaBB)/sum(T2)
  phi_22 <- ((PAaBB)^2)/sum(T2)
  phi_23 <- (2*PAABB*PaaBB)/sum(T2)
  
  T3 <- c(2*PAABB*PAaBb,2*PAaBB*PAABb , 2*PAABb*PAaBb,2*PAaBB*PAAbb,2*PAabb*PAABB , 2*PAAbb*PAaBb,2*PAabb*PAABb,
          2*PAaBB*PAaBb,2*PAABB*PaaBb,2*PAABb*PaaBB, (PAaBb)^2,2*PAABb*PaaBb,2*PAABB*Paabb,2*PAAbb*PaaBB,2*PAabb*PAaBB,
          2*PAabb*PAaBb,2*PAABb*Paabb,2*PAAbb*PaaBb,2*PaaBB*PAaBb,2*PaaBb*PAaBB,2*PaaBb*PAaBb,2*Paabb*PAaBB,2*PaaBB*PAabb,
          2*Paabb*PAaBb,2*PaaBb*PAabb)
  
  phi3 <- T3/sum(T3)
  
  phi3_AABB <- sum(phi3[c(1,5,9,13)])
  phi3_AABb <- sum(phi3[c(2,3,7,10,12,17)])
  phi3_AAbb <- sum(phi3[c(4,6,14,18)])
  phi3_AaBB <- sum(phi3[c(2,4,8,15,20,22)])
  phi3_AaBb <- sum(phi3[c(1,3,6,8,16,19,21,24)])+2*phi3[11]
  phi3_Aabb <- sum(phi3[c(5,7,15,16,23,25)])
  phi3_aaBB <- sum(phi3[c(10,14,19,23)])
  phi3_aaBb <- sum(phi3[c(9,12,18,20,21,25)])
  phi3_aabb <- sum(phi3[c(13,17,22,24)])
  
  
  
  T4 <- c(2*PAabb*PAAbb,(PAabb)^2,2*PAAbb*Paabb,2*Paabb*PAabb)
  phi_41 <- (2*PAabb*PAAbb)/sum(T4)
  phi_42 <- ((PAabb)^2)/sum(T4)
  phi_43 <- (2*PAAbb*Paabb)/sum(T4)
  
  T5 <- c(2*PaaBB*PaaBb + (PaaBb)^2+2*PaaBB*Paabb + 2*PaaBb*Paabb)
  phi_51 <- (2*PaaBB*PaaBb)/sum(T5)
  phi_52 <- ((PaaBb)^2)/sum(T5)
  phi_53 <- (2*PaaBB*Paabb)/sum(T5)
  
  epAABB <- c(2,phi_11+phi_12,0,phi_21+phi_23,phi3_AABB,0,0,0,0)
  epAABb <- c(0,1-phi_12+phi_13,0,0,phi3_AABb,0,0,0,0)
  epAAbb <- c(0,1-phi_11-phi_13,2,0,phi3_AAbb,phi_41+phi_43,0,0,0)
  epAaBB <- c(0,0,0,1+phi_22-phi_23,phi3_AaBB,0,0,0,0)
  epAaBb <- c(0,0,0,0,phi3_AaBb,0,0,0,0)
  epAabb <- c(0,0,0,0,phi3_Aabb,1+phi_42-phi_43,0,0,0)
  epaaBB <- c(0,0,0,1-phi_21-phi_22,phi3_aaBB,0,2,phi_51+phi_53,0)
  epaaBb <- c(0,0,0,0,phi3_aaBb,0,0,1+phi_52-phi_53,0)
  epaabb <- c(0,0,0,0,phi3_aabb,1-phi_41-phi_42,0,1-phi_51-phi_52,2)
  
  ep <- rbind(epAABB,epAABb,epAAbb,epAaBB,epAaBb,epAabb,epaaBB,epaaBb,epaabb)
  
  #ep[which(is.nan(ep))] <- 0
  return(ep)
}





work_test1 <- function(M,mn){
  
  
  n <- dim(M)[2]
  nc <- combn(n,2)
  nci2 <- dim(nc)[2]
  LD <- matrix(NA,nrow=nci2,8)
  LDn <- rep(NA,nci2)
  for(i in 1:nci2){
    
    index <- paste(M[,nc[1,i]],M[,nc[2,i]],sep="")
    index[which(index=="22")] <- 1
    index[which(index=="21")] <- 2
    index[which(index=="20")] <- 3
    index[which(index=="12")] <- 4
    index[which(index=="11")] <- 5
    index[which(index=="10")] <- 6
    index[which(index=="02")] <- 7
    index[which(index=="01")] <- 8
    index[which(index=="00")] <- 9
    
    LDn[i] <- paste(mn[nc[1,i]],mn[nc[2,i]],sep="-")
    tmp <- try(geno_est_autoP(geno=as.numeric(index)),TRUE)
    if(class(tmp)=="try-error"){
      LD[i,] <- rep(NA,8)
    }else{
      LD[i,] <- tmp
    }
    
    LD[i,] <- tmp
    if(i%%1000==0)
      cat("comb=",i,"\n")
  }
  
  colnames(LD) <- c("pA","pB","DA","DB","detal_ab","D_Ab","D_aB","D_AB")
  rownames(LD) <- LDn
  return(LD)
}


get_con_param<-function(parm.id)
{
  for (e in commandArgs())
  {
    ta = strsplit(e,"=", fixed=TRUE);
    if(! is.na( ta[[1]][2]))
    {
      temp = ta[[1]][2];
      if( ta[[1]][1] == parm.id) {
        return (as.character(temp));
      }
    }
  }
  
  return(NA);
}




work_test1_1 <- function(M,nc,mn,interval=c(1,2)){
  
  
  n <- dim(M)[2]
  nci2 <- interval[2]-interval[1]+1
  LD <- matrix(NA,nrow=nci2,8)
  LDn <- rep(NA,nci2)
  k <- 1
  for(i in interval[1]:interval[2]){
    
    index <- paste(M[,nc[1,i]],M[,nc[2,i]],sep="")
    index[which(index=="22")] <- 1
    index[which(index=="21")] <- 2
    index[which(index=="20")] <- 3
    index[which(index=="12")] <- 4
    index[which(index=="11")] <- 5
    index[which(index=="10")] <- 6
    index[which(index=="02")] <- 7
    index[which(index=="01")] <- 8
    index[which(index=="00")] <- 9
    
    LDn[k] <- paste(mn[nc[1,i]],mn[nc[2,i]],sep="-")
    tmp <- try(geno_est_autoP(geno=as.numeric(index)),TRUE)
    if(class(tmp)=="try-error"){
      LD[i,] <- rep(NA,8)
    }else{
      LD[i,] <- tmp
    }
    
    LD[k,] <- tmp
    if(k%%10000==0)
      cat("comb=",k,"\n")
    
    k <- k + 1
  }
  
  colnames(LD) <- c("pA","pB","DA","DB","detal_ab","D_Ab","D_aB","DAB")
  rownames(LD) <- LDn
  return(LD)
}
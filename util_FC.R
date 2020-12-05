


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

auto_DF10 <- function(pall){
  
  
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
  zy2 <- 2*PAABB*PAABb
  zy3 <- 2*PAABB*PAAbb+(PAABb)^2
  zy4 <- 2*PAABb*PAAbb
  zy5 <- (PAAbb)^2
  zy6 <- 2*PAABB*PAaBB
  zy7 <- 2*PAABB*PAaBb+2*PAaBB*PAABb
  zy8 <- 2*PAABb*PAaBb+2*PAaBB*PAAbb+2*PAabb*PAABB
  zy9 <- 2*PAAbb*PAaBb+2*PAabb*PAABb
  zy10 <- 2*PAabb*PAAbb
  zy11 <- (PAaBB)^2+2*PAABB*PaaBB
  zy12 <- 2*PAaBB*PAaBb+2*PAABB*PaaBb+2*PAABb*PaaBB
  zy13 <- (PAaBb)^2+2*PAABb*PaaBb+2*PAABB*Paabb+2*PAAbb*PaaBB+2*PAabb*PAaBB
  zy14 <- 2*PAabb*PAaBb+2*PAABb*Paabb+2*PAAbb*PaaBb
  zy15 <- (PAabb)^2+2*PAAbb*Paabb
  zy16 <- 2*PAaBB*PaaBB
  zy17 <- 2*PaaBB*PAaBb+2*PaaBb*PAaBB
  zy18 <- 2*PaaBb*PAaBb+2*Paabb*PAaBB+2*PaaBB*PAabb
  zy19 <- 2*Paabb*PAaBb+2*PaaBb*PAabb
  zy20 <- 2*Paabb*PAabb
  zy21 <- (PaaBB)^2
  zy22 <- 2*PaaBB*PaaBb
  zy23 <- (PaaBb)^2+2*PaaBB*Paabb
  zy24 <- 2*PaaBb*Paabb
  zy25 <- (Paabb)^2
  
  tmp <- c(zy1,zy2,zy3,zy4,zy5,zy6,zy7,zy8,zy9,zy10,
           zy11,zy12,zy13,zy14,zy15,zy16,zy17,zy18,zy19,
           zy20,zy21,zy22,zy23,zy24,zy25)
  #sum(tmp)
  return(tmp)
}



sim_geno_auto  <- function(allpar,n=1000){
  
  
  conp <- auto_DF10(pall=allpar)
  
  id <- sample(1:25,n,replace = T,prob=conp)
  
  return(id)
}




geno_est_auto <- function(geno){
  
  nts <- rep(0,25)
  nt <- table(geno)
  nts[as.numeric(names(nt))] <- as.numeric(nt)
  nn <- sum(nts)
  
  res4 <- EM_C(nts=nts)
  
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

  #DaB
  DaB <- res4[1]+res4[4]/2-pB*Deab-pA*DB-pA*pB^2

  #DAB
  DAB <- res4[1]-(pA^2*pB^2+pB^2*DA+pA^2*DB+2*pB*DAb+2*pA*DaB+2*pA*pB*Deab)-DA*DB-Deab^2
  estp <- as.numeric(c(pA,pB,DA,DB,Deab,DAb,DaB,DAB))
  return(estp)
}



EM_D <- function(ntsAA){
  
  p4 <- rep(1/3,3)
  nn <- sum(ntsAA)
  iter <- 0
  
  while(1){
    
    p41 <- p4
    phi <- 2*(p4[1]*p4[3])/(2*p4[1]*p4[3]+p4[2]^2)
    p1 <- (2*ntsAA[1]+ntsAA[2]+phi*ntsAA[3])/(2*nn)
    p2 <- (ntsAA[2]+ntsAA[4]+2*(1-phi)*ntsAA[3])/(2*nn)
    p3 <- (2*ntsAA[5]+ntsAA[4]+phi*ntsAA[3])/(2*nn)
    
    
    p4 <- c(p1,p2,p3)
    iter <- iter + 1
    if(max(abs(p4-p41))<1e-5)
      break
  }
  return(p4)
  
}



EM_ab <- function(nts){
  
  
  
  p4 <- rep(0.25,4)
  nn <- sum(nts)
  iter <- 0
  
  while(1){
    
    p41 <- p4
    phi <- EM_E(p=p4)
    emm <- EM_M(phi=phi)
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

EM_M <- function(phi){
  
  phi1 <- phi[1]
  phi2 <- phi[2]
  phi3 <- phi[3]
  phi4 <- phi[4]
  phi51 <- phi[5]
  phi52 <- phi[6]
  phi6 <- phi[7]
  phi7 <- phi[8]
  phi8 <- phi[9]
  phi9 <- phi[10]
  
  epAB <- c(4,3,2,1,0,3,phi1+2,1+phi2,phi3,0,2,1+phi4,phi51-phi52+1,phi6,0,1,phi7,phi8,phi9,0,0,0,0,0,0)
  
  epAb <- c(0,1,2,3,4,0,1-phi1,2-phi2,3-phi3,3,0,1-phi4,phi52-phi51+1,2-phi6,2,0,1-phi7,1-phi8,1-phi9,1,0,0,0,0,0)
  
  epaB<- c(0,0,0,0,0,1,1-phi1,1-phi2,1-phi3,0,2,2-phi4,phi52-phi51+1,1-phi6,0,3,3-phi7,2-phi8,1-phi9,0,4,3,2,1,0)
  
  epab<- c(0,0,0,0,0,0,phi1,phi2,phi3,1,0,phi4,phi51-phi52+1,1+phi6,2,0,phi7,1+phi8,2+phi9,3,0,1,2,3,4)
  
  ep <- rbind(epAB,epAb,epaB,epab)
  return(ep)
}

EM_E <- function(p){
  
  PAB <- p[1]
  PAb <- p[2]
  PaB <- p[3]
  Pab <- p[4]
  
  
  phi_1 <- (4*PAB^3*Pab)/(4*PAB^3*Pab+12*PAB^2*PAb*PaB)
  phi_2 <- (12*PAB^2*PAb*Pab)/(12*PAB^2*PAb*Pab+12*PAB*PAb^2*PaB)
  phi_3 <- (12*PAB*PAb^2*Pab)/(12*PAB*PAb^2*Pab+4*PAb^3*PaB)
  phi_4 <- (12*PAB^2*PaB*Pab)/(12*(PAB^2*PaB*Pab+PAB*PAb*PaB^2))
  phi_51 <- (6*PAB^2*Pab^2)/(6*PAB^2*Pab^2+24*PAB*PAb*PaB*Pab+6*PAb^2*PaB^2)
  phi_52 <- (6*PAb^2*PaB^2)/(6*PAB^2*Pab^2+24*PAB*PAb*PaB*Pab+6*PAb^2*PaB^2)
  phi_6 <- (12*PAB*PAb*Pab^2)/(12*(PAB*PAb*Pab^2+PAb^2*PaB*Pab))
  phi_7 <- (12*PAB*PaB^2*Pab)/(12*PAB*PaB^2*Pab+4*PAb*PaB^3)
  phi_8 <- (12*PAB*PaB*Pab^2)/(12*(PAB*PaB*Pab^2+PAb*PaB^2*Pab))
  phi_9 <- (4*PAB*Pab^3)/(4*PAB*Pab^3+12*PAb*PaB*Pab^2)
  
  phi_all <- c(phi_1,phi_2,phi_3,phi_4,phi_51,phi_52,phi_6,phi_7,phi_8,phi_9)
  return(phi_all)
}



EM_C <- function(nts){
  
  
  
  p4 <- rep(1/9,9)
  nn <- sum(nts)
  iter <- 0
  nts11 <- matrix(rep(nts,9),nrow=9,byrow=T)
  while(1){
    
    p41 <- p4
    phi <- CEM_E(p=p4)
    emm <- CEM_M(phi=phi)
    p4 <- rowSums(emm*nts11)/(2*nn)
    
    iter <- iter + 1
    if(max(abs(p4-p41))<1e-5)
      break
  }
  return(p4)
}




CEM_M <- function(phi){
  
  phi1 <- phi[1]
  phi2 <- phi[2]
  phi31 <- phi[3]
  phi32 <- phi[4]
  phi4 <- phi[5]
  phi51 <- phi[6]
  phi52 <- phi[7]
  phi61 <- phi[8]
  phi62 <- phi[9]
  phi63 <- phi[10]
  phi64 <- phi[11]
  phi7 <- phi[12]
  phi81 <- phi[13]
  phi82 <- phi[14]
  phi9 <- phi[15]
  phi10 <- phi[16]
  phi111 <- phi[17]
  phi112 <- phi[18]
  phi12 <- phi[19]
  phi13 <- phi[20]
  
  epAABB <- c(2,1,phi1,0,0,1,phi2,phi31,0,0,phi4,phi51,phi61,0,0,0,0,0,0,0,0,0,0,0,0)
  
  epAABb <- c(0,1,2*(1-phi1),1,0,0,1-phi2,phi32,phi7,0,0,phi52,phi62,phi81,0,0,0,0,0,0,0,0,0,0,0)
  
  epAAbb <- c(0,0,phi1,1,2,0,0,1-phi31-phi32,1-phi7,1,0,0,phi63,phi82,phi9,0,0,0,0,0,0,0,0,0,0)
  
  epAaBB <- c(0,0,0,0,0,1,1-phi2,1-phi31-phi32,0,0,2*(1-phi4),1-phi51-phi52,phi64,0,0,1,phi10,phi111,0,0,0,0,0,0,0)
  
  epAaBb <- c(0,0,0,0,0,0,phi2,phi32,1-phi7,0,0,1-phi51-phi52,2*(1-phi61-phi62-phi63-phi64),1-phi81-phi82,0,0,1-phi10,phi112,phi12,0,0,0,0,0,0)
  
  epAabb <- c(0,0,0,0,0,0,0,phi31,phi7,1,0,0,phi64,1-phi81-phi82,2*(1-phi9),0,0,1-phi111-phi112,1-phi12,1,0,0,0,0,0)
  
  epaaBB <- c(0,0,0,0,0,0,0,0,0,0,phi4,phi52,phi63,0,0,1,1-phi10,1-phi111-phi112,0,0,2,1,phi13,0,0)
  
  epaaBb <- c(0,0,0,0,0,0,0,0,0,0,0,phi51,phi62,phi82,0,0,phi10,phi112,1-phi12,0,0,1,2*(1-phi13),1,0)
  
  epaabb <- c(0,0,0,0,0,0,0,0,0,0,0,0,phi61,phi81,phi9,0,0,phi111,phi12,1,0,0,phi13,1,2)
  
  ep <- rbind(epAABB,epAABb,epAAbb,epAaBB,epAaBb,epAabb,epaaBB,epaaBb,epaabb)
  return(ep)
}

CEM_E <- function(p){
  
  PAABB <- p[1]
  PAABb <- p[2]
  PAAbb <- p[3]
  PAaBB <- p[4]
  PAaBb <- p[5]
  PAabb <- p[6]
  PaaBB <- p[7]
  PaaBb <- p[8]
  Paabb <- p[9]
  
  phi_1 <- (2*PAABB*PAAbb)/(2*PAABB*PAAbb+(PAABb)^2)
  phi_2 <- (2*PAABB*PAaBb)/(2*PAABB*PAaBb+2*PAaBB*PAABb)
  phi_31 <- (2*PAabb*PAABB)/(2*PAABb*PAaBb+2*PAaBB*PAAbb+2*PAabb*PAABB)
  phi_32 <- (2*PAABb*PAaBb)/(2*PAABb*PAaBb+2*PAaBB*PAAbb+2*PAabb*PAABB)
  phi_4 <- (2*PAABB*PaaBB)/((PAaBB)^2+2*PAABB*PaaBB)
  phi_51 <- (2*PAABB*PaaBb)/(2*PAaBB*PAaBb+2*PAABB*PaaBb+2*PAABb*PaaBB)
  phi_52 <- (2*PAABb*PaaBB)/(2*PAaBB*PAaBb+2*PAABB*PaaBb+2*PAABb*PaaBB)
  phi_61 <- (2*PAABB*Paabb)/((PAaBb)^2+2*PAABb*PaaBb+2*PAABB*Paabb+2*PAAbb*PaaBB+2*PAabb*PAaBB)
  phi_62 <- (2*PAABb*PaaBb)/((PAaBb)^2+2*PAABb*PaaBb+2*PAABB*Paabb+2*PAAbb*PaaBB+2*PAabb*PAaBB)
  phi_63 <- (2*PAAbb*PaaBB)/((PAaBb)^2+2*PAABb*PaaBb+2*PAABB*Paabb+2*PAAbb*PaaBB+2*PAabb*PAaBB)
  phi_64 <- (2*PAabb*PAaBB)/((PAaBb)^2+2*PAABb*PaaBb+2*PAABB*Paabb+2*PAAbb*PaaBB+2*PAabb*PAaBB)
  phi_7 <- (2*PAabb*PAABb)/(2*PAAbb*PAaBb+2*PAabb*PAABb)
  phi_81 <- (2*PAABb*Paabb)/(2*PAabb*PAaBb+2*PAABb*Paabb+2*PAAbb*PaaBb)
  phi_82 <- (2*PAAbb*PaaBb)/(2*PAabb*PAaBb+2*PAABb*Paabb+2*PAAbb*PaaBb)
  phi_9 <- (2*PAAbb*Paabb)/((PAabb)^2+2*PAAbb*Paabb)
  phi_10 <- (2*PaaBb*PAaBB)/(2*PaaBB*PAaBb+2*PaaBb*PAaBB)
  phi_111 <- (2*Paabb*PAaBB)/(2*PaaBb*PAaBb+2*Paabb*PAaBB+2*PaaBB*PAabb)
  phi_112 <- (2*PaaBb*PAaBb)/(2*PaaBb*PAaBb+2*Paabb*PAaBB+2*PaaBB*PAabb)
  phi_12 <- (2*Paabb*PAaBb)/(2*Paabb*PAaBb+2*PaaBb*PAabb)
  phi_13 <- (2*PaaBB*Paabb)/((PaaBb)^2+2*PaaBB*Paabb)
  
  phi_all <- c(phi_1,phi_2,phi_31,phi_32,phi_4,phi_51,phi_52,phi_61,phi_62,phi_63,phi_64,phi_7,phi_81,phi_82,phi_9,
               phi_10,phi_111,phi_112,phi_12,phi_13)
  return(phi_all)
}











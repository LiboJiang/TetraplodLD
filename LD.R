

Lcs <- function(l1,l2,l3){
  
  deln <- c(rownames(l1),rownames(l2))
  
  delni <- match(deln,rownames(l3))
  
  fd <- l3[-delni,]
  
  fd
}

LD_chr <- function(dn,mi){
  
  pos <- mi$pos
  nc <- combn(length(pos),2)
  nci2 <- dim(nc)[2]
  dist <-rep(NA,nci2)
  for(i in 1:nci2){
    
    dist[i] <- pos[nc[2,i]]-pos[nc[1,i]]
  }
  
  ii <- order(dist)
  rsq <- cbind(dist[ii]/10^6,dn[ii,])
  ii <- which(is.na(rowSums(rsq)))
  rsq[-ii,]
}



LD_chr_inv <- function(dn,mi){
  
  
  chr1 <- as.character(unique(mi$chr))
  scafold_len <- c();
  for(sc in chr1){
    ids <- which(mi$chr == sc);
    scafold_len <- c(scafold_len, max((mi$pos[ids])));
  }
  
  leiji_len <- cumsum( c(0,scafold_len) );
  
  position <- rep(0, length(mi$chr));
  for(scafold in 1:length(chr1)){
    ids <- which(mi$chr == chr1[scafold]);
    position[ids] <- as.numeric(mi$pos[ids])/(10^6) + leiji_len[scafold]/(10^6)
  }
  
  pos <- position
  nc <- combn(length(pos),2)
  nci2 <- dim(nc)[2]
  dist <-rep(NA,nci2)
  for(i in 1:nci2){
    
    dist[i] <- pos[nc[2,i]]-pos[nc[1,i]]
  }
  
  ii <- order(dist)
  rsq <- cbind(dist[ii],dn[ii,])
  ii <- which(is.na(rowSums(rsq)))
  rsq[-ii,]
}

LD_points_chr <- function(res,xll=seq(0,50,10),yll=seq(0,1,0.2),
                          h=5,w=7,filen="LD_point_chr03a.pdf"){
  
  res <- as.data.frame(res)
  xl <- c(-2,(max(xll)[1])*1.01);yl <- c(-0.01,max(yll)*1.05)
  
  pdf(filen,height=h,width=w)
  par(oma=c(4,4.6,0.5,1),mar=c(0,0,0,0))
  plot(NA,NA,pch=16,type="n",xlab=" ",ylab=" ",xlim=xl,ylim=yl,xaxt="n",yaxt="n",xaxs="i", yaxs="i")
  points(res[,1],res[,2],col="#436EEE50",pch=20,cex=0.2)
  axis(2,at=yll,labels=yll,cex.axis=1.6,lwd=1)
  axis(1,at=xll,labels=xll,cex.axis=1.6,lwd=1)
  mtext("Distance (Mb)",1,cex=1.8,line=2.7)
  mtext(expression(r^2),2,cex=1.8,line=2.7)
  dev.off()
}



LD_pointsN_chr <- function(res,xll=seq(0,50,10),yll=seq(0,1,0.2),
                          h=5,w=7,filen="LD_point_chr03a.pdf"){
  
  res <- as.data.frame(res)
  xl <- c(-2,(max(xll)[1])*1.01);yl <- c(-0.01,max(yll)*1.05)
  
  pdf(filen,height=h,width=w)
  par(oma=c(4,4.6,0.5,1),mar=c(0,0,0,0))
  plot(NA,NA,pch=16,type="n",xlab=" ",ylab=" ",xlim=xl,ylim=yl,xaxt="n",yaxt="n",xaxs="i", yaxs="i")
  points(res[,1],res[,2],col="#436EEE50",pch=20,cex=0.2)
  axis(2,at=yll,labels=yll,cex.axis=1.6,lwd=1)
  axis(1,at=xll,labels=xll,cex.axis=1.6,lwd=1)
  mtext("Distance (Mb)",1,cex=1.8,line=2.7)
  mtext(expression(D^"'"),2,cex=1.8,line=2.7)
  dev.off()
}

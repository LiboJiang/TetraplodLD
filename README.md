# A theoretical framework to model linkage disequilibria for allotetraploid populations


# 1 Introduction 

At present, the TetraplodLD package is used to linkage disequilibrium analysis in the nature population of allotetraploid. This guide gives some brief instructions on how to perform the tasks of linkage disequilibrium analysis by this package. The outline of this guide is as follows:


# 2 Data format

ID         snp1  snp2  snp3  snp4  snp5  snp6

IND1        2     0     1     0     1      1

IND2        0     2     2     1     1      1

IND3        1     1     2     1     2      1

IND4        0     2     1     1     1      1

IND5        0     2     1     1     1      1

IND6        0     2     1     1     1      1

IND7        1     0     0     1     0      1

IND8        2     2     2     2     2      2

IND9        1     2     2     2     2      1

This table indicates dosage-unknown markers, each marker contain three genotypes (aaaa=0, A___=1, AAAA=2). 


# 3 Work example

#Dosage-unknown marker

#read genotype file

chr03a_m <- read.csv("./data/chr03am.csv")[,-1]

#extract snp ID

snpn <- colnames(chr03a_m)

#load functions

source("util_PC_m.R")

source("util_PC.R")

#Linkage disequilibrium analysis

LD <- work_test1_mix(M=chr03a_m,mn=snpn)

#M is a matrix; mn is a character vector giving the snp ID.

#work_test1_mix produces a matrix with some or all of the following elements in :

LR        The log-likelihood ratio between haplotype and diplotype models.

Pv         P-value is calculated through chi-square distribution based on LR.

m1_pA     The estimated allele frequency of A by haplotype model.

m1_pB     The estimated allele frequency of B by haplotype model.

m1_D      The estimate of the LD coefficient by haplotype model.

m2_pA     The estimated allele frequency of A by diplotype model.

m1_pB     The estimated allele frequency of B by diplotype model.

m1_DA     The estimate of the LD coefficient at the locus A.

m1_DB     The estimate of the LD coefficient at the locus B.

m1_Deab    The sum of the estimate of the LD coefficient between two nonalleles at different loci on the same haplotype and the estimate of the LD coefficient between two nonalleles on different haplotypes.

m1_DAb     The estimate of the LD coefficient between two alleles from SNP A and one allele from SNP B.

m1_DaB     The estimate of the LD coefficient between two alleles from SNP A and one allele from SNP B.

m1_DAB     between two alleles from SNP A and two alleles from SNP B.

m2_DA_n, m2_DB_n, m2_Deab_n, m2_DAb_n, m2_DaB_n, m2_DAB_n  The estimate of the standardized LD coefficient.



# 4 Computer simulation

#Simulation dosage-unknown markers based on diplotype model

#load functions

source("util_FC.R")

#The true value of LD coefficient in computer simulation 

pall <- c(0.6,0.45,0.12,0.12,0.02,0.02,0.02,0.01)

#1000 simulation replicates

allres <- c()

for(i in 1:1000){

    geno <- sim_geno_auto(allpar=pall,n=500) #n represents the sample size
    
    res <- geno_est_auto(geno=geno)
    
    allres <- rbind(allres,res)
 }
 
#Calculate the mean value of the estimation paraeters

colMeans(allres)

#Calculate the standard deviation of the estimation paraeters

apply(allres,2,sd)

#Simulation dosage-known markers based on diplotype model

#load functions

source("util_PC.R")

#The true value of LD coefficient in computer simulation 

pall <- c(0.6,0.45,0.12,0.12,0.02,0.02,0.02,0.01)

#1000 simulation replicates

allres <- c()

for(i in 1:1000){

    geno <- sim_geno_autoP(allpar=pall,n=500) #n represents the sample size

    res <- geno_est_autoP(geno=geno)

    allres <- rbind(allres,res)

 }
 
#Calculate the mean value of the estimation paraeters

colMeans(allres)

#Calculate the standard deviation of the estimation paraeters

apply(allres,2,sd)


#experimental aligner.
#72 504
#445 897
sam1<-gzreader("../HLA_asya/HLA-MRD-July/L4_S17_L001_R1_001.fastq.gz","../HLA_asya/HLA-MRD-July/L4_S17_L001_R1_001.fastq.gz")
#sam1<-gzreader("../HLA_asya/HLA-MRD-July/1_S12_L001_R1_001.fastq.gz","../HLA_asya/HLA-MRD-July/1_S12_L001_R2_001.fastq.gz")
sam1<-list(Iamp1=sam1[grepl("CCCTGACC[GC]AGACCTG",substr(sam1$read1,1,20)),],
           Iamp2=sam1[grepl("CGACGGCAA[AG]GATTAC",substr(sam1$read1,1,20)),],
           Iamp1alt=sam1[grepl("TC[CT]CACTCCATGAGGTATTTC|TCCCACTCCATGAAGTATTTC",substr(sam1$read1,1,22)),],
           Iamp2alt=sam1[grepl("GGCAA[AG]GATTACATCGCC|GGCAAGGATTACATCGCT",substr(sam1$read1,1,20)),],
           Iamp1_inv=sam1[grepl("GGGCCGCCTCC[AC]ACTTG|GGGCCGTCTCCCACTTG|GGACCGCCTCCCACTTG",substr(sam1$read1,1,20)),],
           Iamp2_inv=sam1[grepl("[CT]GGTGG[AG]CTGGGAAGA",substr(sam1$read1,1,20)),],
           Iamp1alt_inv=sam1[grepl("GAGC[GC]ACTCCACGCAC|GAGCCCGTCCACGCAC",substr(sam1$read1,1,22)),],
           Iamp2alt_inv=sam1[grepl("TCAGGGTGAGGGGCT|TCAGGGTGCAGGGCT",substr(sam1$read1,1,20)),])

r1mat<-do.call(cbind,strsplit(substr(sam1$Iamp1$read1,20,250),""))
r1mat2<-do.call(cbind,strsplit(substr(sam1$Iamp2$read1,20,250),""))

r1invmat<-do.call(cbind,strsplit(substr(sam1$Iamp1$read2,20,250),""))
r1mat_alt<-do.call(cbind,strsplit(substr(sam1$Iamp1alt$read1,20,250),""))
r1mat2_alt<-do.call(cbind,strsplit(substr(sam1$Iamp2alt$read1,20,250),""))

r1invaltmat<-do.call(cbind,strsplit(substr(sam1$Iamp1$read2,20,250),""))

basemat<-do.call(cbind,strsplit(substr(AlignHLAspf$C$Alignment,72,72+230),""))
basemat_amp2<-do.call(cbind,strsplit(substr(AlignHLAspf$C$Alignment,445,445+230),""))

basemat_amp2_alt<-do.call(cbind,strsplit(substr(AlignHLAspf$C$Alignment,449,449+230),""))
basemat_alt<-do.call(cbind,strsplit(substr(AlignHLAspf$C$Alignment,95,95+230),""))

r2mat<-do.call(cbind,strsplit(substr(revcomp(sam1$Iamp1$read2),20,250),""))
basemat2<-do.call(cbind,strsplit(substr(AlignHLAspf$C$Alignment,274,504),""))
score1<-matrix(0,ncol = ncol(r1mat),nrow=ncol(basemat));
score2<-matrix(0,ncol = ncol(r1mat2),nrow=ncol(basemat_amp2));
score1alt<-matrix(0,ncol = ncol(r1mat_alt),nrow=ncol(basemat_alt));
score2alt<-matrix(0,ncol = ncol(r1mat2_alt),nrow=ncol(basemat_amp2_alt));

#score1<-numeric(ncol(basemat))
for (i in 1:ncol(basemat))
{
basemat_ref<-matrix(basemat[,i],nrow=nrow(r1mat),ncol=ncol(r1mat))
basemat_ref2<-matrix(basemat_amp2[,i],nrow=nrow(r1mat2),ncol=ncol(r1mat2))
basemat_ref_alt<-matrix(basemat_alt[,i],nrow=nrow(r1mat_alt),ncol=ncol(r1mat_alt))
basemat_ref2_alt<-matrix(basemat_amp2_alt[,i],nrow=nrow(r1mat2_alt),ncol=ncol(r1mat2_alt))

#logic_mat<-(basemat_ref!=r1mat)
#score1[i,]<-colSums(logic_mat)
score1[i,]<-((colSums((((basemat_ref!=r1mat))))))
score2[i,]<-((colSums((((basemat_ref2!=r1mat2))))))
score1alt[i,]<-((colSums((((basemat_ref_alt!=r1mat_alt))))))
score2alt[i,]<-((colSums((((basemat_ref2_alt!=r1mat2_alt))))))
#print(i)
}

rownames(score1)<-AlignHLAspf$C$Allele
rownames(score2)<-AlignHLAspf$C$Allele
rownames(score1alt)<-AlignHLAspf$C$Allele
rownames(score2alt)<-AlignHLAspf$C$Allele

#svd or something!
# d <- diag(t$d)
# 
# u <- t$u
# v <- t$v
# u1 <- as.matrix(u[-1, 1])
# v1 <- as.matrix(v[-1, 1])
# d1 <- as.matrix(d[1, 1])
# l1 <- u1 %*% d1 %*% t(v1)
# 
# depth <- 2
# us <- as.matrix(u[, 1:depth])
# vs <- as.matrix(v[, 1:depth])
# ds <- as.matrix(d[1:depth, 1:depth])
# ls <- us %*% ds %*% t(vs)

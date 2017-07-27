#experimental aligner.
#72 504
#445 897
sam1<-gzreader("../HLA_asya/HLA-MRD-July/1-1_S1_L001_R1_001.fastq.gz","../HLA_asya/HLA-MRD-July/1-1_S1_L001_R2_001.fastq.gz")
sam1<-list(Iamp1=sam1[grepl("CCCTGACC[GC]AGACCTG",substr(sam1$read1,1,20)),],
           Iamp2=sam1[grepl("CGACGGCAA[AG]GATTAC",substr(sam1$read1,1,20)),],)

r1mat<-do.call(cbind,strsplit(substr(sam1$Iamp1$read1,20,250),""))
basemat<-do.call(cbind,strsplit(substr(AlignHLAspf$C$Alignment,72,72+230),""))
r2mat<-do.call(cbind,strsplit(substr(revcomp(sam1$Iamp1$read2,20,250)),""))
basemat2<-do.call(cbind,strsplit(substr(AlignHLAspf$C$Alignment,274,504),""))
score1<-matrix(0,ncol = ncol(r1mat),nrow=ncol(basemat));
#score1<-numeric(ncol(basemat))
for (i in 1:3116)
{
basemat_ref<-matrix(basemat[,i],nrow=nrow(r1mat),ncol=ncol(r1mat))
#basemat_ref2<-matrix(basemat2[,i],nrow=nrow(r1mat),ncol=ncol(r1mat))

#logic_mat<-(basemat_ref!=r1mat)
#score1[i,]<-colSums(logic_mat)
score1[i,]<-(exp(colSums(log(abs((basemat_ref==r1mat)-0.01)))))*
print(i)
}

#svd or something!
d <- diag(t$d)

u <- t$u
v <- t$v
u1 <- as.matrix(u[-1, 1])
v1 <- as.matrix(v[-1, 1])
d1 <- as.matrix(d[1, 1])
l1 <- u1 %*% d1 %*% t(v1)

depth <- 2
us <- as.matrix(u[, 1:depth])
vs <- as.matrix(v[, 1:depth])
ds <- as.matrix(d[1:depth, 1:depth])
ls <- us %*% ds %*% t(vs)

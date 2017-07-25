#HLA analysis from safety2 files

library(dplyr)
library(data.table)
library(igraph)
library(Biostrings)

#safety3<-get_safety3(amplicones_miseq_new_lists_all$MP$safety2)



get_HLA_result<-function (donor) {  print("START!")
  print(format(Sys.time(), "%a %b %d %X %Y"))
  if (sum(sapply(donor$safety2, nrow)[c(1,2,6,7,12,13,17,18)])>0) {
  safety3<-get_safety3(donor$safety2[(donor$safety2$parents==1&donor$safety2$freq>0.005)|(donor$safety2$parents==2&donor$safety2$freq>0.05), ])
  safety4<-get_safety4(safety3)
  safety5<-get_safety5(DT=safety4[[1]])
  list(safety3=safety3, safety4=safety4, safety5=safety5)
  }
}
get_safety3 <- function(safety2) {
  HLA_allele_list<-HLA_base
  #View(HLA_allele_list)
  HLA_allele_Iclass<-select(HLA_allele_list, Allele, Sequence, HLA_class)%>%filter(HLA_class=="A"|HLA_class=="B"|HLA_class=="C")%>%select(Allele, Sequence)
  HLA_allele_Iclass$Iamp1<-0
  HLA_allele_Iclass$Iamp1_inv<-0
  HLA_allele_Iclass$Iamp2<-0
  HLA_allele_Iclass$Iamp2_inv<-0
  HLA_allele_Iclass$Iamp1alt<-0
  HLA_allele_Iclass$Iamp1alt_inv<-0
  HLA_allele_Iclass$Iamp2alt<-0
  HLA_allele_Iclass$Iamp2alt_inv<-0
  #View(HLA_allele_Iclass)
  HLA_Iclass_amps<-select(HLA_allele_Iclass, Allele)
  HLA_Iclass_amps$Iamp1<-0
  HLA_Iclass_amps$Iamp1_inv<-0
  HLA_Iclass_amps$Iamp2<-0
  HLA_Iclass_amps$Iamp2_inv<-0
  HLA_Iclass_amps$Iamp1alt<-0
  HLA_Iclass_amps$Iamp1alt_inv<-0
  HLA_Iclass_amps$Iamp2alt<-0
  HLA_Iclass_amps$Iamp2alt_inv<-0
  #View(HLA_Iclass_amps)
  
  allele_freq_in_ampsI <- function(what,where, indices) {
    
    what$assembled<-gsub(what$assembled, pattern = "N", replacement = ".")
    
    for (i in 1:nrow(what)) {
      where[grepl(x=HLA_allele_Iclass$Sequence,
                  pattern = what$assembled[i], 
                  fixed=F)]<-where[grepl(x=HLA_allele_Iclass$Sequence, 
                                         pattern = what$assembled[i])]+what$freq[i]
      
      indices[grepl(x=HLA_allele_Iclass$Sequence,
                    pattern = what$assembled[i], 
                    fixed=F)] <- paste0(indices[grepl(x=HLA_allele_Iclass$Sequence,
                                                      pattern = what$assembled[i], 
                                                      fixed=F)],",", i)
    }
    list(where,indices)
  }
  Ns<-rep(0,times=8)
  names(Ns)<-colnames(HLA_allele_Iclass)[3:10]
  for (colname in colnames(HLA_allele_Iclass)[3:10]) {
    if (nrow(safety2[[colname]])>0) {
    tempres<-allele_freq_in_ampsI(safety2[[colname]], HLA_allele_Iclass[[colname]], HLA_Iclass_amps[[colname]])
    HLA_allele_Iclass[[colname]]<-tempres[[1]]
    HLA_Iclass_amps[[colname]]<- tempres[[2]]
    Ns[colname]<-sum(safety2[[colname]]$readnumber)
    #Ns[colname]<-0
    }
  }
  
 
  #View(HLA_Iclass_amps)
  safety3<-list(HLA_allele_Iclass=HLA_allele_Iclass, HLA_Iclass_amps=HLA_Iclass_amps, Ns=Ns)
  safety3
}


get_safety4<-function(safety3) {
 safety3$HLA_allele_Iclass$amps<-apply(safety3$HLA_Iclass_amps[2:9], MARGIN = 1, paste0, collapse = "_")
 #safety3$HLA_allele_Iclass[3:10]<-apply(safety3$HLA_allele_Iclass[3:10], MARGIN = 2, prop.table)
 DT<-as.data.table(safety3$HLA_allele_Iclass)
                              DT<-DT[amps!="0_0_0_0_0_0_0_0", .(Allele=paste0(Allele, collapse=" "),
                                                                                Iamp1=unique(Iamp1),
                                                                                Iamp1_inv=unique(Iamp1_inv),
                                                                                Iamp2=unique(Iamp2),
                                                                                Iamp2_inv=unique(Iamp2_inv),
                                                                                Iamp1alt=unique(Iamp1alt),
                                                                                Iamp1alt_inv=unique(Iamp1alt_inv),
                                                                                Iamp2alt=unique(Iamp2alt),
                                                                                Iamp2alt_inv=unique(Iamp2alt_inv)), by=amps]

DT<-as.data.frame(DT)
Ps<-rep(0.05, times=8)
 for (i in 3:10) {
   DT[,i][DT[,i]>0]<- (1-Ps[i-2])#(1-exp(-safety3$Ns[i-2]*Ps[i-2]))
   DT[,i][DT[,i]==0]<- Ps[i-2]#exp(-safety3$Ns[i-2]*Ps[i-2])
 }


DT$Score<-apply(DT[3:10], 1, prod)

# DT[grepl(DT$Allele, pattern = "A*", fixed = T), ]$Score<-prop.table(DT[grepl(DT$Allele, pattern = "A*", fixed = T), ]$Score)
# DT[grepl(DT$Allele, pattern = "B*",fixed = T), ]$Score<-prop.table(DT[grepl(DT$Allele, pattern = "B*", fixed = T), ]$Score)
# DT[grepl(DT$Allele, pattern = "C*", fixed = T), ]$Score<-prop.table(DT[grepl(DT$Allele, pattern = "C*", fixed = T), ]$Score)
#View(DT)

amps_list<-str_split(DT$amps, pattern = "_")
str_split(amps_list[[1]], pattern = ",")
for (i in 1:length(amps_list)) {
  amps_list[[i]]<-str_split(amps_list[[i]], pattern = "," )
}
isdaddy <- function(list1, list2) {
  a<-length(setdiff(list2[[1]], list1[[1]]))==0#(identical(list1[[1]], list2[[1]])+
  b<-length(setdiff(list2[[2]], list1[[2]]))==0
  c<-length(setdiff(list2[[3]], list1[[3]]))==0
  d<-length(setdiff(list2[[4]], list1[[4]]))==0
  e<-length(setdiff(list2[[5]], list1[[5]]))==0
  f<-length(setdiff(list2[[6]], list1[[6]]))==0
  g<-length(setdiff(list2[[7]], list1[[7]]))==0
  h<-length(setdiff(list2[[8]], list1[[8]]))==0
  daddy<-sum(c(a,b,c,d,e,f,g,h))==8
  daddy
}
res<-matrix(0, ncol=nrow(DT), nrow=nrow(DT))
for (i in 1:length(amps_list))
  for (j in 1:length(amps_list)) {
    res[i,j]<-isdaddy(amps_list[[i]], amps_list[[j]])
  } 
gr<-simplify(graph_from_adjacency_matrix(res))
#DT$degree<-degree(gr, mode="out")
#DT$cluster_id<-clusters(gr)$membership
#View(DT)
#DT$rank<-ave(-DT$degree, DT$cluster_id, FUN = rank)
#View(DT)
DT$degree<-degree(gr, mode="in")
list(DT=DT, gr=gr)
}
get_safety5<-function(DT) {HLA_typing_resultsIclass<-filter(DT, DT$degree==0)%>%select(Allele, Score)
HLA_typing_resultsIclass[grepl(HLA_typing_resultsIclass$Allele, pattern = "A*", fixed = T), ]$Score<-prop.table(HLA_typing_resultsIclass[grepl(HLA_typing_resultsIclass$Allele, pattern = "A*", fixed = T), ]$Score)
HLA_typing_resultsIclass[grepl(HLA_typing_resultsIclass$Allele, pattern = "B*",fixed = T), ]$Score<-prop.table(HLA_typing_resultsIclass[grepl(HLA_typing_resultsIclass$Allele, pattern = "B*", fixed = T), ]$Score)
HLA_typing_resultsIclass[grepl(HLA_typing_resultsIclass$Allele, pattern = "C*", fixed = T), ]$Score<-prop.table(HLA_typing_resultsIclass[grepl(HLA_typing_resultsIclass$Allele, pattern = "C*", fixed = T), ]$Score)
#get_safety5<-function(DT) {HLA_typing_resultsIclass<-filter(DT, DT$rank<=1.0)%>%select(Allele, Score)
#View(HLA_typing_resultsIclass)
HLA_typing_resultsIclass
}

#MP<-get_HLA_result(amplicones_miseq_new_lists_all$MP)

TypingResults2<-lapply(amplicones_miseq_new_lists_all, get_HLA_result)

#KochHla_prob<-get_HLA_result(amplicones_miseq_new_lists_all$Koch)

#ggplot(TypingResults$MP$safety5, aes(x=Allele, y=log(1-Score), fill=Score))+geom_bar(stat="identity")+coord_flip()+scale_fill_gradient()+theme_light()




# for (i in KochHla2$safety4[3:10]) {
#   Ps<-rep(0.01, 8)
#   KochHla2$safety4[[1]][,i][KochHla2$safety4[[1]][,i]==0]<-exp(-KochHla2$safety3[[3]][i]*Ps[i])
#   KochHla2$safety4[[1]][,i][KochHla2$safety4[[1]][,i]>0]<-(1-exp(-KochHla2$safety3[[3]][i]*Ps[i]))
# }

Mega_table_hlares<-do.call(rbind, lapply(TypingResults2, function(x){x$safety5}))
Mega_table_hlares<-rownames_to_column(df = Mega_table_hlares)
Mega_table_hlares$donor<-sapply(str_split(Mega_table_hlares$rowname, fixed(".")), function(x) {x[[1]]})
Mega_table_hlares<-select(Mega_table_hlares, donor, Allele, Score)
Mega_table_hlares<-filter(Mega_table_hlares, Score>())
ggplot(Mega_table_hlares[Mega_table_hlares$donor=="1085", ], aes(x=Allele, y=Score, fill=Score))+geom_bar(stat="identity")+coord_flip()+scale_fill_gradient()+theme_light()
#test<-get_HLA_result(amplicones_miseq_new_lists_all$`1739_3`)

#Download all nessesary packages

library(dplyr)
library(data.table)
library(cowsay)
library(Biostrings)
library(stringr)
library(igraph)
library(tibble)
library(parallel)
library(tidyr)

#Function to get one table for all people with tidy alleles and scores. Takes list of results (output of get_HLA_result)
#as an input
get_mega_table <- function(Typing_results) {
  Mega_table_first<-do.call(rbind, lapply(Typing_results, function(x){x$safety5}))
  Mega_table_first<-rownames_to_column(Mega_table_first)
  Mega_table_first$donor<-sapply(str_split(Mega_table_first$rowname, fixed(".")), function(x) {x[[1]]})
  Mega_table_first<-select(Mega_table_first, donor, Allele, Score, no_amps, sumreads)
  Mega_table_first$tidyAllele<-sapply(Mega_table_first$Allele, gettidyHLA)
  Mega_table_first$tidyAllele<-gsub(Mega_table_first$tidyAllele, pattern = ":NA", replacement = "", fixed = T)
  #Mega_table_100717<-Mega_table_100717[Mega_table_100717$Score>1e-5, ]
  Mega_table_first$n_noamps<-sapply(strsplit(Mega_table_first$no_amps, split = "   "), length)
  Mega_table_first<-select(Mega_table_first, donor, tidyAllele, Allele, Score, n_noamps, sumreads, no_amps)
  Mega_table_first
}

#Function to get list of results. Input - safety 2(output of HLA_amplicones_full)
get_HLA_result<-function (donor) {  say("START!")

    print(format(Sys.time(), "%a %b %d %X %Y"))
  
  if (sum(sapply(donor$safety2, nrow)[c(1,2,6,7,12,13,17,18)])>0)
    #if (donor$safety2!="fail")
    {#paretn was freq>0.001
    safety3<-get_safety3(lapply(donor$safety2, function(x){x[((x$parents==1&x$freq>0.01)|(x$parents==2&x$freq>0.01)|(x$parents==3&x$freq>0.01)|(x$parents==4&x$freq>0.01))&!grepl(pattern = "NULL",x$assembled), ]}))
    safety4<-get_safety4(safety3)
    safety5<-get_safety5(DT1=safety4$DT1, DQB=safety4$DQB, DRB=safety4$DRB, DPB=safety4$DPB, DQA=safety4$DQA)
    list(safety3=safety3, safety4=safety4, safety5=safety5)
  }
}

# get_safety3 is a function to obtain safety 3. Safety 3 is a list with variety of objects.
#For each of HLAIclass, DRB, DQB we heve two tables: 1)table with frequencies of reads that match a particular allele
#2) table with indices of reads that match a particular allele 
#Both tables contain whole list of alleles and 8 columns with amps freqs and indices
#Safety 3 also contains number of reads for each amp in HLAI, DRB, DQB
#Safety 3 contains merged safety 2 with information obtained from both inverse and non inverse systems

get_safety3 <- function(safety2) {
  
#HLAbase filter to mask new allele (is used to understand how new allele will work with our pipeline)
  #HLA_base<-HLA_base[!grepl(x = HLA_base$Allele, pattern="A*02:01", fixed = T), ]

#Modify the HLA base to improve the result for short alleles (150N is added to both sides of sequence in base)    
  HLA_allele_list<-HLA_base
  N150<-paste0(rep(x = "N", times=150), collapse="")
  HLA_allele_list$Sequence<-paste0(N150, HLA_allele_list$Sequence, N150)

#Get 3 different bases for the first class, DRB, DQB separetly  
  HLA_allele_Iclass<-select(HLA_allele_list, Allele, Sequence, HLA_class)%>%filter(HLA_class=="A"|HLA_class=="B"|HLA_class=="C")%>%select(Allele, Sequence)
  HLA_allele_DRB<-select(HLA_allele_list, Allele, Sequence, HLA_class)%>%filter(grepl(HLA_allele_list$HLA_class, pattern="DRB*", fixed=F))%>%select(Allele, Sequence)
  HLA_allele_DQB<-select(HLA_allele_list, Allele, Sequence, HLA_class)%>%filter(grepl(HLA_allele_list$HLA_class, pattern="DQB*", fixed=F))%>%select(Allele, Sequence)
  HLA_allele_DPB<-select(HLA_allele_list, Allele, Sequence, HLA_class)%>%filter(grepl(HLA_allele_list$HLA_class, pattern="DPB*", fixed=F))%>%select(Allele, Sequence)
  HLA_allele_DQA<-select(HLA_allele_list, Allele, Sequence, HLA_class)%>%filter(grepl(HLA_allele_list$HLA_class, pattern="DQA*", fixed=F))%>%select(Allele, Sequence)
  
#Initiate safety 3 tables, should add 8 columns with appropriate names
  HLA_allele_Iclass$Iamp1<-0
  HLA_allele_Iclass$Iamp1_inv<-0
  HLA_allele_Iclass$Iamp2<-0
  HLA_allele_Iclass$Iamp2_inv<-0
  HLA_allele_Iclass$Iamp1alt<-0
  HLA_allele_Iclass$Iamp1alt_inv<-0
  HLA_allele_Iclass$Iamp2alt<-0
  HLA_allele_Iclass$Iamp2alt_inv<-0
   
  HLA_allele_DQB$IIamp1_DQB<-0
  HLA_allele_DQB$IIamp1_DQB_inv<-0
  HLA_allele_DQB$IIamp2<-0
  HLA_allele_DQB$IIamp2_inv<-0
  HLA_allele_DQB$IIamp1_DQBalt<-0
  HLA_allele_DQB$IIamp1_DQBalt_inv<-0
  HLA_allele_DQB$IIamp2_DQBalt<-0
  HLA_allele_DQB$IIamp2_DQBalt_inv<-0
  
  HLA_allele_DRB$IIamp1_Others<-0
  HLA_allele_DRB$IIamp1_Others_inv<-0
  HLA_allele_DRB$IIamp2<-0
  HLA_allele_DRB$IIamp2_inv<-0
  HLA_allele_DRB$IIamp1_DRBalt<-0
  HLA_allele_DRB$IIamp1_DRBalt_inv<-0
  HLA_allele_DRB$IIamp2_DRBalt<-0
  HLA_allele_DRB$IIamp2_DRBalt_inv<-0
  
  HLA_allele_DPB$IIamp1_DPB<-0
  HLA_allele_DPB$IIamp1_DPB_inv<-0
  HLA_allele_DPB$IIamp2_DPB<-0
  HLA_allele_DPB$IIamp2_DPB_inv<-0

  
  HLA_allele_DQA$IIamp_DQA<-0
  HLA_allele_DQA$IIamp_DQA_inv<-0

  
  HLA_Iclass_amps<-select(HLA_allele_Iclass, -2)
  HLA_DQB_amps<-select(HLA_allele_DQB, -2)
  HLA_DRB_amps<-select(HLA_allele_DRB, -2)
  
  HLA_DPB_amps<-select(HLA_allele_DPB, -2)
  HLA_DQA_amps<-select(HLA_allele_DQA, -2)
  
#function to obtain list with indices, frequencies, also puts names of matched alleles in col Exact_new of safety2
  allele_freq_in_amps <- function(what,where,indices, base) {
    for (i in 1:nrow(what)) {
      VC<-vcountPattern(pattern = what$assembled[i], subject = DNAStringSet(base$Sequence),
                        max.mismatch = 0, with.indels = F, fixed=F)
      what$Exact_new[i]<-paste0(base$Allele[VC!=0], collapse = " ")
      where[VC!=0]<-where[VC!=0]+what$freq[i]
      indices[VC!=0] <- paste0(indices[VC!=0],",", i)
    }
    list(where,indices, what)
  }

#Initiate vector with names to count readnumber in amps 
  n_readsI<-rep(0,times=8)
  n_readsDQB<-rep(0,times=8)
  n_readsDRB<- rep(0,times=8)
  
  n_readsDPB<-rep(0,times=4)
  n_readsDQA<-rep(0,times=2)
  
  names(n_readsI)<-colnames(HLA_allele_Iclass)[3:10]
  names(n_readsDQB)<-colnames(HLA_allele_DQB)[3:10]
  names(n_readsDRB)<-colnames(HLA_allele_DRB)[3:10]
  
  names(n_readsDPB)<-colnames(HLA_allele_DPB)[3:6]
  names(n_readsDQA)<-colnames(HLA_allele_DQA)[3:4]

#Function to get properly named tables of frequencies and indices
freqs_indices_tables <- function(HLA_allele, colnums, HLA_amps, n_reads, safety2) {
  for (colname in intersect(colnames(HLA_allele)[colnums], names(safety2))) {
    if (nrow(safety2[[colname]])>0)
      {
      tempres<-allele_freq_in_amps(safety2[[colname]], HLA_allele[[colname]], HLA_amps[[colname]], HLA_allele)
      HLA_allele[[colname]]<-tempres[[1]]
      HLA_amps[[colname]]<- tempres[[2]]
      safety2[[colname]]<- tempres[[3]]
      n_reads[colname]<-sum(safety2[[colname]]$readnumber)
      }
  }
  list(HLA_allele, HLA_amps, n_reads, safety2=safety2)
}
HLAI<-freqs_indices_tables(HLA_allele_Iclass, c(3:10), HLA_Iclass_amps, n_readsI, safety2)
DQB<-freqs_indices_tables(HLA_allele_DQB, c(3:10), HLA_DQB_amps, n_readsDQB, safety2=HLAI[[4]])
DRB<-freqs_indices_tables(HLA_allele_DRB, c(3:10), HLA_DRB_amps, n_readsDRB, safety2=DQB[[4]])

DPB<-freqs_indices_tables(HLA_allele_DPB, c(3:6), HLA_DPB_amps, n_readsDPB, safety2=DRB[[4]])
DQA<-freqs_indices_tables(HLA_allele_DQA, c(3:4), HLA_DQA_amps, n_readsDQA, safety2=DPB[[4]])
safety2<-DQA[[4]]
 
#Get merged inv and non inv lists to search for possible new allele variants  
  safety2_merged<-list()
#   HLA_base_str<-DNAStringSet(HLA_allele_list$Sequence)
#   for (name in  names(safety2)[!grepl(pattern = "inv", x = names(safety2), fixed=T)]) 
#   if(nrow(safety2[[name]])!=0&nrow(safety2[[paste0(name, "_inv", sep="")]])!=0)
#   {
#   name2<-paste0(name, "_inv", sep="")
#   safety2_merged[[name]]<-full_join(x = safety2[[name]], y= safety2[[name2]], by=c("assembled", "Exact_new"))
#   safety2_merged[[name]]<-safety2_merged[[name]]%>%select("assembled",  nreads="readnumber.x", nreads_inv="readnumber.y", "Exact_new")
#   safety2_merged[[name]]$nreads[is.na(safety2_merged[[name]]$nreads)]<-0
#   safety2_merged[[name]]$nreads_inv[is.na(safety2_merged[[name]]$nreads_inv)]<-0
#   
# #Match suspicious(inv+non_inv+noExact match) seqs on the whole base
#    for (i in which(safety2_merged[[name]]$nreads_inv!=0&safety2_merged[[name]]$nreads!=0&safety2_merged[[name]]$Exact_new=="", )) {
#   VC1<-vcountPattern(pattern = safety2_merged[[name]]$assembled[i], subject = HLA_base_str,
#                      max.mismatch = 0, with.indels = F, fixed=F)
#   safety2_merged[[name]]$Exact_new[i]<-paste0(HLA_allele_list$Allele[VC1!=0], collapse=" ")
# 
#    }
#   }

  safety3<-list(HLA_allele_Iclass=HLAI[[1]], HLA_Iclass_amps=HLAI[[2]],
                HLA_allele_DQB=DQB[[1]], HLA_DQB_amps=DQB[[2]],
                HLA_allele_DRB=DRB[[1]], HLA_DRB_amps=DRB[[2]],
                HLA_allele_DPB=DPB[[1]], HLA_DPB_amps=DPB[[2]],
                HLA_allele_DQA=DQA[[1]], HLA_DQA_amps=DQA[[2]],
                n_readsI=HLAI[[3]], n_readsDQB=DQB[[3]], n_readsDRB=DRB[[3]],
                n_readsDPB=DPB[[3]], n_readsDQA=DQA[[3]],
                safety2_merged=safety2_merged)

  safety3
}


get_safety4<-function(safety3) {
  safety3$HLA_allele_Iclass$amps<-apply(safety3$HLA_Iclass_amps[2:9], MARGIN = 1, paste0, collapse = "_")
  safety3$HLA_allele_DQB$amps<-apply(safety3$HLA_DQB_amps[2:9], MARGIN = 1, paste0, collapse = "_")
  safety3$HLA_allele_DRB$amps<-apply(safety3$HLA_DRB_amps[2:9], MARGIN = 1, paste0, collapse = "_")
  
  safety3$HLA_allele_DPB$amps<-apply(safety3$HLA_DPB_amps[2:5], MARGIN = 1, paste0, collapse = "_")
  safety3$HLA_allele_DQA$amps<-apply(safety3$HLA_DQA_amps[2:3], MARGIN = 1, paste0, collapse = "_")
  
  for (i in 1:nrow(safety3$HLA_allele_Iclass)) {
    safety3$HLA_allele_Iclass$no_amps[i]<-paste0(colnames(safety3$HLA_Iclass_amps[which(safety3$HLA_Iclass_amps[i, ]==0)]), collapse = "   ")
  }
  for (i in 1:nrow(safety3$HLA_allele_DQB)) {
    safety3$HLA_allele_DQB$no_amps[i]<-paste0(colnames(safety3$HLA_DQB_amps[which(safety3$HLA_DQB_amps[i, ]==0)]), collapse = "   ")
  }
  for (i in 1:nrow(safety3$HLA_allele_DRB)) {
    safety3$HLA_allele_DRB$no_amps[i]<-paste0(colnames(safety3$HLA_DRB_amps[which(safety3$HLA_DRB_amps[i, ]==0)]), collapse = "   ")
  }
  
  for (i in 1:nrow(safety3$HLA_allele_DPB)) {
    safety3$HLA_allele_DPB$no_amps[i]<-paste0(colnames(safety3$HLA_DPB_amps[which(safety3$HLA_DPB_amps[i, ]==0)]), collapse = "   ")
  }
for (i in 1:nrow(safety3$HLA_allele_DQA)) {
    safety3$HLA_allele_DQA$no_amps[i]<-paste0(colnames(safety3$HLA_DQA_amps[which(safety3$HLA_DQA_amps[i, ]==0)]), collapse = "   ")
  }
  
  #safety3$HLA_allele_Iclass[3:10]<-apply(safety3$HLA_allele_Iclass[3:10], MARGIN = 2, prop.table)
  DT1<-as.data.table(safety3$HLA_allele_Iclass)
  DQB<-as.data.table(safety3$HLA_allele_DQB)
  DRB<-as.data.table(safety3$HLA_allele_DRB)
  
  DPB<-as.data.table(safety3$HLA_allele_DPB)
  DQA<-as.data.table(safety3$HLA_allele_DQA)
  
  DT1<-DT1[amps!="0_0_0_0_0_0_0_0", .(Allele=paste0(Allele, collapse=" "),
                                    Iamp1=unique(Iamp1),
                                    Iamp1_inv=unique(Iamp1_inv),
                                    Iamp2=unique(Iamp2),
                                    Iamp2_inv=unique(Iamp2_inv),
                                    Iamp1alt=unique(Iamp1alt),
                                    Iamp1alt_inv=unique(Iamp1alt_inv),
                                    Iamp2alt=unique(Iamp2alt),
                                    Iamp2alt_inv=unique(Iamp2alt_inv), 
                                    no_amps=unique(no_amps)),by=amps]
  
  DQB<-DQB[amps!="0_0_0_0_0_0_0_0", .(Allele=paste0(Allele, collapse=" "),
                                      IIamp1_DQB=unique(IIamp1_DQB),
                                      IIamp1_DQB_inv=unique(IIamp1_DQB_inv),
                                      IIamp2=unique(IIamp2),
                                      IIamp2_inv=unique(IIamp2_inv),
                                      IIamp1_DQBalt=unique(IIamp1_DQBalt),
                                      IIamp1_DQBalt_inv=unique(IIamp1_DQBalt_inv),
                                      IIamp2_DQBalt=unique(IIamp2_DQBalt),
                                      IIamp2_DQBalt_inv=unique(IIamp2_DQBalt_inv),
                                      no_amps=unique(no_amps)), by=amps]
  
  DRB<-DRB[amps!="0_0_0_0_0_0_0_0", .(Allele=paste0(Allele, collapse=" "),
                                      IIamp1_Others=unique(IIamp1_Others),
                                      IIamp1_Others_inv=unique(IIamp1_Others_inv),
                                      IIamp2=unique(IIamp2),
                                      IIamp2_inv=unique(IIamp2_inv),
                                      IIamp1_DRBalt=unique(IIamp1_DRBalt),
                                      IIamp1_DRBalt_inv=unique(IIamp1_DRBalt_inv),
                                      IIamp2_DRBalt=unique(IIamp2_DRBalt),
                                      IIamp2_DRBalt_inv=unique(IIamp2_DRBalt_inv),
                                      no_amps=unique(no_amps)), by=amps]
  
  DPB<-DPB[amps!="0_0_0_0", .(Allele=paste0(Allele, collapse=" "),
                              IIamp1_DPB=unique(IIamp1_DPB),
                              IIamp1_DPB_inv=unique(IIamp1_DPB_inv),
                              IIamp2_DPB=unique(IIamp2_DPB),
                              IIamp2_DPB_inv=unique(IIamp2_DPB_inv),
                              no_amps=unique(no_amps)), by=amps]
  DQA<-DQA[amps!="0_0", .(Allele=paste0(Allele, collapse=" "),
                          IIamp_DQA=unique(IIamp_DQA),
                          IIamp_DQA_inv=unique(IIamp_DQA_inv),
                          no_amps=unique(no_amps)), by=amps]
  
  DT1<-as.data.frame(DT1)
  DQB<-as.data.frame(DQB)
  DRB<-as.data.frame(DRB)
  DPB<-as.data.frame(DPB)
  DQA<-as.data.frame(DQA)
  
  DT1$sumreads<-apply(DT1[3:10], 1, function(x) {sum(x*safety3$n_readsI)})
  DQB$sumreads<-apply(DQB[3:10], 1, function(x) {sum(x*safety3$n_readsDQB)})
  DRB$sumreads<-apply(DRB[3:10], 1, function(x) {sum(x*safety3$n_readsDRB)})
  DPB$sumreads<-apply(DPB[3:6], 1, function(x) {sum(x*safety3$n_readsDPB)})
  DQA$sumreads<-apply(DQA[3:4], 1, function(x) {sum(x*safety3$n_readsDQA)})
  
  # DT1$meanreads<-apply(DT1[3:10], 1, function(x) {mean(x*safety3$n_readsI)})
  # DQB$meanreads<-apply(DQB[3:10], 1, function(x) {mean(x*safety3$n_readsDQB)})
  # DRB$meanreads<-apply(DRB[3:10], 1, function(x) {mean(x*safety3$n_readsDRB)})
  # 
  # DT1$medianreads<-apply(DT1[3:10], 1, function(x) {median(x*safety3$n_readsI)})
  # DQB$medianreads<-apply(DQB[3:10], 1, function(x) {median(x*safety3$n_readsDQB)})
  # DRB$medianreads<-apply(DRB[3:10], 1, function(x) {median(x*safety3$n_readsDRB)})
  
  Ps<-rep(0.05, times=8)
  for (i in 3:10) {
    DT1[,i][DT1[,i]>0]<- (1-Ps[i-2])#(1-exp(-safety3$Ns[i-2]*Ps[i-2]))
    DT1[,i][DT1[,i]==0]<- Ps[i-2]#exp(-safety3$Ns[i-2]*Ps[i-2])
    DQB[,i][DQB[,i]>0]<- (1-Ps[i-2])
    DQB[,i][DQB[,i]==0]<- Ps[i-2]
    DRB[,i][DRB[,i]>0]<- (1-Ps[i-2])
    DRB[,i][DRB[,i]==0]<- Ps[i-2]
  }
  
  DT1$Score<-apply(DT1[3:10], 1, prod)
  DQB$Score<-apply(DQB[3:10], 1, prod)
  DRB$Score<-apply(DRB[3:10], 1, prod)
  
  Ps<-rep(0.05, times=4)
  for (i in 3:6) {
    DPB[,i][DPB[,i]>0]<- (1-Ps[i-2])#(1-exp(-safety3$Ns[i-2]*Ps[i-2]))
    DPB[,i][DPB[,i]==0]<- Ps[i-2]#exp(-safety3$Ns[i-2]*Ps[i-2])
   
  }
  
  DPB$Score<-apply(DPB[3:6], 1, prod)
  
  Ps<-rep(0.05, times=2)
  for (i in 3:4) {
    DQA[,i][DQA[,i]>0]<- (1-Ps[i-2])#(1-exp(-safety3$Ns[i-2]*Ps[i-2]))
    DQA[,i][DQA[,i]==0]<- Ps[i-2]#exp(-safety3$Ns[i-2]*Ps[i-2])
  }
  
  DQA$Score<-apply(DQA[3:4], 1, prod)


  amps_list1<-str_split(DT1$amps, pattern = "_")
  
  for (i in 1:length(amps_list1)) {
    amps_list1[[i]]<-str_split(amps_list1[[i]], pattern = "," )
  }
  amps_listDQB<-str_split(DQB$amps, pattern = "_")

  for (i in 1:length(amps_listDQB)) {
    amps_listDQB[[i]]<-str_split(amps_listDQB[[i]], pattern = "," )
  }
  amps_listDRB<-str_split(DRB$amps, pattern = "_")

  for (i in 1:length(amps_listDRB)) {
    amps_listDRB[[i]]<-str_split(amps_listDRB[[i]], pattern = "," )
  }
  amps_listDPB<-str_split(DPB$amps, pattern = "_")
  for (i in 1:length(amps_listDPB)) {
    amps_listDPB[[i]]<-str_split(amps_listDPB[[i]], pattern = "," )
  }
  amps_listDQA<-str_split(DQA$amps, pattern = "_")
  if(length(amps_listDQA)!=0)#zatychka
  for (i in 1:length(amps_listDQA)) {
    amps_listDQA[[i]]<-str_split(amps_listDQA[[i]], pattern = "," )
  }
  
  
  isdaddy <- function(list1, list2) {
    a<-length(setdiff(list2[[1]], list1[[1]]))==0
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
  
  isdaddy4 <- function(list1, list2) {
    a<-length(setdiff(list2[[1]], list1[[1]]))==0
    b<-length(setdiff(list2[[2]], list1[[2]]))==0
    c<-length(setdiff(list2[[3]], list1[[3]]))==0
    d<-length(setdiff(list2[[4]], list1[[4]]))==0
    daddy<-sum(c(a,b,c,d))==4
    daddy
  }
  isdaddy2 <- function(list1, list2) {
    a<-length(setdiff(list2[[1]], list1[[1]]))==0
    b<-length(setdiff(list2[[2]], list1[[2]]))==0
    daddy<-sum(c(a,b))==2
    daddy
  }
  
  res1<-matrix(0, ncol=nrow(DT1), nrow=nrow(DT1))
  resDQB<-matrix(0, ncol=nrow(DQB), nrow=nrow(DQB))
  resDRB<-matrix(0, ncol=nrow(DRB), nrow=nrow(DRB))
  resDPB<-matrix(0, ncol=nrow(DPB), nrow=nrow(DPB))
  resDQA<-matrix(0, ncol=nrow(DQA), nrow=nrow(DQA))
 
   for (i in 1:length(amps_list1))
    for (j in 1:length(amps_list1)) {
      res1[i,j]<-isdaddy(amps_list1[[i]], amps_list1[[j]])
    } 
  for (i in 1:length(amps_listDQB))
    for (j in 1:length(amps_listDQB)) {
      resDQB[i,j]<-isdaddy(amps_listDQB[[i]], amps_listDQB[[j]])
    } 
  for (i in 1:length(amps_listDRB))
    for (j in 1:length(amps_listDRB)) {
      resDRB[i,j]<-isdaddy(amps_listDRB[[i]], amps_listDRB[[j]])
    } 
  
  for (i in 1:length(amps_listDPB))
    for (j in 1:length(amps_listDPB)) {
      resDPB[i,j]<-isdaddy4(amps_listDPB[[i]], amps_listDPB[[j]])
    } 
  if(length(amps_listDQA)!=0)
  for (i in 1:length(amps_listDQA))
    for (j in 1:length(amps_listDQA)) {
      resDQA[i,j]<-isdaddy2(amps_listDQA[[i]], amps_listDQA[[j]])
    } 
  
  gr1<-igraph::simplify(graph_from_adjacency_matrix(res1))
  grDQB<-igraph::simplify(graph_from_adjacency_matrix(resDQB))
  grDRB<-igraph::simplify(graph_from_adjacency_matrix(resDRB))
  grDPB<-igraph::simplify(graph_from_adjacency_matrix(resDPB))
  grDQA<-igraph::simplify(graph_from_adjacency_matrix(resDQA))
  
  
  DT1$degree<-degree(gr1, mode="in")
  DQB$degree<-degree(grDQB, mode="in")
  DRB$degree<-degree(grDRB, mode="in")
  DPB$degree<-degree(grDPB, mode="in")
  DQA$degree<-degree(grDQA, mode="in")
  
  list(DT1=DT1, gr1=gr1, DQB=DQB, grDQB=grDQB, DRB=DRB, grDRB=grDRB, DPB=DPB, grDPB=grDPB, DQA=DQA, grDQA=grDQA)
}

get_safety5<-function(DT1, DQB, DRB, DPB, DQA) {HLA_typing_resultsIclass<-filter(DT1, (DT1$degree==0)|(DT1$no_amps==""))%>%select(Allele, Score, no_amps, sumreads)

HLA_typing_resultsDQB<-filter(DQB, (DQB$degree==0)|(DQB$no_amps==""))%>%select(Allele, Score, no_amps, sumreads)
HLA_typing_resultsDRB<-filter(DRB, (DRB$degree==0)|(DRB$no_amps==""))%>%select(Allele, Score, no_amps, sumreads)
HLA_typing_resultsDPB<-filter(DPB, (DPB$degree==0)|(DPB$no_amps==""))%>%select(Allele, Score, no_amps, sumreads)
HLA_typing_resultsDQA<-filter(DQA, (DQA$degree==0)|(DQA$no_amps==""))%>%select(Allele, Score, no_amps, sumreads)

HLA_all<-list(HLA_typing_resultsIclass, HLA_typing_resultsDQB, HLA_typing_resultsDRB, HLA_typing_resultsDPB, HLA_typing_resultsDQA)
HLA_typing_results<-do.call(rbind, HLA_all)

HLA_typing_results[grepl(HLA_typing_results$Allele, pattern = "A*", fixed = T), ]$Score<-prop.table(HLA_typing_results[grepl(HLA_typing_results$Allele, pattern = "A*", fixed = T), ]$Score)
HLA_typing_results[grepl(HLA_typing_results$Allele, pattern = "B*",fixed = T), ]$Score<-prop.table(HLA_typing_results[grepl(HLA_typing_results$Allele, pattern = "B*", fixed = T), ]$Score)
HLA_typing_results[grepl(HLA_typing_results$Allele, pattern = "C*", fixed = T), ]$Score<-prop.table(HLA_typing_results[grepl(HLA_typing_results$Allele, pattern = "C*", fixed = T), ]$Score)
HLA_typing_results[grepl(HLA_typing_results$Allele, pattern = "DQB", fixed = T), ]$Score<-prop.table(HLA_typing_results[grepl(HLA_typing_results$Allele, pattern = "DQB", fixed = T), ]$Score)
HLA_typing_results[grepl(HLA_typing_results$Allele, pattern = "DRB1", fixed = T), ]$Score<-prop.table(HLA_typing_results[grepl(HLA_typing_results$Allele, pattern = "DRB1", fixed = T), ]$Score)
HLA_typing_results[grepl(HLA_typing_results$Allele, pattern = "DRB[2-9]", fixed = F), ]$Score<-prop.table(HLA_typing_results[grepl(HLA_typing_results$Allele, pattern = "DRB[2-9]", fixed = F), ]$Score)

HLA_typing_results[grepl(HLA_typing_results$Allele, pattern = "DPB", fixed = F), ]$Score<-prop.table(HLA_typing_results[grepl(HLA_typing_results$Allele, pattern = "DPB", fixed = F), ]$Score)
HLA_typing_results[grepl(HLA_typing_results$Allele, pattern = "DQA", fixed = F), ]$Score<-prop.table(HLA_typing_results[grepl(HLA_typing_results$Allele, pattern = "DQA", fixed = F), ]$Score)


HLA_typing_results<-select(HLA_typing_results, Allele, Score, no_amps, sumreads)
}


signs3<-function(string) {
  paste(str_split(string, pattern = fixed(":"), simplify = T )[1:3], sep=":", collapse = ":") 
  
}

gettidyHLA<-function (longstr) {
  paste(unique(sapply(str_split(longstr, pattern=fixed(" "), simplify = T), signs3)), collapse = " ")
}


hla_wide <- function(x) {
  d3c<-read.csv2(x,stringsAsFactors = F)
  
  d3c$Locus<-strsplit(x = d3c$tidyAllele, split = "*", fixed=T)%>%sapply("[[", 1)
  d3c$FirstAllele<-strsplit(x = d3c$tidyAllele, split = " ", fixed=T)%>%sapply("[[", 1)
  d3c$SimpleAllele<-str_sub(string = d3c$FirstAllele, start = 1, end = -4)
  d3c_long<-d3c%>%
    select(donor, Locus, FirstAllele, SimpleAllele, tidyAllele)
  
  d3c_wide<-d3c%>%
    select(donor, Locus, FirstAllele)
  d3c_wide$FirstAllele<-gsub(pattern = "DPB1*107:01", replacement = "DPB1*13:01:01",
                             x = d3c_wide$FirstAllele, fixed = T)
  d3c_wide$FirstAllele<-gsub(pattern = "B*07:100 ", replacement = "B*15:01:01",
                             x = d3c_wide$FirstAllele, fixed = T)
  d3c_wide$FirstAllele<-gsub(pattern = "DPB1*138:01", replacement = "DQB1*23:01:01",
                             x = d3c_wide$FirstAllele, fixed = T)
  
  d3c_wide<-d3c_wide%>%
    group_by(donor)%>%
    mutate(n=as.integer(duplicated(Locus))+1)
  
  d3c_wide$Locus<-paste0(d3c_wide$Locus, "_", d3c_wide$n)
  d3c_wide<-d3c_wide%>%
    select(1:3)
  d3c_wide<-spread(d3c_wide, Locus, FirstAllele)
  
  d3c_wide$genotype<-apply(d3c_wide[,c(-1)], MARGIN=1, paste0, collapse="_")
  d3c_wide$genotype2<-apply(d3c_wide[,c(2:7,12:15)], MARGIN=1, paste0, collapse="_")
  DRB<-grepl(pattern="DRB1[*]03:01",
             x=d3c_wide$genotype)
  DQB<-grepl(pattern="DQB1[*]02:01",
             x=d3c_wide$genotype)
  DQA<-grepl(pattern="DQA1[*]05:01",
             x=d3c_wide$genotype)
  d3c_wide$DR3<-DRB==T&DQB==T&DQA==T
  
  DRB1<-grepl(pattern="DRB1[*]04:[01|02|04|05|08]",
              x=d3c_wide$genotype)
  DQB2<-grepl(pattern="DQB1[*]03:[01|05]",
              x=d3c_wide$genotype)
  DQB1<-grepl(pattern="DQB1[*]02",
              x=d3c_wide$genotype)
  DQA3<-grepl(pattern="DQA1[*]03:01",
              x=d3c_wide$genotype)
  
  DQB3<-(DQB2==T|DQB1==T)
  
  d3c_wide$DR4<-as.logical(DRB1*DQB3*DQA3)
  
  d3c_wide$same_genotype<-0
  d3c_wide$same_genotype2<-0
  
  for (i in (1:nrow(d3c_wide))) {
    d3c_wide$same_genotype[i]<-paste(d3c_wide$donor[grepl(pattern = d3c_wide$genotype[i], 
                                                          x = d3c_wide$genotype, fixed = T)], collapse="_")
    d3c_wide$same_genotype2[i]<-paste(d3c_wide$donor[grepl(pattern = d3c_wide$genotype2[i], 
                                                           x = d3c_wide$genotype2, fixed = T)], collapse="_")
    
  }
  
  list(wide=d3c_wide, long=d3c_long)
}


#Megatable filtered by hand as an input (no more than 2 alleles of each type: A, B,C..etc.)
#diab<-hla_wide("diab_all.csv")

#write.csv2(diab$wide, file="hla_diab_all_wide_final.csv")
#write.csv2(diab$long, file="hla_diab_all_long.csv")

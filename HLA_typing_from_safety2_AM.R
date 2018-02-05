#Download all nessesary packages

library(dplyr)
library(data.table)
library(cowsay)
library(Biostrings)
library(stringr)
library(igraph)
library(tibble)
library(parallel)


#Function to get one table for all people with tidy alleles and scores. Takes list of results (output of get_HLA_result)
#as an input
get_mega_table <- function(Typing_results) {
  Mega_table_first<-do.call(rbind, lapply(Typing_results, function(x){x$safety5}))
  Mega_table_first<-rownames_to_column(df = Mega_table_first)
  Mega_table_first$donor<-sapply(str_split(Mega_table_first$rowname, fixed(".")), function(x) {x[[1]]})
  Mega_table_first<-select(Mega_table_first, donor, Allele, Score, no_amps, sumreads, meanreads, medianreads)
  Mega_table_first$tidyAllele<-sapply(Mega_table_first$Allele, gettidyHLA)
  Mega_table_first$tidyAllele<-gsub(Mega_table_first$tidyAllele, pattern = ":NA", replacement = "", fixed = T)
  #Mega_table_100717<-Mega_table_100717[Mega_table_100717$Score>1e-5, ]
  Mega_table_first$n_noamps<-sapply(strsplit(Mega_table_first$no_amps, split = "   "), length)
  Mega_table_first<-select(Mega_table_first, donor, tidyAllele, Allele, Score, n_noamps, sumreads, meanreads, medianreads, no_amps)
  Mega_table_first
}

#Function to get list of results. Input - safety 2(output of HLA_amplicones_full)
get_HLA_result<-function (donor) {  say("START!")

    print(format(Sys.time(), "%a %b %d %X %Y"))
  
  if (sum(sapply(donor$safety2, nrow)[c(1,2,6,7,12,13,17,18)])>0)
    #if (donor$safety2!="fail")
    {
    safety3<-get_safety3(lapply(donor$safety2, function(x){x[((x$parents==1&x$freq>0.001)|(x$parents==2&x$freq>0.05))&!grepl(pattern = "NULL",x$assembled), ]}))
    safety4<-get_safety4(safety3)
    safety5<-get_safety5(DT1=safety4$DT1, DQB=safety4$DQB, DRB=safety4$DRB)
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
  HLA_allele_DRB<-select(HLA_allele_list, Allele, Sequence, HLA_class)%>%filter(grepl(HLA_allele_list$HLA_class, pattern="D.B*", fixed=F))%>%select(Allele, Sequence)
  HLA_allele_DQB<-select(HLA_allele_list, Allele, Sequence, HLA_class)%>%filter(grepl(HLA_allele_list$HLA_class, pattern="D.B*", fixed=F))%>%select(Allele, Sequence)
  
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
  
  HLA_Iclass_amps<-select(HLA_allele_Iclass, -2)
  HLA_DQB_amps<-select(HLA_allele_DQB, -2)
  HLA_DRB_amps<-select(HLA_allele_DRB, -2)
  
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
  names(n_readsI)<-colnames(HLA_allele_Iclass)[3:10]
  names(n_readsDQB)<-colnames(HLA_allele_DQB)[3:10]
  names(n_readsDRB)<-colnames(HLA_allele_DRB)[3:10]

#Function to get properly named tables of frequencies and indices
freqs_indices_tables <- function(HLA_allele, HLA_amps, n_reads, safety2) {
  for (colname in intersect(colnames(HLA_allele)[3:10], names(safety2))) {
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
HLAI<-freqs_indices_tables(HLA_allele_Iclass, HLA_Iclass_amps, n_readsI, safety2)
DQB<-freqs_indices_tables(HLA_allele_DQB, HLA_DQB_amps, n_readsDQB, safety2=HLAI[[4]])
DRB<-freqs_indices_tables(HLA_allele_DRB, HLA_DRB_amps, n_readsDRB, safety2=DQB[[4]])
safety2<-DRB[[4]]
 
#Get merged inv and non inv lists to search for possible new allele variants  
  safety2_merged<-list()
  HLA_base_str<-DNAStringSet(HLA_allele_list$Sequence)
  for (name in  names(safety2)[!grepl(pattern = "inv", x = names(safety2), fixed=T)]) 
  if(nrow(safety2[[name]])!=0&nrow(safety2[[paste0(name, "_inv", sep="")]])!=0)
  {
  name2<-paste0(name, "_inv", sep="")
  safety2_merged[[name]]<-full_join(x = safety2[[name]], y= safety2[[name2]], by=c("assembled", "Exact_new"))
  safety2_merged[[name]]<-safety2_merged[[name]]%>%select("assembled",  nreads="readnumber.x", nreads_inv="readnumber.y", "Exact_new")
  safety2_merged[[name]]$nreads[is.na(safety2_merged[[name]]$nreads)]<-0
  safety2_merged[[name]]$nreads_inv[is.na(safety2_merged[[name]]$nreads_inv)]<-0
  
#Match suspicious(inv+non_inv+noExact match) seqs on the whole base
   for (i in which(safety2_merged[[name]]$nreads_inv!=0&safety2_merged[[name]]$nreads!=0&safety2_merged[[name]]$Exact_new=="", )) {
  VC1<-vcountPattern(pattern = safety2_merged[[name]]$assembled[i], subject = HLA_base_str,
                     max.mismatch = 0, with.indels = F, fixed=F)
  safety2_merged[[name]]$Exact_new[i]<-paste0(HLA_allele_list$Allele[VC1!=0], collapse=" ")

   }
  }

  safety3<-list(HLA_allele_Iclass=HLAI[[1]], HLA_Iclass_amps=HLAI[[2]],
                HLA_allele_DQB=DQB[[1]], HLA_DQB_amps=DQB[[2]],
                HLA_allele_DRB=DRB[[1]], HLA_DRB_amps=DRB[[2]],
                n_readsI=HLAI[[3]], n_readsDQB=DQB[[3]], n_readsDRB=DRB[[3]],
                safety2_merged=safety2_merged)

  safety3
}


get_safety4<-function(safety3) {
  safety3$HLA_allele_Iclass$amps<-apply(safety3$HLA_Iclass_amps[2:9], MARGIN = 1, paste0, collapse = "_")
  safety3$HLA_allele_DQB$amps<-apply(safety3$HLA_DQB_amps[2:9], MARGIN = 1, paste0, collapse = "_")
  safety3$HLA_allele_DRB$amps<-apply(safety3$HLA_DRB_amps[2:9], MARGIN = 1, paste0, collapse = "_")
  
  for (i in 1:nrow(safety3$HLA_allele_Iclass)) {
    safety3$HLA_allele_Iclass$no_amps[i]<-paste0(colnames(safety3$HLA_Iclass_amps[which(safety3$HLA_Iclass_amps[i, ]==0)]), collapse = "   ")
  }
  for (i in 1:nrow(safety3$HLA_allele_DQB)) {
    safety3$HLA_allele_DQB$no_amps[i]<-paste0(colnames(safety3$HLA_DQB_amps[which(safety3$HLA_DQB_amps[i, ]==0)]), collapse = "   ")
  }
  for (i in 1:nrow(safety3$HLA_allele_DRB)) {
    safety3$HLA_allele_DRB$no_amps[i]<-paste0(colnames(safety3$HLA_DRB_amps[which(safety3$HLA_DRB_amps[i, ]==0)]), collapse = "   ")
  }
  
  #safety3$HLA_allele_Iclass[3:10]<-apply(safety3$HLA_allele_Iclass[3:10], MARGIN = 2, prop.table)
  DT1<-as.data.table(safety3$HLA_allele_Iclass)
  DQB<-as.data.table(safety3$HLA_allele_DQB)
  DRB<-as.data.table(safety3$HLA_allele_DRB)
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
  DT1<-as.data.frame(DT1)
  DQB<-as.data.frame(DQB)
  DRB<-as.data.frame(DRB)
  
  DT1$sumreads<-apply(DT1[3:10], 1, function(x) {sum(x*safety3$n_readsI)})
  DQB$sumreads<-apply(DQB[3:10], 1, function(x) {sum(x*safety3$n_readsDQB)})
  DRB$sumreads<-apply(DRB[3:10], 1, function(x) {sum(x*safety3$n_readsDRB)})
  
  DT1$meanreads<-apply(DT1[3:10], 1, function(x) {mean(x*safety3$n_readsI)})
  DQB$meanreads<-apply(DQB[3:10], 1, function(x) {mean(x*safety3$n_readsDQB)})
  DRB$meanreads<-apply(DRB[3:10], 1, function(x) {mean(x*safety3$n_readsDRB)})
  
  DT1$medianreads<-apply(DT1[3:10], 1, function(x) {median(x*safety3$n_readsI)})
  DQB$medianreads<-apply(DQB[3:10], 1, function(x) {median(x*safety3$n_readsDQB)})
  DRB$medianreads<-apply(DRB[3:10], 1, function(x) {median(x*safety3$n_readsDRB)})
  
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
  res1<-matrix(0, ncol=nrow(DT1), nrow=nrow(DT1))
  resDQB<-matrix(0, ncol=nrow(DQB), nrow=nrow(DQB))
  resDRB<-matrix(0, ncol=nrow(DRB), nrow=nrow(DRB))
 
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
  
  gr1<-igraph::simplify(graph_from_adjacency_matrix(res1))
  grDQB<-igraph::simplify(graph_from_adjacency_matrix(resDQB))
  grDRB<-igraph::simplify(graph_from_adjacency_matrix(resDRB))
  
  DT1$degree<-degree(gr1, mode="in")
  DQB$degree<-degree(grDQB, mode="in")
  DRB$degree<-degree(grDRB, mode="in")
  
  list(DT1=DT1, gr1=gr1, DQB=DQB, grDQB=grDQB, DRB=DRB, grDRB=grDRB)
}

get_safety5<-function(DT1, DQB, DRB) {HLA_typing_resultsIclass<-filter(DT1, DT1$degree==0)%>%select(Allele, Score, no_amps, sumreads, meanreads, medianreads)
HLA_typing_resultsDQB<-filter(DQB, DQB$degree==0)%>%select(Allele, Score, no_amps, sumreads, meanreads, medianreads)
HLA_typing_resultsDRB<-filter(DRB, DRB$degree==0)%>%select(Allele, Score, no_amps, sumreads, meanreads, medianreads)
HLA_all<-list(HLA_typing_resultsIclass, HLA_typing_resultsDQB, HLA_typing_resultsDRB)
HLA_typing_results<-do.call(rbind, HLA_all)

HLA_typing_results[grepl(HLA_typing_results$Allele, pattern = "A*", fixed = T), ]$Score<-prop.table(HLA_typing_results[grepl(HLA_typing_results$Allele, pattern = "A*", fixed = T), ]$Score)
HLA_typing_results[grepl(HLA_typing_results$Allele, pattern = "B*",fixed = T), ]$Score<-prop.table(HLA_typing_results[grepl(HLA_typing_results$Allele, pattern = "B*", fixed = T), ]$Score)
HLA_typing_results[grepl(HLA_typing_results$Allele, pattern = "C*", fixed = T), ]$Score<-prop.table(HLA_typing_results[grepl(HLA_typing_results$Allele, pattern = "C*", fixed = T), ]$Score)
HLA_typing_results[grepl(HLA_typing_results$Allele, pattern = "DQB", fixed = T), ]$Score<-prop.table(HLA_typing_results[grepl(HLA_typing_results$Allele, pattern = "DQB", fixed = T), ]$Score)
HLA_typing_results[grepl(HLA_typing_results$Allele, pattern = "DRB1", fixed = T), ]$Score<-prop.table(HLA_typing_results[grepl(HLA_typing_results$Allele, pattern = "DRB1", fixed = T), ]$Score)
HLA_typing_results[grepl(HLA_typing_results$Allele, pattern = "DRB[2-9]", fixed = F), ]$Score<-prop.table(HLA_typing_results[grepl(HLA_typing_results$Allele, pattern = "DRB[2-9]", fixed = F), ]$Score)


HLA_typing_results<-select(HLA_typing_results, Allele, Score, no_amps, sumreads, meanreads, medianreads)
}


signs3<-function(string) {
  paste(str_split(string, pattern = fixed(":"), simplify = T )[1:3], sep=":", collapse = ":") 
  
}

gettidyHLA<-function (longstr) {
  paste(unique(sapply(str_split(longstr, pattern=fixed(" "), simplify = T), signs3)), collapse = " ")
}



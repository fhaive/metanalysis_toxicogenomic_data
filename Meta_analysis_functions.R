#meta analysis functions

#needed libraries
library(esc)
library(metafor)
library(metap)
library(TopKLists)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("RankProd")
library(RankProd)
library(matrixStats)

#compute effect sizes for meta analysis
calc_effect_size_rank<- function(meta_data) {
  index_pval<-which(grepl('P.Value',colnames(meta_data)))
  meta_data_pval<-meta_data[,index_pval]
  colnames(meta_data_pval)<-colnames(meta_data)[index_pval]
  #heatmap(as.matrix(meta_data_pval))
  tmp_data<-meta_data_pval
  row.names(tmp_data)<-NULL
  pval<-as.vector(rowMeans(tmp_data))
  tmp <- data.frame(
    pval = pval,
    n =rep(ncol(meta_data_pval),length(pval)),
    studyname = rownames(meta_data_pval)
  )
  effect_sizes_values<-effect_sizes(tmp, p = pval, totaln = n, study = studyname, fun = "chisq")
  #check is study label is kept, in case add colnames(meta_data_log)
  effect_size_rank<-effect_sizes_values[,c('study','es')]
  rownames(effect_size_rank)<-effect_size_rank$study
  effect_size_rank$study<-NULL
  #ranked final list 
  effect_size_rank<-effect_size_rank[order(-effect_size_rank$es), , drop = FALSE]
  return(effect_size_rank)
}

#compute pvalue based rank for meta analysis 
calc_pvalue_based_rank<-function(meta_data){
  #retrieve data
  index_pval<-which(grepl('P.Value',colnames(meta_data)))
  meta_data_pval<-as.data.frame(meta_data[,index_pval])
  colnames(meta_data_pval)<-colnames(meta_data)[index_pval]
  #for each gene combine the p-values by the sum of logs method
  fisher_based_res<-list()
  for(i in 1:length(rownames(meta_data_pval))){
    metap<-sumlog(meta_data_pval[i,])
    fisher_based_res[[i]]<-metap$p
  } 
  fisher_based_pvalues<-as.data.frame(fisher_based_res)
  colnames(fisher_based_pvalues)<-rownames(meta_data_pval)
  fisher_based_pvalues<-t(fisher_based_pvalues)
  fisher_based_pvalues<-as.data.frame(fisher_based_pvalues)
  colnames(fisher_based_pvalues)<-"pValue"
  #ranked final list 
  fisher_based_pvalues<- fisher_based_pvalues[order(fisher_based_pvalues$pValue), , drop = FALSE]
  return(fisher_based_pvalues)
}

#compute rank based rank for meta analysis 
calc_rank_base_rank<-function(meta_data){
  index_pval_adj<-which(grepl('_adj.P.Val',colnames(meta_data)))
  meta_data_adpval<-as.data.frame(meta_data[,index_pval_adj])
  colnames(meta_data_adpval)<-colnames(meta_data)[index_pval_adj]
  cl<-rep.int(1,times = length(colnames(meta_data_adpval)))
  rp.advance.input<-meta_data_adpval
  colnames(rp.advance.input)<-NULL
  rownames(rp.advance.input)<-NULL
  rp.advance.input<-as.matrix(rp.advance.input)
  #origin contains the  labels for different studies
  origin<-gsub(pattern = "_adj.P.Val",colnames(meta_data_adpval),replacement = "")
  origin<-gsub("\\_.*","",origin)
  o<-rep(1,dim(meta_data_adpval)[2])
  RP_advance_out<-RP.advance(data = meta_data_adpval, cl = o, origin = o, calculateProduct =TRUE)
  #ranked  final list
  ranks_based_pvalues<-as.data.frame(RP_advance_out$pval)
  rownames(ranks_based_pvalues)<-rownames(meta_data_adpval)
# ranks_based_pvalues<-ranks_based_pvalues[order(abs(ranks_based_pvalues$`class1 < class2`),abs(ranks_based_pvalues$`class1 > class2`)), , drop = FALSE]
  return(ranks_based_pvalues)
}


calc_rank_base_rank_subsets<-function(meta_data){
  index_pval_adj<-which(grepl('_adj.P.Val',colnames(meta_data)))
  meta_data_adpval<-as.data.frame(meta_data[,index_pval_adj])
  colnames(meta_data_adpval)<-colnames(meta_data)[index_pval_adj]
  ranks<-list()
  for(i in 1:10){
    data<-meta_data_adpval
    index_col<-sample(colnames(data),size = 75)
    data_partition<-data[,index_col]
    ranks_product<-calc_rank_base_rank(data_partition)
    ranks[[i]]<-ranks_product
  }
  class1<-data.frame()
  class2<-data.frame()
  for(j in 1:length(ranks)){
    if(plyr::empty(class1)){
      class1<-as.data.frame(ranks[[j]][,1])
    } else{
      class1<-qpcR:::cbind.na(class1,ranks[[j]][,1])
    } 
    if(plyr::empty(class2)){
      class2<-as.data.frame(ranks[[j]][,2])
    } 
    else{
      class2<-qpcR:::cbind.na(class2,ranks[[j]][,2])
    } 
  }
  ranks_based_pvalues<-data.frame("class1 < class2"= rowMeans(class1),
                                  "class1 > class2" = rowMeans(class2))
  rownames(ranks_based_pvalues)<-rownames(meta_data_adpval)
  ranks_based_pvalues<-ranks_based_pvalues[order(abs(ranks_based_pvalues$class1...class2),abs(ranks_based_pvalues$class1...class2.1)), , drop = FALSE]
  return(ranks_based_pvalues)
}


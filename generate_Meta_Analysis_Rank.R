##metaanalysis functions 
source(file = './Meta_analysis_functions.R')
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ReactomePA")

##retrieve data
library(readxl)
library(openxlsx)
library(biomaRt)
library(dplyr)
library(ReactomePA)

read_excel_allsheets <- function(filename, tibble = FALSE) {
  # I prefer straight data.frames
  # but if you like tidyverse tibbles (the default with read_excel)
  # then just pass tibble = TRUE
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}


directory<-'./ENM_public_data/'
datasets_type_folders<-c('I','II','III','IV')
FC_pVal_tables<-list()
for(i in 1:length(datasets_type_folders)){
  path<-paste0(directory,datasets_type_folders[i])
  files<-list.files(path)
  FC_pVal_table<-list()
  for(j in 1:length(files)){
    GSE<-gsub("\\_.*","",files[j])
    file<-paste0(path,'/',files[j],'/',GSE,'_Unfiltered_DEG.xlsx')
    if(files[j]=="EMTAB6396_A549_files"){
      file<-'./ENM_public_data/IV/EMTAB6396_A549_files/EMTAB6396_A549_Unfiltered_DEG.xlsx'
      GSE<-"EMTAB6396_A549" ###correct
    } 
    if(files[j]=="EMTAB6396_Beas2B_files"){
      file<-'./ENM_public_data/IV/EMTAB6396_Beas2B_files/EMTAB6396_Beas2B_Unfiltered_DEG.xlsx'
      GSE<-"EMTAB6396_Beas2B"
    } 
    if(files[j]=="EMTAB6396_MRC9_files"){
      file<-'./ENM_public_data/IV/EMTAB6396_MRC9_files/EMTAB6396_MRC9_Unfiltered_DEG.xlsx'
      GSE<-"EMTAB6396_MRC9"
    } 
    if(files[j]=="EMTAB6396_THP1_files"){
      file<-'./ENM_public_data/IV/EMTAB6396_THP1_files/EMTAB6396_THP1_Unfiltered_DEG.xlsx'
      GSE<-"EMTAB6396_THP1"
    } 
    if(files[j]=="GSE19487_liver_files"){
      file<-'./ENM_public_data/IV/GSE19487_liver_files/GSE19487_liver_Unfiltered_DEG.xlsx'
      GSE<-"GSE19487_liver"
    } 
    if(files[j]=="GSE19487_lung_files"){
      file<-'./ENM_public_data/IV/GSE19487_lung_files/GSE19487_lung_Unfiltered_DEG.xlsx'
      GSE<-"GSE19487_lung"
    } 
    if(files[j]=="GSE99929_Jurkat_files"){
      file<-'./ENM_public_data/IV/GSE99929_Jurkat_files/GSE99929_Jurkat_Unfiltered_DEG.xlsx'
      GSE<-"GSE99929_Jurkat"
    } 
    if(files[j]=="GSE99929_THP1_files"){
      file<-'./ENM_public_data/IV/GSE99929_THP1_files/GSE99929_THP1_Unfiltered_DEG.xlsx'
      GSE<-"GSE99929_THP1"
    } 
    if(files[j]=="GSE39330_donor_files"){
      file<-'./ENM_public_data/III/GSE39330_donor_files/GSE39330_donor_Unfiltered_DEG.xlsx'
      GSE<-"GSE39330_donor"
    } 
    if(files[j]=="GSE39330_jurkat_files"){
      file<-'./ENM_public_data/III/GSE39330_jurkat_files/GSE39330_jurkat_Unfiltered_DEG.xlsx'
      GSE<-"GSE39330_jurkat"
    } 
    if(files[j]=="GSE84982_Caco2_files"){
      file<-'./ENM_public_data/III/GSE84982_Caco2_files/GSE84982_Caco2_Unfiltered_DEG.xlsx'
      GSE<-"GSE84982_Caco2"
    } 
    if(files[j]=="GSE84982_MCF7_files"){
      file<-"./ENM_public_data/III/GSE84982_MCF7_files/GSE84982_MCF7_Unfiltered_DEG.xlsx"
      GSE<-"GSE84982_MCF7"
    } 
    if(files[j]=="GSE148705_THP1_files"){
      file<-"./ENM_public_data/IV/GSE148705_THP1_files/GSE148705_THP1_Unfiltered_DEG.xlsx"
      GSE<-"GSE148705_THP1"
    } 
    if(files[j]=="GSE148705_BEAS2B_files"){
      file<-"./ENM_public_data/IV/GSE148705_BEAS2B_files/GSE148705_BEAS2B_Unfiltered_DEG.xlsx"
      GSE<-"GSE148705_BEAS2B"
    } 
    limma.outputs<-read_excel_allsheets(filename = file)
    limma.data<-list()
    for(k in 1:length(limma.outputs)){
      data<-limma.outputs[[k]]
      data<-data[,c('logFC','P.Value','adj.P.Val')]
      rownames(data)<-limma.outputs[[k]]$ID
      colnames(data)<-c(paste0(GSE,'_logFC'),paste0(GSE,'_P.Value'),paste0(GSE,'_adj.P.Val'))
      limma.data[[k]]<-data
      names(limma.data)[k]<-names(limma.outputs)[k]
    }
    FC_pVal_table[[j]]<-limma.data
    names(FC_pVal_table)[j]<-GSE
  }
  FC_pVal_tables[[i]]<-FC_pVal_table
  names(FC_pVal_tables)[i]<-datasets_type_folders[i]
}

#merge all the data in a single list
FC_pVal_complete_list<-list()
for (n in 1:length(FC_pVal_tables)) FC_pVal_complete_list<-append(FC_pVal_complete_list,FC_pVal_tables[[n]])
rm(FC_pVal_tables)

#retrieving selected GSEs and common genes
metaAnalysis_data<-read.table("./common_genes_dataset.txt")
selected_GSEs_complete<-gsub(pattern = '_logFC',replacement = '',x = colnames(metaAnalysis_data))
selected_GSEs_complete<-gsub(pattern = '_P.Value',replacement = '',x = selected_GSEs_complete)
selected_GSEs_complete<-gsub(pattern = '_adj.P.Val',replacement = '',x = selected_GSEs_complete)
selected_GSEs_complete<-unique(selected_GSEs_complete)
to_be_discarded<-c('GSE14452','GSE7010','GSE92899','GSE30213','GSE30214','GSE30215','GSE30180',
                  'GSE30200','GSE100500','GSE44294')
selected_GSEs<-setdiff(selected_GSEs_complete,to_be_discarded)
selected_genes<-rownames(metaAnalysis_data)
rm(selected_GSEs_complete)

##discard GSE which are not considered in our pool
GSE_list_tmp<-list()
for(i in 1:length(FC_pVal_complete_list)){
  GSE<-names(FC_pVal_complete_list)[i]
  if((GSE %in% selected_GSEs)){
    GSE_list_tmp[[i]]<-FC_pVal_complete_list[[i]]
    names(GSE_list_tmp)[i]<-names(FC_pVal_complete_list)[i]
    print(GSE)
  } 
}
GSE_list<-GSE_list_tmp[-which(sapply(GSE_list_tmp, is.null))]
rm(GSE_list_tmp)

##convert mouse ensemblID to human
mouse = useMart("ensembl", dataset="mmusculus_gene_ensembl") 
human =  useMart("ensembl", dataset="hsapiens_gene_ensembl")
GSE_list_converted<-list()
for(j in 1:length(GSE_list)){
  if(grepl('ENSMUSG',rownames(GSE_list[[j]][[1]])[1])){
    genes<-rownames(GSE_list[[j]][[1]])
    mouse_human<-getLDS(attributes=c("ensembl_gene_id"),
                        filters="ensembl_gene_id", values=genes, mart=mouse,
                        attributesL=c("ensembl_gene_id"), martL=human, uniqueRows = TRUE, bmHeader=FALSE)
    colnames(mouse_human)<-c('mouse_ensembl_ID','human_ensembl_ID')
    mouse_human<-mouse_human[!duplicated(mouse_human$mouse_ensembl_ID), ]
    mouse_human<-mouse_human[!duplicated(mouse_human$human_ensembl_ID), ]
    mouse_dataset<-list()
    for(k in 1:length(GSE_list[[j]])){
      mouse_dataset[[k]]<-GSE_list[[j]][[k]]
      mouse_dataset[[k]]<-mouse_dataset[[k]][mouse_human$mouse_ensembl_ID,]
      rownames(mouse_dataset[[k]])<-mouse_human$human_ensembl_ID
      names(mouse_dataset)[k]<-names(GSE_list[[j]])[k]
    }
    GSE_list_converted[[j]]<-mouse_dataset
    names(GSE_list_converted)[j]<-names(GSE_list)[j]
  }else{
    GSE_list_converted[[j]]<-GSE_list[[j]]
    names(GSE_list_converted[[j]])<-names(GSE_list[[j]])
    names(GSE_list_converted)[j]<-names(GSE_list)[j]
  }
}


##select only intersection genes
GSE_list_converted_cut<-list()
for(i in 1:length(GSE_list_converted)){
  pairwise_comparisons<-GSE_list_converted[[i]]
  for(j in 1:length(pairwise_comparisons)){
    pairwise_comparisons[[j]]<-pairwise_comparisons[[j]][selected_genes,]
    name<-paste(names(GSE_list_converted)[i],names(GSE_list_converted[[i]])[j],sep = '_')
    names(pairwise_comparisons)[j]<-name
  }
  GSE_list_converted_cut[[i]]<-pairwise_comparisons
  names(GSE_list_converted_cut)[i]<-names(GSE_list_converted)[i]
}

experiments<-list()
for(i in 1:length(GSE_list_converted_cut)) experiments<-append(experiments,GSE_list_converted_cut[[i]])

metanalysis_dataframe<-data.frame()
for(j in 1:length(experiments)){
  experiment<-experiments[[j]]
  experiment_name<-names(experiments)[j]
  colname_1<-paste(experiment_name,'logFC',sep ='_')
  colname_2<-paste(experiment_name,'P.Value',sep ='_')
  colname_3<-paste(experiment_name,'adj.P.Val',sep ='_')
  colnames(experiment)<-c(colname_1,colname_2,colname_3)
  if(plyr::empty(metanalysis_dataframe)) metanalysis_dataframe<-experiment
  else metanalysis_dataframe<-qpcR:::cbind.na(metanalysis_dataframe,as.data.frame(experiment))
}
write.table(x = metanalysis_dataframe, file = './data_space.txt', sep ='\t')


##run the metanalysis
effect_size_rank<-calc_effect_size_rank(metanalysis_dataframe)
fisher_based_pvalues<-calc_pvalue_based_rank(metanalysis_dataframe)
ranks_based_pvalues<-calc_rank_base_rank_subsets(metanalysis_dataframe)
data<-list()
data[[1]]<-rownames(effect_size_rank)
data[[2]]<-rownames(fisher_based_pvalues)
data[[3]]<-rownames(ranks_based_pvalues)
names(data)<-c('Effect_size','Fisher_test','Rank_based')
outputBorda<-Borda(data)
final_ranked_genes_mean<-as.data.frame(outputBorda$TopK$mean)
final_ranked_genes_median<-as.data.frame(outputBorda$TopK$median)
write.table(x = final_ranked_genes_mean,file = "./final_ranked_genes_mean.txt")
write.table(x = final_ranked_genes_mean,file = "./final_ranked_genes_mean.csv")
write.table(x = final_ranked_genes_median,file = "./final_ranked_genes_median.txt")


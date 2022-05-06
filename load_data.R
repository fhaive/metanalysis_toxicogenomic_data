load_foldchange_dataset = function(){

  metanalysis_dataframe<-read.table("/Users/meanser/onedrive/correlation_between_gene_and_molecular_descriptors/data/data_space_original.txt")
  expression_matrix<-metanalysis_dataframe[,which(grepl('_logFC',colnames(metanalysis_dataframe)))]
  colnames(expression_matrix)<-gsub(x = colnames(expression_matrix),pattern = '_logFC',replacement = '')
  colnames(expression_matrix)<-gsub(pattern = '.control',replacement = '-control',x = colnames(expression_matrix))
  colnames(expression_matrix)<-gsub(pattern = '.cont',replacement = '-cont',x = colnames(expression_matrix))
  colnames(expression_matrix)<-gsub(pattern = '.Ctrl',replacement = '-Ctrl',x = colnames(expression_matrix))
  colnames(expression_matrix)[30]<-"GSE122197_CuO_COOH_PBS_high-control_PB"
  colnames(expression_matrix)[18]<-'GSE122197_CuO_COOH_OVA_high-control_OVA'
  colnames(expression_matrix)[545]<-"GSE17676_Nanotubes30nm-flat_titanium"
  colnames(expression_matrix)[546]<-"GSE17676_Nanotubes100nm-flat_titanium"
  colnames(expression_matrix)[541]<-"GSE157266_CdTE_PEG_D-Ctrl_D"
  colnames(expression_matrix)[581]<-"GSE99929_Jurkat_Jurkat_GO-Jurkat_control"
  colnames(expression_matrix)[582]<-"GSE99929_Jurkat_Jurkat_GONH2-Jurkat_control"
  colnames(expression_matrix)[583]<-"GSE99929_THP1_THP1_GO-THP1_control"
  colnames(expression_matrix)[584]<-"GSE99929_THP1_THP1_GONH2-THP1_control"
  
  expression_matrix = as.matrix(expression_matrix)
  expression_matrix = t(apply(expression_matrix,1,DescTools::Winsorize,probs = c(0.01, 0.99)))
  
  return(expression_matrix) 
  
}


#### LOADING Descriptors dataset

read_molecular_descriptors_datasets = function(min_unique_number = 10){
  library(readxl)
  library(ggplot2)
  library(dplyr)
  library(hrbrthemes)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(rcompanion)
  
  cube_root_f = function(x){
    sign(x) * abs(x)^(1/3)
  }
  

  load("data_input/mds_num.RData")
  
  numeric_descriptors = mds_num
  rownames(numeric_descriptors) = numeric_descriptors[,1]
  numeric_descriptors = numeric_descriptors[,-1]
  
  options(warn = 2)
  
  for(i in 1:ncol(numeric_descriptors)){
    for(j in 1:nrow(numeric_descriptors)){
      numeric_descriptors[j,i]= as.numeric(numeric_descriptors[j,i])
    }
  }
  
  load("data_input/mds_cat.RData")
  
  text_descriptors = mds_cat
  rownames(text_descriptors) = text_descriptors$Sample
  text_descriptors = text_descriptors[,-1]
  
  for(i in 1:ncol(text_descriptors)) text_descriptors[,i] = as.factor(text_descriptors[,i])
  
  idx_samples = which(text_descriptors$`Core material` %in% c("UFP", "CoCl2"))
  
  # <3 samples
  idx_samples = c(idx_samples, which(text_descriptors$Chemistry 
                                     %in% c("aluminum (948), copper (150), iron (879), lead (181), silicon (1098) (nanograms per milligrams)",
                                            "Cobalt chloride")))
  
  # only one sample
  idx_samples = c(idx_samples, which(text_descriptors$Material %in% c("pollutant")))
  
  if(length(idx_samples)>0) text_descriptors = text_descriptors[-idx_samples,]
  
  # ensure same order
  numeric_descriptors = numeric_descriptors[rownames(text_descriptors),]
  
  desc = apply(numeric_descriptors, 2,function(colonna) length(table(colonna)))
  index_correlation = which(desc > min_unique_number)
  index_anova = which(desc <= min_unique_number)
  
  # Identify continous descriptors (with more than 5 different numbers, with less than 5 different numbers)
  continous_descriptors = numeric_descriptors[,index_correlation]
  continous_descriptors = cube_root_f(continous_descriptors)
  factor_descriptors = numeric_descriptors[,index_anova]
  factor_descriptors = cbind(factor_descriptors, text_descriptors[rownames(factor_descriptors),])
  
  md_datasets = list("continous_descriptors"=continous_descriptors,"factor_descriptors"=factor_descriptors)
  
  return(md_datasets)
}


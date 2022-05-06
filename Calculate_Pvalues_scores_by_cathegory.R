##select genes from groupings
library(readxl)
library(openxlsx)
library(vcd)
rank_file <- read.table("./final_ranked_genes_mean.txt")
top_rank<-as.vector(rank_file[0:1872,])


meta_analysis_groups<-openxlsx::read.xlsx('./Meta-analysis samples grouping.xlsx')
meta_analysis_groups$`Aspect.Ratio.(length/diameter)`<-NULL
meta_analysis_groups$`Dose.(repeated.or.single)`<-NULL
meta_analysis_groups$other.notes<-NULL
meta_analysis_groups$Additional.control<-NULL

metanalysis_dataframe<-read.table("./data_space.txt")

##prepare the contingency table 
biological_systems<-table(meta_analysis_groups$Biological.system)
organism<-table(meta_analysis_groups$Organism)
Biological_system_grouping<-table(meta_analysis_groups$Biological.system_grouping)
invitro_invivo<-table(meta_analysis_groups$`in.vivo-in.vitro`)
Core_material<-table(meta_analysis_groups$Core.material)
Functionalisation<-table(meta_analysis_groups$Functionalisation)
Chemistry<-table(meta_analysis_groups$Chemistry)
Metal<-table(meta_analysis_groups$Metal)
Oxide<-table(meta_analysis_groups$Oxide)
Geometry<-table(meta_analysis_groups$Geometry)
Material<-table(meta_analysis_groups$Material)
Time_h<-table(meta_analysis_groups$Time_h)
Time_period<-table(meta_analysis_groups$Time_period)

#Create a vector containing all the possible categories 
group_names<-c()
for(j in 2:length(colnames(meta_analysis_groups))){
  column<-meta_analysis_groups[,j]
  column_group<-paste(colnames(meta_analysis_groups)[j],column,sep = '__')
  group_names<-append(group_names,unique(column_group))
  
}


#Create a list in which each group contains the samples which correspond to that group
group_samples<-list()
for(i in 1:length(group_names)){
  group<-gsub(pattern = '\\__.*',x = group_names[i],replacement = '')
  type<-gsub(pattern = '.*\\__',x = group_names[i],replacement = '')
  column_index<-which(colnames(meta_analysis_groups)==group)
  row_index<-which(meta_analysis_groups[,column_index] == type)
  samples<-meta_analysis_groups[row_index,1]
  group_samples[[i]]<-samples
  names(group_samples)[i]<-group_names[i]
}
##Create a table in which each group has its corresponding number of samples inside
numb_samples_group<-matrix(ncol = length(group_samples), nrow=1)
colnames(numb_samples_group)<-names(group_samples)
for(g in 1:length(group_samples)) numb_samples_group[1,g]<-length(group_samples[[g]])
numb_samples_group<-t(numb_samples_group)
colnames(numb_samples_group)<-'Number of samples'

##Get ignificant genes for each sample pvalue
pval_metanalysis_dataframe<- metanalysis_dataframe[,which(grepl('_P.Val',colnames(metanalysis_dataframe)))]
is_significant_pval<-pval_metanalysis_dataframe
for(i in 1:dim(is_significant_pval)[2]){
  for(e in 1:dim(is_significant_pval)[1]){
    if(is_significant_pval[e,i] <= 0.05) is_significant_pval[e,i]<- 'Yes'
    else is_significant_pval[e,i]<- 'No'
  }
}
colnames(is_significant_pval)<-gsub(x = colnames(is_significant_pval),pattern = '_P.Value',replacement = '')
colnames(is_significant_pval)<-gsub(pattern = '.control',replacement = '-control',x = colnames(is_significant_pval))
colnames(is_significant_pval)<-gsub(pattern = '.cont',replacement = '-cont',x = colnames(is_significant_pval))
colnames(is_significant_pval)<-gsub(pattern = '.Ctrl',replacement = '-Ctrl',x = colnames(is_significant_pval))
colnames(is_significant_pval)[30]<-"GSE122197_CuO_COOH_PBS_high-control_PB"
colnames(is_significant_pval)[18]<-'GSE122197_CuO_COOH_OVA_high-control_OVA'
colnames(is_significant_pval)[545]<-"GSE17676_Nanotubes30nm-flat_titanium"
colnames(is_significant_pval)[546]<-"GSE17676_Nanotubes100nm-flat_titanium"
colnames(is_significant_pval)[541]<-"GSE157266_CdTE_PEG_D-Ctrl_D"
colnames(is_significant_pval)[581]<-"GSE99929_Jurkat_Jurkat_GO-Jurkat_control"
colnames(is_significant_pval)[582]<-"GSE99929_Jurkat_Jurkat_GONH2-Jurkat_control"
colnames(is_significant_pval)[583]<-"GSE99929_THP1_THP1_GO-THP1_control"
colnames(is_significant_pval)[584]<-"GSE99929_THP1_THP1_GONH2-THP1_control"


##calculate the score for each gene for each group

genes_score_bygroup<-matrix(nrow=dim(is_significant_pval)[1], ncol= length(group_samples))
for(g in 1:length(group_samples)){
  samples<-group_samples[[g]]
  n_samples<-length(samples)
  Pval_samples<-as.matrix(is_significant_pval[,samples])
  for(e in 1:dim(is_significant_pval)[1]){
    if(!is.null(dim(Pval_samples))){
      t<-table(Pval_samples[e,])
      if(is.na(t['Yes'])) t['Yes']<-0
    }
    else if(Pval_samples[e]=='Yes') t['Yes']<-1
    else t['Yes']<-0
    genes_score_bygroup[e,g]<-t['Yes']/n_samples
  }
}
colnames(genes_score_bygroup)<-names(group_samples)
rownames(genes_score_bygroup)<-rownames(is_significant_pval)

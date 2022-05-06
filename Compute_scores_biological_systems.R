##Compute scores for biological systems 
##needs variables computed in the R code Calculate Pvalues by cathegories (group names and is_significant)
##select genes from groupings
library(readxl)
library(openxlsx)
library(vcd)

rank_file <- read.table("./final_ranked_genes_mean.txt")
top_rank<-rank_file[0:1872,]
meta_analysis_groups<-openxlsx::read.xlsx('./Meta-analysis samples grouping.xlsx')
meta_analysis_groups$`Aspect.Ratio.(length/diameter)`<-NULL
meta_analysis_groups$`Dose.(repeated.or.single)`<-NULL
meta_analysis_groups$other.notes<-NULL
meta_analysis_groups$Additional.control<-NULL

###___________________________________________________________________________
##select biological system specific samples only
key_names<-c('in.vivo-in.vitro__in vitro', 'Biological.system_grouping__cancer cell line',
             'in.vivo-in.vitro__in vivo','Biological.system_grouping__primary cell line')
restricted_groups<-list()
for(k in 1:length(key_names)){
  restricted_group<-list()
  length(restricted_group)<-length(group_samples)
  biological_system<-key_names[k]
  for(i in 1:length(group_samples)){
    cut<-intersect(group_samples[[i]],group_samples[[biological_system]])
    restricted_group[[i]]<-cut
    names(restricted_group)[i]<-names(group_samples)[i]
    print(names(restricted_group)[i])
    
  }
  restricted_group_no0<-restricted_group[-which(sapply(restricted_group, function(x) length(x)==0))]
  restricted_groups[[k]]<-restricted_group_no0
  names(restricted_groups)[k]<-biological_system
}

#select the directory to store the generated files
path='./'
genes_score_pvalbased_by_groups<-list()
for(b in 1:length(restricted_groups)){
  restricted_groups_no0<-restricted_groups[[b]]
  genes_score_pvalbased<-matrix(nrow=dim(is_significant_pval)[1], ncol= length(restricted_groups_no0))
  for(g in 1:length(restricted_groups_no0)){
    samples<-restricted_groups_no0[[g]]
    n_samples<-length(samples)
    Pval_samples<-is_significant_pval[,which(colnames(is_significant_pval) %in% samples)] 
    for(e in 1:dim(is_significant_pval)[1]){
      if(!is.null(dim(Pval_samples))){
        t<-table(as.factor(Pval_samples[e,]))
        if(is.na(t['Yes'])) t['Yes']<-0
      }
      else if(Pval_samples[e]=='Yes') t['Yes']<-1
      else t['Yes']<-0
      genes_score_pvalbased[e,g]<-t['Yes']/n_samples
    }
  }
  colnames(genes_score_pvalbased)<-names(restricted_groups_no0)
  rownames(genes_score_pvalbased)<-rownames(is_significant_pval)
  genes_score_pvalbased_by_groups[[b]]<-genes_score_pvalbased
  names(genes_score_pvalbased_by_groups)[b]<-names(restricted_groups)[b]
  file_name<-paste0('single_gene_significance_score_restricted_groups',names(restricted_groups)[b],sep='_')
  file_exstension<-paste0(file_name,'.txt',sep='')
  path_complete<-paste0(path,file_exstension,sep='')
  write.table(x = genes_score_pvalbased, file = path_complete,sep = '\t')
}

##PLOT HEATMAP SMALLER IN VIVO  CLASSES
row_indeces<-top_rank
##TEMP
in_vivo_data_total<-as.data.frame(genes_score_pvalbased_by_groups[['in.vivo-in.vitro__in vivo']])
core_material_indeces<-which(grepl(x=colnames(in_vivo_data_total),pattern = 'Core.material'))
time_indeces<-c(154,155,156)
in_vivo_data<-as.matrix(in_vivo_data_total[row_indeces,time_indeces])

pdf(file = "./Heatmap_in_vivo_timeperiod.pdf", width=10)
out<-gplots::heatmap.2(in_vivo_data,
                        hclustfun = function(x) hclust(x,method = 'ward'),
                              distfun = function(x) dist(x,method = 'euclidean'),
                              trace = 'none', margins = c(10,5),col=brewer.pal(9,"Blues"),
                              srtCol=45,xlab="", ylab= "",
                              main='Time period in vivo top rank genes scores')

dev.off()

clustered_in_vivo<-rownames(in_vivo_data)[out$rowInd]
long_term_genes<-c()
for(i in 1:length(clustered_in_vivo)){
  index<-which(rownames(in_vivo_data_total)==clustered_in_vivo[i])
  value<-in_vivo_data_total[index,156]
  long_term_genes[i]<-value
  names(long_term_genes)[i]<-clustered_in_vivo[i]
}

first_element_cluster<-which(names(long_term_genes)=='ENSG00000166851')
last_element_cluster<-which(names(long_term_genes)=='ENSG00000142634')

in_vivo_long_term_genes<-names(long_term_genes)[first_element_cluster:last_element_cluster]
write.table(x = in_vivo_long_term_genes, file =  './in_vivo_long_term.txt',sep='\t')
write.table(x = in_vivo_long_term_genes, file =  './in_vivo_long_term_july.txt',sep='\t')






##PLOT HEATMAP SMALLER IN VITRO  CLASSES
row_indeces<-top_rank
##TEMP
in_vitro_data_total<-as.data.frame(genes_score_pvalbased_by_groups[['in.vivo-in.vitro__in vitro']])
core_material_indeces<-which(grepl(x=colnames(in_vitro_data_total),pattern = 'Core.material'))
time_indeces<-c(187,188,189)
in_vitro_data<-as.matrix(in_vitro_data_total[row_indeces,time_indeces])
pdf(file = "./Heatmap_in_vitro_timeperiod.pdf", width=10)
out<-gplots::heatmap.2(in_vitro_data,
                        hclustfun = function(x) hclust(x,method = 'ward'),
                        distfun = function(x) dist(x,method = 'euclidean'),
                        trace = 'none', margins = c(10,5),col=brewer.pal(9,"Blues"),
                        srtCol=45,xlab="", ylab= "",
                        main='Time period in vitro top rank genes scores')
dev.off()

clustered_in_vitro<-rownames(in_vitro_data)[out$rowInd]
long_term_genes<-c()
for(i in 1:length(clustered_in_vitro)){
  index<-which(rownames(in_vitro_data_total)==clustered_in_vitro[i])
  value<-in_vitro_data_total[index,189]
  long_term_genes[i]<-value
  names(long_term_genes)[i]<-clustered_in_vitro[i]
}
first_index_cluster<-which(names(long_term_genes)=='ENSG00000149476')
last_index_cluster<-which(names(long_term_genes)=='ENSG00000163875')
in_vitro_long_term_genes<-names(long_term_genes)[first_index_cluster:last_index_cluster]
write.table(x = in_vitro_long_term_genes, file =  './in_vitro_long_term_july.txt',sep='\t')


always_long_term<-intersect(in_vitro_long_term_genes,in_vivo_long_term_genes)
write.table(x = always_long_term, file =  './always_long_term_july.txt',sep='\t')

library(gplots)
long_term_genes<-list(in_vivo=in_vivo_long_term_genes,in_vitro=in_vitro_long_term_genes)
v.table <- venn(long_term_genes)
print(v.table)

write.table(x = in_vitro_long_term_genes, file =  './in_vitro_long_term.txt',sep='\t')

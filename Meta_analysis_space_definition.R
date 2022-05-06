###################################### meta analysis data space definition ##################################################
####################### selection of the shared genes across platforms and building the data space ###############
install.packages('xlsx')
library(openxlsx)
library(biomaRt)
library(dplyr)
library(xlsx)

#retrieve limma output FC and pvalues from Zenodo data

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
    limma.output<-openxlsx::read.xlsx(file)
    data<-limma.output[,c('logFC','P.Value','adj.P.Val')]
    rownames(data)<-limma.output$ID
    colnames(data)<-c(paste0(GSE,'_logFC'),paste0(GSE,'_P.Value'),paste0(GSE,'_adj.P.Val'))
    FC_pVal_table[[j]]<-data
    names(FC_pVal_table)[j]<-GSE
  }
  FC_pVal_tables[[i]]<-FC_pVal_table
  names(FC_pVal_tables)[i]<-datasets_type_folders[i]
}

#merge all the data in a single list
FC_pVal_complete_list<-list()
for (n in 1:length(FC_pVal_tables)) FC_pVal_complete_list<-append(FC_pVal_complete_list,FC_pVal_tables[[n]])

#--------------------------------------------Data exploration-------------------------------------------------
#identify datasets based on old platforms that share little genes with the other more recent platforms

upset_plot_input_data<-list()

for(dataset in 1:length(FC_pVal_complete_list)){
  upset_plot_input_data[[dataset]]<-as.data.frame(rep.int(1,length(rownames(FC_pVal_complete_list[[dataset]]))))
  colnames(upset_plot_input_data[[dataset]])<-names(FC_pVal_complete_list)[dataset]
  rownames(upset_plot_input_data[[dataset]])<-rownames(FC_pVal_complete_list[[dataset]])
  names(upset_plot_input_data)[dataset]<-colnames(upset_plot_input_data[[dataset]])
  if(!grepl('ENSG',rownames(FC_pVal_complete_list[[dataset]])[1])) upset_plot_input_data[[dataset]]<-NULL
}
small_datasets<-c('GSE14452','GSE7010','GSE92899','GSE30213','GSE30214','GSE30215','GSE30180',
                  'GSE30200')
upset_plot_input_data[small_datasets]<-NULL
upset_plot_input<-upset_plot_input_data[-which(sapply(upset_plot_input_data, is.null))]


#___________________________plots array length___________________________________________________________
length_array<-vector('logical', length = length(upset_plot_input))
for(i in 1:length(upset_plot_input)){
  length_array[i]<-dim(upset_plot_input[[i]])[1]
}
y_values<-as.vector(length_array)
hist(y_values,plot = TRUE, breaks = 100, xlim = c(5000,35000), xlab = "Array length")
x_labels<-names(upset_plot_input_data)[which(names(upset_plot_input_data)!="")]
all_genes<-c()
for(i in 1:length(upset_plot_input)){
  if(length(all_genes)==0) all_genes<-rownames(upset_plot_input[[i]])
  else all_genes<-union(all_genes,rownames(upset_plot_input[[i]]))
}
all_genes_length<-length(unique(all_genes))
entry_spaces<-grepl(pattern = " ",x = all_genes)
for(i in 1:length(entry_spaces)) if(entry_spaces[i]==TRUE) print(entry_spaces[i])
entry_dot<-all_genes[which(grepl(pattern = ".",x = all_genes)==TRUE)]
for(i in 1:length(entry_dot)) if(entry_dot[i]==TRUE) print(entry_dot[i])
entry_dash<-all_genes[which(grepl(pattern = "_",x = all_genes)==TRUE)]

dataframe_input<-data.frame()
for(i in 1:length(upset_plot_input)){
  if(plyr::empty(dataframe_input)) dataframe_input<-upset_plot_input[[i]]
  else dataframe_input<-qpcR:::cbind.na(dataframe_input,as.data.frame(upset_plot_input[[i]]))
}
dataframe_input[is.na(dataframe_input)] <- 0
UpSetR::upset(dataframe_input,nintersects = NA, nsets = ncol(dataframe_input))

#heatmap 
common_genes_heatmap<-matrix(nrow = length(upset_plot_input),ncol=length(upset_plot_input))
for(j in 1:length(upset_plot_input)){
  row_dataset<-rownames(upset_plot_input[[j]])
  for(m in 1:length(upset_plot_input)){
    col_dataset<-rownames(upset_plot_input[[m]])
    common_genes_heatmap[j,m]<-length(intersect(row_dataset,col_dataset))
  }
}
common_genes_heatmap<-as.data.frame(common_genes_heatmap)
rownames(common_genes_heatmap)<-names(upset_plot_input)
colnames(common_genes_heatmap)<-rownames(common_genes_heatmap)
h<-heatmap(as.matrix(common_genes_heatmap))
indeces<-h$colInd[1:20]
cluster1<-upset_plot_input[indeces]
common_cluster1<-c()
for(j in 1:length(cluster1)){
  genes<-rownames(cluster1[[j]])
  if(length(common_cluster1)==0) common_cluster1<-genes
  else common_cluster1<-intersect(common_cluster1,genes)
  print(length(common_cluster1))
}
#1906
cluster2<-upset_plot_input[-indeces]
common_cluster2<-c()
for(j in 1:length(cluster2)){
  genes<-rownames(cluster2[[j]])
  if(length(common_cluster2)==0) common_cluster2<-genes
  else common_cluster2<-intersect(common_cluster2,genes)
  print(length(common_cluster2))
}
#9141
#platform intersections
metadata<-openxlsx::read.xlsx("/Users/giusy/Downloads/Online-only_Table1.xlsx")
affy_datasets<-metadata[which(metadata$platform=="Affymetrix" & metadata$organism=="Homo sapiens"),]
affy_datasets<-affy_datasets$ID
agilent_datasets<-metadata[which(metadata$platform=="Agilent" & metadata$organism=="Homo sapiens"),]
agilent_datasets<-agilent_datasets$ID
agilent_datasets[10]<-"EMTAB6396_A549r"
illumina_datasets<-metadata[which(metadata$platform=="Illumina" & metadata$organism=="Homo sapiens"),]
illumina_datasets<-illumina_datasets$ID
affy_common_data<-upset_plot_input_data[affy_datasets]
t<-upset_plot_input_data[agilent_datasets]
#agilent_common_data<-list()
#for(n in 1:length(agilent_datasets)) agilent_common_data[n]<-upset_plot_input_data[agilent_datasets[n]]
illumina_common_data<-upset_plot_input_data[illumina_datasets]

#calculate affy genes
common_affy<-list()
for(j in 1:length(affy_common_data)){
  genes<-rownames(affy_common_data[[j]])
  if(length(common_affy)==0) common_affy<-genes
  else common_affy<-intersect(common_affy,genes)
  print(length(common_affy))
}
#calculate affy genes
agilent_common_data<-t
common_agilent<-list()
for(j in 1:length(agilent_common_data)){
  genes<-rownames(agilent_common_data[[j]])
  if(length(common_agilent)==0) common_agilent<-genes
  else common_agilent<-intersect(common_agilent,genes)
  print(length(common_agilent))
}
#calculate illumina genes
common_illumina<-list()
for(j in 1:length(illumina_common_data)){
  genes<-rownames(illumina_common_data[[j]])
  if(length(common_illumina)==0) common_illumina<-genes
  else common_illumina<-intersect(common_illumina,genes)
  print(length(common_illumina))
}

affy_illumina<-intersect(common_illumina,common_affy) #2577
affy_agilent<-intersect(common_agilent,common_affy) #3241
agilent_illumina<-intersect(common_agilent,common_illumina) #4456
complete<-intersect(affy_agilent,common_illumina)#1798
venn_data_platforms= list("Illumina"= common_illumina,
                "Agilent"=common_agilent,
                "Affymetrix"= common_affy)
venn.plot.platforms<-VennDiagram::venn.diagram(x = venn_data_platforms, filename = NULL)


#unique genes 
all_genes_platforms<-union(common_affy,common_agilent)
all_genes_platforms<-union(common_illumina,all_genes_platforms)
affy_only<-setdiff(all_genes_platforms,common_agilent) #2658
affy_only<-setdiff(affy_only,common_illumina) #2658

agilent_only<-setdiff(all_genes_platforms,common_affy) #2658
agilent_only<-setdiff(agilent_only,common_illumina) #2658

illumina_only<-setdiff(all_genes_platforms,common_agilent) #2658
illumina_only<-setdiff(illumina_only,common_affy) #2658

length(illumina_only)
write(x=illumina_only,file = './illumina_unique_genes.xls')
length(agilent_only)
write(x=agilent_only,file = './agilent_unique_genes.xls')
length(affy_only)
write(x=affy_only,file = './affymetrix_unique_genes.xls')
####_____________________________________________________________________________________________________________________________________----

#calculate common genes
common_genes_human<-list()
for(j in 1:length(upset_plot_input)){
    genes<-rownames(upset_plot_input[[j]])
    if(length(common_genes_human)==0) common_genes_human<-genes
    else common_genes_human<-intersect(common_genes_human,genes)
    print(length(intersect(common_genes_human,genes)))
    print(j)
  }
print(length(common_genes_human))


#retrieve from the original tables only the information about the common genes,
human_datasets<-list()
for(j in 1:length(FC_pVal_complete_list)){
  if(grepl('ENSG',rownames(FC_pVal_complete_list[[j]])[1])){
    human_dataset<-FC_pVal_complete_list[[j]]
    human_dataset<-human_dataset[common_genes_human,]
    rownames(human_dataset)<-common_genes_human
    human_datasets[[j]]<-human_dataset
    names(human_datasets)[j]<-names(FC_pVal_complete_list)[j]
  }
  else human_datasets[[j]]<-NULL
}
human_datasets_NN<-human_datasets[-which(sapply(human_datasets, is.null))]


#merge all datasets in a single dataframe
human_dataframe<-data.frame()
for(i in 1:length(human_datasets_NN)){
  if(plyr::empty(human_dataframe)) human_dataframe<-human_datasets_NN[[i]]
  else human_dataframe<-qpcR:::cbind.na(human_dataframe,as.data.frame(human_datasets_NN[[i]]))
}



#--------------------------------common mouse genes
#similarly to human, old platforms have been discarded

upset_plot_input_data_mouse<-list()

for(dataset in 1:length(FC_pVal_complete_list)){
  upset_plot_input_data_mouse[[dataset]]<-as.data.frame(rep.int(1,length(rownames(FC_pVal_complete_list[[dataset]]))))
  colnames(upset_plot_input_data_mouse[[dataset]])<-names(FC_pVal_complete_list)[dataset]
  rownames(upset_plot_input_data_mouse[[dataset]])<-rownames(FC_pVal_complete_list[[dataset]])
  names(upset_plot_input_data_mouse)[dataset]<-colnames(upset_plot_input_data_mouse[[dataset]])
  if(!grepl('ENSMUSG',rownames(FC_pVal_complete_list[[dataset]])[1])) upset_plot_input_data_mouse[[dataset]]<-NULL
}
old_platforms<-c('GSE100500','GSE44294')
upset_plot_input_data_mouse[old_platforms]<-NULL
upset_plot_input_mouse<-upset_plot_input_data_mouse[-which(sapply(upset_plot_input_data_mouse, is.null))]


#comupute the genes shared across all the mouse datasets
common_genes_mouse<-list()
mouse_GSE<-c()
mouse_ind<-c()
complete_list_mouse_pre<-FC_pVal_complete_list
old_platforms<-c('GSE100500','GSE44294','GSE29042','GSE35193','GSE19487_liver','GSE19487_lung')
complete_list_mouse_pre[old_platforms]<-NULL
complete_list_mouse<-complete_list_mouse_pre[-which(sapply(complete_list_mouse_pre, is.null))]

for(j in 1:length(upset_plot_input_mouse)){
  if(grepl('ENSMUSG',rownames(upset_plot_input_mouse[[j]])[1])){
    print(j)
    print(rownames(upset_plot_input_mouse[[j]])[1])
    genes<-rownames(upset_plot_input_mouse[[j]])
    if(length(common_genes_mouse)==0) common_genes_mouse<-genes
    else common_genes_mouse<-intersect(common_genes_mouse,genes)
    mouse_GSE[j]<-names(upset_plot_input_mouse)[j]
  }
  print(length(common_genes_mouse))
}



#convert shared mouse genes to human

mouse_GSE<-unique(mouse_GSE)
mouse_GSE<-mouse_GSE[!is.na(mouse_GSE)]
#convert
human = useMart("ensembl", dataset="hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset="mmusculus_gene_ensembl") 
mouse_human<-getLDS(attributes=c("ensembl_gene_id"),
                    filters="ensembl_gene_id", values=common_genes_mouse, mart=mouse,
                    attributesL=c("ensembl_gene_id"), martL=human, uniqueRows = TRUE, bmHeader=FALSE)
colnames(mouse_human)<-c('mouse_ensembl_ID','human_ensembl_ID')
mouse_human<-mouse_human[!duplicated(mouse_human$mouse_ensembl_ID), ]
mouse_human<-mouse_human[!duplicated(mouse_human$human_ensembl_ID), ]


#retrieve from the original tables only the information about the common genes,
#and convert the rownames to the ENSG IDs
common_genes_mouse_converted<-mouse_human$human_ensembl_ID
mouse_datasets<-list()
for(j in 1:length(FC_pVal_complete_list)){
  if(grepl('ENSMUSG',rownames(FC_pVal_complete_list[[j]])[1])){
    mouse_dataset<-FC_pVal_complete_list[[j]]
    mouse_dataset<-mouse_dataset[mouse_human$mouse_ensembl_ID,]
    rownames(mouse_dataset)<-mouse_human$human_ensembl_ID
  }
  else mouse_dataset<-data.frame()
  mouse_datasets[[j]]<-mouse_dataset
  names(mouse_datasets)[j]<-names(FC_pVal_complete_list)[j]
}

#merge all datasets in a single dataframe
mouse_dataframe<-data.frame()
for(n in 1:length(mouse_datasets)){
  if(plyr::empty(mouse_dataframe)){
    mouse_dataframe<-mouse_datasets[[n]]
  } 
  else if(length(mouse_datasets[[n]])==0){}
  else{
    last_index<-ncol(mouse_dataframe)+1
    mouse_dataframe<-cbind(mouse_dataframe,mouse_datasets[[n]])
    new_final_index<-ncol(mouse_dataframe)
    colnames(mouse_dataframe)[last_index:new_final_index]<-colnames(mouse_datasets[[n]])
  }
}


###_______________________ merge all datasets ____________________________
mouse_genes<-rownames(mouse_dataframe)
human_genes<-rownames(human_dataframe)
length(intersect(mouse_genes,human_genes))

venn_data= list("Mouse"= mouse_genes,
            "Human"= human_genes)
venn.plot<-VennDiagram::venn.diagram(x = venn_data, filename = NULL)

v.table<-gplots::venn(venn_data)

final_common_genes<-intersect(mouse_genes,human_genes)

common_mouse<-mouse_dataframe[final_common_genes,]
common_human<-human_dataframe[final_common_genes,]
common_genes_dataset<-common_human
common_genes_dataset<-qpcR:::cbind.na(common_genes_dataset,common_mouse)

write.table(x = common_genes_dataset,file = "./common_genes_dataset.txt")


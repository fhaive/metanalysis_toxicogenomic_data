library(tidyverse)
library(clusterProfiler)
library(DOSE)
library(biomaRt)
library(dplyr)
#upload the meta-analysis rank of genes and the dataframe
metanalysis_dataframe <- read.table('./data_space.txt', sep ='\t')
mainfile <- read.table("./final_ranked_genes_mean.txt")

#convert the gene to gene symbols
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
mainfile_gene_symbol<-getBM(attributes=c("hgnc_symbol", "ensembl_gene_id"),
                    filters="ensembl_gene_id", values=mainfile$outputBorda.TopK.mean, mart=mart,
                    uniqueRows = TRUE, bmHeader=FALSE)

######GSEA needs the pvalue for each gene, the minimum value across experiments is selected
pvalue_table<-metanalysis_dataframe[,which(grepl('P.Value',colnames(metanalysis_dataframe)))]
expression_value_table<-metanalysis_dataframe[,which(grepl('_logFC',colnames(metanalysis_dataframe)))]
min_pvalues<-as.data.frame(apply(pvalue_table, 1, FUN=min))
colnames(min_pvalues)<-'min.Pvalue'
gene_list<-vector(mode = "logical", length = 3676)
for(i in 1:3676){
  print(i)
  ENSG<-as.character(mainfile[i,1])
  if(ENSG %in% rownames(min_pvalues)) gene_list[i]<-min_pvalues[which(rownames(min_pvalues)==ENSG),1]
  gene_symbol_index<-which(mainfile_gene_symbol$ensembl_gene_id==ENSG)
  names(gene_list)[i]<-mainfile_gene_symbol$hgnc_symbol[gene_symbol_index]
}
#### FGSEA for each database ####
GO <- "/data_input/Curateed gene sets from MSigDB/c5.go.v7.2.symbols.gmt"
REACTOME <- "/data_input/Curateed gene sets from MSigDB/c2.cp.reactome.v7.2.symbols.gmt"
KEGG <- "/data_input/Curateed gene sets from MSigDB/c2.cp.kegg.v7.2.symbols.gmt"
WIKI <- "/data_input/Curateed gene sets from MSigDB/c2.cp.wikipathways.v7.2.symbols.gmt"
MSIGDB <- "/data_input/Curateed gene sets from MSigDB/msigdb.v7.2.symbols.gmt"
TRANSCRIPTION_FACTORS<-"/data_input/Curateed gene sets from MSigDB/c3.tft.v7.2.symbols.gmt"
myGO <- fgsea::gmtPathways(GO)
myREACTOME <- fgsea::gmtPathways(REACTOME)
myKEGG <- fgsea::gmtPathways(KEGG)
myWIKI <- fgsea::gmtPathways(WIKI)
myMSIGDB <- fgsea::gmtPathways(MSIGDB)
myTS <- fgsea::gmtPathways(TRANSCRIPTION_FACTORS)

pathlength <- c()
for(i in 1:length(myMSIGDB)){
  pathlength <- c(pathlength, length(myMSIGDB[[i]]))
}

summary(pathlength)

fgResGO <- fgsea::fgsea(pathways = myGO, 
                      stats = gene_list,
                      minSize=10,
                      nperm=10000)
# library
library(ggridges)
library(ggplot2)

ordered_fgRes<-fgResGO[order(fgResGO$pval)]
expression_values<-as.data.frame(matrix(0, ncol = 3, nrow = 0))
colnames(expression_values)<-c('Pathway','Gene','Expression_value')
for(i in 1:21){
  enriching_genes<-ordered_fgRes$leadingEdge[[i]]
  indeces<-which(mainfile_gene_symbol$hgnc_symbol %in% enriching_genes)
  enriching_genes_ENS<-mainfile_gene_symbol$ensembl_gene_id[indeces]
  expression_values_mean<-rowMeans(expression_value_table[enriching_genes_ENS,],na.rm = TRUE)
  expression_value_df_tmp<-as.data.frame(matrix(0, ncol = 3, nrow = length(enriching_genes_ENS)))
  colnames(expression_value_df_tmp)<-c('Pathway','Gene','Expression_value')
  expression_value_df_tmp$Pathway<-replicate(nrow(expression_value_df_tmp),ordered_fgRes$pathway[i])
  expression_value_df_tmp$Gene<-enriching_genes_ENS
  expression_value_df_tmp$Expression_value<-expression_values_mean
  expression_values<-rbind(expression_values,expression_value_df_tmp)
}
# basic example
p<-ggplot(expression_values, aes(x = Expression_value, y = Pathway, fill = Pathway)) +
  stat_density_ridges()+
  theme_ridges() + 
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(1,1,1,1), "cm"))
p  

topPathwaysUp_GO <- fgResGO[ES > 0][head(order(pval), n=20), pathway]
topPathwaysDown_GO <- fgResGO[ES < 0][head(order(pval), n=20), pathway]
topPathways_GO <- c(topPathwaysUp_GO, topPathwaysDown_GO)
fgsea::plotGseaTable(myGO[topPathwaysUp_GO], gene_list, fgResGO, 
              gseaParam=0.05)


fgResREACTOME <- fgsea::fgsea(pathways = myREACTOME, 
                      stats = gene_list,
                      minSize=10,
                      nperm=10000)
ordered_fgRes<-fgResREACTOME[order(fgResREACTOME$pval)]
expression_values<-as.data.frame(matrix(0, ncol = 3, nrow = 0))
colnames(expression_values)<-c('Pathway','Gene','Expression_value')
for(i in 1:21){
  enriching_genes<-ordered_fgRes$leadingEdge[[i]]
  indeces<-which(mainfile_gene_symbol$hgnc_symbol %in% enriching_genes)
  enriching_genes_ENS<-mainfile_gene_symbol$ensembl_gene_id[indeces]
  expression_values_mean<-rowMeans(expression_value_table[enriching_genes_ENS,],na.rm = TRUE)
  expression_value_df_tmp<-as.data.frame(matrix(0, ncol = 3, nrow = length(enriching_genes_ENS)))
  colnames(expression_value_df_tmp)<-c('Pathway','Gene','Expression_value')
  expression_value_df_tmp$Pathway<-replicate(nrow(expression_value_df_tmp),substr(x = ordered_fgRes$pathway[i],0,20))
  expression_value_df_tmp$Gene<-enriching_genes_ENS
  expression_value_df_tmp$Expression_value<-expression_values_mean
  expression_values<-rbind(expression_values,expression_value_df_tmp)
}
# basic example
p<-ggplot(expression_values, aes(x = Expression_value, y = Pathway, fill = Pathway)) +
  stat_density_ridges()+
  theme_ridges() + 
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0.01,0.01,0.01,0.01), "cm"))
p  


topPathwaysUp_REACTOME <- fgResREACTOME[ES > 0][head(order(pval), n=20), pathway]
topPathwaysDown_REACTOME <- fgResREACTOME[ES < 0][head(order(pval), n=20), pathway]
topPathways_REACTOME <- c(topPathwaysUp_REACTOME, rev(topPathwaysDown_REACTOME))
fgsea::plotGseaTable(myREACTOME[topPathwaysUp_REACTOME], gene_list, fgResREACTOME, 
                     gseaParam=0.05)

fgResKEGG <- fgsea::fgsea(pathways = myKEGG, 
                      stats = gene_list,
                      minSize=10,
                      nperm=10000)
ordered_fgRes<-fgResKEGG[order(fgResKEGG$pval)]
expression_values<-as.data.frame(matrix(0, ncol = 3, nrow = 0))
colnames(expression_values)<-c('Pathway','Gene','Expression_value')
for(i in 1:21){
  enriching_genes<-ordered_fgRes$leadingEdge[[i]]
  indeces<-which(mainfile_gene_symbol$hgnc_symbol %in% enriching_genes)
  enriching_genes_ENS<-mainfile_gene_symbol$ensembl_gene_id[indeces]
  expression_values_mean<-rowMeans(expression_value_table[enriching_genes_ENS,],na.rm = TRUE)
  expression_value_df_tmp<-as.data.frame(matrix(0, ncol = 3, nrow = length(enriching_genes_ENS)))
  colnames(expression_value_df_tmp)<-c('Pathway','Gene','Expression_value')
  expression_value_df_tmp$Pathway<-replicate(nrow(expression_value_df_tmp),substr(x = ordered_fgRes$pathway[i],0,20))
  expression_value_df_tmp$Gene<-enriching_genes_ENS
  expression_value_df_tmp$Expression_value<-expression_values_mean
  expression_values<-rbind(expression_values,expression_value_df_tmp)
}
# basic example
p<-ggplot(expression_values, aes(x = Expression_value, y = Pathway, fill = Pathway)) +
  stat_density_ridges()+
  theme_ridges() + 
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0.01,0.01,0.01,0.01), "cm"))
p  


topPathwaysUp_KEGG <- fgResKEGG[ES > 0][head(order(pval), n=20), pathway]
topPathwaysDown_KEGG <- fgResKEGG[ES < 0][head(order(pval), n=20), pathway]
topPathways_KEGG <- c(topPathwaysUp_KEGG, rev(topPathwaysDown_KEGG))
fgsea::plotGseaTable(myKEGG[topPathwaysUp_KEGG], gene_list, fgResKEGG, 
                     gseaParam=0.05)

fgResWIKI <- fgsea::fgsea(pathways = myWIKI, 
                      stats = gene_list,
                      minSize=10,
                      nperm=10000)

topPathwaysUp_WIKI <- fgResWIKI[ES > 0][head(order(pval), n=20), pathway]
topPathwaysDown_WIKI <- fgResWIKI[ES < 0][head(order(pval), n=20), pathway]
topPathways_WIKI <- c(topPathwaysUp_WIKI, rev(topPathwaysDown_WIKI))
fgsea::plotGseaTable(myWIKI[topPathwaysUp_WIKI], gene_list, fgResWIKI, 
                     gseaParam=0.05)

fgResMSIGDB <- fgsea::fgsea(pathways = myMSIGDB, 
                      stats = gene_list,
                      minSize=10,
                      nperm=10000)

topPathwaysUp_MSIGDB <- fgResMSIGDB[ES > 0][head(order(pval), n=20), pathway]
topPathwaysDown_MSIGDB <- fgResMSIGDB[ES < 0][head(order(pval), n=20), pathway]
topPathways_MSIGDB <- c(topPathwaysUp_MSIGDB, rev(topPathwaysDown_MSIGDB))
fgsea::plotGseaTable(myMSIGDB[topPathwaysUp_MSIGDB], gene_list, fgResMSIGDB, 
                     gseaParam=0.05)


#########################compute the rank threshold based on the single GSEA results##############################

compute_gsea_thresh <- function(geneList, fgsea_res, myGO=myGO){
  gseaParam=1
  stats <- geneList
  fgsea_res <- fgsea_res[order(fgsea_res$pval, decreasing = FALSE),]
  max_vec <- c()
  sigpath <- which(fgsea_res$pval<0.05)
  for(i in 1:length(sigpath)){
  
    pathway <- myGO[[fgsea_res$pathway[i]]]
    rnk <- rank(-stats)
    ord <- order(rnk)
    statsAdj <- stats[ord]
    statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
    statsAdj <- statsAdj/max(abs(statsAdj))
    pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
    pathway <- sort(pathway)
    gseaRes <- fgsea::calcGseaStat(statsAdj, selectedStats = pathway, 
                          returnAllExtremes = TRUE)
    bottoms <- gseaRes$bottoms
    tops <- gseaRes$tops
    n <- length(statsAdj)
    xs <- as.vector(rbind(pathway - 1, pathway))
    ys <- as.vector(rbind(bottoms, tops))
    toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))

    max_vec <- c(max_vec, which(names(geneList) %in% names(gseaRes$tops)[which(gseaRes$tops==max(gseaRes$tops))]))
}

  return(max_vec)
}

max_vec_GO <- compute_gsea_thresh(geneList = gene_list, fgsea_res = fgResGO, myGO = myGO)
max_vec_REACTOME <- compute_gsea_thresh(geneList = gene_list, fgsea_res = fgResREACTOME, myGO = myREACTOME)
max_vec_KEGG <- compute_gsea_thresh(geneList = gene_list, fgsea_res = fgResKEGG, myGO = myKEGG)
max_vec_WIKI <- compute_gsea_thresh(geneList = gene_list, fgsea_res = fgResWIKI, myGO = myWIKI)
max_vec_MSIGDB <- compute_gsea_thresh(geneList = gene_list, fgsea_res = fgResMSIGDB, myGO = myMSIGDB)

total_peaks_vec <- c(max_vec_GO, max_vec_REACTOME, max_vec_KEGG, max_vec_WIKI, max_vec_MSIGDB)
summary(total_peaks_vec)

quantile(total_peaks_vec, 0.1)

##cut Rank at 1872


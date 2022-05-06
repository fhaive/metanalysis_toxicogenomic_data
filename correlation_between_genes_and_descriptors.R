source('load_data.R')
source('compute_correlation.R')

md_datasets = read_molecular_descriptors_datasets(min_unique_number = 10)
gene_data = load_foldchange_dataset()
data_space = gene_data[,rownames(md_datasets$continous_descriptors)]

MDC = md_datasets$continous_descriptors
MDC = MDC[,-c(16:39)]
cor_res = compute_correlation(fold_changes=data_space,
                              descriptors = MDC,
                              nSamples = 10,method = "pearson")

res = threshold_correlation(CorMat = cor_res$CorMat,PVal = cor_res$PVal,top_perc = 10, th_pval = 0.05,th_corr = 0)


ranked_list = read.table("data_input/final_ranked_genes_mean.csv", sep = ";")
ranked_list = ranked_list$V2

gene_lists = res$gene_list
for(i in 1:length(gene_lists)){
  gene_lists[[i]] = rownames(gene_lists[[i]])
} 

stats = length(ranked_list):1
names(stats) = ranked_list

library(fgsea)
fgseaRes = fgsea(gene_lists, stats, scoreType = "pos", eps = 0)
fgseaRes2 = fgseaRes[fgseaRes$padj<0.01,]
fgseaRes2 = fgseaRes2[order(fgseaRes2$padj, decreasing = T),]
fgseaRes2$pathway = factor(fgseaRes2$pathway, levels = fgseaRes2$pathway)
colnames(fgseaRes2)[1] = "descriptor"

ggplot(fgseaRes2, aes(y=descriptor, x=-log(padj))) + 
  geom_point(size=3) + 
  geom_segment(aes(y=descriptor, 
                   yend=descriptor, 
                   x=0, 
                   xend=-log(padj))) + theme_bw()+
  theme(axis.text.y = element_text(angle=0, vjust=0.6,size = 20),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle=0,size = 20),
        axis.title.x = element_text(size = 20)) 

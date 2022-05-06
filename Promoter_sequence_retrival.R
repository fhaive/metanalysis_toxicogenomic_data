if (!requireNamespace("remotes", quietly=TRUE))
  install.packages("remotes")
remotes::install_github("snystrom/memes")
library(meme)
library(openxlsx)
library("biomaRt", lib.loc="/Library/Frameworks/R.framework/Versions/4.1/Resources/library/")


rank_file <- read.table("./final_ranked_genes_mean.txt")
top_rank<-as.vector(rank_file[0:1872,])
in_vitro_long_term_genes<-read.delim("./in_vitro_long_term_july.txt")
in_vivo_long_term_genes<-read.delim("./in_vivo_long_term_july.txt")

ensembl <- useMart("ensembl")
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
searchAttributes(mart = ensembl, pattern = "hgnc")
promoter_sequences_upstream <- biomaRt::getSequence(id=in_vivo_long_term_genes,
                                                    type="ensembl_gene_id",
                                                    seqType="coding_gene_flank",
                                                    upstream=500,
                                                    mart=ensembl) 
colnames(promoter_sequences_upstream)<-c('Upstream_sequence','ensembl_gene_id')
promoter_sequences_downstream <- biomaRt::getSequence(id=in_vivo_long_term_genes,
                                                    type="ensembl_gene_id",
                                                    seqType="coding_gene_flank",
                                                    downstream =100,
                                                    mart=ensembl) 
colnames(promoter_sequences_downstream)<-c('Downstream_sequence','ensembl_gene_id')
promoter_sequences<-cbind(promoter_sequences_upstream,promoter_sequences_downstream$Downstream_sequence)
rownames(promoter_sequences)<-promoter_sequences$ensembl_gene_id
promoter_sequences$ensembl_gene_id<-NULL
for(i in 1:nrow(promoter_sequences)){
  complete_sequence<-paste0(promoter_sequences$Upstream_sequence[i],promoter_sequences$`promoter_sequences_downstream$Downstream_sequence`[i])
  promoter_sequences$complete_sequence[i]<-complete_sequence
}
colnames(promoter_sequences)<-c('Upstream_sequence', 'Downstream_sequence', 'Complete_sequence')
input_exportFasta<-as.data.frame(promoter_sequences[,3])
input_exportFasta$ensembl_gene_id<-rownames(promoter_sequences)
exportFASTA(input_exportFasta, file='./promoters_sequences_long_term_biomarkers.fasta')



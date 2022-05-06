##RUN THE META_ANALYSIS on tggagte data
source(file = 'Documents/Github/Metananalysis/Meta_analysis_functions.R')

#input_data need to be a data frame with genes on the rows and samples as columns.
drug_data<-read.table(file = 'data_input/opentggate_meta_analysis_input.txt')

######################################run the metanalysis

#The effect_size function selects from the original dataframe columns that have the 'P.Valuel in the pattern of their colnmaes. Make sure that the colnames are correct before running. 
#It returns the ranked list of genes with the effect "es" calculated
effect_size_rank<-calc_effect_size_rank(drug_data)

#The pvalue_based function selects from the original dataframe columns that have the 'P.Value in the pattern of their colnmaes. Make sure that the colnames are correct before running. 
#It returns the ranked list of genes with the pval calculated
fisher_based_pvalues<-calc_pvalue_based_rank(drug_data)

#The rank_base function selects from the original dataframe columns that have the _adj.P.Val in the pattern of their colnmaes. Make sure that the colnames are correct before running. 
#It returns the ranked list of genes with the pval calculated
ranks_based_pvalues<-calc_rank_base_rank_subsets(drug_data)

######prepare data for Borda
data<-list()
data[[1]]<-rownames(effect_size_rank)
data[[2]]<-rownames(fisher_based_pvalues)
data[[3]]<-rownames(ranks_based_pvalues)
names(data)<-c('Effect_size','Fisher_test','Rank_based')
outputBorda<-Borda(data)

###select median or mean 
final_ranked_genes_mean<-as.data.frame(outputBorda$TopK$mean)
final_ranked_genes_median<-as.data.frame(outputBorda$TopK$median)
View(final_ranked_genes_median)

#Write file
# write.table(x = final_ranked_genes_mean,file = "/Users/giusy/Documents/Papers/Metanalysis/Meta_analysis_with_Nanosolutions_data/final_ranked_drug_genes_mean.txt")
# write.table(x = final_ranked_genes_mean,file = "/Users/giusy/Documents/Papers/Metanalysis/Meta_analysis_with_Nanosolutions_data/final_ranked_drug_genes_mean.csv")
# write.table(x = final_ranked_genes_median,file = "/Users/giusy/Documents/Papers/Metanalysis/Meta_analysis_with_Nanosolutions_data/final_ranked_drug_genes_median.txt")


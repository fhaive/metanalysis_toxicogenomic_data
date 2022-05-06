# Analyse continous descriptors and their correlation with genes expressed as FC
# PARAMETERS
# nSamples = 5 # PERFORM CORRELATION IF AT LEAST nSamples HAVE VALUES
# th_pval = 0.001 # FOR EACH PAIR OF GENE&MD WE RUN THE COR.TEST FUNCTION. WE STORE THE CORRELATION ONLY IF IT'S PVALUE < th_corr
# th_corr = 0.8 # FOR EACH PAIR OF GENE&MD with PVALUE < th_corr WE consider relevant only abs(correlation values) > th_corr
library(tidyinftheo)

compute_correlation = function(fold_changes,descriptors, nSamples = 10, method = "pearson"){
  
  CorMat = matrix(0, nrow = nrow(fold_changes), ncol = ncol(descriptors))
  rownames(CorMat) = rownames(fold_changes)
  colnames(CorMat) = colnames(descriptors)
  
  PVal = matrix(0, nrow = nrow(fold_changes), ncol = ncol(descriptors))
  rownames(PVal) = rownames(fold_changes)
  colnames(PVal) = colnames(descriptors)
  
  pd = md_datasets$factor_descriptors
 
  pb = txtProgressBar(min =1, max = nrow(fold_changes), style=3)
  for(i in 1:nrow(fold_changes)){
    gi = fold_changes[i,]
    
    for(j in 1:ncol(descriptors)){
      good_idx = which(is.na(descriptors[,j])==FALSE)
      
      cd = descriptors[good_idx,j]
      if(length(good_idx)>nSamples){ 
        
        res = cor.test(gi[good_idx], cd,method = method)
        CorMat[i,j] = res$estimate
        PVal[i,j]=res$p.value
        
      }else{
        print("less than 1 sample")
      }
    }
    setTxtProgressBar(pb,i)
  }
  close(pb)
  
  return(list("CorMat" = CorMat,"PVal"=PVal))
}


threshold_correlation = function(CorMat,PVal,th_pval = 0.001,top_perc = 10, th_corr = 0.8){
  CorMat2 = CorMat
  CorMat2[PVal>th_pval] = 0
  
  n_genes = round((nrow(CorMat) * top_perc)/100)
  
  gene_list = list()
  for(i in 1:ncol(CorMat2)){
    values = CorMat2[,i]
    ranked = order(abs(values), decreasing = T)
    to_rem = ranked[-(1:n_genes)]  
    values[to_rem] = 0
    CorMat2[,i] = values
    gene_list[[i]] = data.frame(genes = rownames(CorMat2)[ranked[1:n_genes]],
                                 correlation = CorMat2[ranked[1:n_genes],i])
  }
  names(gene_list) = colnames(CorMat2)
  non_binary_cor = CorMat2
  
  CorMat2[CorMat2!=0]=1
  res = which(CorMat2==1, arr.ind = T)
  
  actual_cor = rep(0, nrow(res))
  for(i in 1:nrow(res))actual_cor[i] = CorMat[res[i,1],res[i,2]]
  
  res = cbind(res, actual_cor)
  res[,2] = colnames(CorMat)[as.numeric(res[,2])]
  return(list("gene_list"=gene_list,"Thresholded"=non_binary_cor,"BinaryCor" = CorMat2, "Genes_MD_Pairs"=res))
}



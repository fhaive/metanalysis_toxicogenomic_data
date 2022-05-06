<div id="top"></div>

<br />
<div align="center">

<h3 align="center">A gene regulation model reveals an ancestral adaptation response to particulate exposure triggered by nanomaterials</h3>

</div>



<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
      <ul>
        <li><a href="#built-with">Built With</a></li>
      </ul>
    </li>
    <li>
      <a href="#getting-started">Prerequisites</a>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#contact">Contact</a></li>
  </ol>
</details>



<!-- ABOUT THE PROJECT -->
## About The Project

Toxicogenomics aims at characterising the mechanism of action (MOA) of environmental exposures, and often relies on transcriptomics to investigate the responses of exposed biological systems. However, the identification of shared toxicogenomics-derived signatures across exposed biological systems is hampered by the complexity and heterogeneity of transcriptomics data. Given the lack of a clear transcriptomic signature of exposure, we hypothesise that common patterns of gene regulation could explain the response to engineered nanomaterials (ENM) across biological systems, disentangling the complexity of their MOA. We performed meta-analysis of a large collection of transcriptomics data from various ENM exposure studies and identified deregulation of immune functions as a prominent response across different ENM exposures. This pattern of transcriptional deregulation differed significantly from exposure to drugs. By investigating the promoter regions of genes frequently altered both in vitro and in vivo following exposure to ENM, we identified a set of binding sites for zinc finger transcription factors C2H2 involved in chromatin remodelling and immunomodulation. We further demonstrate that this gene regulatory model also underlies the transcriptomic MOA in non-mammal species of ecotoxicological interest exposed to ENMs, suggesting that it may be part of the innate immune system conserved by natural selection.

<p align="right">(<a href="#top">back to top</a>)</p>



### Built With

* [RStudio](https://www.rstudio.com)

<p align="right">(<a href="#top">back to top</a>)</p>



<!-- GETTING STARTED -->

### Prerequisites

To be able to run the code, you'll first need to have a working version of both R and RStudio.

The project has been developed with the following R version and environment:

R version 4.1.0 (2021-05-18)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS Mojave 10.14.6

For running the code you will need to install external libraries.
You can install them using the following commands:

* Libraries
  ```sh
  install.packages("xlsx")
  install.packages("readxl")
  install.packages("openxlsx")
  install.packages("vcd")
  install.packages("tidyinftheo")
  install.packages("ggplot2")
  install.packages("fgsea")
  install.packages("tidyverse")
  install.packages("ClusterProfiler")
  install.packages("DOSE")
  install.packages("biomaRt")
  install.packages("dplyr")
  install.packages("ggridges")
  install.packages("ReactomePA")
  install.packages("hrbrthemes")
  install.packages("stringr")
  install.packages("tidyr")
  install.packages("rcompanion")
  install.packages("esc")
  install.packages("metap")
  install.packages("metafor")
  install.packages("TopKLists")
  install.packages("RankProd")
  install.packages("matrixStats")
  install.packages("meme")
  ```
<p align="right">(<a href="#top">back to top</a>)</p>



<!-- USAGE EXAMPLES -->
## Usage

When downloading the code, the folder will contain the following scripts and a folder of input data (/data_input)



The analysis performed in this study can be reproduced as follows:

- [1] Meta_analysis_space_definition.R

	This script generated the meta_analysis data space. Input required for this script is the data downloaded from  https://zenodo.org/record/3949890#.YlPUri0RqH0.  


- [2]  generate_Meta_Analysis_Rank.R 

	This script generated the meta_analysis rank. The script uses function in Meta_analysis_functions.R 

- [3] feature_selection_based_on_GSEA.R
	
	This script prioritises the top of the meta analysis rank based on GSEA. The script requires the files in the folder Curateed gene sets from MSigDB.  


- [4] Calculate_Pvalues_scores_by_cathegory.R

	This script computes the score for each gene. It requires the meta analysis data space. 


- [6] Compute_scores_biological_systems.R

	This script computes the score for each gene in specific biological contexts. It requires the output of the previous script Calculate_Pvalues_scores_by_cathegory.R.

- [7] Promoter_sequence_retrival.R

	This script outputs DNA sequences of the promoters around the transcription starting site. It requires as input a list of genes of interest.

- [8] Meta_analysis_tggate_data.R

	This script generated the meta_analysis rank for the TG-Gates data. The script uses function in Meta_analysis_functions.R 

- [9] correlation_between_genes_and_descriptors.R 
	
	This script computes the correlation between gene and molecular descriptors and performs the GSEA analysis to identify the molecular descriptors correlated to the genes at the top of the rank. It requires in input the following data files: mds_num.RData, mds_cat.RData, and final_ranked_genes_mean.csv. 

	The script requires load_data.R and compute_correlation.R


<p align="right">(<a href="#top">back to top</a>)</p>



<!-- CONTACT -->
## Contact

Giusy del Giudice - [@_GiusydG](https://twitter.com/_GiusydG)
email: giusy.delgiudice@tuni.fi

Project Link: [https://github.com/fhaive/metanalysis_toxicogenomic_data](https://github.com/fhaive/metanalysis_toxicogenomic_data)

<p align="right">(<a href="#top">back to top</a>)</p>


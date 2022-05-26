# ZFP36_L1_CD4_code
Code used for analysis in Matheson et al manuscript on the role of ZFP36/ZFP36L1 in activated CD4 T cells
<br>
<br>This repository contains the code used for analysis and figure generation for the manuscript entitled "ZFP36 and ZFP36L1 integrate transcriptome and metabolic reprograming of activated mouse CD4+ T cells" (Matheson et al). 
<br>
<br>Analysis workflow and code for plotting figures and generation of supplementary tables:
<br>final_code_for_all_figures.Rmd
<br>HTML R notebook for this is also provided.
<br>
<br>A few other scripts are used as indicated in the workflow:
<br>
<br>R code used for calculation of z-scores for RNA stability:
<br>intensityDifference.R (modified from https://github.com/s-andrews/intensitydiff)
<br>
<br>R scripts used to generate matched control gene sets and assess enrichment for CLIP targets and AREs within gene sets of interest:
<br>expression_match_iCLIP.R
<br>TATT_analysis.R
<br>
<br>Initial annotation of CLIP data also utilises the iCLIP_Genialis_feature_annotation.R available here:
<br>https://github.com/LouiseMatheson/Process_CLIP_data
<br>
<br>Scripts used to download phyloP scores, convert to BED format and then filter for regions of interest:
<br>download_mm10_phyloP60_chr_wig2bed.sh
<br>pull_out_phyloP.R 
<br>
<br>All analysis was performed using R v4.1.2 and packages as shown in the sessionInfo() section at the end of the R notebook. 


Latest release (v1.1):

[![DOI](https://zenodo.org/badge/489508294.svg)](https://zenodo.org/badge/latestdoi/489508294)


# ZFP36_L1_CD4_code
R code used for analysis in Matheson et al manuscript on the role of ZFP36/ZFP36L1 in activated CD4 T cells


This repository contains the code used for analysis and figure generation for the manuscript entitled "ZFP36 and ZFP36L1 integrate transcriptome and metabolic reprograming of activated mouse CD4+ T cells" (Matheson et al). 

final_code_for_all_figures.Rmd - analysis workflow and code for plotting figures and generation of supplementary tables.
HTML R notebook for this is also provided.

A few other scripts are used as indicated in the workflow:
intensityDifference.R (modified from https://github.com/s-andrews/intensitydiff) - used for calculation of z-scores for RNA stability

expression_match_iCLIP.R and TATT_analysis.R - used to generate matched control gene sets and assess enrichment for CLIP targets and AREs within gene sets of interest.

Initial annotation of CLIP data also utilises the iCLIP_Genialis_feature_annotation.R available here:
https://github.com/LouiseMatheson/Process_CLIP_data



All analysis was performed using R v4.1.2 and packages as shown in the sessionInfo() section at the end of the R notebook. 

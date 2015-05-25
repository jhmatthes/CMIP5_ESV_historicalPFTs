# CMIP5_ESV_historicalPFTs
This set of code was developed to analyze differences between the CMIP5 piControl (pre-industrial control) modeled distributions of PFTs and the Euro-American settlement vegetation dataset for the upper Midwest through the northeastern United States.

format_pls_pls.R was executed first to format the dataset for the Euro-American settlement period. 

esv_cmip5_comparison_workflow.R contains the overall workflow to execute the comparison analysis between the set of CMIP5 models that contributed 'landCoverFrac' (PFT fraction) for the pre-industrial control period and the ESV dataset. All other functions are called from this workflow script.

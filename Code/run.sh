#1.Calculate the depletion score for multiplexed CRISPRi screen
Rscript DepletionScore_Calculator.R
#2.Calculate the epistasis score for sgRNA pairs targeting enhancers
python3 gi_score_calculation_ForGithub.py D30vsD0_Norm_LogFC_matrix_MeanCombined.txt D30vsD0_Norm_LogFC_matrix_MeanCombined

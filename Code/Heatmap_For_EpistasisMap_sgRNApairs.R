library(pheatmap)

gi_m = read.table("D30vsD0_Norm_LogFC_matrix_MeanCombined_emap_ave_abba_l.txt",header=T,row.names = 1)
neg_index = grepl("negative_control",colnames(gi_m))
e5_index = grep("Enhancer_5", rownames(gi_m))
pr_index = grep("promoter", rownames(gi_m))
e1234_index = grep("Enhancer_1|Enhancer_2|Enhancer_3|Enhancer_4", rownames(gi_m))
e67_index = grep("Enhancer_6|Enhancer_7", rownames(gi_m))
new_index = c(pr_index,e1234_index,e5_index,e67_index)
nop_new_index = c(e1234_index,e5_index,e67_index)

loc_annot = c()
for (i in seq(1,length(rownames(gi_m)[!neg_index][new_index]))){
  loc_annot = c(loc_annot, strsplit(rownames(gi_m)[!neg_index][new_index][i],";")[[1]][1])
}
loc_annot_rdf = data.frame(Loc=loc_annot,row.names=rownames(gi_m)[!neg_index][new_index])
loc_annot_cdf = data.frame(Loc=loc_annot,row.names=colnames(gi_m)[!neg_index][new_index])
paletteLength=5
myBreaks <- c(-3,-2,-1,1,2,3)
myColor <- colorRampPalette(c("#b6b3b3","#d6d4d4", "white", "pink","red"))(paletteLength)
loc_color = c("#CE2525","#F2A470","#D98068","#BC3331","#681615","#00DDF9","#0087CA","#003B64")
names(loc_color) = c("Myc_promoter","Myc_Enhancer_1","Myc_Enhancer_2",
                     "Myc_Enhancer_3", "Myc_Enhancer_4", "Myc_Enhancer_5",
                     "Myc_Enhancer_6", "Myc_Enhancer_7")
loc_color <- list(Loc = loc_color)

pdf("D30vsD0_Norm_LogFC_matrix_MeanCombined_GIscore_emap_ave_abba_l_sgRNA_WGreyRed.pdf", 10,8)
pheatmap(as.matrix(gi_m[!neg_index,!neg_index][nop_new_index,nop_new_index]), color = myColor, breaks = myBreaks,margins = c(8, 8),
         annotation_row = loc_annot_rdf,annotation_colors = loc_color,border_color=F,
         labels_col="",labels_row = "",annotation_col = loc_annot_cdf,
         cluster_rows=FALSE,cluster_cols=FALSE,na.color = "grey",main="GIscore_sgRNA")
while (!is.null(dev.list()))  dev.off()


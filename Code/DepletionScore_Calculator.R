library('limma')
library(edgeR) ###why good turingproportions?
library(ggplot2)
library(dplyr)
library(UsingR)
library(reshape2)
library(pheatmap)
library(stringr)
###replicate 1
d0_r1 = read.table("dl_mapping_d0_trimmed_rep1.txt",header=T,sep="\t")
rn_r1 = paste(d0_r1[,1],d0_r1[,2],sep="_")
d0m_r1 = d0_r1[3:ncol(d0_r1)]
row.names(d0m_r1) <- rn_r1
d30_r1 = read.table("dl_mapping_d30_trimmed_rep1.txt",header=T,sep="\t")
d30m_r1 = d30_r1[3:ncol(d30_r1)]
row.names(d30m_r1) <- rn_r1

###replicate 2
d0_r2 = read.table("dl_mapping_d0_trimmed_rep2.txt",header=T,sep="\t")
rn_r2 = paste(d0_r2[,1],d0_r2[,2],sep="_")
d0m_r2 = d0_r2[3:ncol(d0_r2)]
row.names(d0m_r2) <- rn_r2
d30_r2 = read.table("dl_mapping_d30_trimmed_rep2.txt",header=T,sep="\t")
d30m_r2 = d30_r2[3:ncol(d30_r2)]
row.names(d30m_r2) <- rn_r2

sgrna_marker = read.table("Library_Info.txt",header = F, sep="\t")
sgrna_seq_marker_list = list()
for (i in seq(1,nrow(d0m_r1))){
  message(sgrna_marker[i,4])
  each_list = list()
  sgrna_seq = as.character(sgrna_marker[i,4])
  if (length(sgrna_seq) > 20){
    each_list[[toupper(substring(sgrna_seq,3,10))]] = sgrna_marker[i,7]
  }
  else{
    each_list[[toupper(substring(sgrna_seq,2,9))]] = sgrna_marker[i,7]
  }
  sgrna_seq_marker_list = append(sgrna_seq_marker_list,each_list )
}

###remove median count < 30
rowmedian_r1 = apply(d0m_r1,1,median)
colmedian_r1 = apply(d0m_r1,2,median)
rowmedian_r2 = apply(d0m_r2,1,median)
colmedian_r2 = apply(d0m_r2,2,median)
low_readcount_r1_index = unique(c(rownames(d0m_r1)[rowmedian_r1< 30],rownames(d0m_r1)[colmedian_r1<30]))
low_readcount_r2_index = unique(c(rownames(d0m_r2)[rowmedian_r2< 30],rownames(d0m_r2)[colmedian_r2<30]))
low_readcount_all = unique(c(low_readcount_r1_index,low_readcount_r2_index))
keep_sgrna_index = !rownames(d0m_r1) %in% low_readcount_all


sgrnap_list = c()
genep_list = c()
same_sgrnap_list = c()
same_gene_list = c()
locp_list = c()
for(i in rownames(d0m_r1)[keep_sgrna_index])
{
  message(i)
  for(j in rownames(d0m_r1)[keep_sgrna_index])
  {
    sgrnap_list=c(sgrnap_list,paste(i,j,sep="|"))
    gene1 = paste(strsplit(i,"_")[[1]][2],strsplit(i,"_")[[1]][3],strsplit(i,"_")[[1]][4],sep="_")
    gene2 = paste(strsplit(j,"_")[[1]][2],strsplit(j,"_")[[1]][3],strsplit(j,"_")[[1]][4],sep="_")
    genep_list=c(genep_list,paste(gene1,gene2,sep="|"))
    gene1_m = sgrna_seq_marker_list[[strsplit(i,"_")[[1]][1]]]
    gene2_m = sgrna_seq_marker_list[[strsplit(j,"_")[[1]][1]]]
    locp_list = c(locp_list, paste(gene1_m,gene2_m,sep="|"))
    if (i==j){
      same_sgrnap_list=c(same_sgrnap_list,paste(i,j,sep="|"))
      same_gene_list=c(same_gene_list,paste(gene1,gene2,sep="|"))
    }
  }
}
####Add pseudocount 10
d0_d30_r1_df = data.frame(sgrnap=sgrnap_list,genep=genep_list,locp=locp_list,
                          d0=c(as.matrix(d0m_r1[keep_sgrna_index,keep_sgrna_index]))+10,
                          d30=c(as.matrix(d30m_r1[keep_sgrna_index,keep_sgrna_index]))+10)
d0_prop_r1 = goodTuringProportions(d0_d30_r1_df$d0)
d30_prop_r1 = goodTuringProportions(d0_d30_r1_df$d30)
logfc_r1_df = data.frame(sgrnap=d0_d30_r1_df$sgrnap,genep=d0_d30_r1_df$genep,
                         d30VSd0=log(d30_prop_r1/d0_prop_r1,2),locp=locp_list,
                         groups=rep("sgRNA|sgRNA",nrow(d0_d30_r1_df)))

d0_d30_r2_df = data.frame(sgrnap=sgrnap_list,genep=genep_list,locp=locp_list,
                          d0=c(as.matrix(d0m_r2[keep_sgrna_index,keep_sgrna_index]))+10,
                          d30=c(as.matrix(d30m_r2[keep_sgrna_index,keep_sgrna_index]))+10)
d0_prop_r2 = goodTuringProportions(d0_d30_r2_df$d0)
d30_prop_r2 = goodTuringProportions(d0_d30_r2_df$d30)
logfc_r2_df = data.frame(sgrnap=d0_d30_r2_df$sgrnap,genep=d0_d30_r2_df$genep,
                         d30VSd0=log(d30_prop_r2/d0_prop_r2,2),locp=locp_list,
                         groups=rep("sgRNA|sgRNA",nrow(d0_d30_r2_df)))

###normalized the logfc based on neg-neg sgRNA pairs 
###(log2FoldChange - mean(Log2FoldChange of Neg_Neg pairs))/sd(Log2FoldChange of Neg_Neg pairs)
new_groups= as.character(logfc_r1_df$groups)
filtered_rownames = rownames(d0m_r1)[keep_sgrna_index]
new_groups[grepl("ve.*\\|",logfc_r1_df$locp) & grepl("\\|.*Negat",logfc_r1_df$locp)] <- "Neg_Neg"
logfc_r1_df_neg_sub = logfc_r1_df[new_groups=="Neg_Neg",]
logfc_norm_r1_df = data.frame(sgrnap=logfc_r1_df$sgrnap,genep=logfc_r1_df$genep,
                                d30VSd0=(logfc_r1_df$d30VSd0-mean(logfc_r1_df_neg_sub$d30VSd0))/sd(logfc_r1_df_neg_sub$d30VSd0),
                                group=new_groups,locp=logfc_r1_df$locp)
logfc_norm_r1_df_matrix <- matrix(logfc_norm_r1_df$d30VSd0,nrow=length(filtered_rownames),byrow=F)

logfc_r2_df_neg_sub = logfc_r2_df[new_groups=="Neg_Neg",]
logfc_norm_r2_df = data.frame(sgrnap=logfc_r2_df$sgrnap,genep=logfc_r2_df$genep,
                                d30VSd0=(logfc_r2_df$d30VSd0-mean(logfc_r2_df_neg_sub$d30VSd0))/sd(logfc_r2_df_neg_sub$d30VSd0),
                                group=new_groups,locp=logfc_r2_df$locp)
logfc_norm_r2_df_matrix <- matrix(logfc_norm_r2_df$d30VSd0,nrow=length(filtered_rownames),byrow=F)

single_sgrna_list = c()
single_sgrna_loc_list = c()
for(i in seq(1,length(logfc_norm_r2_df$sgrnap))){
  message(i)
  seq=strsplit(as.character(logfc_norm_r2_df$sgrnap[i]),"_")[[1]][1]
  marker=strsplit(as.character(logfc_norm_r2_df$locp[i]),"\\|")[[1]][1]
  if (length(grep("Neg", marker))>0) { marker = "negative_control"}
  single_sgrna_list = c(single_sgrna_list, paste(marker, seq, sep=";"))
}
rownames(logfc_norm_r1_df_matrix)=unique(single_sgrna_list);colnames(logfc_norm_r1_df_matrix)=unique(single_sgrna_list)
rownames(logfc_norm_r2_df_matrix)=unique(single_sgrna_list);colnames(logfc_norm_r2_df_matrix)=unique(single_sgrna_list)
logfc_norm_combined = (logfc_norm_r1_df$d30VSd0+logfc_norm_r2_df$d30VSd0)/2.0
logfc_norm_combined_matrix <- matrix(logfc_norm_combined,nrow=length(filtered_rownames),byrow=F)
rownames(logfc_norm_combined_matrix)=unique(single_sgrna_list)
colnames(logfc_norm_combined_matrix)=unique(single_sgrna_list)
###generate the depletion score for sgRNA pairs
write.table(logfc_norm_combined_matrix,file="D30vsD0_Norm_LogFC_matrix_MeanCombined.txt",sep="\t",quote = F,row.names=T)

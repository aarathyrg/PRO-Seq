#Hierarchical clustering_heatmap with deseq2-selected genes.
#Updated 21_09_22
#Aarathy
#################
#packages
#################
library("DESeq2")
library("pheatmap")
##################
getwd()
setwd("D:\\vm_data_cache_hdd/Katrin_F_Aarathy/bw_proseq_correct_index/Proseq_paper/")
# count matrix from PRO-Seq data
countdata <- read.table("counts_proseq_correct_index_11_07_22_new.txt",header = T)
countdata <- as.matrix(countdata)
head(countdata)
cts <-countdata
cts <- cts[,55:81]
cts<- as.data.frame(cts)
# creating sample info coldata
treatment <- factor(c(rep("ut",3), rep(c(rep("IFNb",3),rep("IFNg",3)),4)))
time <- factor(c(rep("ut",3),rep("1h30min",6),rep("4h",6),rep("24h",6),rep("48h",6)))
treatment_time <- factor(paste(time,treatment,sep = "_"))
replicate  <- factor(rep(c("R1","R2","R3"),9))
names<- paste("WT",time,treatment,replicate, sep = "_")
coldata <-data.frame(row.names = names,treatment_time)
cts <-as.matrix(cts)
cts <- cts[ ,rownames(coldata)]
all(rownames(coldata)==colnames(cts))
is.na(cts) %>% table()
head(cts)
#########################################################################
dds <- DESeqDataSetFromMatrix(countData = cts,colData = coldata, design = ~treatment_time)
##############################################
library(pheatmap)
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}
#######################
#variance stabilization
#######################
vsd<-vst(dds,blind = F)
vsd<-assay(vsd)
head(vsd)
############################
#first 1000 genes selected based on padj values from deseq2(treatment to untreated) 
#genes pooled from IFN beta and IFN gamma treatment
###########################
genes<-read.table("before_clust_candidate_padj_0.01_log2fc1_first1000each_basedonpadj_total3126_genes.txt",header = T)
nrow(genes)
gene_list <- c(genes$genes)
# select candidate genes for plotting
data_plot<-vsd[rownames(vsd) %in% gene_list, ]
nrow(data_plot)
head(data_plot)
####################################################
breaksList<-seq(-3,3,by =0.5)
col<- hcl.colors(11, "YlGnBu")
col1<-rev(col)# reverse the order of selected colours
#z-score normalization
data_subset_norm <- t(apply(data_plot, 1, cal_z_score))

#####################################################
pdf("Results_DESeq2/Wald_beta_gamma/Plots/Annotated_heatmap_from_wald_clust_11_first1000each_basedonpadj.pdf")

out <- pheatmap(data_subset_norm, color = col1,breaks = breaksList,
                cluster_rows = T, cluster_cols = F, 
                clustering_method = "ward.D2",
                show_rownames = F, 
                show_colnames = T, 
                cutree_rows = 11)
dev.off()
clusters <- factor(cutree(out$tree_row, k=11))[out$tree_row[["order"]]]

                          cluster=as.factor(clusters))

#####################################################
#Obtain_clusters
#####################################################
LRT.clust <- cbind(data_subset_norm, 
                   cluster = cutree(res$tree_row, 
                                    k = 11))
LRT.clust <-as.data.frame(LRT.clust)
LRT.clust["genes"]<-rownames(LRT.clust)
write.table(LRT.clust,"TABLES/candidate_first1000_basedonpadj_clust11_genes.txt",row.names = T, col.names = T,sep = "\t")
################################################################################################
*************************************
#################################################################################################


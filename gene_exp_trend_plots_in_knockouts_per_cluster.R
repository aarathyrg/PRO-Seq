#updated 11-09-22
#Aarathy
#To plot trend plots of wildtype and knockous across timepoints separately for IFNb and IFNg
#Out of the vst reads, select candidate genes and calculate zscore across treatment conditions
#Sections-
#Create plot for all genes -IFN beta
#Create plot for all genes -IFN gamma
#Create plot for individual genes
############################################
#Packages
############################################
library("DESeq2")
library("tidyverse")
library("dplyr")
############################################
#DESeq2 again to obtain normalized reads this time together with KO
############################################

setwd("D:\\vm_data_cache_hdd/Katrin_F_Aarathy/bw_proseq_correct_index/Proseq_paper/")
countdata <- read.table("counts_proseq_correct_index_11_07_22_new.txt",header = T)
countdata <- as.matrix(countdata)
head(countdata)
cts <-countdata

######################
# creating coldata, sample names matching columns of count table
######################
genotype <- factor(c(rep("WT", 27), rep("IRF9", 27),rep("IRF1",27)))
genotype<-factor(genotype,levels = c("WT","IRF9","IRF1"))
treatment <- factor(c(rep(c(rep("ut",3),rep(c(rep("IFNb", 3),rep("IFNg", 3)),4)),3)))
time <- factor(rep(c(rep("ut",3),rep("1h30min",6),rep("4h",6),rep("24h",6),rep("48h",6)),3))
treatment_time <- factor(paste(treatment,time,sep = "_"))

replicate <- factor(rep(c("R1","R2","R3"),27))
names <- factor(paste(genotype,time,treatment,replicate,sep = "_"))
coldata <-data.frame(row.names = names, genotype, treatment, time, treatment_time)
rownames(coldata)
colnames(cts)

#just for labelling
timepoints <- factor(rep(c(rep("untreated",3),rep("1h30min",6),rep("4h",6),rep("24h",6),rep("48h",6)),3))
timepoints<- factor(timepoints,levels = c("untreated","1h30min","4h", "24h", "48h"))
sample_info <-data.frame(row.names = names, genotype, treatment, timepoints)


#select columns of countdata based on sample names in coldata 
cts_all <- cts[ ,rownames(coldata)]
all(rownames(coldata)==colnames(cts_all))
is.na(cts_all) %>% table()
head(cts_all)
#DESeq2 function
dds <- DESeqDataSetFromMatrix(countData = cts_all,colData = coldata, design = ~treatment_time)
##############################################
#varaiance stabilized transformation
##############################################

vsd_all<-vst(dds,blind = F)
vsd_all<- assay(vsd_all)
vsd_all<-as.data.frame(vsd_all)
vsd_all["genes"]<-rownames(vsd_all)
head(vsd_all)
######################Create plot for all genes -IFN beta##############################################
#select IFN beta samples specifically
##############################################
vsd_all_IFNb<- vsd_all[,c(grep("ut_ut",colnames(vsd_all)),grep("IFNb",colnames(vsd_all)))]
head(vsd_all_IFNb)
vsd_all_IFNb["genes"]<-rownames(vsd_all_IFNb)
#############################################
#covert to tibble
############################################

vsd_all_IFNb<-as_tibble(vsd_all_IFNb)

############################################
#selecting required sample information (rows) from coldata
############################################
sample_info_IFNb<-sample_info[c(grep("ut_ut",colnames(vsd_all)),grep("IFNb",colnames(vsd_all))),]
sample_info_IFNb["sample"]<-rownames(sample_info_IFNb)
head(vsd_all_IFNb)
############################################
#convert table to long format
############################################
vsd_long_IFNb <- vsd_all_IFNb %>% pivot_longer(cols =c(WT_ut_ut_R1:IRF1_48h_IFNb_R3), 
                                               names_to = "sample", 
                                               values_to = "cts",
                                               values_transform = list(cts=as.numeric))
vsd_long_IFNb <- full_join(vsd_long_IFNb, sample_info_IFNb, by = "sample")
head(vsd_long_IFNb)

##############################################################
#candidate genes are selected based on DESeq2 -untreated and each treatment condition,
#select first 1000 genes based on Padj--> pool all conditions
##############################################################
candidate_genes <- read.table("/Users/aarathyrg/Dropbox/Aarathy_katrin_proseq/Results_DESeq2/For_paper/Wald_beta_gamma/IFNb_IFNg_wald_combined_unique_with_clust_11_genes.txt",header = T)
candidate_genes<-candidate_genes$gene
candidate_genes <- candidate_genes%>%unique()
head(candidate_genes)
##############################################################
#Select for candidate genes , group by genes and calculate z score 
##############################################################
#to plot in the right order of treatment
vsd_long_IFNb$treatment_time<-factor(vsd_long_IFNb$treatment_time, levels=c("ut_ut","IFNb_1h30min","IFNb_4h","IFNb_24h","IFNb_48h"))

data_wt_ko_long_IFNb<- vsd_long_IFNb %>% dplyr::filter(genes %in% candidate_genes)%>%
  group_by(genes)%>%
  mutate(cts_zscore=(cts-mean(cts))/sd(cts))%>%
  group_by(genes,genotype,timepoints)%>%
  #take the mean within one genotype_timepoint across the three replicate for a specific gene)
  summarise(mean_zscore_of_replicates = mean(cts_zscore),
            nrep=n())%>%ungroup()

head(data_wt_ko_long_IFNb)
#table with gene-cluster info
gene_cluster<- read.table("Results_DESeq2/Wald_beta_gamma/candidate_first1000_basedonpadj_clust11_genes.txt",sep ="\t",header = T)
#merge with data_wt_ko_long_IFNb
vsd_cluster_IFNb <- data_wt_ko_long_IFNb %>% 
  inner_join(gene_cluster, by = "genes")

head(vsd_cluster_IFNb)

#######################
#plot trends of IFNb induced genes across timepoints, 
#seperated by clusters(columns) 
#genotype(rows)
#######################

png("/FIGURES/PRO-Seq/IFNb_trend_plot.png")
vsd_cluster_IFNb %>% 
  ggplot(aes(timepoints, mean_zscore_of_replicates)) +
  geom_line(aes(group = genes), alpha = 0.1, colour ="gray") +
  geom_line(stat = "summary", fun = "median", colour ="#420B55", size = 0.5, 
            aes(group = 1)) +
  facet_grid(rows = vars(genotype) , cols =vars(cluster))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12,face="bold"))+
  theme(
    axis.text.x = element_text(size = 9, angle = 90,hjust = 1,vjust = 0.5)
  )+
  ggtitle("Genes in clusters in IFNb treatment")
dev.off()


######################*Create plot for all genes -IFN gamma*##############################################
#select IFN gamma samples specifically
##############################################
vsd_all_IFNg<- vsd_all[,c(grep("ut_ut",colnames(vsd_all)),grep("IFNg",colnames(vsd_all)))]
head(vsd_all_IFNg)
vsd_all_IFNg["genes"]<-rownames(vsd_all_IFNg)
#############################################
#covert to tibble
############################################

vsd_all_IFNg<-as_tibble(vsd_all_IFNg)

############################################
#selecting required sample information (rows) from coldata
############################################
sample_info_IFNg<-sample_info[c(grep("ut_ut",colnames(vsd_all)),grep("IFNg",colnames(vsd_all))),]
sample_info_IFNg["sample"]<-rownames(sample_info_IFNg)
head(vsd_all_IFNg)
############################################
#convert table to long format
############################################
vsd_long_IFNg <- vsd_all_IFNg %>% pivot_longer(cols =c(WT_ut_ut_R1:IRF1_48h_IFNg_R3), 
                                               names_to = "sample", 
                                               values_to = "cts",
                                               values_transform = list(cts=as.numeric))
vsd_long_IFNg <- full_join(vsd_long_IFNg, sample_info_IFNg, by = "sample")
head(vsd_long_IFNg)

##############################################################
#Select for candidate genes , group by genes and calculate z score 
##############################################################

data_wt_ko_long_IFNg<- vsd_long_IFNg %>% dplyr::filter(genes %in% candidate_genes)%>%
  group_by(genes)%>%
  mutate(cts_zscore=(cts-mean(cts))/sd(cts))%>%
  group_by(genes,genotype,timepoints)%>%
  #take the mean within one genotype_timepoint across the three replicate for a specific gene)
  summarise(mean_zscore_of_replicates = mean(cts_zscore),
            nrep=n())%>%ungroup()

head(data_wt_ko_long_IFNg)
#merge with data_wt_ko_long_IFNg
vsd_cluster_IFNg <- data_wt_ko_long_IFNg %>% 
  inner_join(gene_cluster, by = "genes")

#######################
#plot trends of IFNg induced genes across timepoints, 
#seperated by clusters(columns) &
#genotype(rows)
#######################
######################*Create plot for individual genes*##############################################

png("FIGURES/PRO-Seq/IFNg_trend_plot.png")

vsd_cluster_IFNg %>% 
  ggplot(aes(timepoints, mean_zscore_of_replicates)) +
  geom_line(aes(group = genes), alpha = 0.1, colour ="gray") +
  geom_line(stat = "summary", fun = "median", colour ="#065287", size = 0.5, 
            aes(group = 1)) +
  facet_grid(rows = vars(genotype) , cols =vars(cluster))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12,face="bold"))+
  theme(
    axis.text.x = element_text(size = 10, angle = 90,hjust = 1,vjust = 0.5)
  )+
  ggtitle("Genes in clusters in IFNg treatment")
dev.off()
###########################################################################
#Genes_of_interest
######################################################################################################
#cd72 IFNb
Cd72<-vsd_cluster_IFNb%>% filter(genes=="Cd72")
png("FIGURES/Genes_in_clusters/Cd72_cluster_3_IFNb.png")
Cd72%>%
  ggplot(aes(timepoints, mean_zscore_of_replicates)) +
  geom_line(aes(group = genes), alpha = 0.1, colour ="gray") +
  geom_line(stat = "summary", fun = "median", colour ="#420B55", size = 0.5, 
            aes(group = 1)) +
  facet_grid(rows = vars(genotype) , cols =vars(cluster))+
  
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=10,face="bold"))+
  theme(
    axis.text.x = element_text(size = 7, angle = 90,hjust = 1,vjust = 0.5)
  )+
  ggtitle("Cd72")

dev.off()
######################
#cd72 ifng
Cd72<-vsd_cluster_IFNg%>% filter(genes=="Cd72")
png("FIGURES/Genes_in_clusters/Cd72_cluster_3_IFNg.pnf")
Cd72%>%
  ggplot(aes(timepoints, mean_zscore_of_replicates)) +
  geom_line(aes(group = genes), alpha = 0.1, colour ="gray") +
  geom_line(stat = "summary", fun = "median", colour ="#065287", size = 0.5, 
            aes(group = 1)) +
  facet_grid(rows = vars(genotype) , cols =vars(cluster))+
  
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=10,face="bold"))+
  theme(
    axis.text.x = element_text(size = 7, angle = 90,hjust = 1,vjust = 0.5)
  )+
  ggtitle("Cd72_IFNg")

dev.off()

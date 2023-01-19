#updated 16-02-22
#Aarathy
# to obtain pausing index: ratio of read density of promoter region to gene body
#plot log(pausing indices at treatment/pausing indices at homeostasis) in clusters of
#genes defined by hierarchical clustering
####################
#packages
###################
library("dplyr")
library("gplots")
library("tidyverse")
####################
#set work directory
####################
getwd()
setwd("/Users/aarathyrg/Dropbox/Aarathy_katrin_proseq/PRO-seq/FIGURES/")
######################
######################
#refgene annotation file 
#######################
#annotated genes
refGene <- read.table("../../PRO-seq/gencode.mm10.annotation_11_07_22.bed")[,1:6]
refGene <- refGene[grep("random|Un|hap", refGene$V1, invert=TRUE),]
refGene.excluded <- refGene[(refGene$V3-refGene$V2)<=1000,]
refGene <- refGene[(refGene$V3-refGene$V2)>1000,]
bodies<-refGene

#regions excluding first 500bp used for estimating reads in gene body
bodies$V2[bodies$V6 == "+"] <- bodies$V2[bodies$V6 == "+"]+500
bodies$V3[bodies$V6 == "+"] <- bodies$V3[bodies$V6 == "+"]-500
bodies$V3[bodies$V6 == "-"] <- bodies$V3[bodies$V6 == "-"]-500
bodies$V2[bodies$V6 == "-"] <- bodies$V2[bodies$V6 == "-"]+500
bodies <-unique(bodies)
head(bodies)
bodies<-bodies[!duplicated(bodies$V5),]

# info about which gene body-chromosome position used in the analysis
colnames(bodies)<-c("chr","start","end","ensemble","genes","strand")

# similarly first 500 bp from TSS was used for estimating reads in promoter region
####################
#loading countdata file
####################

#countdata from gene body-excluding first 500bp from annotated gene
countdata_gene_body <- read.table("counts_proseq_TSS_plus_TSS_minus_500_for_polII_pausing.txt",header = T)
countdata_gene_body<-as.data.frame(countdata_gene_body[55:81])

#countdata  first 500bp of annotated gene
countdata_TSS <- read.table("counts_on_promoters_TSS500_proseq_correct.txt",header = T)
countdata_TSS<-as.data.frame(countdata_TSS[55:81])
colnames(countdata_TSS)<-paste("TSS",colnames(countdata_gene_body),sep = "_")
countdata_TSS["genes"]<-rownames(countdata_TSS)
countdata_gene_body["genes"]<-rownames(countdata_gene_body)
TSS_gene_body<-merge(countdata_TSS,countdata_gene_body,by="genes")


##############################################
#sample_info
###############################################

genotype <- factor(c(rep("WT", 27)))

treatment <- factor((c(rep("ut",3),rep(c(rep("IFNb", 3),rep("IFNg", 3)),4))))
time <- factor(c(rep("ut",3),rep("1h30min",6),rep("4h",6),rep("24h",6),rep("48h",6)))
treatment_time <- factor(paste(treatment,time,sep = "_"))

replicate <- factor(rep(c("R1","R2","R3"),9))
names <- factor(paste(genotype,time,treatment,replicate,sep = "_"))
sample_info <-data.frame(row.names = names, genotype, treatment, time, treatment_time)

############################################
#merge tables to get genes and chromosome position info
############################################

#combine gene and chromosome position info
polII_combined_table<-merge(TSS_gene_body,bodies,by="genes")

rownames(polII_combined_table)<-polII_combined_table$genes
########################################################

#merge tables to get info about which cluster each gene belongs
gene_cluster<- read.table("../../Results_DESeq2//Wald_beta_gamma/candidate_first1000_basedonpadj_clust11_genes.txt",sep ="\t", header = T)
genes_in_cluster<-gene_cluster[c("cluster","genes")]
##########################################################

polII_combined_genes_cluster<-polII_combined_table[genes_in_cluster$genes,]
polII_combined_genes_cluster<-merge(polII_combined_genes_cluster,genes_in_cluster,by="genes")
rownames(polII_combined_genes_cluster)<- polII_combined_genes_cluster$genes
########################################################

polII_combined_genes_cluster["gene_length"]<-polII_combined_genes_cluster$end-polII_combined_genes_cluster$start
head(polII_combined_genes_cluster[2:55])
#############################################################
#normalizing to cpm both tss and genebody
#############################################################
polII_combined_genes_cluster[2:55]<-polII_combined_genes_cluster[2:55]/1000000
#########################################################

polII_wt_cluster<-polII_combined_genes_cluster[2:55]#rename as polII_wt_cluster excluding first column

polII_wt_cluster[1:27]<-polII_wt_cluster[1:27]/500
polII_wt_cluster[28:54]<-polII_wt_cluster[28:54]/polII_combined_genes_cluster$gene_length
polII_wt_cluster[28:54]

#############################################
ratio_polII<-data.frame()
ratio_polII<-polII_wt_cluster[1:27]/polII_wt_cluster[28:54]

colnames(ratio_polII)<-gsub("TSS_","",colnames(ratio_polII))
ratio_polII_IFNb<-ratio_polII[c(grep("ut_ut",colnames(ratio_polII)),grep("IFNb",colnames(ratio_polII)))]
#filter for infinity and 0
ratio_polII_IFNb<- filter_if(ratio_polII_IFNb, is.numeric, all_vars((.) != "inf"))
ratio_polII_IFNb<- filter_if(ratio_polII_IFNb, is.numeric, all_vars((.) != 0))
ratio_polII_IFNb["genes"]<-rownames(ratio_polII_IFNb)
###################################
#combine cluster info by merging table
ratio_polII_IFNb<-full_join(ratio_polII_IFNb,genes_in_cluster,by="genes")

#########################################################################################
#longer version of table
polII_index_IFNb<-ratio_polII_IFNb%>%pivot_longer(
  cols = c(WT_ut_ut_R1:WT_4h_IFNb_R3),
  names_to = "sample",
  values_to = "polII_index",
  values_transform = list(polII_index=as.numeric))


sample_info["sample"]<-rownames(sample_info)
polII_index_IFNb<-full_join(polII_index_IFNb,sample_info,by ="sample")
# arranging timepoints in correct order
polII_index_IFNb$treatment_time<-factor(polII_index_IFNb$treatment_time,levels = c("ut_ut","IFNb_1h30min","IFNb_4h","IFNb_24h","IFNb_48h"))

#group together and take mean
polII_index_IFNb<-polII_index_IFNb %>% group_by(cluster,genes,treatment_time)%>%
  summarise(mean_polII_index= mean(polII_index))
########################################################################################

polII_index_IFNb<-polII_index_IFNb%>%drop_na()
#wide table format
polII_index_IFNb_wide <- polII_index_IFNb %>% pivot_wider(names_from = treatment_time,values_from = mean_polII_index)
head(polII_index_IFNb_wide)
#add columns with ratios to untreated
polII_index_IFNb_wide<-polII_index_IFNb_wide%>% mutate_at(vars(ut_ut:IFNb_48h), list(norm_to_ut=~./ut_ut))
#select only the ratio columns
polII_index_IFNb_wide<-polII_index_IFNb_wide[,c(1,2,8:12)]

colnames(polII_index_IFNb_wide)<-c("cluster", "genes","ut", "IFNb_1h30min", "IFNb_4h","IFNb_24h", "IFNb_48h")
polII_index_IFNb_wide[3:7]<-sapply(polII_index_IFNb_wide[3:7], function(x)log(x))# applying log function
head(polII_index_IFNb_wide)
polII_log_ratio<-polII_index_IFNb_wide%>% pivot_longer(cols = c(ut:IFNb_48h),
                                                      names_to = "treatment_time",
                                                      values_to = "polII_index_ratio")
#order correctly
polII_log_ratio$treatment_time<-factor(polII_log_ratio$treatment_time,
                                       levels =c("ut","IFNb_1h30min","IFNb_4h","IFNb_24h","IFNb_48h"))


polII_log_ratio<-polII_log_ratio%>%filter(polII_index_ratio!="-inf")# remove inf values

#violiin plots showing log ratio of untreated to each treatment condition
png("../../../../Desktop/CODES/IFNb_log_ratio_pol_II_index_to_ut_violin.png")
polII_log_ratio%>%#group_by(treatment_time)%>%
  ggplot(aes(treatment_time, polII_index_ratio)) +
  geom_violin()+
  coord_cartesian(ylim=c(-3.5,3.5))+stat_summary(fun = median,geom = "point")+
  #geom_jitter()+
  facet_grid(cols=vars(cluster))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=10,face="bold"))+
  theme(
    axis.text.x = element_text(size = 7, angle = 90,hjust = 1,vjust = 0.5))
dev.off()

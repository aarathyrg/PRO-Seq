#To obtain putative TREs (enhancer regions defined by dREG(Dankolab) likely belonging to selected genes
#here- genes in cluster1-typical ISG cluster
# Te TRE regions are already filtered for those having(DESeq2 based) upregulated transcripts and
#in the range of 50kb from TSS of the selected genes(first 1000 genes based on padj)
#updated 09-11-22
#Aarathy

#packages
#####################
library("DESeq2")
library("dplyr")
library("TxDb.Mmusculus.UCSC.mm10.knownGene")

#######################################
#function
#######################################
center.bed<-function(df,number){
  df[3]<-round(((df[3]+df[4])/2)-number)
  df[4]<-round(df[3]+2*number)
  return(df)}

# selected TRE regions columns:chr, start, end,  respective log2fc,  padj (for each condition), genes
TRE_dreg_genes<-read.table("TRE_up_50kb_range_from_genes_of_interest.txt",header = T)
#genes and respective clusters
gene_cluster<-read.table("genes_in_cluster.txt",header = T)
gene_cluster<-gene_cluster[,c("cluster","genes")]
head(gene_cluster$genes[gene_cluster$cluster == 1])

TRE_dreg_genes_cluster<-merge(TRE_dreg_genes,gene_cluster[,c("genes","cluster")],by="genes")

######################################################################
#selecting _TRES_BELONGING TO GENES IN CLUSTER1
######################################################################
cluster<-paste(rep("cluster",11),1:11,sep = "")

clst<-c(1,2,3,9,10)

datf<-data.frame()

for (i in clst){
 datf<-as.data.frame(TRE_dreg_genes_cluster[TRE_dreg_genes_cluster$cluster == i,])
 datf<-center.bed(datf,200)
 if (i==1) { 
   datf1<-datf %>%filter(IFNb_1h30min_l2fc>=1) %>% filter(IFNb_1h30min_padj<=0.01)
   datf1<-rbind(datf1,datf %>%filter(IFNg_1h30min_l2fc>=1) %>% filter(IFNg_1h30min_padj<=0.01))%>%distinct()
   #removing_duplicated region resulting from selecting the same dreg region due to multiple genes
   #first three clolumns-chr start end
   datf_f<-datf1<-datf1[!duplicated(datf[1:3]),] 
    } else if (i==2) {
   datf1<-datf %>%filter(IFNb_1h30min_l2fc>=1) %>% filter(IFNb_1h30min_padj<=0.01)
   datf2%>%filter(IFNg_1h30min_l2fc>=1) %>% filter(IFNg_1h30min_padj<=0.01)
   datf_f<-rbind(datf1,datf2) %>%distinct()
   datf_f<-datf_f<-datf_f[!duplicated(datf_f[1:3]),] 
    } else if (i==3) {
   datf1<-datf %>%filter(IFNb_4h_l2fc>=1) %>% filter(IFNb_4h_padj<=0.01)
   datf2<-datf %>%filter(IFNb_24h_l2fc>=1) %>% filter(IFNb_24h_padj<=0.01)
   datf3<-datf %>%filter(IFNb_48h_l2fc>=1) %>% filter(IFNb_48h_padj<=0.01)
   datf_f<-rbind(datf1,datf2,datf3)%>%distinct()
   datf_f<-datf_f<-datf_f[!duplicated(datf_f[1:3]),] 
    } else if (i==9) {
   datf1<-datf%>%filter(IFNg_24h_l2fc>=1) %>% filter(IFNg_24h_padj<=0.01)
   datf2<-datf%>%filter(IFNg_4h_l2fc>=1) %>% filter(IFNg_4h_padj<=0.01)
   datf3<-datf%>%filter(IFNg_48h_l2fc>=1) %>% filter(IFNg_48h_padj<=0.01)
   datf4<-datf%>%filter(IFNg_1h30min_l2fc>=1) %>% filter(IFNg_1h30min_padj<=0.01)
   datf_f<-rbind(datf1,datf2,datf3,datf4)%>%distinct()
   datf_f<-datf_f<-datf_f[!duplicated(datf_f[1:3]),] 
    } else if (i==10) {
   datf1<-datf%>%filter(IFNg_24h_l2fc>=1) %>% filter(IFNg_24h_padj<=0.01)
   datf2<-datf%>%filter(IFNg_48h_l2fc>=1) %>% filter(IFNg_48h_padj<=0.01)
   datf3<-datf%>%filter(IFNg_4h_l2fc>=1) %>% filter(IFNg_4h_padj<=0.01)
   datf_f<-rbind(datf1,datf2,datf3)%>%distinct()
    }

  write.table(datf_f,file = paste("/Users/aarathyrg/Desktop/CODES/Enhancers/",
                               cluster[i],".txt",sep = ""),sep = "\t",col.names = T,row.names = F)
  }



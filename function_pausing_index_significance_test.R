#updated 2023_01_17
# to perform wilcoxon rank sum test on pausing indices of untreated vs treated samples
#Aarathy
#set work directory
setwd("/Users/aarathyrg/Desktop/CODES/")
# 
#input table containing pausing indices of each gene(row) per condition(column)

PI<-read.table("pol_II_pausing_index_per_gene.txt",header = T)
PI<-as.data.frame(PI)

# remove non finite values from colums containing indices
PI<-pi[is.finite(rowSums(pi[3:6])),]

#df: data frame; ut_col_index: column index of untreated condition
#tr_col_index:column index of untreated condition to be tested with wilcoxon test
PI_significance<-function(df,ut_col_index,tr_col_index){
  significance<-wilcox.test(df[,ut_col_index],df[,tr_col_index])
  return(significance)
}



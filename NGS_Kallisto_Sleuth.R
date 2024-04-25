#Set working directory
setwd("D:/NGS_Kai_Jenny")

# Loading libraries
library(sleuth)
library(biomaRt)
library(dplyr)

#read the design matrix for protein coding Genes or ncRNA
s2c = read.table('./DesignMatrix/DesignTablecDNATAZ_NTC.csv', sep=',', header=TRUE, stringsAsFactors = FALSE)
#s2c = read.table('./DesignMatrix/DesignTablecDNAYAP_NTC.csv', sep=',', header=TRUE, stringsAsFactors = FALSE)
#s2c = read.table('./DesignMatrix/DesignTablencRNATAZ_NTC.csv', sep=',', header=TRUE, stringsAsFactors = FALSE)
#s2c = read.table('./DesignMatrix/DesignTablencRNAYAP_NTC.csv', sep=',', header=TRUE, stringsAsFactors = FALSE)


#Annotations to map with the reference genome for mouse data
AnnotationReference <- biomaRt::useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
annotations <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "transcript_version", "ensembl_gene_id", "external_gene_name", "description", "transcript_biotype", "chromosome_name"), mart = AnnotationReference)
annotations <- dplyr::rename(annotations, target_id = ensembl_transcript_id, ens_gene = ensembl_gene_id, ext_gene = external_gene_name, biotype = transcript_biotype, location=chromosome_name)
annotations <- dplyr::select(annotations, c('target_id', 'ens_gene', 'ext_gene', 'biotype', 'location'))
head(annotations)

#Annotations to map with the reference Genome for Human data
#AnnotationReference <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
#annotations <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "transcript_version", "ensembl_gene_id", "external_gene_name", "description", "transcript_biotype", "chromosome_name"), mart = AnnotationReference)
#annotations <- dplyr::rename(annotations, target_id = ensembl_transcript_id, ens_gene = ensembl_gene_id, ext_gene = external_gene_name, biotype = transcript_biotype, location=chromosome_name)
#annotations <- dplyr::select(annotations, c('target_id', 'ens_gene', 'ext_gene', 'biotype', 'location'))
#head(annotations)


#Prepare the sleuth model and fit the data to get information on all genes all biotypes
so <- sleuth_prep(s2c, extra_bootstrap_summary = TRUE, target_mapping=annotations, gene_mode=T, aggregation_column='ens_gene', transformation_function = function(x) log2(x + 0.5))
so <- sleuth_fit(so, ~grouping, 'full')
so <- sleuth_fit(so, ~1, 'reduced')

#Perform the likelihood ratio and wald test
so <- sleuth_lrt(so, 'reduced', 'full')
so <- sleuth_wt(so, which_beta = 'groupingTAZ')

#Output Results of lrt and wald test into dataframes, then merge by target ids
res_lrt <- sleuth_results(so, 'reduced:full', 'lrt', show_all = TRUE)
res_wt <- sleuth_results(so, test = 'groupingTAZ', show_all = TRUE)
res <- merge(res_lrt, res_wt[, c('target_id', 'b', 'se_b', 'mean_obs')], on = 'target_id', sort = FALSE)

#Rename the default sleuth columns to something more understandable
res <-dplyr::rename(res,LRTS=test_stat,log2FC=b, se_log2FC=se_b)

#Extract est expression values and tpm transcripts per million counts for each sample
expressionvalues<-sleuth_to_matrix(so,'obs_norm','scaled_reads_per_base')
expressiondata=as.data.frame(expressionvalues)
expressiondata=tibble::rownames_to_column(expressiondata,'target_id')

transcriptspermillion<-sleuth_to_matrix(so,'obs_norm','tpm')
transcriptpermilliondata=as.data.frame(transcriptspermillion)
transcriptspermilliondata=tibble::rownames_to_column(transcriptpermilliondata,'target_id')

Results_est <-dplyr::inner_join(res,expressiondata,by='target_id')
Results_est_tpm <-dplyr::inner_join(Results_est,transcriptspermilliondata,by='target_id',suffix = c('.exp','.tpm'))


# Saving Unfiltered list of GENES with counts, tpm, Fold changes FC, lrts
write.csv(Results_est_tpm,'./Results/Unfiltered/MouseTAZvsNTC.csv')

#Filter by FDR which is represented as q.val
res_significant <- dplyr::filter(Results_est_tpm, qval <= 0.05)

#Save significant by FDR
write.csv(res_significant,'./Results/Filtered/Significant_TAZvsNTC.csv')

#Load the significantly expressed list 
res_significant= read.csv('./Results/Filtered/Significant_TAZvsNTC.csv', row.names = 1)

#Filter by genetype for protein Coding
proteincoding= res_significant[res_significant$biotype=='protein_coding',]
#noncodingRNA= res_significant[res_significant$biotype=='lncRNA',]

#Remove duplicated entries
proteincoding= dplyr::distinct(proteincoding,target_id,.keep_all = TRUE)
#noncodingRNA= dplyr::distinct(noncodingRNA,target_id,.keep_all = TRUE)

#Filter by Fold change threshold
FinalDEGsP <- proteincoding %>% dplyr::filter(!between(log2FC,log2(1/1.5), log2(1.5)))

#FinalDEGslncRNA <- noncodingRNA %>% dplyr::filter(!between(log2FC,log2(1/1.5), log2(1.5)))

#save as filtered final list for differentially expression Protein coding List
write.csv(FinalDEGsP,'./Results/Filtered/Filtered_DEGsP_TAZvsNTC_FDR_FC.csv')

#write.csv(FinalDEGslncRNA,'./Results_ncRNA/Filtered_DEGslncRNA_YAP2vsNTC_1_5FC.csv')'


#Heat map plotting,Read the filtered file list

#Load the list of significantly expressed genes and define it as dataframe for heatmap
InputHeat= read.csv('./Results/Filtered/Filtered_DEGsP_TAZvsNTC_FDR_FC.csv', row.names = 1)

#Filtered to get the Top 100 list, you can change this number to your choice e.g. n=50
Top=dplyr::slice_min(InputHeat,order_by=pval, n=100)

#Filter by direction
Upregulated <- Top %>% dplyr::filter(log2FC>0)
Downregulated <- Top %>% dplyr::filter(log2FC<0)
CombinedUD=rbind(Upregulated,Downregulated)

#Use the following commands if you want a rank ordered list
#Top50=dplyr::slice_max(InputHeat,order_by=log2FC, n=50)
#Bot50=dplyr::slice_min(InputHeat,order_by=log2FC, n=50)
#Top100UD=rbind(Top50,Bot50)

BiocManager::install("ComplexHeatmap")
library(circlize)
library(ComplexHeatmap)
library(ggplot2)

#Column number depends on how many samples you have tpm
mat= as.matrix(CombinedUD[26:33])

rownames(mat)=CombinedUD[,3]
matplusone=mat+1
logmat=log(matplusone,10)

#center and scale each column (Z-score) then transpose
mat.scaled <- t(apply(logmat, 1, scale)) 
colnames(mat.scaled)<-colnames(logmat)

#getting log2 value for each gene we are keeping
l2_val <- as.matrix(CombinedUD$log2FC) 
colnames(l2_val)<-"log2FC"

mycols <- colorRamp2(breaks = c(0,1.5,3 ), 
                     colors = c("blue", "white", "red"))

#show row name False when you have too many candidates for removing genenames
h2 <-Heatmap(mat.scaled, width=1,height=8,show_row_names = TRUE, name = "Z-score",row_names_gp = gpar(fontsize = 7, fontface ="bold"))
h1 <-Heatmap(l2_val, cluster_rows=FALSE,width=0.1,height=8, show_column_names=TRUE, name='log2FC',
             cell_fun = function(j, i, x, y, w, h, col) {grid.text(round(l2_val[i, j],2), x, y,gp = gpar(fontsize = 4, fontface = "bold"))})
h <-h1+h2
h


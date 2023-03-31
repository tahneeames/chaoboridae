
##data wrangling and data merging from trinotate + DE analysis 
### read in a trinotate files
chaoTrinotate=read.csv2("~/Desktop/Chaoborid/resultsFiles/chaoborus/trimmed/chaoborus.trimmed.Trinotate.txt",sep="\t")
mochTrinotate=read.csv2("~/Desktop/Chaoborid/resultsFiles/mochlonyx/trimmed/mochlonyx.trimmed.Trinotate.txt",sep="\t")
eucTrinotate=read.csv2("~/Desktop/Chaoborid/resultsFiles/eucorethra/trimmed/eucorethra.trimmed.Trinotate.txt",sep="\t")
#filter trinotate output by DE gene list (see .sh file for trinotate generation)
#first rename the rows to merge in both dataframes for both organisms 
names(chaoTrinotate)[1] <- "gene_id"
names(mochTrinotate)[1] <- "gene_id"
names(eucTrinotate)[1] <- "gene_id"

#merge full pairwise dataframes from DESeq 
#Chaoborus 
chao_mtPVAL_lfc=filter(chao_deMT5_2repDF,pvalue <0.05&log2FoldChange > 2)
chao_giPVAL_lfc=filter(chao_deGI5_2repDF,pvalue <0.05&log2FoldChange > 2)
chao_de_merge_pVal_lfc=merge(chao_mtPVAL_lfc, chao_giPVAL_lfc, by=0)
#Mochlonyx
moch_mtPVAL_lfc=filter(moch_deMT5_2repDF,pvalue <0.05&log2FoldChange > 2)
moch_giPVAL_lfc=filter(moch_deGI5_2repDF,pvalue <0.05&log2FoldChange > 2)
moch_de_merge_pVal_lfc=merge(moch_mtPVAL_lfc, moch_giPVAL_lfc, by=0)
#Eucorethra  
euc_mtPVAL_lfc=filter(euc_deMT5_DF,pvalue <0.05&log2FoldChange > 2)
euc_giPVAL_lfc=filter(euc_deGI5_DF,pvalue <0.05&log2FoldChange > 2)
euc_de_merge_pVal_lfc=merge(euc_mtPVAL_lfc, euc_giPVAL_lfc, by=0)


#then merge with trinotate file 
#rename these ones too just like with the other df
names(chao_de_merge_pVal_lfc)[1] <- "gene_id"
names(moch_de_merge_pVal_lfc)[1] <- "gene_id"
names(euc_de_merge_pVal_lfc)[1] <- "gene_id"
#merge to have an annotated dataframe with lfc values and gene homologies/ontologies/pathways
chao_all_trinotate=merge(x = chaoTrinotate, y = chao_de_merge_pVal_lfc, by = "gene_id",all = FALSE)
moch_all_trinotate=merge(x = mochTrinotate, y = moch_de_merge_pVal_lfc, by = "gene_id",all = FALSE)
euc_all_trinotate=merge(x = eucTrinotate, y = euc_de_merge_pVal_lfc, by = "gene_id",all = FALSE)

#adding averages for TPM values of tissue replicates 
#chaoborus
chao_geneCounts_2reps$as_TPM <- rowMeans(chao_geneCounts_2reps[ , c(1,2)], na.rm=TRUE)
chao_geneCounts_2reps$mt_TPM <- rowMeans(chao_geneCounts_2reps[ , c(5,6)], na.rm=TRUE)
chao_geneCounts_2reps$gi_TPM <- rowMeans(chao_geneCounts_2reps[ , c(3,4)], na.rm=TRUE)
chao_geneCounts_2reps$gene_id<-rownames(chao_geneCounts_2reps)

#mochlonyx
moch_geneCounts_2reps$as_TPM <- rowMeans(moch_geneCounts_2reps[ , c(1,2)], na.rm=TRUE)
moch_geneCounts_2reps$mt_TPM <- rowMeans(moch_geneCounts_2reps[ , c(5,6)], na.rm=TRUE)
moch_geneCounts_2reps$gi_TPM <- rowMeans(moch_geneCounts_2reps[ , c(3,4)], na.rm=TRUE)
moch_geneCounts_2reps$gene_id<-rownames(moch_geneCounts_2reps)

#eucorethra 
euc_geneCounts$tr_TPM <- rowMeans(euc_geneCounts[ , c(7,9)], na.rm=TRUE)
euc_geneCounts$mt_TPM <- rowMeans(euc_geneCounts[ , c(4,6)], na.rm=TRUE)
euc_geneCounts$gi_TPM <- rowMeans(euc_geneCounts[ , c(1,3)], na.rm=TRUE)
euc_geneCounts$gene_id<-rownames(euc_geneCounts)

#mean lfc in merged dataset 
chao_all_trinotate$meanLFC <- rowMeans(chao_all_trinotate[ , c(19,27)], na.rm=TRUE)
moch_all_trinotate$meanLFC <- rowMeans(moch_all_trinotate[ , c(19,27)], na.rm=TRUE)
euc_all_trinotate$meanLFC <- rowMeans(euc_all_trinotate[ , c(19,27)], na.rm=TRUE)

#merge 
chao_trinotate_geneTPM=merge(x = chao_all_trinotate, y = chao_geneCounts_2reps, by = "gene_id")
moch_trinotate_geneTPM=merge(x = moch_all_trinotate, y = moch_geneCounts, by = "gene_id")
euc_trinotate_geneTPM=merge(x = euc_all_trinotate, y = euc_geneCounts, by = "gene_id")
#reorder by mean lfc 
chao_LFC=arrange(chao_trinotate_geneTPM,-meanLFC)
moch_LFC=arrange(moch_trinotate_geneTPM,-meanLFC)
euc_LFC=arrange(euc_trinotate_geneTPM,-meanLFC)
#change . to na
is.na(chao_LFC) <- chao_LFC == "."
is.na(moch_LFC) <- moch_LFC == "."
is.na(euc_LFC) <- euc_LFC == "."
#remove rows that don't have blast or pfam hits 
chao_LFC %>% filter(!if_all(c(sprot_Top_BLASTX_hit, sprot_Top_BLASTP_hit, Pfam), is.na)) -> chao_LFC_pop
moch_LFC %>% filter(!if_all(c(sprot_Top_BLASTX_hit, sprot_Top_BLASTP_hit, Pfam), is.na)) -> moch_LFC_pop
euc_LFC %>% filter(!if_all(c(sprot_Top_BLASTX_hit, sprot_Top_BLASTP_hit, Pfam), is.na)) -> euc_LFC_pop

#remove rows of isoforms for the same gene  
chao_LFC_genes=distinct(chao_LFC_pop, gene_id, .keep_all = TRUE)
moch_LFC_genes=distinct(moch_LFC_pop, gene_id, .keep_all = TRUE)
euc_LFC_genes=distinct(euc_LFC_pop, gene_id, .keep_all = TRUE)
#save
write.csv(chao_LFC_genes,file="chaoborus_asUp_trinotate_gene&TPM_pfamBlAST_lfc2.csv")
write.csv(moch_LFC_genes,file="mochlonyx_trUp_trinotate_gene&TPM_pfamBlAST_lfc2.csv")
write.csv(euc_LFC_genes,file="eucorethra_trUp_trinotate_gene&TPM_pfamBlAST_lfc2.csv")


#make gene TPM matrix for full trinotate output 
chao_trinotate_geneTPM_all=merge(x = chaoTrinotate, y = chao_geneCounts_2reps, by = "gene_id")
moch_trinotate_geneTPM_all=merge(x = mochTrinotate, y = moch_geneCounts, by = "gene_id")
euc_trinotate_geneTPM_all=merge(x = eucTrinotate, y = euc_geneCounts, by = "gene_id")
#remove duplicates/isoforms
chao_trinotate_geneTPM_all=distinct(chao_trinotate_geneTPM_all, gene_id, .keep_all = TRUE)
moch_trinotate_geneTPM_all=distinct(moch_trinotate_geneTPM_all, gene_id, .keep_all = TRUE)
euc_trinotate_geneTPM_all=distinct(euc_trinotate_geneTPM_all, gene_id, .keep_all = TRUE)

#get frequencies of unique pfam identities in this table 
#first remove information to exract the Pfam identities only (here I'm removing everything following the first ^ from the trinotate output to only keep the pfam id)
chao_LFC_genes$Pfam_id <- sub('(\\^).*$', '', chao_LFC_genes$Pfam, perl=TRUE)
#frequency table of each pfam hit from this new column 
chao_pfamFreq=as.data.frame(table(chao_LFC_genes$Pfam_id))

#add counts to the original dataframe too! 
chao_LFC_genes%>%
  group_by(Pfam_id) %>%
  mutate(pfamCount = n()) -> chao_LFC_genes
#save the new sheets with frequencies added 
write.csv(chao_LFC_genes,file="chaoborus_asUp_trinotate_gene&TPM_pfamBlAST_lfc2_PfamFeq.csv")


#search function for pfam terms 
#chaoborus
c_pf <- function(TERM){
  #search for a string in the pfam column 
  chaoFun <- chao_trinotate_geneTPM_all[grepl(TERM, chao_trinotate_geneTPM_all$Pfam), ]
  #gather hits into data frame 
  chaoFun <<- chaoFun
  #write pfam domain table onto directory 
  write.table(chaoFun,paste0('chao_',TERM),sep = '\t',row.names = FALSE )
  #print number of hits 
  return(dim(chaoFun))
}

#mochlonyx same function 
m_pf <- function(TERM){
  mochFun <- moch_trinotate_geneTPM_all[grepl(TERM, moch_trinotate_geneTPM_all$Pfam), ]
  mochFun <<- mochFun
  write.table(mochFun,paste0('moch_',TERM),sep = '\t',row.names = FALSE )
  return(dim(mochFun))
}

#eucorethra same function 
e_pf <- function(TERM){
  eucFun <- euc_trinotate_geneTPM_all[grepl(TERM, euc_trinotate_geneTPM_all$Pfam), ]
  eucFun <<- eucFun
  write.table(eucFun,paste0('euc_',TERM),sep = '\t',row.names = FALSE )
  return(dim(eucFun))
}

#heatmap for searched genes 
#chaoborus
chaoFun.matrix=chaoFun[,18:20]
rownames(chaoFun.matrix)<- chaoFun$gene_id
chaoFun.matrix=as.matrix(chaoFun.matrix)
chaoFun.matrix=chaoFun.matrix[rowSums(chaoFun.matrix[, -1])>0, ]
pheatmap(chaoFun.matrix,scale="row",color=rev(natparks.pals("WindCave",n=150,type="continuous")),labels_col=c("as_TPM"="air-sac","mt_TPM"="malpighian \ntubules","gi_TPM"="midgut")) 

#mochlonyx
mochFun.matrix=mochFun[,18:20]
rownames(mochFun.matrix)<- mochFun$gene_id
mochFun.matrix=as.matrix(mochFun.matrix)
mochFun.matrix=mochFun.matrix[rowSums(mochFun.matrix[, -1])>0, ]
pheatmap(mochFun.matrix,scale="row",color=rev(natparks.pals("WindCave",n=150,type="continuous")),labels_col=c("tr_TPM"="trachea","mt_TPM"="malpighian \ntubules","gi_TPM"="midgut")) 

#eucorethra
eucFun.matrix=eucFun[,18:20]
rownames(eucFun.matrix)<- eucFun$gene_id
eucFun.matrix=as.matrix(eucFun.matrix)
eucFun.matrix=eucFun.matrix[rowSums(eucFun.matrix[, -1])>0, ]
pheatmap(eucFun.matrix,scale="row",color=rev(natparks.pals("WindCave",n=150,type="continuous")),labels_col=c("tr_TPM"="trachea","mt_TPM"="malpighian \ntubules","gi_TPM"="midgut")) 


#function to build a table from pfam domains and number of hits from upregulated gene list 
extract_PF_substrings <- function(string_list) {
  # Define regular expression to match desired substrings
  regex <- "PF.*?\\^.*?\\^.*?\\^"
  # Apply regex to each string in list and extract matching substrings
  substrings <- lapply(string_list, function(x) {
    matches <- regmatches(x, gregexpr(regex, x))
    unlist(matches)
  })
  # Return list of extracted substrings
  return(substrings)
}

chao.pfam.hits <- unlist(extract_PF_substrings(chao_LFC_genes$Pfam))
moch.pfam.hits <- unlist(extract_PF_substrings(moch_LFC_genes$Pfam))
euc.pfam.hits <- unlist(extract_PF_substrings(euc_LFC_genes$Pfam))

write.csv(sort(table(chao.pfam.hits),decreasing=TRUE),'chao.pfam.csv')
write.csv(sort(table(moch.pfam.hits),decreasing=TRUE),'moch.pfam.csv')
write.csv(sort(table(euc.pfam.hits),decreasing=TRUE),'euc.pfam.csv')



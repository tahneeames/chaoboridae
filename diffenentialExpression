###differential expression anaylsis
library(DESeq2)
library(ggplot2)
library(ggthemes)
library(scales)
library(patchwork)
library(gridExtra)
library(MetBrewer)
library(dplyr)
#write salmon gene count matrix file into a csv or tab delimited txt file
#write sample info file 
#upload both into R 
#output genetated in Unix using salmon - code for generating count matrix file is available in bash script file  
chao_geneCounts=read.table("~/Desktop/Chaoborid/resultsFiles/chaoborus/trimmed/chaoborus-trimmed.salmon.gene.TPM.not_cross_norm",header=T,sep="\t",row.names=1)
chao_metadata=read.table("~/Desktop/Chaoborid/resultsFiles/chaoborus/chao_metadata.txt",header=T,sep="\t", row.names=1)
View(chao_geneCounts)
View(chao_metadata)

#convert gene count dataframe into matrix
chaoCount_matrix=as.matrix(chao_geneCounts)
View(chaoCount_matrix)
#double checking! 
is.matrix(chaoCount_matrix)
#convert replicate number from numeric
chao_metadata$replicate=as.factor(chao_metadata$replicate)
#reordering dataframes by tissue type
chaoCount_matrix <- chaoCount_matrix[, c(1,2,3,7,8,9,4,5,6)]
chao_metadata=chao_metadata[c(2,5,9,1,3,7,4,6,8),]


#make deseq2 dataframe from these dataframes
chao_dds2 <- DESeqDataSetFromMatrix(countData=round(chaoCount_matrix), colData=chao_metadata, design= ~ tissue)
#remove values less than 5
chao_dds5 <- chao_dds2[rowSums(counts(chao_dds2)) >= 5,]
#run deseq2 on both 
chao_deSeq5=DESeq(chao_dds5)
#extract results
chao_deMT5 = results(chao_deSeq5,contrast=c("tissue","airSac","malTube"))
chao_deGI5=results(chao_deSeq5,contrast=c("tissue","airSac","gut"))
chao_deMG5=results(chao_deSeq5,contrast=c("tissue","gut","malTube"))

#PCA
chao_rld<-rlog(chao_deSeq5,blind=TRUE)
#plot 
chao_3pca=plotPCA(chao_rld,intgroup=c("tissue","replicate"),returnData=TRUE)
chao_pca=ggplot(chao_3pca,aes(PC1,PC2,col=tissue,shape=replicate,fill=tissue))+
geom_point(size=3)+
scale_color_manual(name="Tissue",labels=c("Air sac","Midgut","Malpighian tubule"),values=c("#8b6e28","#d4ba66","#b7751f"),limits=c("airSac","malTube","gut"))+
scale_fill_manual(name="Tissue",labels=c("Air sac","Midgut","Malpighian tubule"),values=c("#8b6e28","#d4ba66","#b7751f"),limits=c("airSac","malTube","gut"))+
scale_shape_manual(name="Replicate",labels=c("1","2","3"),values=c(21,22,24),limits=c("1","2","3"))+
theme_few()+
theme(plot.title = element_text(size = 12, family = "Tahoma", face = "bold"),
  text = element_text(size = 10, family = "Tahoma"),
  axis.title = element_text(face="bold",size=12),
  axis.text.x=element_text(size = 10),
  axis.text.y=element_text(size = 10),
  strip.text = element_text(face="bold"))+
theme(legend.title = element_text(size = 10, face = "bold"),text=element_text(size = 10, family = "Tahoma"))

#remaking with two reps to see what the pca looks like without these wild outlliers 
chao_geneCounts_2reps=read.table("~/Desktop/Chaoborid/resultsFiles/chaoborus/trimmed/chaoborus-trimmed.salmon.gene.TPM.not_cross_norm_2reps.txt",header=T,sep="\t",row.names=1)
chao_metadata_2reps=read.table("~/Desktop/Chaoborid/resultsFiles/chaoborus/chao_metadata_2reps.txt",header=T,sep="\t", row.names=1)
#match the orders 
chao_metadata_2reps=chao_metadata_2reps[c(1,2,5,6,3,4),]
#convert gene count dataframe into matrix
chao_geneCounts_2reps_matrix=as.matrix(chao_geneCounts_2reps)
#convert replicate number from numeric
chao_metadata_2reps$rep=as.character(chao_metadata_2reps$rep)
#make deseq2 dataframe from these dataframes
chao_dds_2rep <- DESeqDataSetFromMatrix(countData=round(chao_geneCounts_2reps_matrix), colData=chao_metadata_2reps, design= ~ tissue)
#remove values less than 5
chao_dds5_2reps <- chao_dds_2rep[rowSums(counts(chao_dds_2rep)) >= 5,]
#run deseq2
chao_deSeq_2reps=DESeq(chao_dds5_2reps)
chao_DE=data.frame(results(chao_deSeq_2reps)) 
chao_rld_2reps<-rlog(chao_deSeq_2reps,blind=TRUE)

#checking pca again 
chao_2pca=plotPCA(chao_vst_2reps,intgroup=c("tissue","rep"),returnData=TRUE)
chao_pca_2reps=ggplot(chao_2pca,aes(PC1,PC2,col=tissue,shape=factor(rep),fill=tissue))+
geom_point(size=3)+
scale_color_manual(name="Tissue",labels=c("Air sac","Midgut","Malpighian tubule"),values=c("#8b6e28","#d4ba66","#b7751f"),limits=c("airSac","malTube","gut"))+
scale_fill_manual(name="Tissue",labels=c("Air sac","Midgut","Malpighian tubule"),values=c("#8b6e28","#d4ba66","#b7751f"),limits=c("airSac","malTube","gut"))+
scale_shape_manual(name="Replicate",labels=c("1","3"),values=c(21,24),limits=c("1","2"))+
theme_few()+
theme(plot.title = element_text(size = 12, family = "Tahoma", face = "bold"),
  text = element_text(size = 10, family = "Tahoma"),
  axis.title = element_text(face="bold",size=12),
  axis.text.x=element_text(size = 10),
  axis.text.y=element_text(size = 10),
  strip.text = element_text(face="bold"))+
theme(legend.title = element_text(size = 10, face = "bold"),text=element_text(size = 10, family = "Tahoma"))

ggsave("chao_pca_2reps.png",plot=last_plot(),bg="transparent")

#building pairwise contrasts
chao_deMT5_2reps = results(chao_deSeq_2reps,contrast=c("tissue","airSac","malTube"))
chao_deGI5_2reps=results(chao_deSeq_2reps,contrast=c("tissue","airSac","gut"))
chao_deMT5_2repDF =as.data.frame(chao_deMT5_2reps)
chao_deGI5_2repDF = as.data.frame(chao_deGI5_2reps)
#differentual expression column based on p-values and log2-fold change values 
chao_deMT5_2repDF$diffexpressed10 <- "NO"
chao_deMT5_2repDF$diffexpressed10[chao_deMT5_2repDF$log2FoldChange > 10 & chao_deMT5_2repDF$pvalue < 0.05] <- "UP"
chao_deMT5_2repDF$diffexpressed10[chao_deMT5_2repDF$log2FoldChange < -10 & chao_deMT5_2repDF$pvalue < 0.05] <- "DOWN"
chao_deGI5_2repDF$diffexpressed10 <- "NO"
chao_deGI5_2repDF$diffexpressed10[chao_deGI5_2repDF$log2FoldChange > 10 & chao_deGI5_2repDF$pvalue < 0.05] <- "UP"
chao_deGI5_2repDF$diffexpressed10[chao_deGI5_2repDF$log2FoldChange < -10 & chao_deGI5_2repDF$pvalue < 0.05] <- "DOWN"

chao_deMT5_2repDF$diffexpressed5 <- "NO"
chao_deMT5_2repDF$diffexpressed5[chao_deMT5_2repDF$log2FoldChange > 5 & chao_deMT5_2repDF$pvalue < 0.05] <- "UP"
chao_deMT5_2repDF$diffexpressed5[chao_deMT5_2repDF$log2FoldChange < -5 & chao_deMT5_2repDF$pvalue < 0.05] <- "DOWN"
chao_deGI5_2repDF$diffexpressed5 <- "NO"
chao_deGI5_2repDF$diffexpressed5[chao_deGI5_2repDF$log2FoldChange > 5 & chao_deGI5_2repDF$pvalue < 0.05] <- "UP"
chao_deGI5_2repDF$diffexpressed5[chao_deGI5_2repDF$log2FoldChange < -5 & chao_deGI5_2repDF$pvalue < 0.05] <- "DOWN"
#trying new stuff in the middle of this randomly feb 2023 
chao_mt_asUp_10=subset(chao_deMT5_2repDF, diffexpressed10 == 'UP')
chao_gi_asUp_10=subset(chao_deGI5_2repDF, diffexpressed10 == 'UP')
chao_as_up_Join10lfc=merge(chao_mt_asUp_10, chao_gi_asUp_10, by=0)
chao_mt_asUp_5=subset(chao_deMT5_2repDF, diffexpressed5 == 'UP')
chao_gi_asUp_5=subset(chao_deGI5_2repDF, diffexpressed5 == 'UP')
chao_as_up_Join5lfc=merge(chao_mt_asUp_5, chao_gi_asUp_5, by=0)
#print 
write.csv(chao_as_up_Join10lfc,"chao_airSac_upregulatedGene_list_LFC10.csv")
write.csv(chao_as_up_Join5lfc,"chao_airSac_upregulatedGene_list_LFC5.csv")

#mochlonyx 
moch_geneCounts=read.table("~/Desktop/Chaoborid/resultsFiles/mochlonyx/trimmed/mochlonyx-trimmed.salmon.gene.TPM.not_cross_norm",header=T,sep='\t',row.names=1)
moch_metadata=read.table("~/Desktop/Chaoborid/resultsFiles/mochlonyx/moch_metadata.txt",header=T,sep="\t", row.names=1)
#reordering dataframes by tissue type 
#note! might not acrually have to do this if the order makes sense! :) 
moch_metadata=moch_metadata[c(8,5,3,7,2,9,6,4,1),]
#convert gene count file into matrix
moch_geneCounts_matrix=as.matrix(euc_geneCounts)
#convert replicate number from numeric 
moch_metadata$rep=as.factor(moch_metadata$rep)

#make deseq2 dataframe from these dataframes
moch_dds <- DESeqDataSetFromMatrix(countData=round(moch_geneCounts_matrix), colData=moch_metadata, design= ~ tissue)
#remove values less than 5
moch_dds5 <- moch_dds[rowSums(counts(moch_dds)) >= 5,]
#run deseq2
moch_deSeq5=DESeq(moch_dds5)
moch_DE=data.frame(results(moch_deSeq5))
moch_deMT5 = results(moch_deSeq5,contrast=c("tissue","trach","malTube"))
moch_deGI5=results(moch_deSeq5,contrast=c("tissue","trach","gut"))
moch_deMT5_DF =as.data.frame(moch_deMT5)
moch_deGI5_DF = as.data.frame(moch_deGI5)
#differentual expression column based on p-values and log2-fold change values 
moch_deMT5_DF$diffexpressed10 <- "NO"
moch_deMT5_DF$diffexpressed10[moch_deMT5_DF$log2FoldChange > 10 & moch_deMT5_DF$pvalue < 0.05] <- "UP"
moch_deMT5_DF$diffexpressed10[moch_deMT5_DF$log2FoldChange < -10 & moch_deMT5_DF$pvalue < 0.05] <- "DOWN"
moch_deGI5_DF$diffexpressed10 <- "NO"
moch_deGI5_DF$diffexpressed10[moch_deGI5_DF$log2FoldChange > 10 & moch_deGI5_DF$pvalue < 0.05] <- "UP"
moch_deGI5_DF$diffexpressed10[moch_deGI5_DF$log2FoldChange < -10 & moch_deGI5_DF$pvalue < 0.05] <- "DOWN"

moch_deMT5_DF$diffexpressed5 <- "NO"
moch_deMT5_DF$diffexpressed5[moch_deMT5_DF$log2FoldChange > 5 & moch_deMT5_DF$pvalue < 0.05] <- "UP"
moch_deMT5_DF$diffexpressed5[moch_deMT5_DF$log2FoldChange < -5 & moch_deMT5_DF$pvalue < 0.05] <- "DOWN"
moch_deGI5_DF$diffexpressed5 <- "NO"
moch_deGI5_DF$diffexpressed5[moch_deGI5_DF$log2FoldChange > 5 & moch_deGI5_DF$pvalue < 0.05] <- "UP"
moch_deGI5_DF$diffexpressed5[moch_deGI5_DF$log2FoldChange < -5 & moch_deGI5_DF$pvalue < 0.05] <- "DOWN"

moch_mt_trUp_10=subset(moch_deMT5_DF, diffexpressed10 == 'UP')
moch_gi_trUp_10=subset(moch_deGI5_DF, diffexpressed10 == 'UP')
moch_tr_up_Join10lfc=merge(moch_mt_trUp_10, moch_gi_trUp_10, by=0)
moch_mt_trUp_5=subset(moch_deMT5_DF, diffexpressed5 == 'UP')
moch_gi_trUp_5=subset(moch_deGI5_DF, diffexpressed5 == 'UP')
moch_tr_up_Join5lfc=merge(moch_mt_trUp_5, moch_gi_trUp_5, by=0)

write.csv(moch_tr_up_Join10lfc,"moch_trachea_upregulatedGene_list_LFC10.csv")
write.csv(moch_tr_up_Join5lfc,"moch_trachea_upregulatedGene_list_LFC5.csv")

moch_rld<-rlog(moch_deSeq5,blind=TRUE)
moPCA=plotPCA(moch_rld,intgroup=c("tissue","rep"),returnData=TRUE)
#plot in ggplot2
moPCAplot=ggplot(moPCA,aes(PC1,PC2,col=tissue,shape=rep,fill=tissue))+
geom_point(size=3)+
scale_color_manual(name="Tissue",labels=c("Trachea","Midgut","Malpighian tubule"),values=c("#8b6e28","#d4ba66","#b7751f"),limits=c("trach","malTube","gut"))+
scale_fill_manual(name="Tissue",labels=c("Trachea","Midgut","Malpighian tubule"),values=c("#8b6e28","#d4ba66","#b7751f"),limits=c("trach","malTube","gut"))+
scale_shape_manual(name="Replicate",labels=c("1","2","3"),values=c(21,22,24),limits=c("1","2","3"))+
theme_few()+
theme(plot.title = element_text(size = 12, family = "Tahoma", face = "bold"),
  text = element_text(size = 10, family = "Tahoma"),
  axis.title = element_text(face="bold",size=12),
  axis.text.x=element_text(size = 10),
  axis.text.y=element_text(size = 10),
  strip.text = element_text(face="bold"))+
theme(legend.title = element_text(size = 10, face = "bold"),text=element_text(size = 10, family = "Tahoma"))

ggsave("moch_PCAplot.png",plot=last_plot(),bg="transparent")


#now with eucorethra
euc_geneCounts=read.table("~/Desktop/Chaoborid/resultsFiles/eucorethra/trimmed/eucorethra-trimmed.salmon.gene.TPM.not_cross_norm",header=T,sep='\t',row.names=1)
euc_metadata=read.table("~/Desktop/Chaoborid/resultsFiles/eucorethra/euc_metadata.txt",header=T,sep="\t", row.names=1)
#reordering dataframes by tissue type
euc_metadata=euc_metadata[c(8,5,3,7,2,9,6,4,1),]
#convert gene count file into matrix
euc_geneCounts_matrix=as.matrix(euc_geneCounts)
#convert replicate number from numeric 
euc_metadata$rep=as.factor(euc_metadata$rep)

#make deseq2 dataframe from these dataframes
euc_dds <- DESeqDataSetFromMatrix(countData=round(euc_geneCounts_matrix), colData=euc_metadata, design= ~ tissue)
#remove values less than 5
euc_dds5 <- euc_dds[rowSums(counts(euc_dds)) >= 5,]
#run deseq2
euc_deSeq5=DESeq(euc_dds5)
euc_DE=data.frame(results(euc_deSeq5))
euc_deMT5 = results(euc_deSeq5,contrast=c("tissue","trach","malTube"))
euc_deGI5=results(euc_deSeq5,contrast=c("tissue","trach","gut"))
euc_deMT5_DF =as.data.frame(euc_deMT5)
euc_deGI5_DF = as.data.frame(euc_deGI5)
#differentual expression column based on p-values and log2-fold change values 
euc_deMT5_DF$diffexpressed10 <- "NO"
euc_deMT5_DF$diffexpressed10[euc_deMT5_DF$log2FoldChange > 10 & euc_deMT5_DF$pvalue < 0.05] <- "UP"
euc_deMT5_DF$diffexpressed10[euc_deMT5_DF$log2FoldChange < -10 & euc_deMT5_DF$pvalue < 0.05] <- "DOWN"
euc_deGI5_DF$diffexpressed10 <- "NO"
euc_deGI5_DF$diffexpressed10[euc_deGI5_DF$log2FoldChange > 10 & euc_deGI5_DF$pvalue < 0.05] <- "UP"
euc_deGI5_DF$diffexpressed10[euc_deGI5_DF$log2FoldChange < -10 & euc_deGI5_DF$pvalue < 0.05] <- "DOWN"

euc_deMT5_DF$diffexpressed5 <- "NO"
euc_deMT5_DF$diffexpressed5[euc_deMT5_DF$log2FoldChange > 5 & euc_deMT5_DF$pvalue < 0.05] <- "UP"
euc_deMT5_DF$diffexpressed5[euc_deMT5_DF$log2FoldChange < -5 & euc_deMT5_DF$pvalue < 0.05] <- "DOWN"
euc_deGI5_DF$diffexpressed5 <- "NO"
euc_deGI5_DF$diffexpressed5[euc_deGI5_DF$log2FoldChange > 5 & euc_deGI5_DF$pvalue < 0.05] <- "UP"
euc_deGI5_DF$diffexpressed5[euc_deGI5_DF$log2FoldChange < -5 & euc_deGI5_DF$pvalue < 0.05] <- "DOWN"

euc_mt_trUp_10=subset(euc_deMT5_DF, diffexpressed10 == 'UP')
euc_gi_trUp_10=subset(euc_deGI5_DF, diffexpressed10 == 'UP')
euc_tr_up_Join10lfc=merge(euc_mt_trUp_10, euc_gi_trUp_10, by=0)
euc_mt_trUp_5=subset(euc_deMT5_DF, diffexpressed5 == 'UP')
euc_gi_trUp_5=subset(euc_deGI5_DF, diffexpressed5 == 'UP')
euc_tr_up_Join5lfc=merge(euc_mt_trUp_5, euc_gi_trUp_5, by=0)

write.csv(euc_tr_up_Join10lfc,"euc_trachea_upregulatedGene_list_LFC10.csv")
write.csv(euc_tr_up_Join5lfc,"euc_trachea_upregulatedGene_list_LFC5.csv")

euc_rld<-rlog(euc_deSeq5,blind=TRUE)
euPCA=plotPCA(euc_rld,intgroup=c("tissue","rep"),returnData=TRUE)
#plot in ggplot2
euPCAplot=ggplot(euPCA,aes(PC1,PC2,col=tissue,shape=rep,fill=tissue))+
geom_point(size=3)+
scale_color_manual(name="Tissue",labels=c("Trachea","Midgut","Malpighian tubule"),values=c("#8b6e28","#d4ba66","#b7751f"),limits=c("trach","malTube","gut"))+
scale_fill_manual(name="Tissue",labels=c("Trachea","Midgut","Malpighian tubule"),values=c("#8b6e28","#d4ba66","#b7751f"),limits=c("trach","malTube","gut"))+
scale_shape_manual(name="Replicate",labels=c("1","2","3"),values=c(21,22,24),limits=c("1","2","3"))+
theme_few()+
theme(plot.title = element_text(size = 12, family = "Tahoma", face = "bold"),
  text = element_text(size = 10, family = "Tahoma"),
  axis.title = element_text(face="bold",size=12),
  axis.text.x=element_text(size = 10),
  axis.text.y=element_text(size = 10),
  strip.text = element_text(face="bold"))+
theme(legend.title = element_text(size = 10, face = "bold"),text=element_text(size = 10, family = "Tahoma"))

ggsave("euc_PCAplot.png",plot=last_plot(),bg="transparent")

#make some plots before moving on to more data wrangling and viz 
##Chaoborus
#MA plot 
#plot in ggplot 
chao_MA_MT=ggplot(chao_deMT5_DF, aes(baseMean, log2FoldChange, color=padj)) + 
geom_point(size=1,alpha=0.5) + 
scale_y_continuous(oob=squish) + 
scale_x_log10() + 
labs(title="Air-sac v Malpighian Tubules",x= element_blank(),y= element_blank()) + 
scale_color_gradientn(colors = met.brewer("VanGogh3"))+ 
theme_few()

chao_MA_GI=ggplot(chao_deGI5_DF, aes(baseMean, log2FoldChange, colour=padj)) + 
geom_point(size=1,alpha=0.5) + 
scale_y_continuous(oob=squish) + 
scale_x_log10()+ 
labs(title="Air-sac v Midgut",x= element_blank(),y= element_blank()) + 
scale_color_gradientn(colors = met.brewer("VanGogh3")) + 
theme_few()

chao_MA_patch <- chao_MA_MT + chao_MA_GI + plot_layout(guides = "collect")
chao_MA <- patchwork::patchworkGrob(chao_MA_patch)
gridExtra::grid.arrange(chao_MA, left = "Log Fold Change", bottom = "Mean of Normalized Counts")
chao_MA
ggsave("chao_MAplot.png",plot=last_plot(),bg="transparent")

#volcano plot 
#plot in ggplot 
chao_volcano_MT <- ggplot(data=chao_deMT5_2repDF, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed5)) + geom_point(alpha=0.75) + theme_few() + 
geom_vline(xintercept=c(-5, 5), col="#8d642c") +
geom_hline(yintercept=-log10(0.05), col="#8d642c")+
labs(title="Air-sac v Malpighian Tubules",x= element_blank(),y= element_blank()) + 
theme(plot.title = element_text(size = 10, family = "Tahoma", face = "bold"),
  text = element_text(size = 10, family = "Tahoma"),
  axis.title = element_text(face="bold"),
  axis.text.x=element_text(size = 10),
  axis.text.y=element_text(size = 10),
  strip.text = element_text(face="bold"))+
scale_color_manual(name="",labels=c("Up","Down","No"), values=c("#8b6e28","#b7751f","#C9D5FD"),limits=c("UP","DOWN","NO"))+
ylim(NA,130)+
xlim(-30,30)+
theme(legend.position=c(.89,.92))+
theme(legend.background = element_rect(fill = "transparent"),legend.key = element_rect(fill = "transparent", color = NA))

chao_volcano_GI <- ggplot(data=chao_deGI5_2repDF, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed5)) + geom_point(alpha=0.75) + theme_few() + 
geom_vline(xintercept=c(-5, 5), col="#8d642c") +
geom_hline(yintercept=-log10(0.05), col="#8d642c")+
labs(title="Air-sac v Midut",x= element_blank(),y= element_blank()) + 
theme(plot.title = element_text(size = 10, family = "Tahoma", face = "bold"),
  text = element_text(size = 10, family = "Tahoma"),
  axis.title = element_text(face="bold"),
  axis.text.x=element_text(size = 10),
  axis.text.y=element_text(size = 10),
  strip.text = element_text(face="bold"))+
scale_color_manual(name="",labels=c("Up","Down","No"), values=c("#8b6e28","#d4ba66","#C9D5FD"),limits=c("UP","DOWN","NO"))+
ylim(NA,130)+
xlim(-30,30)+
theme(legend.position=c(.89,.92))+
theme(legend.background = element_rect(fill = "transparent"),legend.key = element_rect(fill = "transparent", color = NA))

chao_volcano <- chao_volcano_MT + chao_volcano_GI 
ggsave("chaoborus_volcano_scaleMatch_legendMv.png",plot=last_plot(),bg="transparent")

##Mochlonyx
#MA plot 
#plot in ggplot 
moch_MA_MT=ggplot(moch_deMT5_DF, aes(baseMean, log2FoldChange, color=padj)) + 
geom_point(size=1,alpha=0.5) + 
scale_y_continuous(oob=squish) + 
scale_x_log10() + 
labs(title="Air-sac v Malpighian Tubules",x= element_blank(),y= element_blank()) + 
scale_color_gradientn(colors = met.brewer("VanGogh3"))+ 
theme_few()

moch_MA_GI=ggplot(moch_deGI5_DF, aes(baseMean, log2FoldChange, colour=padj)) + 
geom_point(size=1,alpha=0.5) + 
scale_y_continuous(oob=squish) + 
scale_x_log10()+ 
labs(title="Air-sac v Midgut",x= element_blank(),y= element_blank()) + 
scale_color_gradientn(colors = met.brewer("VanGogh3")) + 
theme_few()

moch_MA_patch <- moch_MA_MT + moch_MA_GI + plot_layout(guides = "collect")
moch_MA <- patchwork::patchworkGrob(moch_MA_patch)
gridExtra::grid.arrange(moch_MA, left = "Log Fold Change", bottom = "Mean of Normalized Counts")
moch_MA
ggsave("chao_MAplot.png",plot=last_plot(),bg="transparent")

#volcano plot 
#plot in ggplot 
moch_volcano_MT <- ggplot(data=moch_deMT5_2repDF, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed5)) + geom_point(alpha=0.75) + theme_few() + 
geom_vline(xintercept=c(-5, 5), col="#8d642c") +
geom_hline(yintercept=-log10(0.05), col="#8d642c")+
labs(title="Air-sac v Malpighian Tubules",x= element_blank(),y= element_blank()) + 
theme(plot.title = element_text(size = 10, family = "Tahoma", face = "bold"),
  text = element_text(size = 10, family = "Tahoma"),
  axis.title = element_text(face="bold"),
  axis.text.x=element_text(size = 10),
  axis.text.y=element_text(size = 10),
  strip.text = element_text(face="bold"))+
scale_color_manual(name="",labels=c("Up","Down","No"), values=c("#8b6e28","#b7751f","#C9D5FD"),limits=c("UP","DOWN","NO"))+
ylim(NA,130)+
xlim(-30,30)+
theme(legend.position=c(.89,.92))+
theme(legend.background = element_rect(fill = "transparent"),legend.key = element_rect(fill = "transparent", color = NA))

moch_volcano_GI <- ggplot(data=moch_deGI5_2repDF, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed5)) + geom_point(alpha=0.75) + theme_few() + 
geom_vline(xintercept=c(-5, 5), col="#8d642c") +
geom_hline(yintercept=-log10(0.05), col="#8d642c")+
labs(title="Air-sac v Midut",x= element_blank(),y= element_blank()) + 
theme(plot.title = element_text(size = 10, family = "Tahoma", face = "bold"),
  text = element_text(size = 10, family = "Tahoma"),
  axis.title = element_text(face="bold"),
  axis.text.x=element_text(size = 10),
  axis.text.y=element_text(size = 10),
  strip.text = element_text(face="bold"))+
scale_color_manual(name="",labels=c("Up","Down","No"), values=c("#8b6e28","#d4ba66","#C9D5FD"),limits=c("UP","DOWN","NO"))+
ylim(NA,130)+
xlim(-30,30)+
theme(legend.position=c(.89,.92))+
theme(legend.background = element_rect(fill = "transparent"),legend.key = element_rect(fill = "transparent", color = NA))

moch_volcano <- moch_volcano_MT + moch_volcano_GI 
ggsave("mochlonyx_volcano_scaleMatch_legendMv.png",plot=last_plot(),bg="transparent")


##Eucorethra 
#MA plot 
#plot in ggplot 
euc_MA_MT=ggplot(euc_deMT5_DF, aes(baseMean, log2FoldChange, color=padj)) + 
geom_point(size=1,alpha=0.5) + 
scale_y_continuous(oob=squish) + 
scale_x_log10() + 
labs(title="Trachea v Malpighian Tubules",x= element_blank(),y= element_blank()) + 
scale_color_gradientn(colors = met.brewer("VanGogh3"))+ 
theme_few()

euc_MA_GI=ggplot(euc_deGI5_DF, aes(baseMean, log2FoldChange, colour=padj)) + 
geom_point(size=1,alpha=0.5) + 
scale_y_continuous(oob=squish) + 
scale_x_log10()+ 
labs(title="Trachea v Midgut",x= element_blank(),y= element_blank()) + 
scale_color_gradientn(colors = met.brewer("VanGogh3")) + 
theme_few()

euc_MA_patch <- euc_MA_MT + euc_MA_GI + plot_layout(guides = "collect")
euc_MA <- patchwork::patchworkGrob(euc_MA_patch)
gridExtra::grid.arrange(euc_MA, left = "Log Fold Change", bottom = "Mean of Normalized Counts")
euc_MA
ggsave("eucMAplot.png",plot=last_plot(),bg="transparent")

#volcano plot 
#plot in ggplot 
euc_volcano_MT <- ggplot(data=euc_deMT5_DF, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed5)) + geom_point(alpha=0.75) + theme_few() + 
geom_vline(xintercept=c(-5, 5), col="#8d642c") +
geom_hline(yintercept=-log10(0.05), col="#8d642c")+
labs(title="Trachea v Malpighian Tubules",x= element_blank(),y= element_blank()) + 
theme(plot.title = element_text(size = 10, family = "Tahoma", face = "bold"),
  text = element_text(size = 10, family = "Tahoma"),
  axis.title = element_text(face="bold"),
  axis.text.x=element_text(size = 10),
  axis.text.y=element_text(size = 10),
  strip.text = element_text(face="bold"))+
scale_color_manual(name="",labels=c("Up","Down","No"), values=c("#8b6e28","#b7751f","#C9D5FD"),limits=c("UP","DOWN","NO"))+
ylim(NA,130)+
xlim(-30,30)+
theme(legend.position=c(.89,.92))+
theme(legend.background = element_rect(fill = "transparent"),legend.key = element_rect(fill = "transparent", color = NA))

euc_volcano_GI <- ggplot(data=euc_deGI5_DF, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed5)) + geom_point(alpha=0.75) + theme_few() + 
geom_vline(xintercept=c(-5, 5), col="#8d642c") +
geom_hline(yintercept=-log10(0.05), col="#8d642c")+
labs(title="Trachea v Midut",x= element_blank(),y= element_blank()) + 
theme(plot.title = element_text(size = 10, family = "Tahoma", face = "bold"),
  text = element_text(size = 10, family = "Tahoma"),
  axis.title = element_text(face="bold"),
  axis.text.x=element_text(size = 10),
  axis.text.y=element_text(size = 10),
  strip.text = element_text(face="bold"))+
scale_color_manual(name="",labels=c("Up","Down","No"), values=c("#8b6e28","#d4ba66","#C9D5FD"),limits=c("UP","DOWN","NO"))+
ylim(NA,130)+
xlim(-30,30)+
theme(legend.position=c(.89,.92))+
theme(legend.background = element_rect(fill = "transparent"),legend.key = element_rect(fill = "transparent", color = NA))

euc_volcano <- euc_volcano_MT + euc_volcano_GI 
ggsave("eucorethra_volcano_scaleMatch_legendMv.png",plot=last_plot(),bg="transparent")

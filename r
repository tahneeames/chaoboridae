#R script for Chaoboridae 'omics project
#multiple steps throughout this script work in tandem with the .sh script file, but all necessary files are uploaded onto data dryad (at some point)
setwd("~/Desktop/rDirectory/chaoboridae")
#graphing contig lengths from trinity stats output generated in Unix, comparing an assembly done using trimmomatic to an untrimmed assembly (Trinity)
library(ggplot2)
library(ggthemes)
#reading in the files 
chaoUntrimmed=read.table("~/Desktop/Chaoborid/resultsFiles/chaoborus/untrimmed/chaoborus.untrimmed_contigLength.txt") #nrow 662288
chaoTrimmed=read.table("~/Desktop/Chaoborid/resultsFiles/chaoborus/trimmed/chaoborus.trimmed_contigLength.txt") #nrow 602688 
#combine the data frames
chao_contigLength=data.frame(Trimmed=c(chaoTrimmed,rep(NA,662288-length(chaoTrimmed))), Untrimmed=c(chaoUntrimmed,rep(NA,662288-length(chaoUntrimmed))))
chao_contigLength=data.frame(Trimmed=c(chaoTrimmed,rep(NA,59600)), Untrimmed=c(chaoUntrimmed))
chaoTrimmed=cbind(rownames(chaoTrimmed),chaoTrimmed)
colnames(chaoTrimmed)<-c('contig','Trimmed')
chaoUntrimmed=cbind(rownames(chaoUntrimmed),chaoUntrimmed)
colnames(chaoUntrimmed)<-c('contig','Untrimmed')
chao_contigLength<-qpcR:::cbind.na(chaoTrimmed,chaoUntrimmed)
run=c(rep("trimmed",602688),rep("untrimmed",662288))
#merge into new data frame with contig length and assembly 
colnames(chao_contig)=c("readLength","run")
chao_contig=cbind(chao_contig,run)

#plot 
chao_contLen=ggplot(data=chao_contig,aes(x = run, y = readLength,fill=run))+
geom_violin(color="#ab8060",trim=FALSE)+ 
labs(x = "Run", y = "Contig Length")+
scale_fill_manual(values=c("#cd9e50","#d5a069"))+
theme(legend.position ="none")+
theme(legend.title = element_text(size = 4, face = "bold"))+
theme_tufte()+
theme(plot.title = element_text(size = 15, family = "Tahoma", face = "bold"),text = element_text(size = 14, family = "Tahoma",color="#BF8B67",),axis.title = element_text(color="#BF8B67",face="bold"),axis.text.x=element_text(color="#BF8B67",face = "bold",size = 10),axis.text.y=element_text(color="#BF8B67",size = 10),strip.text = element_text(color="#BF8B67",face="bold"), axis.ticks = element_blank())+
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill = "#FFF4CC",colour = "#BF8B67", size = 1.5, linetype = "solid"))+ 
theme(plot.background = element_rect(fill = "#FFF4CC"))+
theme(legend.position ="none")
#save
ggsave("contLeng_boldBG.png",plot=last_plot())

#trimmomatic assembly seems to be better - going with that method for eucorethra 
#contig length chaoborus v eucorethra 
eucorethraAssembly=read.table("~/Desktop/Chaoborid/resultsFiles/eucoretrha/trimmed/eucorethra.trimmed_contigLength.txt") #259551
chaoborusAssembly=read.table("~/Desktop/Chaoborid/resultsFiles/chaoborus/trimmed/chaoborus.trimmed_contigLength.txt") #602688
#add column with length identity
sp=rep(c("Eucorethra underwoodi"),times=c(259551))
eucorethraAssembly= cbind(eucorethraAssembly,sp)
sp=rep(c("Chaoborus trivittatus"),times=c(602688))
chaoborusAssembly= cbind(chaoborusAssembly,sp)
#make new data frame with both of these contig lengths
assemblyStat=rbind(eucorethraAssembly,chaoborusAssembly)

#violin plots in ggplot2
assemb_vi=ggplot(data=assemblyStat,aes(x = sp, y = log10(V1),fill=sp))+
geom_violin(color="#ab8060",trim=FALSE)+ 
labs(x = "Species", y = "Contig Length")+
ylim(NA,1000)+
scale_fill_manual(values=c("#cd9e50","#d5a069"))+
theme(legend.position ="none")+
theme(legend.title = element_text(size = 4, face = "bold"))+
theme_tufte()+
theme(plot.title = element_text(size = 15, family = "Tahoma", face = "bold"),text = element_text(size = 14, family = "Tahoma",color="#BF8B67",),axis.title = element_text(color="#BF8B67",face="bold"),axis.text.x=element_text(color="#BF8B67",face = "bold",size = 10),axis.text.y=element_text(color="#BF8B67",size = 10),strip.text = element_text(color="#BF8B67",face="bold"), axis.ticks = element_blank())+
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill = "#FFF4CC",colour = "#BF8B67", size = 1.5, linetype = "solid"))+ 
theme(plot.background = element_rect(fill = "#FFF4CC"))+
theme(legend.position ="none")


##data for targeted query
library(ggplot2)
library(ggthemes)
library(gridExtra)
library(patchwork)
library(MetBrewer)
library(tidyverse)

# load BLAST summary results
geneContent=read.table("~/Desktop/Chaoborid/resultsFiles/chaoborus/trimmed/chaoborus.trimmed_targetedQuery.txt",header=T,sep="\t")  

# plot identity of BLAST hit (x-axis) and length of hit divided by length of query (y-axis)
id_by_length <- ggplot(geneContent,aes(x=percIdent,y=percSeq,fill=class,color=class)) + 
geom_point(shape=21,alpha=0.65,size=2.3) + 
labs(x="Percent identity of BLAST hit",y="BLAST hit length / query length") + 
scale_color_manual(values=met.brewer('Pillement',14),
  name ="Class",
  breaks=c("carbonicAnhydrase","channel_Cl","channel_K","channel_water","cuticularProtein","DEGENaC_pickpocket","exchanger_CPA","exchanger_other","NaKatpase_subunit","receptor_CAPA","receptor_leukokinin","receptor_serotonin","Resilin","vatpaseSubunit"),
  labels=c("Carbonic anhydrase",expression(paste("Cl"^"-","channel")),expression(paste("K"^"+","channel")),"Water channel","Cuticular protein","Pickpocket gene","CPA exchanger","Other exchanger",expression(paste("Na"^"+","K"^"+","ATPase subunit")),"CAPA receptor","Leukokinin receptor","Serotonin receptor","Resilin",expression(paste("Vacuolar H"^"+","ATPase subunit"))))+
scale_fill_manual(values=met.brewer('Pillement',14),
  name ="Class",
  breaks=c("carbonicAnhydrase","channel_Cl","channel_K","channel_water","cuticularProtein","DEGENaC_pickpocket","exchanger_CPA","exchanger_other","NaKatpase_subunit","receptor_CAPA","receptor_leukokinin","receptor_serotonin","Resilin","vatpaseSubunit"),
  labels=c("Carbonic anhydrase",expression(paste("Cl"^"-","channel")),expression(paste("K"^"+","channel")),"Water channel","Cuticular protein","Pickpocket gene","CPA exchanger","Other exchanger",expression(paste("Na"^"+","K"^"+","ATPase subunit")),"CAPA receptor","Leukokinin receptor","Serotonin receptor","Resilin",expression(paste("Vacuolar H"^"+","ATPase subunit"))))+
theme_few() +
  theme(
  text = element_text(size = 10, family = "Tahoma"),
  axis.title = element_text(face="bold"),
  axis.text.x=element_text(size = 10),
  axis.text.y=element_text(size = 10),
  strip.text = element_text(face="bold"))+ 
 theme(legend.title = element_text(
      family = "Tahoma",size = 9),legend.text=element_text(family = "Tahoma",size = 8))+
 theme(legend.text.align = 0)

# generate histograms of percent identity and of length of hit divided by length of query
id_hist <- ggplot(data=geneContent,aes(x=percIdent)) + 
geom_histogram(binwidth=5,color="#A1A18A",fill="#A1A18A",alpha=0.7) + 
labs(x="Percent identity of BLAST hit",y=element_blank())+
theme_few()+
theme(plot.title = element_text(size = 12, family = "Tahoma", face = "bold"),
  text = element_text(size = 10, family = "Tahoma"),
  axis.title = element_text(face="bold"),
  axis.text.x=element_text(size = 10),
  axis.text.y=element_text(size = 10),
  strip.text = element_text(face="bold"))

length_hist <- ggplot(data=geneContent,aes(x=percSeq)) + 
geom_histogram(binwidth=.05,color="#e3b577",fill="#e3b577",alpha=0.7) + 
labs(x="BLAST hit length / query length",y=element_blank())+
theme_few()+
theme(plot.title = element_text(size = 12, family = "Tahoma", face = "bold"),
  text = element_text(size = 10, family = "Tahoma"),
  axis.title = element_text(face="bold"),
  axis.text.x=element_text(size = 10),
  axis.text.y=element_text(size = 10),
  strip.text = element_text(face="bold"))

# generate multi-panel plot
(id_hist| length_hist) /
      id_by_length  + plot_layout(guides = 'collect') & theme(legend.position = c(1,0))

ggsave("chaob_targetQuery.png",plot=last_plot(),bg="transparent")

(id_hist| length_hist) /
      id_by_length 

ggsave("chaob_targetQuery_legendDown.png",plot=last_plot(),bg="transparent")

#wrapping them with the legend only on the scatterplot
library(cowplot)

id_by_lengthCP <- ggplot(geneContent,aes(x=percIdent,y=percSeq,fill=class,color=class)) + 
geom_point(shape=21,alpha=0.65,size=2.3) + 
labs(x="Percent identity of BLAST hit",y="BLAST hit length / query length") + 
scale_color_manual(values=met.brewer('Pillement',14),
  name ="Class",
  breaks=c("carbonicAnhydrase","channel_Cl","channel_K","channel_water","cuticularProtein","DEGENaC_pickpocket","exchanger_CPA","exchanger_other","NaKatpase_subunit","receptor_CAPA","receptor_leukokinin","receptor_serotonin","Resilin","vatpaseSubunit"),
  labels=c("Carbonic anhydrase",expression(paste("Cl"^"-","channel")),expression(paste("K"^"+","channel")),"Water channel","Cuticular protein","Pickpocket gene","CPA exchanger","Other exchanger",expression(paste("Na"^"+","K"^"+","ATPase subunit")),"CAPA receptor","Leukokinin receptor","Serotonin receptor","Resilin",expression(paste("Vacuolar H"^"+","ATPase subunit"))))+
scale_fill_manual(values=met.brewer('Pillement',14),
  name ="Class",
  breaks=c("carbonicAnhydrase","channel_Cl","channel_K","channel_water","cuticularProtein","DEGENaC_pickpocket","exchanger_CPA","exchanger_other","NaKatpase_subunit","receptor_CAPA","receptor_leukokinin","receptor_serotonin","Resilin","vatpaseSubunit"),
  labels=c("Carbonic anhydrase",expression(paste("Cl"^"-","channel")),expression(paste("K"^"+","channel")),"Water channel","Cuticular protein","Pickpocket gene","CPA exchanger","Other exchanger",expression(paste("Na"^"+","K"^"+","ATPase subunit")),"CAPA receptor","Leukokinin receptor","Serotonin receptor","Resilin",expression(paste("Vacuolar H"^"+","ATPase subunit"))))+
theme_few() +
  theme(
  text = element_text(size = 10, family = "Tahoma"),
  axis.title = element_text(face="bold"),
  axis.text.x=element_text(size = 10),
  axis.text.y=element_text(size = 10),
  strip.text = element_text(face="bold"))+ 
 theme(legend.title = element_text(
      family = "Tahoma",size = 4),legend.text=element_text(family = "Tahoma",size = 3))+
 theme(legend.text.align = 0)

id_grid <- plot_grid(id_hist,length_hist)
all_grid <- plot_grid(id_grid,id_by_length,ncol=1)

ggsave("targetdGene_query_cowplot.png",plot=last_plot(),bg="transparent")

#same plots but for eucorethra 
geneContent=read.table("~/Desktop/Chaoborid/resultsFiles/eucorethra/trimmed/eucorethra.trimmed_targetedQuery.txt",header=T,sep="\t")  

#plot identity of BLAST hit (x-axis) and length of hit divided by length of query (y-axis)
id_by_length <- ggplot(geneContent,aes(x=percIdent,y=percSeq,fill=class,color=class)) + 
geom_point(shape=21,alpha=0.65,size=2.3) + 
labs(x="Percent identity of BLAST hit",y="BLAST hit length / query length") + 
scale_color_manual(values=met.brewer('Pillement',16),
  name ="Class",
  breaks=c("carbonicAnhydrase","channel_Cl","channel_K","channel_water","cuticularProtein","DEGENaC_pickpocket","exchanger_CPA","exchanger_other","NaKatpase_subunit","receptor_CAPA","receptor_leukokinin","receptor_serotonin","Resilin","vatpaseSubunit","V0_domain"),
  labels=c("Carbonic anhydrase",expression(paste("Cl"^"-","channel")),expression(paste("K"^"+","channel")),"Water channel","Cuticular protein","Pickpocket gene","CPA exchanger","Other exchanger",expression(paste("Na"^"+","K"^"+","ATPase subunit")),"CAPA receptor","Leukokinin receptor","Serotonin receptor","Resilin",expression(paste("Vacuolar H"^"+","ATPase subunit")),"V0 domain"))+
scale_fill_manual(values=met.brewer('Pillement',16),
  name ="Class",
  breaks=c("carbonicAnhydrase","channel_Cl","channel_K","channel_water","cuticularProtein","DEGENaC_pickpocket","exchanger_CPA","exchanger_other","NaKatpase_subunit","receptor_CAPA","receptor_leukokinin","receptor_serotonin","Resilin","vatpaseSubunit"),
  labels=c("Carbonic anhydrase",expression(paste("Cl"^"-","channel")),expression(paste("K"^"+","channel")),"Water channel","Cuticular protein","Pickpocket gene","CPA exchanger","Other exchanger",expression(paste("Na"^"+","K"^"+","ATPase subunit")),"CAPA receptor","Leukokinin receptor","Serotonin receptor","Resilin",expression(paste("Vacuolar H"^"+","ATPase subunit"))))+
theme_few() +
  theme(
  text = element_text(size = 10, family = "Tahoma"),
  axis.title = element_text(face="bold"),
  axis.text.x=element_text(size = 10),
  axis.text.y=element_text(size = 10),
  strip.text = element_text(face="bold"))+ 
 theme(legend.title = element_text(
      family = "Tahoma",size = 4),legend.text=element_text(family = "Tahoma",size = 3))+
 theme(legend.text.align = 0)

#generate histograms of percent identity and of length of hit divided by length of query
id_hist <- ggplot(data=geneContent,aes(x=percIdent)) + 
geom_histogram(binwidth=5,color="#A1A18A",fill="#A1A18A",alpha=0.7) + 
labs(x="Percent identity of BLAST hit",y=element_blank())+
theme_few()+
theme(plot.title = element_text(size = 12, family = "Tahoma", face = "bold"),
  text = element_text(size = 10, family = "Tahoma"),
  axis.title = element_text(face="bold"),
  axis.text.x=element_text(size = 10),
  axis.text.y=element_text(size = 10),
  strip.text = element_text(face="bold"))

length_hist <- ggplot(data=geneContent,aes(x=percSeq)) + 
geom_histogram(binwidth=.05,color="#e3b577",fill="#e3b577",alpha=0.7) + 
labs(x="BLAST hit length / query length",y=element_blank())+
theme_few()+
theme(plot.title = element_text(size = 12, family = "Tahoma", face = "bold"),
  text = element_text(size = 10, family = "Tahoma"),
  axis.title = element_text(face="bold"),
  axis.text.x=element_text(size = 10),
  axis.text.y=element_text(size = 10),
  strip.text = element_text(face="bold"))

#compile into one plot
id_grid <- plot_grid(id_hist,length_hist)
all_grid <- plot_grid(id_grid,id_by_length,ncol=1)
#save
ggsave("targetdGene_query_cowplot_eucoretrha.png",plot=last_plot(),bg="transparent")

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
chaoCount_matrix <- chaoCount_matrix[, c(2,5,9,1,3,7,4,6,8)]
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
chao_3pca=plotPCA(chao_rld_2reps,intgroup=c("tissue","replicate"),returnData=TRUE)
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
chao_metadata_2reps$rep=as.factor(chao_metadata_2reps$rep)
#make deseq2 dataframe from these dataframes
chao_dds_2rep <- DESeqDataSetFromMatrix(countData=round(chao_geneCounts_2reps_matrix), colData=chao_metadata_2reps, design= ~ tissue)
#remove values less than 5
chao_dds5_2reps <- chao_dds_2rep[rowSums(counts(chao_dds_2rep)) >= 5,]
#run deseq2
chao_deSeq_2reps=DESeq(chao_dds5_2reps)

chao_rld_2reps<-rlog(chao_deSeq_2reps,blind=TRUE)

#checking pca again 
chao_2pca=plotPCA(chao_rld_2reps,intgroup=c("tissue","replicate"),returnData=TRUE)
chao_pca_2reps=ggplot(chao_2pca,aes(PC1,PC2,col=tissue,shape=replicate,fill=tissue))+
geom_point(size=3)+
scale_color_manual(name="Tissue",labels=c("Air sac","Midgut","Malpighian tubule"),values=c("#8b6e28","#d4ba66","#b7751f"),limits=c("airSac","malTube","gut"))+
scale_fill_manual(name="Tissue",labels=c("Air sac","Midgut","Malpighian tubule"),values=c("#8b6e28","#d4ba66","#b7751f"),limits=c("airSac","malTube","gut"))+
scale_shape_manual(name="Replicate",labels=c("1","3"),values=c(21,22,24),limits=c("1","3"))+
theme_few()+
theme(plot.title = element_text(size = 12, family = "Tahoma", face = "bold"),
  text = element_text(size = 10, family = "Tahoma"),
  axis.title = element_text(face="bold",size=12),
  axis.text.x=element_text(size = 10),
  axis.text.y=element_text(size = 10),
  strip.text = element_text(face="bold"))+
theme(legend.title = element_text(size = 10, face = "bold"),text=element_text(size = 10, family = "Tahoma"))


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

#make some plots before moving on to more data wrangling and viz 
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
chao_volcano_MT <- ggplot(data=chao_deMT5_DF, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed5)) + geom_point(alpha=0.75) + theme_few() + 
geom_vline(xintercept=c(-2, 2), col="#8d642c") +
geom_hline(yintercept=-log10(0.05), col="#8d642c")+
labs(title="Air-sac v Malpighian Tubules",x= element_blank(),y= element_blank()) + 
theme(plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
  text = element_text(size = 10, family = "Tahoma"),
  axis.title = element_text(face="bold"),
  axis.text.x=element_text(size = 10),
  axis.text.y=element_text(size = 10),
  strip.text = element_text(face="bold"))+
scale_color_manual(name="Differently \nExpressed",labels=c("Up","Down","No"), values=c("#E2BE62","#DFD1A6","#C9D5FD"),limits=c("UP","DOWN","NO"))+
ylim(NA,130)+
theme(legend.position="none")

chao_volcano_GI <- ggplot(data=chao_deGI5_DF, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed5)) + geom_point(alpha=0.75) + theme_few() + 
geom_vline(xintercept=c(-2, 2), col="#8d642c") +
geom_hline(yintercept=-log10(0.05), col="#8d642c")+
labs(title="Air-sac v Midut",x= element_blank(),y= element_blank()) + 
theme(plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
  text = element_text(size = 10, family = "Tahoma"),
  axis.title = element_text(face="bold"),
  axis.text.x=element_text(size = 10),
  axis.text.y=element_text(size = 10),
  strip.text = element_text(face="bold"))+
scale_color_manual(name="Differently \nExpressed",labels=c("Up","Down","No"), values=c("#E2BE62","#DFD1A6","#C9D5FD"),limits=c("UP","DOWN","NO"))+
ylim(NA,130)

chao_volcano <- chao_volcano_MT + chao_volcano_GI + plot_layout(guides = "collect")
ggsave("chaoborus_volcano.png",plot=last_plot(),bg="transparent")


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
geom_vline(xintercept=c(-2, 2), col="#8d642c") +
geom_hline(yintercept=-log10(0.05), col="#8d642c")+
labs(title="Trachea v Malpighian Tubules",x= element_blank(),y= element_blank()) + 
theme(plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
  text = element_text(size = 10, family = "Tahoma"),
  axis.title = element_text(face="bold"),
  axis.text.x=element_text(size = 10),
  axis.text.y=element_text(size = 10),
  strip.text = element_text(face="bold"))+
scale_color_manual(name="Differently \nExpressed",labels=c("Up","Down","No"), values=c("#E2BE62","#DFD1A6","#C9D5FD"),limits=c("UP","DOWN","NO"))+
ylim(NA,130)+
theme(legend.position="none")

euc_volcano_GI <- ggplot(data=euc_deGI5_DF, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed5)) + geom_point(alpha=0.75) + theme_few() + 
geom_vline(xintercept=c(-2, 2), col="#8d642c") +
geom_hline(yintercept=-log10(0.05), col="#8d642c")+
labs(title="Trachea v Midut",x= element_blank(),y= element_blank()) + 
theme(plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
  text = element_text(size = 10, family = "Tahoma"),
  axis.title = element_text(face="bold"),
  axis.text.x=element_text(size = 10),
  axis.text.y=element_text(size = 10),
  strip.text = element_text(face="bold"))+
scale_color_manual(name="Differently \nExpressed",labels=c("Up","Down","No"), values=c("#E2BE62","#DFD1A6","#C9D5FD"),limits=c("UP","DOWN","NO"))+
ylim(NA,130)

euc_volcano <- euc_volcano_MT + euc_volcano_GI + plot_layout(guides = "collect")
ggsave("eucorethra_volcano.png",plot=last_plot(),bg="transparent")


##data wrangling and data merging 
### read in a trinotate files
chaoTrinotate=read.csv2("~/Desktop/Chaoborid/resultsFiles/chaoborus/trimmed/chaoborus.trimmed.Trinotate.txt",sep="\t")
eucTrinotate=read.csv2("~/Desktop/Chaoborid/resultsFiles/eucorethra/trimmed/eucorethra.trimmed.Trinotate.txt",sep="\t")

#filter trinotate output by DE gene list (see .sh file for trinotate generation)
#first rename the rows to merge in both dataframes for both organisms 
names(chaoTrinotate)[1] <- "gene_id"
names(chao_as_up_Join5lfc)[1] <- "gene_id"
names(eucTrinotate)[1] <- "gene_id"
names(euc_tr_up_Join5lfc)[1] <- "gene_id"

#filter by that row (gene identity in this case)
filter(chaoTrinotate,gene_id %in% c(chao_as_up_Join5lfc$gene_id)) -> chaoborusUpreg_annotate
filter(eucTrinotate,gene_id %in% c(euc_tr_up_Join5lfc$gene_id)) -> eucorethraUpreg_annotate
#export dataframes
write.csv(chaoborusUpreg_annotate,file="trivittatus_asUp_trinotate.csv")
write.csv(eucorethraUpreg_annotate,file="eucorethra_trUp_trinotate.csv")

#preserve both dataframes so I can sort by lfc values 
chao_UP_5lfc_trinotate=merge(x = chao_as_up_Join5lfc, y = chaoTrinotate, by = "gene_id")
euc_UP_5lfc_trinotate=merge(x = euc_tr_up_Join5lfc, y = eucTrinotate, by = "gene_id")
#export
write.csv(chao_UP_5lfc_trinotate,file="trivittatus_asUp_trinotateFULL.csv")
write.csv(euc_UP_5lfc_trinotate,file="eucorethra_trUp_trinotateFULL.csv")

#sort them by lfc value 
chao_orderUP=arrange(chao_UP_5lfc_trinotate,log2FoldChange.x,log2FoldChange.y)
euc_orderUP=arrange(euc_UP_5lfc_trinotate,log2FoldChange.x,log2FoldChange.y)
#export these too 
write.csv(chao_orderUP,file="trivittatus_asUp_trinotateFULL_lfcOrder.csv")
write.csv(euc_orderUP,file="eucorethra_trUp_trinotateFULL_lfcOrder.csv")

#making another one with no cutoff value, just the raw lfc values to see what will happen 
#merge full pairwise dataframes from DESeq 
chao_de_merge=merge(chao_deMT5_2repDF, chao_deGI5_2repDF, by=0)
chao_de_merge %>% filter(if_all(contains("p.value"),~ (.x< 0.5)))  -> chao_de_merge_pVal
#eucorethra too 
euc_de_merge=merge(euc_deMT5_DF, euc_deGI5_DF, by=0)
euc_de_merge %>% filter(if_all(contains("p.value"),~ (.x< 0.5)))  -> euc_de_merge_pVal
#then merge with trinotate file '

#rename these ones too just like with the other df
names(chao_de_merge_pVal)[1] <- "gene_id"
names(euc_de_merge_pVal)[1] <- "gene_id"

#merge to have an annotated dataframe with lfc values and gene homologies/ontologies/pathways
chao_all_trinotate=merge(x = chaoTrinotate, y = chao_de_merge_pVal, by = "gene_id")
euc_all_trinotate=merge(x = eucTrinotate, y = euc_de_merge_pVal, by = "gene_id")
#rearrange by lfc 
chao_orderALL=arrange(chao_all_trinotate,log2FoldChange.x,log2FoldChange.y)
euc_orderALL=arrange(euc_all_trinotate,log2FoldChange.x,log2FoldChange.y)
#export 
write.csv(chao_orderALL,file="trivittatus_trinotateFULL_lfcOrder.csv")
write.csv(euc_orderALL,file="eucorethra_trinotateFULL_lfcOrder.csv")

#okay but a lot of the annotations are blank! Let's do more sifting 
#filter only genes that have a blast result 
#for some reason the trinotate output uses little periods for blank spaces - going to remove them so I can filter the data 
is.na(chao_orderALL) <- chao_orderALL == "."
#filter based on presence of blast result
chao_orderALL %>% filter(!if_all(c(sprot_Top_BLASTX_hit, sprot_Top_BLASTP_hit), is.na)) -> chao_order_BLAST
#another based on pfam 
chao_orderALL %>% filter(!is.na(Pfam)) -> chao_order_Pfam

#eucorethra now 
is.na(euc_orderALL) <- euc_orderALL == "."
euc_orderALL %>% filter(!if_all(c(sprot_Top_BLASTX_hit, sprot_Top_BLASTP_hit), is.na)) -> euc_order_BLAST
euc_orderALL %>% filter(!is.na(Pfam)) -> euc_order_Pfam
#save these! 
write.csv(chao_order_BLAST,file="trivittatus_trinotateBLAST_lfcOrder.csv")
write.csv(chao_order_Pfam,file="trivittatus_trinotatePFAM_lfcOrder.csv")
write.csv(euc_order_BLAST,file="eucorethra_trinotateBLAST_lfcOrder.csv")
write.csv(euc_order_Pfam,file="eucorethra_trinotatePFAM_lfcOrder.csv")


#plotting specific genes functions 
chao_isoCount=read.table("~/Desktop/Chaoborid/resultsFiles/chaoborus/trimmed/chaoborus-trimmed.salmon.isoform.TPM.not_cross_norm_2reps.txt",header=T,sep='\t',row.names=1)
#rename columns to match tissue (no rep)
colnames(chao_isoCount)[1:2] <- "trachea"
colnames(chao_isoCount)[3:4] <- "midgut"
colnames(chao_isoCount)[5:6] <- "malTube"
#calc row means for each tissue type into new df 
chao_isoMean <-t(apply(chao_isoCount,1, function(x) tapply(x,colnames(chao_isoCount),mean)))
chao_isoMean=as.data.frame(chao_isoMean)
#make new column for species name
org=rep("chaoborus",times=602688)
chao_isoMean=cbind(chao_isoMean,org)

#doing the same thing with the eucorethra data set - converting to mean values for easy plotting :) 
euc_isoCount=read.table("~/Desktop/Chaoborid/resultsFiles/eucorethra/trimmed/eucorethra-trimmed.salmon.isoform.TPM.not_cross_norm",header=T,sep='\t',row.names=1)
colnames(euc_isoCount)[1:3] <- "midgut"
colnames(euc_isoCount)[4:6] <- "malTube"
colnames(euc_isoCount)[7:9] <- "trachea"
euc_isoMed<-t(apply(euc_isoCount,1, function(x) tapply(x,colnames(euc_isoCount),median)))
euc_isoMean=as.data.frame(euc_isoMean)
org=rep("eucorethra",times=259551)
euc_isoMean=cbind(euc_isoMean,org)

#combine both dataframes 
isoMeans=rbind(chao_isoMean,euc_isoMean)
#move rownames to a new column 
isoMeans$transcriptID <- rownames(isoMeans)
#pivot tissue columns into a new column by count
isoMeans=isoMeans %>% pivot_longer(cols=c('trachea', 'malTube','midgut'),
  names_to='tissue',
  values_to='count')

#flex ggplot code for easy print values based on rowname ???
#this is for one species, comparing transcription levels of various genes 
isoMean_modular <- filter(isoMeans, transcriptID %in% c("TRINITY_DN23587_c0_g1_i7","TRINITY_DN2775_c0_g1_i8","TRINITY_DN2739_c0_g2_i2","TRINITY_DN7923_c0_g1_i1","TRINITY_DN11599_c0_g1_i6","TRINITY_DN292_c2_g1_i4","TRINITY_DN992_c1_g2_i2","TRINITY_DN3016_c0_g1_i33","TRINITY_DN2528_c0_g1_i8","TRINITY_DN6144_c0_g1_i30","TRINITY_DN2528_c0_g1_i6","TRINITY_DN102613_c5_g1_i3","TRINITY_DN19936_c0_g1_i16","TRINITY_DN8940_c3_g1_i1","TRINITY_DN6335_c0_g1_i1"))


#comparing transcription by species, one gene 
isoMean_modular <- filter(isoMeans, transcriptID %in% c("TRINITY_DN177292_c0_g1_i1","TRINITY_DN7614_c0_g1_i91"))
ggplot(isoMean_modular,aes(x=tissue, y=count, group=org)) +
  geom_line(aes(color=org),size=1)+
  geom_point(aes(color=org),size=1.6)+
scale_color_manual(values=c("#5F6247","#B49355"),
  name ="Species",
  breaks=c("eucorethra","trivittatus"),
  labels=c("Eucorethra \nunderwoodi","Chaoborus \ntrivittatus"))+
scale_x_discrete(name =element_blank(),labels=c("Midgut","Malpighian tubules","Tracheae"))+
ylab("Median TPM")+
theme_classic()+
theme(text = element_text(size = 11, family = "sans",face="bold"),axis.title = element_text(face="bold"),axis.text.x=element_text(size = 9),axis.text.y=element_text(size = 9),legend.text=element_text(size=9,family="sans",face="italic"))+
theme(legend.background = element_rect(fill = "transparent"),legend.key = element_rect(fill = "transparent", color = NA))+
theme(legend.justification=c(1,1), legend.position=c(.24,.94))+
theme(panel.background = element_rect(fill = "transparent",colour = "transparent"))+
theme(plot.background = element_rect(fill = "transparent"))+
 theme(legend.key = element_rect(color = NA, fill = NA),
        legend.key.size = unit(1, "cm")) +
  theme(legend.title.align = 0.5)


#keeping separate by replicate 

#combine both dataframes 
#rename columns to match
chao_isoCount <- rename(chao_isoCount,c("tr_rep1" = "as_rep1", "tr_rep2" = "as_rep2","tr_rep3"="as_rep3"))
isoAll=rbind(chao_isoCount,euc_isoCount)
#move rownames to a new column 
isoAll$transcriptID <- rownames(isoAll)
#pivot tissue columns into a new column by count
isoAll=isoAll %>% pivot_longer(cols=c('tr_rep1','tr_rep2','tr_rep3', 'mt_rep1','mt_rep2','mt_rep3','gi_rep1','gi_rep2','gi_rep3'),
  names_to='sample',
  values_to='count')
#adding columns 
isoAll <- isoAll %>%
  mutate(tissue = case_when(
    startsWith(sample, "tr") ~ "trachea",
    startsWith(sample, "mt") ~ "malTube",  
    startsWith(sample, "gi") ~ "midgut"
    ))

isoAll <- isoAll %>%
  mutate(replicate = case_when(
    endsWith(sample, "1") ~ "one",
    endsWith(sample, "2") ~ "two",  
    endsWith(sample, "3") ~ "three"
    ))


median <- isoAll %>% group_by(transcriptID, tissue) %>% 
mutate(med = median(count))

isoAll=cbind(isoAll,median[,7])



#extract values from dataframe and plot 
#comparing one gene, both species 
iso_modular <- filter(isoAll, transcriptID %in% c("TRINITY_DN177292_c0_g1_i1","TRINITY_DN22179_c0_g1_i1"))
ggplot(iso_modular) +
  geom_line(aes(x=tissue,y=med,group=transcriptID,color=org),size=1)+
  geom_point(aes(x=tissue,y=count,group=sample,color=org),size=1.6)+
scale_color_manual(values=c("#5F6247","#B49355"),
  name ="Species",
  breaks=c("eucorethra","trivittatus"),
  labels=c("Eucorethra \nunderwoodi","Chaoborus \ntrivittatus"))+
scale_x_discrete(name =element_blank(),labels=c("Midgut","Malpighian tubules","Tracheae"))+
ylab("TPM")+
theme_classic()+
theme(text = element_text(size = 12, face="bold"),axis.title = element_text(face="bold"),axis.text.x=element_text(size = 10,face="bold"),axis.text.y=element_text(size = 8),legend.text=element_text(size=10,face="italic"),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
theme(legend.background = element_rect(fill = "transparent"),legend.key = element_rect(fill = "transparent", color = NA))+
theme(legend.justification=c(1,1), legend.position=c(.24,.94))+
theme(panel.background = element_rect(fill = "transparent",colour = "transparent"))+
theme(plot.background = element_rect(fill = "transparent"))+
 theme(legend.key = element_rect(color = NA, fill = NA),
        legend.key.size = unit(1, "cm")) +
  theme(legend.title.align = 0.5)

#comparing one species, multiple genes 
iso_modular <- filter(isoAll, transcriptID %in% c("TRINITY_DN23587_c0_g1_i7","TRINITY_DN2775_c0_g1_i8","TRINITY_DN2739_c0_g2_i2","TRINITY_DN7923_c0_g1_i1","TRINITY_DN11599_c0_g1_i6","TRINITY_DN292_c2_g1_i4","TRINITY_DN992_c1_g2_i2","TRINITY_DN3016_c0_g1_i33","TRINITY_DN2528_c0_g1_i8","TRINITY_DN6144_c0_g1_i30","TRINITY_DN2528_c0_g1_i6","TRINITY_DN102613_c5_g1_i3","TRINITY_DN19936_c0_g1_i16","TRINITY_DN8940_c3_g1_i1","TRINITY_DN6335_c0_g1_i1"))
ggplot(iso_modular) +
geom_point(aes(x=tissue, y=count,group=transcriptID,color=transcriptID),size=1.6,position=position_dodge(0.2))+
geom_line(aes(x=tissue, y=med,group=interaction(replicate, transcriptID),color=transcriptID))+
scale_color_manual(values=met.brewer("Morgenstern", 15))+
scale_x_discrete(name =element_blank(),labels=c("Midgut","Malpighian tubules","Tracheae"))+
ylab("Median TPM")+
theme_classic()+
theme(legend.position="none")+
theme(text = element_text(size = 12, face="bold"),axis.title = element_text(face="bold"),axis.text.x=element_text(size = 10,face="bold"),axis.text.y=element_text(size = 10))+
theme(panel.background = element_rect(fill = "transparent",colour = "transparent"))+
theme(plot.background = element_rect(fill = "transparent"))



#adding all the targeted genes to make a data table
#merge with count table 
chao_targetCount=merge(x = chao_isoCount, y = chao_target, by = "transcriptID",all.y=TRUE)
euc_targetCount=merge(x = euc_isoCount, y = euc_target, by = "transcriptID",all.y=TRUE)
#dumb 
chao_targetCount <- rename(chao_targetCount,c("percIdent" = "X.ident","percSeq"="X.seq"))
#binding these two df
targetAll=rbind(chao_targetCount,euc_targetCount)
#adding columns 
targetAll <- targetAll %>%
  mutate(tissue = case_when(
    startsWith(sample, "tr") ~ "trachea",
    startsWith(sample, "mt") ~ "malTube",  
    startsWith(sample, "gi") ~ "midgut"
    ))


#adding median values 
#make df with medians 
medians <- targetAll %>% group_by(transcriptID, tissue) %>% 
mutate(med = median(count))
#add med column to original df (probably not necessary)
targetAll=cbind(targetAll,medians[,27])

#making the table! 
#removing replicates, who needs em! 
targetFilter <- targetAll %>% distinct(transcriptID, tissue, .keep_all = TRUE)
#taking just the columns we want in the table :) 
targetLess<-targetFilter[,c(2,6,8,25,26)]

#pivoting tissues 
targetLess=targetLess %>% pivot_wider(names_from = tissue, values_from = med)
targetLess_euc<- targetLess %>% filter(org == "eucorethra" ) 
targetLess_triv<- targetLess %>% filter(org == "chaoborus" ) 

targetWide=merge(x=targetLess_triv,y=targetLess_euc,by="geneTag")

targetWideDF<-targetWide[c(1,3,4,5,6,9,10,11)]

targetWideDF  %>% 
  kable(col.names = c("Gene", "Trachea", "Malpighian tubules", "Midgut")) %>% 
  add_header_above(header = c(" " = 1, "C. trivittatus" = 3, "E. underwoodi" = 3))

targetLess %>% gt(groupname_col="category.x") %>% 
   cols_label(geneTag = "Gene", trachea.x = "Tracnea",midgut.x="Midgut",malTube.x="Malpighian tubules",trachea.x = "Tracnea",midgut.x="Midgut",malTube.x="Malpighian tubules") %>% 
   tab_header(title = md("Median TPM of seretory epithelial tissue, C. trivittatus vs E. underwoodi"))

newTable <- targetWideDF %>% 
mutate(category.x = case_when(
  category.x == "neuropeptide"~"Neurotransmitters",
  category.x == "membraneProt"~"Membrane Proteins",
  category.x == "cuticular"~"Cuticular",
  category.x == "vacuolarATPase"~"Vacuolar ATPAse Subunits")) %>% 
 gt(targetWideDF,groupname_col="category.x") %>%
  tab_header(
    title = "Median TPM of seretory epithelial tissue") %>%
  tab_spanner(
    label = "C. trivittatus",
    id="ct",
    columns = c(trachea.x, malTube.x, midgut.x)
  ) %>%
  tab_spanner(
    label = "E. underwoodi",
    id="eu",
    columns = c(trachea.y, malTube.y, midgut.y)
  ) %>%
   tab_style(
    style =list(cell_text(style = "italic"),cell_borders(color="#989898",sides=c("bottom"), weight = px(2))),
    locations = cells_column_spanners(spanners = "eu")
  ) %>%
   tab_style(
    style =list(cell_text(style = "italic"),cell_borders(color="#989898",sides=c("bottom"), weight = px(2))),
    locations = cells_column_spanners(spanners = "ct")
  )  %>% 
  cols_label(
    geneTag = "Gene",
    trachea.x = "Trachea",
    malTube.x = "Malpighian tubule",
    midgut.x = "Midgut",
    trachea.y = "Trachea",
    malTube.y = "Malpighian tubule",
    midgut.y = "Midgut"
   ) %>% 
  tab_options(table.background.color = "#e6dbc8",
    heading.background.color="#c6c3b3",
    row_group.background.color="#cec5bc",
    table_body.hlines.color="#989898",
    table_body.border.top.color="#989898",
    row_group.border.top.color="#989898",
    table_body.border.bottom.color="#989898",
    row_group.border.bottom.color="#989898",
    heading.border.bottom.color="#989898",
    heading.padding = px(25)) %>% 
  tab_style(
    style = list(
      cell_fill(color = "#DED3A6")
    ),
    locations = cells_body(
      columns = vars(trachea.x, malTube.x, midgut.x))) %>%
   tab_style(
    style = list(
      cell_fill(color = "#d1d3be")
    ),
    locations = cells_body(
      columns = vars(trachea.y, malTube.y, midgut.y)))

#merging trinotate file with TPM file to get transcription per gene 
chao_isoCount$transcript_id <- rownames(chao_isoCount)
euc_isoCount$transcript_id <- rownames(euc_isoCount)
#merge into one dataframe :) 
chao_trinotateTPM=merge(x = chao_isoCount, y = chaoTrinotate, by = "transcript_id")
euc_trinotateTPM=merge(x = euc_isoCount, y = eucTrinotate, by = "transcript_id")
#save! 
write.csv(chao_trinotateTPM,"chao_trinotateTPM.csv")
write.csv(euc_trinotateTPM,"euc_trinotateTPM.csv")


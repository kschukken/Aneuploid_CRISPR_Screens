#### UBE2H in TCGA, DepMap, Yeast, Stingele, Durrbaum  ####
# UBE2H in TCGA_CC.R

### 2025.12.12
### Author: Klaske M. Schukken

### Overview: Get TCGA Data: UBE2H expression with tissue type, 7q copy number
##            also get DepMap data. See UBE2H expression per 7q copy number. correlation with aneuploidy score
##            For Supplementary Figure 7 and parts of Figure 6


### Expression and protein data in aneuploid yeast 
## Yeast Data from: Dephoure et al. 2014, eLife; 3:e03023

## Paired human data: Stingele et al. 2012

## Durrbaum et al. 2014 Protein data 


## !! Update these locations to pathway where data was downloaded: 
# Folder with CRISPR screening and proteomics datasets: 
SchukkenData<- "/Volumes/Schukken_SSD/CRISPR SCREEN PAPER/R Code Data Files/Data"
# Folder with dependency datasets: 
Dependency<- "/Volumes/Schukken_SSD/CRISPR SCREEN PAPER/R Code Data Files/Dependency files"
# Folder with results: 
ResultsFile<- "/Volumes/Schukken_SSD/CRISPR SCREEN PAPER/R Code Data Files/Results"



#### Libraries #####

library(ggplot2) 
library(reshape2)
library(tidyr)
library(tidyverse)
library(readxl)
library(ggpubr)
library("cowplot")
library(plyr)
library('gprofiler2')
library(xlsx)
library("readr")
library("reshape2")
library(survival)
library(ggsurvfit)
library(dplyr)
library(survminer) # ggforest hazard ratio plot

#### Get Data #####
setwd(Dependency)

## TCGA
# 1-s2.0-S1535610818301119-mmc2.xlsx
# TCGA human tumor data, from Alison Taylor et al. (2018) 
# Genomic and Functional Approaches to Understanding Cancer Aneuploidy. 
# Volume 33, Issue 4, 9 April 2018, Pages 676-689.e3
# https://www.sciencedirect.com/science/article/pii/S1535610818301119#app2, Table S2
# I removed the top row of Table S2, which was just the name, so the actual column names would be the headers. 
# Downloaded November 10, 2023
TCGA_Aneuploidy<- read.csv("1-s2.0-S1535610818301119-mmc2.csv", header=TRUE)

# TCGA human tumor data, patient data
# Compiled by Ryan Hagenson, November 2023 
# Merges TCGA data, downloaded November 14, 2023
TCGA_Patient <- read.delim("TCGA-patients.tsv", header=TRUE, sep="\t")


# TCGA human tumor data, Sample data
# Compiled by Ryan Hagenson, November 2023
# Merges TCGA data, downloaded November 14, 2023
TCGA_Sample <- read.delim("TCGA-samples.tsv", header=TRUE, sep="\t")



# TCGA human tumor data, mutational data
# made by Ryan Hagenson, November 2023
# Merges TCGA data, downloaded November 14, 2023
TCGA_Mutation <- read.delim("TCGA-mutations.tsv", header=TRUE, sep="\t")

# Get TP53 Mut vs WT: TCGA
ignoreMutations <- c("Silent","Intron","5'UTR","3'Flank","5'Flank","3'UTR","IGR","RNA")
TCGA_TP53Mut<- subset(TCGA_Mutation, Hugo_Symbol == "TP53" & !(Variant_Classification %in% ignoreMutations))$Tumor_Sample_Barcode
TCGA_TP53WT<- setdiff(TCGA_Sample$SAMPLE_ID, TCGA_TP53Mut) # all in list, remove from list

TCGA_Mutation<- NA # No need to keep giant file. we only needed TP53 mutation



# TCGA RNA expression data for UBE2H 
# Downloaded December 12, 2025  from cBioportal
# mRNA Expression, RSEM (Batch normalized from Illumina HiSeq_RNASeqV2)
# https://www.cbioportal.org/results/download?case_set_id=all&gene_list=UBE2H&cancer_study_list=5c8a7d55e4b046111fee2296&plots_horz_selection=%7B%22dataType%22%3A%22clinical_attribute%22%2C%22selectedDataSourceOption%22%3A%22DFS_MONTHS%22%7D&plots_vert_selection=%7B%22selectedGeneOption%22%3A7328%7D&plots_coloring_selection=%7B%22colorByMutationType%22%3A%22false%22%2C%22colorByCopyNumber%22%3A%22true%22%2C%22colorBySv%22%3A%22false%22%7D
TCGA_RNA.UBE2H<- read.delim("TCGA UBE2H mRNA Expression, RSEM (Batch normalized from Illumina HiSeq_RNASeqV2).txt", header= TRUE, sep="\t")


# transcription cluster 76- metabolism #
# downloaded from cBioportal HPA on January 14, 2025
# genes in cluster 76- non tissue specific metabolism 
Cluster76<- read.delim("expressionclustertissue_76_Non-specific.tsv", header= TRUE, sep="\t")




### DEPMAP ##
# Model.csv
# Sample info about DepMap cell lines
# Metadata for all of DepMap’s cancer models/cell lines. A full description of each column is available in the DepMap Release README file.
#Current DepMap Release data, including CRISPR Screens, PRISM Drug Screens, Copy Number, Mutation, Expression, and Fusions
#DepMap, Broad (2024). DepMap 24Q2 Public. Figshare+. Dataset. https://doi.org/10.25452/figshare.plus.25880521.v1

# Downloaded November 8, 2024
# 24Q2
Cell_Info<- read.csv("Model.csv", header=TRUE)


# aneuploidy_scores.csv
# Aneuploidy data per cell line.  Including total aneuploidy score, and cell line mean ploidy and +/- Genome doubling
# Downloaded November 9, 2023, from DepMap
# May 2023 data
Aneuploidy_Scores<- read.csv("aneuploidy_scores.csv", header=TRUE)


# arm_call_scores.csv
# depmap aneuploidy chromosome copy number dataset 
arm_call_scores<- read.csv("arm_call_scores.csv", header= TRUE)


# CRISPRGeneEffect.csv
# CRISPR screen data from DepMap
# Downloaded November 9, 2023
# May 2023 data
CRISPR_Gene_all<- read.csv("CRISPRGeneEffect.csv", header=TRUE)
colnames(CRISPR_Gene_all)[which(names(CRISPR_Gene_all) == "X")] <- "ModelID"


# OmicsExpressionProteinCodingGenesTPMLogp1.csv
# DepMap RNA expression files
# DepMap, Broad (2024). DepMap 24Q4 Public. Figshare+. Dataset.
# downloaded April 9, 2025
# DepMap Public 24Q4
DepMap_Expression<- read.csv("OmicsExpressionProteinCodingGenesTPMLogp1.csv", header=TRUE)




# DepMap Proteomics data
# downloaded August 2020
## Normalized Protein Expression data 
##  https://www.cell.com/cell/fulltext/S0092-8674(19)31385-6#secsectitle0190
## Nusinow et al. Cell, 2020, Quantitative proteomics of the cancer cell line encyclopedia

Protein_Expression.1<-read_excel("mmc2.xlsx", 
                                 sheet= "Normalized Protein Expression")
Protein_Expression.1<-Protein_Expression.1[,c(1:426)] #delete empty collumns

m <- Protein_Expression.1$Protein_Id   
ID_Protein_Expression2<- as.data.frame(t(as.matrix(Protein_Expression.1[,-1] ))) #switch rows and collumns
colnames(ID_Protein_Expression2) <- m

ID_Protein_Expression2$Cell_Lines<-factor(rownames(ID_Protein_Expression2)) 

for (i in 1:(length(ID_Protein_Expression2)-1)){
  ID_Protein_Expression2[,i]<- as.numeric(ID_Protein_Expression2[,i])
}

ID_Protein_Expression2$Cell_Lines<-gsub("_TenPx.*","", ID_Protein_Expression2$Cell_Lines)



CCLE.Depmap.names<-read.delim2("DepMap-2018q3-celllines.csv", 
                               dec=".", header = TRUE, sep=",")
CCLE.Depmap.name<-CCLE.Depmap.names[,1:2] #get CCLE names and depmap ID correlated table

ID_Protein_effect_4<-merge(x= CCLE.Depmap.name, y= ID_Protein_Expression2, 
                           by.x="CCLE_Name", by.y="Cell_Lines", 
                           sort = TRUE)



# Name dataframe with protein_ID and corresponding Protein name, this way we have both gene ID and Protein name
Protein_ProID<-Protein_Expression.1[,c(1:2, 6)] 
#This data is from protein paper, not HGNC
Protein_ProID$Protein_Id<- str_replace_all(Protein_ProID$Protein_Id, "[|]", ".")
Protein_ProID$Protein_Id<- str_replace_all(Protein_ProID$Protein_Id, "[-]", ".")




###
## Compare chromosome loss/gain difference with protein expression changes. 
## Data from Dephoure et al. 2014, eLife; 3:e03023
## Quantitative proteomic analysis reveals posttranslational responses to aneuploidy in yeast
## This is the data Jason M. Sheltzer prepared. 

# Yeast TMT means: 
Yeast.Prot<- read_xlsx("yeast-human aneuploidy.xlsx") 


##  Get Stingele et al. 2012  HCT116 aneuploidy data ##
Stingele.HCT116<- read_xls("Stingele.etal.2012.Genomics.xls", sheet=2)
Stingele.HCT116$Gene<- sub(";.*", "", Stingele.HCT116$`Gene Names`) 

Stingele<- Stingele.HCT116[,c(3,4,9,20, 26, 30, 34,45)]


### 
# Durrbaum et al. 2014 Unique features of the transcriptional response to model aneuploidy in human cells 
# BMC Genomics. 2014 Feb 18;15:139. doi: 10.1186/1471-2164-15-139
# HCT116 aneuploid clones RNA data
H.R.RNA<- read.xlsx("1471-2164-15-139-S2.xlsx", sheetIndex = 1)




#### Pathways #####

# Get pathways: 
Proteosome_Genes<- read.csv("group-690.csv", header=TRUE)

Olfactory_Genes<- read.csv("group-141.csv", skip = 1 , header=TRUE)
Ribosomal_Genes<- read.csv("group-1054.csv", skip = 1 , header=TRUE) #HGNC Ribosomal Proteins
Ribosomal_Genes_noMT<- subset(Ribosomal_Genes, Group %in% c("L ribosomal proteins", "S ribosomal proteins"))

RNA_Processing<- read_tsv("GO_term_summary_20250720_134600.txt")
RNA_Processing1<-unique(toupper(RNA_Processing$Symbol))
RNA_Processing_Spliceosome<- unique(toupper(subset(RNA_Processing, `Annotated Term` %in% 
                                                     c("mRNA splicing, via spliceosome", 
                                                       "regulation of mRNA splicing, via spliceosome") )$Symbol ) )
RNA_Processing_rRNA<- unique(toupper(subset(RNA_Processing, `Annotated Term` %in% 
                                              c("rRNA processing", "regulation of rRNA processing"))$Symbol))


MitoTranslationTranscription <- c("TRMT10C","HSD17B10","PRORP", 
                                  "POLRMT", "TFAM", "TFB2M",
                                  "MRPL58", "MTRFR", "MTRF1", "MTRF1L", 
                                  "MRPL1", "MRPL2", "MRPL3", "MRPL4", "MRPL9", "MRPL10", "MRPL11", "MRPL12", "MRPL13", "MRPL14", "MRPL15", "MRPL16", "MRPL17", 
                                  "MRPL18", "MRPL19", "MRPL20", "MRPL21", "MRPL22", "MRPL23", "MRPL24", "MRPL27", "MRPL28", "MRPL30", "MRPL32", "MRPL33", "MRPL34", "MRPL35", "MRPL36", 
                                  "MRPL37", "MRPL38", "MRPL39", "MRPL40", "MRPL41", "MRPL42", "MRPL43", "MRPL44", "MRPL45", "MRPL46", "MRPL47", "MRPL48", "MRPL49", "MRPL50", 
                                  "MRPL51", "MRPL52", "MRPL53", "MRPL54", "MRPL55", "MRPL58", "MRPL57", "GADD45GIP1", "MRPS30", "MRPS18A", "MRPS2", "MRPS5", "MRPS6", "MRPS7", "MRPS9", "MRPS10", 
                                  "MRPS11", "MRPS12", "MRPS14", "MRPS15", "MRPS16", "MRPS17", "MRPS18C", "MRPS21", "MRPS22", "MRPS23", "MRPS24", "MRPS25", "MRPS26", "MRPS27", 
                                  "MRPS28", "DAP3", "MRPS31", "MRPS33", "MRPS34", "MRPS35", "CHCHD1", "AURKAIP1", "PTCD3", "MRPS18B", 
                                  "FASTK", "FASTKD1", "FASTKD2", "FASTKD3", "TBRG4", "FASTKD5")
MitoElectronChainnAssembly <- c("NDUFS1","NDUFS2","NDUFS3","NDUFS7","NDUFS8","NDUFV1","NDUFV2", "NDUFAB1","NDUFA1",
                                "NDUFA2","NDUFA3","NDUFA5","NDUFA6","NDUFA7","NDUFA8","NDUFA9","NDUFA10","NDUFA11","NDUFA12",
                                "NDUFA13","NDUFB1","NDUFB2","NDUFB3","NDUFB4","NDUFB5","NDUFB6","NDUFB7","NDUFB8","NDUFB9","NDUFB10",
                                "NDUFB11","NDUFC1","NDUFC2","NDUFS4","NDUFS5","NDUFS6","NDUFV3", 
                                "SDHA", "SDHB", "SDHC", "SDHD", 
                                "UQCRB","UQCRQ","UQCRC1","UQCRC2","MT-CYB","CYC1","UQCRFS1","UQCRH","UQCR10","UQCR11", 
                                "COX4I1","COX4I2","COX5A","COX5B","COX6A1", "COX6A2","COX6B1","COX6B2","COX6C", "COX7A1","COX7A2","COX7B","COX7B2","COX7C",
                                "COX8A","COX8C","MT-CO1","MT-CO2","MT-CO3", 
                                "ATP5F1A", "ATP5F1B", "ATP5F1C",  "ATP5F1D", "ATP5F1E", "ATP5MC1", "ATP5MC2", "ATP5MC3", "ATP5ME", "ATP5MF", "ATP5MG", "ATP5MJ", "ATP5MK", 
                                "MT-ATP6", "MT-ATP8", "ATP5PB", "ATP5PD", "ATP5PF", "ATP5PO", "ATP5IF1", 
                                "ATPAF1","ATPAF2","BCS1L","CHCHD7","CMC1","CMC2", "COA1",  "COA3", "COA4","COA5", "COA6",  "COA7", "COA8",
                                "COX7A2L","COX10", "COX11", "COX14","COX15", "COX16",  "COX17", "COX18", "COX19", "COX20",
                                "DMAC1", "DMAC2", "FDXR", "FDX2", "FMC1", "FOXRED1","HCCS","HIGD1A","HIGD2A","LRPPRC", "LYRM7",
                                "NDUFAF1", "NDUFAF2", "NDUFAF3", "NDUFAF4", "NDUFAF5", "NDUFAF6","NDUFAF7",  "NDUFAF8","NUBPL",
                                "OCIAD2", "OXA1L", "PET100","PET117","PNKD", "SCO1", "SCO2",  "SDHAF1","SDHAF2","SDHAF3", "SDHAF4", "SFXN4",
                                "SMIM20", "SURF1","TACO1", "TIMMDC1","TMEM70","TMEM126A", "TMEM177","TMEM186","TMEM223","TMEM242","TTC19", 
                                "UQCC1","UQCC2","UQCC3", "UQCC4","UQCC5","UQCC6", 
                                "ACAD9","COA1","ECSIT","NDUFAF1","TMEM126B","TMEM186", 
                                "BOLA3", "FDXR", "FXN", "GLRX5", "IBA57", "ISCA1", "ISCA2", "ISCU", "LYRM4", "NFS1", "NFU1", "NUBPL") # HGNC groups

### Durrbaum et al. 2014 HCT116 RNA #####

H.R.RNA
# "hCG_18056;tcag7.352;Ube2h;UBE2H"
H.R.RNA$log2HCT116.5_4_2N<- as.numeric(as.character(H.R.RNA$log2HCT116.5_4_2N))
H.R.RNA$p.value_HCT116.5_4<- as.numeric(as.character(H.R.RNA$p.value_HCT116.5_4))

H.R.RNA$log2HCT116.H2B.GFP.5_4_2N<- as.numeric(as.character(H.R.RNA$log2HCT116.H2B.GFP.5_4_2N))
H.R.RNA$p.value_HCT116.H2B.5_4<- as.numeric(as.character(H.R.RNA$p.value_HCT116.H2B.5_4))

H.R.RNA$log2HCT116.3_3_2N<- as.numeric(as.character(H.R.RNA$log2HCT116.3_3_2N))
H.R.RNA$p.value_HCT116.3_3<- as.numeric(as.character(H.R.RNA$p.value_HCT116.3_3))

H.R.RNA$log2HPT1_2N<- as.numeric(as.character(H.R.RNA$log2HPT1_2N))
H.R.RNA$p.value.HPT1<- as.numeric(as.character(H.R.RNA$p.value.HPT1))

H.R.RNA$log2HPT2_2N<- as.numeric(as.character(H.R.RNA$log2HPT2_2N))
H.R.RNA$p.value.HPT2<- as.numeric(as.character(H.R.RNA$p.value.HPT2))

H.R.RNA$log2.RPE1.H2B.GFP.21_3_2N<- as.numeric(as.character(H.R.RNA$log2.RPE1.H2B.GFP.21_3_2N))
H.R.RNA$p.value_RPE1.H2B.21_3<- as.numeric(as.character(H.R.RNA$p.value_RPE1.H2B.21_3))

H.R.RNA$log2.RPE1.5_3.12_3_2N<- as.numeric(as.character(H.R.RNA$log2.RPE1.5_3.12_3_2N))
H.R.RNA$p.value_RPE1.12_3.5_3<- as.numeric(as.character(H.R.RNA$p.value_RPE1.12_3.5_3))




H.R.RNA$meanlog2FC_H<- rowMeans(H.R.RNA[, c(3,6,9)] )
H.R.RNA$RankMeanFC_H<- rank(H.R.RNA$meanlog2FC_H)

ggplot(H.R.RNA, aes(x= RankMeanFC_H, y= meanlog2FC_H))+
  geom_point()+
  xlim(0,13847)+
  xlab("Rank HCT116 Aneu (Durrbaum et al.)")+
  ylab("log2FC RNA")+
  geom_point(data= subset(H.R.RNA, Gene.name == "hCG_18056;tcag7.352;Ube2h;UBE2H"), color="red", size=3)+
  theme_classic()
# Durrbaum et al. 2014 
# plot.Durrbaum.RNA.AneuFC.rank.UBE2H_HCT116Ts5.Ts3.Tetra5.pdf
# Figure S7D
# UBE2H is in top 7.2% of aneuploid dependencies in HCT116. 



### Stingele et al. 2012 HCT116 protein ####

colnames(Stingele)

Stingele[,3]<- as.numeric(Stingele$`mRNA HCT116 5/4`)
Stingele[,4]<- as.numeric(Stingele$`Protein HCT116 5/4`)
Stingele[,5]<- as.numeric(Stingele$`Protein HCT116 H2B-GFP 5/4`)
Stingele[,6]<- as.numeric(Stingele$`Protein HCT116 H2B-GFP 5/3`)
Stingele[,7]<- as.numeric(Stingele$`Protein HCT116 3/3`)

Stingele$MeanHCT116<- rowMeans(Stingele[,c(4,5,6,7)])

Stingele$RankH<- rank(Stingele$MeanHCT116)

ggplot(Stingele, aes(x= RankH, y= MeanHCT116))+
  geom_point()+
  xlim(0,3564)+
  xlab("Rank HCT116 Aneu (Stingele et al.)")+
  ylab("log2FC Protein")+
  geom_point(data= subset(Stingele, Gene == "UBE2H"), color="red", size=3)+
  theme_classic()
# Plot.Stingele.HCT116.Protein.Rank.pdf 
# UBE2H is 3552 out of 3564 top 99.6% of upregulated proteins in aneuploid HCT116

# Figure S7E



### Yeast: Protein & RNA difference in aneuploid yeast  #####

###Plot difference upon trisomy in yeast 
subset(Yeast.Prot, `Yeast Gene`== "UBC8")

ggplot(Yeast.Prot, aes(x=`Yeast RNA`, y=`Yeast TMT - protein`))+
  geom_point(color="black")+
  geom_point(data= subset(Yeast.Prot, `Yeast Gene`== "UBC8"), color="red", size=2)+
  theme_classic()+
  geom_hline(yintercept = 0)+ #add lines at y=0
  geom_vline(xintercept = 0)+ #add lines at x=0
  ylab("Protein difference (Dephoure et al.)\nin aneuploid yeast")+
  xlab("RNA difference (Dephoure et al.)\nin aneuploid yeast")
# size 4x4 
# plot.Dephoure.YeastTrisomyExpression_UBE2H
# Human gene names of proteins upregulated in aneuploid cells: 
#  subset(Yeast.Prot, `Yeast TMT - protein`>1.3 & `Yeast RNA`>1.3)$`Human Gene`
# "PCNA"     "ZMPSTE24" "ATG3"     "ELAC2"    "UBE2H"    "NOM1"     "CYC1"     "PRDX6"    "BLMH"     "TUFM"     "SLC25A35" "SIRT2"    "C1QBP"   
# "IDH3A"    "ALDH5A1"  "SUCLG1"   "VPS18"    "SHPRH"    "ISCA2"    "BLVRB"    "STX16"    "ARL1"     "MDH2" 


Yeast.Prot$RankYeastRNA<- rank(Yeast.Prot$`Yeast RNA`)
ggplot(Yeast.Prot, aes(x= RankYeastRNA, y= `Yeast RNA`))+
  geom_point()+
  xlab("Rank Yeast Aneuploid (Dephoure et al.)")+
  ylab("log2FC RNA")+
  geom_point(data= subset(Yeast.Prot, `Yeast Gene`== "UBC8"), color="red")+
  theme_classic()
# size 4x4 
# Figure S7F
# plot.Dephoure.YeastTrisomyRNA_UBE2H_Rank
# 665/738= top 9.9% of RNA


Yeast.Prot$RankYeast<- rank(Yeast.Prot$`Yeast TMT - protein`)
ggplot(Yeast.Prot, aes(x= RankYeast, y= `Yeast TMT - protein`))+
  geom_point()+
  xlab("Rank Yeast Aneuploid (Dephoure et al.)")+
  ylab("log2FC Protein")+
  geom_point(data= subset(Yeast.Prot, `Yeast Gene`== "UBC8"), color="red")+
  theme_classic()
# size 4x4 
# plot.Dephoure.YeastTrisomyProtein_UBE2H_Rank
# 700/738= top 5.2% of proteins
# Figure S7G





### TCGA UBE2H RNA expression #####
# Merge data to get more info and about patients
TCGA_Patient_info<- merge(TCGA_RNA.UBE2H, TCGA_Sample[, c(1,2,4, 13)], 
                                by.x= "SAMPLE_ID", by.y= "SAMPLE_ID")
TCGA_Patient_info<- merge(x= TCGA_Patient_info, y=TCGA_Patient, 
                                by.y= "PATIENT_ID", by.x= "PATIENT_ID")


ggplot(TCGA_Patient_info, aes(y=log2(UBE2H), x=CANCER_TYPE_ACRONYM))+
  geom_boxplot(outlier.shape=NA)+  
  #ylim(0,5000)+
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  theme_classic()
# TCGA.UBE2H.RNA.perTissue.pdf

#Note Leukemias and Mature B-cell Neoplasms have much lower UBE2H expression than other cell lines 
# DLBC and LAML


# Aneuploidy correlation:
TCGA_Patient_info_Aneu<- merge(x= TCGA_Patient_info, y=TCGA_Aneuploidy, 
                          by.x= "SAMPLE_ID", by.y= "Sample")
TCGA_Patient_info_Aneu$PloidEst<- TCGA_Patient_info_Aneu$Genome_doublings
TCGA_Patient_info_Aneu <- TCGA_Patient_info_Aneu %>%
  mutate(PloidEst = recode(PloidEst, # make estimated ploidy: 0 doublings = 2, 1 doubling = tetraploid 4, 2 doubling is octoploid
                    `0` = 2,
                    `1` = 4,
                    `2` = 8))
TCGA_Patient_info_Aneu$AneuScore.ploidy<- TCGA_Patient_info_Aneu$AneuploidyScore.AS. / TCGA_Patient_info_Aneu$PloidEst

ggplot(TCGA_Patient_info_Aneu, aes(y=log2(UBE2H), x=AneuScore.ploidy))+
  geom_point()+
  geom_smooth(method="lm")+
  xlab("Aneuploidy Score")+ 
  ylab("UBE2H expression")+
  theme_classic()
# TCGA.UBE2H.RNA.Aneuploidy.ploidy.pdf
# Figure S7J

x<- cor.test(log2(TCGA_Patient_info_Aneu$UBE2H), TCGA_Patient_info_Aneu$AneuScore.ploidy)
x$p.value
x
# cor= 0.1028087
# p= 4.149232e-24


# UBE2H expression by 7q status 
ggplot(data = subset(TCGA_Patient_info_Aneu, ! X7q == "NA"), # remove samples where 7q CN status is NA
       aes(y=log2(UBE2H), x=as.factor(X7q)))+
  geom_boxplot()+
  xlab("7q Copy Number")+ 
  ylab("UBE2H expression")+
  theme_classic()
# TCGA.UBE2H.RNA.7qCNboxplot.pdf
# Figure S7L
mean(subset(TCGA_Patient_info_Aneu, X7q == "0" & ! UBE2H %in% c(NA))$UBE2H) - mean(subset(TCGA_Patient_info_Aneu, X7q == "1" & ! UBE2H %in% c(NA))$UBE2H)
x<- t.test(subset(TCGA_Patient_info_Aneu, X7q == "0")$UBE2H, subset(TCGA_Patient_info_Aneu, X7q == "1")$UBE2H)
x$p.value # neutral gain p= 1.894639e-39
mean(subset(TCGA_Patient_info_Aneu, X7q == "0" & ! UBE2H %in% c(NA))$UBE2H) - mean(subset(TCGA_Patient_info_Aneu, X7q == "-1" & ! UBE2H %in% c(NA))$UBE2H)
x<- t.test(subset(TCGA_Patient_info_Aneu, X7q == "0")$UBE2H, subset(TCGA_Patient_info_Aneu, X7q == "-1")$UBE2H)
x$p.value # neutral loss 5.852346e-09


# correlation of UBE2H expression and aneuploidy score in 7q neutral cells only
ggplot(subset(TCGA_Patient_info_Aneu, X7q==0), aes(y=log2(UBE2H), x=AneuScore.ploidy))+
  geom_point()+
  geom_smooth(method="lm")+
  xlab("Aneuploidy Score (neutral 7q only)")+ 
  ylab("UBE2H expression")+
  theme_classic()
# TCGA.UBE2H.RNA.Aneuploidy.ploidy_Neutral7q.pdf
# Figure S7O
x<- cor.test(log2(subset(TCGA_Patient_info_Aneu, X7q==0)$UBE2H), subset(TCGA_Patient_info_Aneu, X7q==0)$AneuScore.ploidy)
x$p.value
x
# p=5.311165e-20 
# corr= 0.1188832




# TCGA  UBE2H RNA correspond with aneuscore per tissue
UBE2H_RNA.Aneu_Tissue<- data.frame (TissueType = character(), 
                                     UBE2H.AnCorr = numeric(), 
                                    UBE2H.An.p = character())

for (i in 1: length(unique(TCGA_Patient_info_Aneu$CANCER_TYPE_ACRONYM))) {
  Tissue<- unique(TCGA_Patient_info_Aneu$CANCER_TYPE_ACRONYM)[i]
  Corr<- cor.test(log2(subset(TCGA_Patient_info_Aneu, CANCER_TYPE_ACRONYM == Tissue)$UBE2H), 
                  subset(TCGA_Patient_info_Aneu, CANCER_TYPE_ACRONYM == Tissue)$AneuScore.ploidy, 
                  method = "pearson")
    
  UBE2H_RNA.Aneu_Tissue<- rbind(UBE2H_RNA.Aneu_Tissue, data.frame(TissueType = Tissue, 
                                                                  UBE2H.AnCorr = Corr$estimate, 
                                                                  UBE2H.An.p = Corr$p.value ))
  
}
UBE2H_RNA.Aneu_Tissue
# Most significant tissue types: 

ggplot(UBE2H_RNA.Aneu_Tissue, aes(x=UBE2H.AnCorr, y= -log2(UBE2H.An.p), label=TissueType))+
  geom_point()+
  geom_hline(yintercept=-log2(0.05))+
  geom_text(hjust=0, vjust=0)+
  theme_classic()
# plot.TCGA.TissueType.corr.pvalue.pdf
# Figure 6H


ggplot(subset(TCGA_Patient_info_Aneu, CANCER_TYPE_ACRONYM == "KIRC"), aes(y=log2(UBE2H), x=AneuScore.ploidy))+
  geom_point()+
  geom_smooth(method="lm")+
  theme_classic()
# TCGA.UBE2H.RNA.Aneuploidy_KIRC.pdf
x<- cor.test(log2(subset(TCGA_Patient_info_Aneu, CANCER_TYPE_ACRONYM == "KIRC")$UBE2H), 
             subset(TCGA_Patient_info_Aneu, CANCER_TYPE_ACRONYM == "KIRC")$AneuScore.ploidy)
x
# 0.2447174
# p= 5.45e-08
# Figure S7L

ggplot(subset(TCGA_Patient_info_Aneu, CANCER_TYPE_ACRONYM == "ACC"), aes(y=log2(UBE2H), x=AneuScore.ploidy))+
  geom_point()+
  geom_smooth(method="lm")+
  theme_classic()
# TCGA.UBE2H.RNA.Aneuploidy_ACC.pdf
x<- cor.test(log2(subset(TCGA_Patient_info_Aneu, CANCER_TYPE_ACRONYM == "ACC")$UBE2H), 
             subset(TCGA_Patient_info_Aneu, CANCER_TYPE_ACRONYM == "ACC")$AneuScore.ploidy, method = "pearson")
x
# 0.4418611
# p= 6.444e-05
# Figure 6M


ggplot(subset(TCGA_Patient_info_Aneu, CANCER_TYPE_ACRONYM == "PAAD"), aes(y=log2(UBE2H), x=AneuScore.ploidy))+
  geom_point()+
  geom_smooth(method="lm")+
  theme_classic()
# TCGA.UBE2H.RNA.Aneuploidy_Pancreatic.pdf
x<- cor.test(log2(subset(TCGA_Patient_info_Aneu, CANCER_TYPE_ACRONYM == "PAAD")$UBE2H), 
             subset(TCGA_Patient_info_Aneu, CANCER_TYPE_ACRONYM == "PAAD")$AneuScore.ploidy)
x
# 0.2145797
# p= 0.006781


ggplot(subset(TCGA_Patient_info_Aneu, CANCER_TYPE_ACRONYM == "SKCM"), aes(y=log2(UBE2H), x=AneuScore.ploidy))+
  geom_point()+
  geom_smooth(method="lm")+
  theme_classic()
# TCGA.UBE2H.RNA.Aneuploidy_SKCM.pdf
x<- cor.test(log2(subset(TCGA_Patient_info_Aneu, CANCER_TYPE_ACRONYM == "SKCM")$UBE2H), 
             subset(TCGA_Patient_info_Aneu, CANCER_TYPE_ACRONYM == "SKCM")$AneuScore.ploidy)
x
# 0.187036
# p= 9.028e-05



ggplot(subset(TCGA_Patient_info_Aneu, CANCER_TYPE_ACRONYM %in% c("READ")), aes(y=log2(UBE2H), x=AneuScore.ploidy))+
  geom_point()+
  geom_smooth(method="lm")+
  theme_classic()
# TCGA.UBE2H.RNA.Aneuploidy_Pancreatic.pdf
x<- cor.test(log2(subset(TCGA_Patient_info_Aneu, CANCER_TYPE_ACRONYM == "READ")$UBE2H), 
             subset(TCGA_Patient_info_Aneu, CANCER_TYPE_ACRONYM == "READ")$AneuScore.ploidy)
x
# 0.2272697
# p= 0.004589






  ##### TCGA Genes in cluster 76 Metabolism  #####

length(Cluster76$Gene[Cluster76$Gene %in% MitoElectronChainnAssembly]) #14 genes
length(Cluster76$Gene[Cluster76$Gene %in% MitoTranslationTranscription]) #17


  ##### TCGA Percent 7q gain #####
# Merge data to get more info and about patients
TCGA_Patient_info<- merge(TCGA_RNA.UBE2H, TCGA_Sample[, c(1,2,4, 13)], 
                          by.x= "SAMPLE_ID", by.y= "SAMPLE_ID")
TCGA_Patient_info<- merge(x= TCGA_Patient_info, y=TCGA_Patient, 
                          by.y= "PATIENT_ID", by.x= "PATIENT_ID")
TCGA_Patient_info_Aneu<- merge(x= TCGA_Patient_info, y=TCGA_Aneuploidy, 
                               by.x= "SAMPLE_ID", by.y= "Sample")


table(TCGA_Aneuploidy$X7q)
#-1    0    1 
#507 6171 2490

# Percent of TCGA tumors with 7q gain: 
PercentGainLoss<- data.frame(Arm= character(),
                             percentGain= as.numeric(), 
                             percentLoss = as.numeric())
nTumor<- length(TCGA_Aneuploidy$X1p)

for (i in 14:52){
  TableCN <- table(TCGA_Aneuploidy[i])
  PercentGainLoss<- rbind(PercentGainLoss, data.frame(Arm= colnames(TCGA_Aneuploidy)[i],
                                                      percentGain = TableCN["1"]/nTumor, 
                                                      percentLoss= TableCN["-1"]/nTumor))
}
PercentGainLoss<- PercentGainLoss[order(PercentGainLoss$percentGain),]
PercentGainLoss
# 7q percent gain: 0.23664703  
# 7q percent loss: 0.04818476

# 7q gain 5th most common in TCGA

  ##### Difference in survival High & low UBE2H Per CANCER TYPE: #####

# No Blood cancer: 
TCGA_Patient_info_NoBlood<- subset(TCGA_Patient_info, ! CANCER_TYPE %in% c("Leukemia", "Mature B-Cell Neoplasms"))
q10 <- quantile(TCGA_Patient_info_NoBlood$UBE2H, 0.10, na.rm = TRUE)
q90 <- quantile(TCGA_Patient_info_NoBlood$UBE2H, 0.90, na.rm = TRUE)
TCGA_Patient_info_NoBlood <- TCGA_Patient_info_NoBlood %>%
  mutate(UBE2Hgroup2 = case_when(
    UBE2H <= q10 ~ "bottom_10%",
    UBE2H >= q90 ~ "top_10%",
    TRUE ~ "middle"
  ))

survfit2(Surv(OS_MONTHS, OS_STATUS== "1:DECEASED") ~ UBE2Hgroup2, 
         data =  TCGA_Patient_info_NoBlood ) %>% 
  ggsurvfit() +
  scale_color_manual(values=c("#00BFC4","black", "#F8766D"))+
  scale_fill_manual(values=c("#00BFC4","black", "#F8766D"))+
  add_confidence_interval() +
  labs(
    title = "TCGA survival:\n High vs low UBE2H Expression", 
    x = "Months",
    y = "Survival probability"
  ) 
# Plot.TCGA.Survival.UBE2HMeanRNA_NoBlood.pdf
survdiff(Surv(OS_MONTHS, OS_STATUS== "1:DECEASED") ~ UBE2Hgroup2 == "top_10%", 
         data = TCGA_Patient_info_NoBlood)$pvalue # Top 10% decrease survival? 
# NS
survdiff(Surv(OS_MONTHS, OS_STATUS== "1:DECEASED") ~ UBE2Hgroup2 == "bottom_10%", 
         data = TCGA_Patient_info_NoBlood)$pvalue # bottom 10% increase survival? 
# p=0.02703888

survdiff(Surv(OS_MONTHS, OS_STATUS== "1:DECEASED") ~ UBE2Hgroup2 == "bottom_10%", 
         data = subset(TCGA_Patient_info_NoBlood, UBE2Hgroup2 %in% c("top_10%", "bottom_10%")) )$pvalue # bottom 10% increase survival?  
# 0.2825731 
# 20% 0.4626602 





# TCGA  Survival High vs low UBE2H: (Pancreatic Cancer) 
TCGA_Patient_info_Pancreatic<- subset(TCGA_Patient_info, CANCER_TYPE_ACRONYM %in% c("PAAD"))
q10 <- quantile(TCGA_Patient_info_Pancreatic$UBE2H, 0.10, na.rm = TRUE)
q90 <- quantile(TCGA_Patient_info_Pancreatic$UBE2H, 0.90, na.rm = TRUE)
TCGA_Patient_info_Pancreatic <- TCGA_Patient_info_Pancreatic %>%
  mutate(UBE2Hgroup2 = case_when(
    UBE2H <= q10 ~ "bottom_10%",
    UBE2H >= q90 ~ "top_10%",
    TRUE ~ "middle"
  ))

survfit2(Surv(OS_MONTHS, OS_STATUS== "1:DECEASED") ~ UBE2Hgroup2, 
         data =  TCGA_Patient_info_Pancreatic ) %>% 
  ggsurvfit() +
  scale_color_manual(values=c("#00BFC4","black", "#F8766D"))+
  scale_fill_manual(values=c("#00BFC4","black", "#F8766D"))+
  add_confidence_interval() +
  labs(
    title = "TCGA PAAD survival:\n High vs low UBE2H Expression", 
    x = "Months",
    y = "Survival probability"
  ) 
# Plot.TCGA.Survival.UBE2HMeanRNA_PAAD.pdf
# Figure 6I
survdiff(Surv(OS_MONTHS, OS_STATUS== "1:DECEASED") ~ UBE2Hgroup2 == "top_10%", 
         data =TCGA_Patient_info_Pancreatic)$pvalue # Top 10% decrease survival? 
# NS
survdiff(Surv(OS_MONTHS, OS_STATUS== "1:DECEASED") ~ UBE2Hgroup2 == "bottom_10%", 
         data = TCGA_Patient_info_Pancreatic)$pvalue # bottom 10% increase survival?  
# 0.0219755 

survdiff(Surv(OS_MONTHS, OS_STATUS== "1:DECEASED") ~ UBE2Hgroup2 == "bottom_10%", 
         data = subset(TCGA_Patient_info_Pancreatic, UBE2Hgroup2 %in% c("top_10%", "bottom_10%")) )$pvalue # bottom 10% increase survival?  
# 0.01182899 
# 20% 0.001983687



# TCGA  Survival High vs low UBE2H: (LUAD) 
TCGA_Patient_info_LUAD<- subset(TCGA_Patient_info, CANCER_TYPE_ACRONYM %in% c("LUAD"))
q10 <- quantile(TCGA_Patient_info_LUAD$UBE2H, 0.10, na.rm = TRUE)
q90 <- quantile(TCGA_Patient_info_LUAD$UBE2H, 0.90, na.rm = TRUE)
TCGA_Patient_info_LUAD <- TCGA_Patient_info_LUAD %>%
  mutate(UBE2Hgroup2 = case_when(
    UBE2H <= q10 ~ "bottom_10%",
    UBE2H >= q90 ~ "top_10%",
    TRUE ~ "middle"))

survfit2(Surv(OS_MONTHS, OS_STATUS== "1:DECEASED") ~ UBE2Hgroup2, 
         data =  TCGA_Patient_info_LUAD ) %>% 
  ggsurvfit() +
  scale_color_manual(values=c("#00BFC4","black", "#F8766D"))+
  scale_fill_manual(values=c("#00BFC4","black", "#F8766D"))+
  add_confidence_interval() +
  labs(
    title = "TCGA LUAD survival:\n High vs low UBE2H Expression", 
    x = "Months",
    y = "Survival probability"
  ) 
# Plot.TCGA.Survival.UBE2HMeanRNA_LUAD.pdf
# Supplementary Figure S7T
survdiff(Surv(OS_MONTHS, OS_STATUS== "1:DECEASED") ~ UBE2Hgroup2 =="top_10%" , 
         data = TCGA_Patient_info_LUAD)$pvalue # Top 10% decrease survival? 
# 0.0345219
survdiff(Surv(OS_MONTHS, OS_STATUS== "1:DECEASED") ~ UBE2Hgroup2 =="bottom_10%", 
         data = TCGA_Patient_info_LUAD)$pvalue # bottom 10% increase survival?  
# 0.2210343

survdiff(Surv(OS_MONTHS, OS_STATUS== "1:DECEASED") ~ UBE2Hgroup2 == "bottom_10%", 
         data = subset(TCGA_Patient_info_LUAD, UBE2Hgroup2 %in% c("top_10%", "bottom_10%")) )$pvalue # bottom 10% increase survival?  
# 0.02051004 
# 20% 0.002575058



# TCGA  Survival High vs low UBE2H: (KIRP) 
TCGA_Patient_info_KIRP<- subset(TCGA_Patient_info, CANCER_TYPE_ACRONYM %in% c("KIRP"))
q10 <- quantile(TCGA_Patient_info_KIRP$UBE2H, 0.10, na.rm = TRUE)
q90 <- quantile(TCGA_Patient_info_KIRP$UBE2H, 0.90, na.rm = TRUE)
TCGA_Patient_info_KIRP <- TCGA_Patient_info_KIRP %>%
  mutate(UBE2Hgroup2 = case_when(
    UBE2H <= q10 ~ "bottom_10%",
    UBE2H >= q90 ~ "top_10%",
    TRUE ~ "middle"))

survfit2(Surv(OS_MONTHS, OS_STATUS== "1:DECEASED") ~ UBE2Hgroup2, 
         data =  TCGA_Patient_info_KIRP ) %>% 
  ggsurvfit() +
  scale_color_manual(values=c("#00BFC4","black", "#F8766D"))+
  scale_fill_manual(values=c("#00BFC4","black", "#F8766D"))+
  add_confidence_interval() +
  labs(
    title = "TCGA KIRP survival:\n High vs low UBE2H Expression", 
    x = "Months",
    y = "Survival probability"
  ) 
# Plot.TCGA.Survival.UBE2HMeanRNA_KIRP.pdf
# Figure 6I
survdiff(Surv(OS_MONTHS, OS_STATUS== "1:DECEASED") ~ UBE2Hgroup2 =="top_10%" , 
         data = TCGA_Patient_info_KIRP)$pvalue # Top 10% decrease survival? 
# 0.03564965
survdiff(Surv(OS_MONTHS, OS_STATUS== "1:DECEASED") ~ UBE2Hgroup2 =="bottom_10%", 
         data = TCGA_Patient_info_KIRP)$pvalue # bottom 10% increase survival?  
# 0.02922772

survdiff(Surv(OS_MONTHS, OS_STATUS== "1:DECEASED") ~ UBE2Hgroup2 == "bottom_10%", 
         data = subset(TCGA_Patient_info_KIRP, UBE2Hgroup2 %in% c("top_10%", "bottom_10%")) )$pvalue # bottom 10% increase survival?  
# 0.00576646 
#20% NS




# TCGA  Survival High vs low UBE2H: (LGG) 
TCGA_Patient_info_LGG<- subset(TCGA_Patient_info, CANCER_TYPE_ACRONYM %in% c("LGG"))
q10 <- quantile(TCGA_Patient_info_LGG$UBE2H, 0.10, na.rm = TRUE)
q90 <- quantile(TCGA_Patient_info_LGG$UBE2H, 0.90, na.rm = TRUE)
TCGA_Patient_info_LGG <- TCGA_Patient_info_LGG %>%
  mutate(UBE2Hgroup2 = case_when(
    UBE2H <= q10 ~ "bottom_10%",
    UBE2H >= q90 ~ "top_10%",
    TRUE ~ "middle"))

survfit2(Surv(OS_MONTHS, OS_STATUS== "1:DECEASED") ~ UBE2Hgroup2, 
         data =  TCGA_Patient_info_LGG ) %>% 
  ggsurvfit() +
  scale_color_manual(values=c("#00BFC4","black", "#F8766D"))+
  scale_fill_manual(values=c("#00BFC4","black", "#F8766D"))+
  add_confidence_interval() +
  labs(
    title = "TCGA LGG survival:\n High vs low UBE2H Expression", 
    x = "Months",
    y = "Survival probability"
  ) 
# Plot.TCGA.Survival.UBE2HMeanRNA_LGG.pdf
# Supplementary Figure S7S
survdiff(Surv(OS_MONTHS, OS_STATUS== "1:DECEASED") ~ UBE2Hgroup2 =="top_10%" , 
         data = TCGA_Patient_info_LGG)$pvalue # Top 10% decrease survival? 
# 0.04049098
survdiff(Surv(OS_MONTHS, OS_STATUS== "1:DECEASED") ~ UBE2Hgroup2 =="bottom_10%", 
         data = TCGA_Patient_info_LGG)$pvalue # bottom 10% increase survival?  
# 0.008912855

survdiff(Surv(OS_MONTHS, OS_STATUS== "1:DECEASED") ~ UBE2Hgroup2 == "bottom_10%", 
         data = subset(TCGA_Patient_info_LGG, UBE2Hgroup2 %in% c("top_10%", "bottom_10%")) )$pvalue # bottom 10% increase survival?  
# 0.001196229 
# 20% NS




  ##### TCGA Hazard ratio UBE2H ####
# Add TP53 mutation data, cancer type, and relative UBE2H expression

TCGA_Patient_info <- TCGA_Patient_info %>%
  separate(col = OS_STATUS, into = c("OS_CENSOR", "OS_STATUS"), sep = ":") 
TCGA_Patient_info$OS_CENSOR<- as.numeric(TCGA_Patient_info$OS_CENSOR)

TCGA_Patient_info$TP53<- TCGA_Patient_info$SAMPLE_ID %in% TCGA_TP53Mut

TCGA_Patient_info$RelativeUBE2H<- TCGA_Patient_info$UBE2H/mean(TCGA_Patient_info$UBE2H, na.rm=TRUE)

Cox_TCGA_UBE2H <- coxph(
  Surv(OS_MONTHS, OS_CENSOR) ~ RelativeUBE2H + TP53 + AGE + CANCER_TYPE,
  data = subset(TCGA_Patient_info, OS_CENSOR %in% c(0,1) ) ) 

summary(Cox_TCGA_UBE2H)
coef(summary(Cox_TCGA_UBE2H))[,5] #this is pvalue for coxph model. Pr(>|z|)


ggforest(Cox_TCGA_UBE2H, data= TCGA_Patient_info) #hazard ratio plot
# plot.hazardratio.TCGA_UBE2HExpression.pdf
# Supplementary Figure S7U


#                                                     coef  exp(coef)  se(coef)      z Pr(>|z|)    
#RelativeUBE2H                                     0.113883  1.120621  0.038347  2.970 0.002980 ** 
#  TP53TRUE                                        0.140326  1.150649  0.043163  3.251 0.001150 ** 
#  AGE                                             0.025318  1.025641  0.001675 15.117  < 2e-16 ***

## UBE2H expression is a significant hazard! 





  ##### TCGA Calculate what percent of tumors have trisomy: ####
# What percent of tumors have Ts1q, Ts7p, Ts13, Ts2p, Ts2q, Ts5p, Ts8p, Ts8q, Ts10p Ts10q

#Make a table with total number and percent of tumors with gain or loss of arm. 
# Assume NA are neutral copy numbers 

PercentArmGainLoss<- data.frame(ChrmArm= character(), 
                                GainCount = numeric(), 
                                NeutralCount = numeric(), 
                                LossCount = numeric(), 
                                GainPercent = numeric(), 
                                LossPercent = numeric())
NumberTumors<- length(TCGA_Aneuploidy$'1q')
ArmCollumns<- c("1p","1q","2p", "2q","3p","3q", "4p",
                "4q","5p","5q", "6p","6q","7p","7q", "8p",
                "8q","9p","9q", "10p","10q","11p","11q", "12p",
                "12q","13q","14q", "15q","16p","16q","17p", "17q",
                "18p","18q","19p", "19q","20p","20q","21q", "22q")
for (i in ArmCollumns){
  table<- table(TCGA_Aneuploidy[i]) #table of number of Loss, Neutral, Gain calls for that arm in TCGA
  
  PercentArmGainLoss<- rbind(PercentArmGainLoss, data.frame(ChrmArm= i, 
                                                            GainCount = table["1"],# count of number of gains
                                                            NeutralCount = table["0"], 
                                                            LossCount = table["-1"], 
                                                            GainPercent = 100*table["1"]/NumberTumors,  #number Gain divided by total number of tumors
                                                            LossPercent = 100*table["-1"]/NumberTumors ) ) #number Gain divided by total number of tumors
}
PercentArmGainLoss$rank<- rank(-PercentArmGainLoss$GainPercent)
PercentArmGainLoss

#    ChrmArm GainCount NeutralCount LossCount GainPercent LossPercent rank  
#1        1p       355         6645      1705    3.373883   16.204144 38.0
#11       1q      2674         6281       411   25.413420    3.906102  4.0   # 25.4% Gain 1q (Rank: 4th)  (more gains than losses)
#12       2p      1075         7847       641   10.216689    6.091998 14.0   # 10.2% Gain 2p (Rank: 14th) (more gains than losses)
#13       2q       658         7627       746    6.253564    7.089907 25.5   # 6.3%  Gain 2q (Rank: 25th)
#14       3p       500         6284      2397    4.751948   22.780840 32.0
#15       3q      1819         6539       720   17.287588    6.842806  8.0
#16       4p       479         7050      2253    4.552366   21.412279 35.0
#17       4q       290         6779      2196    2.756130   20.870557 39.0
#18       5p      2326         6809       618   22.106063    5.873408  6.0  # 22.1% Gain 5p (Rank: 6th). (more gains than losses)
#19       5q       715         6359      2032    6.795286   19.311918 23.0 
#110      6p      1179         7068       938   11.205094    8.914655 12.0
#111      6q       469         6505      2154    4.457328   20.471393 36.0
#112      7p      2960         6336       434   28.131534    4.124691  2.0  # 28.1% Gain 7p (2nd). (more gains than losses)
#113      7q      2490         6171       507   23.664703    4.818476  5.0
#114      8p       979         5403      3118    9.304315   29.633150 18.0  # 9.3%  Gain 8p (18th)
#115      8q      3046         5711       349   28.948869    3.316860  1.0  # 28.9% Gain 8q (1st). (more gains than losses)
#116      9p       658         5896      2659    6.253564   25.270861 25.5 
#117      9q       708         6792      2143    6.728759   20.366850 24.0
#118     10p       984         6789      1922    9.351834   18.266489 17.0  # 9.3% Gain 10p (17th) 
#119     10q       387         6472      2313    3.678008   21.982513 37.0  # 3.7% Gain 10q (37th) 
#120     11p       515         7207      1715    4.894507   16.299183 30.0
#121     11q       586         6712      1611    5.569283   15.310777 29.0
#122     12p      1688         7064       777   16.042577    7.384528  9.0
#123     12q      1032         7449       813    9.808021    7.726668 15.0
#124     13q       867         5735      2288    8.239878   21.744915 19.0  # 8.2% Gain 13 (19th) 
#125     14q       650         6691      1829    6.177533   17.382627 28.0
#126     15q       481         6931      1701    4.571374   16.166128 34.0
#127     16p      1327         7365       990   12.611671    9.408858 10.0
#128     16q       759         6252      2421    7.213458   23.008934 22.0
#129     17p       503         5770      3444    4.780460   32.731420 31.0
#130     17q      1125         6663       959   10.691884    9.114237 13.0
#131     18p      1012         6873      2059    9.617943   19.568523 16.0
#132     18q       499         6272      2693    4.742444   25.593994 33.0
#133     19p       824         6994      1230    7.831211   11.689793 21.0
#134     19q      1193         6840       978   11.338149    9.294811 11.0
#135     20p      2317         6686       718   22.020528    6.823798  7.0
#136     20q      2888         6599       227   27.447253    2.157385  3.0
#137     21q       859         7009      1795    8.163847   17.059494 20.0
#138     22q       653         6430      2329    6.206044   22.134575 27.0




## now find most common gained/lost in colorectal cancers only (for Ts13) 

# first write a function to replace NA values with 0 so loop doesn't crash if CRC has no aneuploidy
clean_value <- function(x) { 
  if (is.na(x)) {
    return(0)
  } else {
    return(x)
  }
} 

# now do loop to get number and percent fo tumors in TCGA with gain
PercentArmGainLoss.CRC<- data.frame(ChrmArm= character(), 
                                    GainCount = numeric(), 
                                    NeutralCount = numeric(), 
                                    LossCount = numeric(), 
                                    GainPercent = numeric(), 
                                    LossPercent = numeric())
NumberCRC<- length(subset(TCGA_Aneuploidy, Type %in% c("COAD", "READ"))$'1q')
ArmCollumns<- c("1p","1q","2p", "2q","3p","3q", "4p",
                "4q","5p","5q", "6p","6q","7p","7q", "8p",
                "8q","9p","9q", "10p","10q","11p","11q", "12p",
                "12q","13q","14q", "15q","16p","16q","17p", "17q",
                "18p","18q","19p", "19q","20p","20q","21q", "22q")
for (i in ArmCollumns){
  table<- table(subset(TCGA_Aneuploidy, Type %in% c("COAD", "READ"))[i]) #table of number of Loss, Neutral, Gain calls for that arm in TCGA
  
  PercentArmGainLoss.CRC<- rbind(PercentArmGainLoss.CRC, data.frame(ChrmArm= i, 
                                                                    GainCount = clean_value(table["1"]),# count of number of gains
                                                                    NeutralCount = clean_value(table["0"]), 
                                                                    LossCount = clean_value(table["-1"]), 
                                                                    GainPercent = 100*table["1"]/NumberCRC,  #number Gain divided by total number of tumors
                                                                    LossPercent = 100*table["-1"]/NumberCRC ) ) #number Gain divided by total number of tumors
}
PercentArmGainLoss.CRC$rank<- rank(-PercentArmGainLoss.CRC$GainPercent)
PercentArmGainLoss.CRC
#    ChrmArm GainCount NeutralCount LossCount GainPercent LossPercent rank
#0        1p         0          331       149          NA  25.3401361 39.0
#1        1q        89          409        53  15.1360544   9.0136054  8.0
#11       2p        56          475        28   9.5238095   4.7619048 20.0
#12       2q        64          470        23  10.8843537   3.9115646 17.5
#13       3p        31          441        75   5.2721088  12.7551020 27.0
#14       3q        53          456        49   9.0136054   8.3333333 21.5
#15       4p        10          367       181   1.7006803  30.7823129 33.5
#16       4q         4          363       168   0.6802721  28.5714286 37.5
#17       5p        80          421        59  13.6054422  10.0340136 11.0
#18       5q        24          385       115   4.0816327  19.5578231 28.0
#19       6p        64          438        39  10.8843537   6.6326531 17.5
#110      6q        48          438        59   8.1632653  10.0340136 23.0
#111      7p       294          274         3  50.0000000   0.5102041  3.0
#112      7q       236          314        10  40.1360544   1.7006803  5.0
#113      8p        76          225       248  12.9251701  42.1768707 14.0
#114      8q       271          262        12  46.0884354   2.0408163  4.0
#115      9p        86          409        70  14.6258503  11.9047619  9.0
#116      9q        63          422        74  10.7142857  12.5850340 19.0
#117     10p        32          446        90   5.4421769  15.3061224 26.0
#118     10q        10          406       105   1.7006803  17.8571429 33.5
#119     11p        45          452        70   7.6530612  11.9047619 24.0
#120     11q        33          429        88   5.6122449  14.9659864 25.0
#121     12p        93          390        61  15.8163265  10.3741497  7.0
#122     12q        78          419        61  13.2653061  10.3741497 13.0
#123     13q       318          216        15  54.0816327   2.5510204  2.0  # 54.1% (rank 2nd) (more gains than losses)
#124     14q        21          344       190   3.5714286  32.3129252 30.0
#125     15q         9          322       197   1.5306122  33.5034014 35.0
#126     16p        80          444        34  13.6054422   5.7823129 11.0
#127     16q        80          426        40  13.6054422   6.8027211 11.0
#128     17p         6          229       324   1.0204082  55.1020408 36.0
#129     17q        65          404        67  11.0544218  11.3945578 16.0
#130     18p        23          207       326   3.9115646  55.4421769 29.0
#131     18q         4          173       386   0.6802721  65.6462585 37.5
#132     19p        53          444        60   9.0136054  10.2040816 21.5
#133     19q        66          431        40  11.2244898   6.8027211 15.0
#134     20p       219          203        90  37.2448980  15.3061224  6.0
#135     20q       426          155         0  72.4489796          NA  1.0 #
#136     21q        19          360       173   3.2312925  29.4217687 31.0
#137     22q        13          362       183   2.2108844  31.1224490 32.0








#### Genes in cluster 76 Metabolism  ####
# transcription cluster 76- metabolism #
# downloaded from cBioportal HPA on January 14, 2025

Cluster76<- read.delim("expressionclustertissue_76_Non-specific.tsv", header= TRUE, sep="\t")



MitoTranslationTranscription <- c("TRMT10C","HSD17B10","PRORP", 
                                  "POLRMT", "TFAM", "TFB2M",
                                  "MRPL58", "MTRFR", "MTRF1", "MTRF1L", 
                                  "MRPL1", "MRPL2", "MRPL3", "MRPL4", "MRPL9", "MRPL10", "MRPL11", "MRPL12", "MRPL13", "MRPL14", "MRPL15", "MRPL16", "MRPL17", 
                                  "MRPL18", "MRPL19", "MRPL20", "MRPL21", "MRPL22", "MRPL23", "MRPL24", "MRPL27", "MRPL28", "MRPL30", "MRPL32", "MRPL33", "MRPL34", "MRPL35", "MRPL36", 
                                  "MRPL37", "MRPL38", "MRPL39", "MRPL40", "MRPL41", "MRPL42", "MRPL43", "MRPL44", "MRPL45", "MRPL46", "MRPL47", "MRPL48", "MRPL49", "MRPL50", 
                                  "MRPL51", "MRPL52", "MRPL53", "MRPL54", "MRPL55", "MRPL58", "MRPL57", "GADD45GIP1", "MRPS30", "MRPS18A", "MRPS2", "MRPS5", "MRPS6", "MRPS7", "MRPS9", "MRPS10", 
                                  "MRPS11", "MRPS12", "MRPS14", "MRPS15", "MRPS16", "MRPS17", "MRPS18C", "MRPS21", "MRPS22", "MRPS23", "MRPS24", "MRPS25", "MRPS26", "MRPS27", 
                                  "MRPS28", "DAP3", "MRPS31", "MRPS33", "MRPS34", "MRPS35", "CHCHD1", "AURKAIP1", "PTCD3", "MRPS18B", 
                                  "FASTK", "FASTKD1", "FASTKD2", "FASTKD3", "TBRG4", "FASTKD5")
MitoElectronChainnAssembly <- c("NDUFS1","NDUFS2","NDUFS3","NDUFS7","NDUFS8","NDUFV1","NDUFV2", "NDUFAB1","NDUFA1",
                                "NDUFA2","NDUFA3","NDUFA5","NDUFA6","NDUFA7","NDUFA8","NDUFA9","NDUFA10","NDUFA11","NDUFA12",
                                "NDUFA13","NDUFB1","NDUFB2","NDUFB3","NDUFB4","NDUFB5","NDUFB6","NDUFB7","NDUFB8","NDUFB9","NDUFB10",
                                "NDUFB11","NDUFC1","NDUFC2","NDUFS4","NDUFS5","NDUFS6","NDUFV3", 
                                "SDHA", "SDHB", "SDHC", "SDHD", 
                                "UQCRB","UQCRQ","UQCRC1","UQCRC2","MT-CYB","CYC1","UQCRFS1","UQCRH","UQCR10","UQCR11", 
                                "COX4I1","COX4I2","COX5A","COX5B","COX6A1", "COX6A2","COX6B1","COX6B2","COX6C", "COX7A1","COX7A2","COX7B","COX7B2","COX7C",
                                "COX8A","COX8C","MT-CO1","MT-CO2","MT-CO3", 
                                "ATP5F1A", "ATP5F1B", "ATP5F1C",  "ATP5F1D", "ATP5F1E", "ATP5MC1", "ATP5MC2", "ATP5MC3", "ATP5ME", "ATP5MF", "ATP5MG", "ATP5MJ", "ATP5MK", 
                                "MT-ATP6", "MT-ATP8", "ATP5PB", "ATP5PD", "ATP5PF", "ATP5PO", "ATP5IF1", 
                                "ATPAF1","ATPAF2","BCS1L","CHCHD7","CMC1","CMC2", "COA1",  "COA3", "COA4","COA5", "COA6",  "COA7", "COA8",
                                "COX7A2L","COX10", "COX11", "COX14","COX15", "COX16",  "COX17", "COX18", "COX19", "COX20",
                                "DMAC1", "DMAC2", "FDXR", "FDX2", "FMC1", "FOXRED1","HCCS","HIGD1A","HIGD2A","LRPPRC", "LYRM7",
                                "NDUFAF1", "NDUFAF2", "NDUFAF3", "NDUFAF4", "NDUFAF5", "NDUFAF6","NDUFAF7",  "NDUFAF8","NUBPL",
                                "OCIAD2", "OXA1L", "PET100","PET117","PNKD", "SCO1", "SCO2",  "SDHAF1","SDHAF2","SDHAF3", "SDHAF4", "SFXN4",
                                "SMIM20", "SURF1","TACO1", "TIMMDC1","TMEM70","TMEM126A", "TMEM177","TMEM186","TMEM223","TMEM242","TTC19", 
                                "UQCC1","UQCC2","UQCC3", "UQCC4","UQCC5","UQCC6", 
                                "ACAD9","COA1","ECSIT","NDUFAF1","TMEM126B","TMEM186", 
                                "BOLA3", "FDXR", "FXN", "GLRX5", "IBA57", "ISCA1", "ISCA2", "ISCU", "LYRM4", "NFS1", "NFU1", "NUBPL") # HGNC groups

Cluster76$Gene[Cluster76$Gene %in% MitoElectronChainnAssembly]
Cluster76$Gene[Cluster76$Gene %in% MitoTranslationTranscription]


### DepMap data: plot UBE2H RNA expression #####

Cell_Info

Aneuploidy_Scores[1:10,]

CRISPR_Gene_all

DepMap_Expression

Aneuploidy_Scores$Aneu.ploidy<- Aneuploidy_Scores$Aneuploidy.score/Aneuploidy_Scores$Ploidy
DepMap_RNA_Aneu<- merge(Aneuploidy_Scores, DepMap_Expression, by.x="DepMap_ID", by.y="X")
DepMap_RNA_Aneu<- merge(DepMap_RNA_Aneu, arm_call_scores, by.x="DepMap_ID", by.y= "X")
DepMap_RNA_Aneu<- merge(DepMap_RNA_Aneu, CRISPR_Gene_all, by.x="DepMap_ID", by.y="ModelID")
# note: gene ID with ".x" is expression, gene id with ".y" is CRISPR effect score

#UBE2H correlation with aneuploidy score
ggplot(DepMap_RNA_Aneu, aes(y=log2(UBE2H..7328..x), x=Aneu.ploidy))+
  geom_point()+
  geom_smooth(method="lm")+
  xlab("DepMap: Aneuploidy Score")+ 
  ylab("UBE2H expression")+
  theme_classic()
# DepMap.UBE2H.RNA.Aneuploidy.ploidy.pdf
# Figure S7K
cor.test(DepMap_RNA_Aneu$UBE2H..7328..x, DepMap_RNA_Aneu$Aneu.ploidy)
# 0.0007549
# 0.1316111
# n= 652


# UBE2H expression by 7q status 
ggplot(data = DepMap_RNA_Aneu, # remove samples where 7q CN status is NA
       aes(y=log2(UBE2H..7328..x), x=as.factor(X7q)))+
  geom_boxplot()+
  xlab("7q Copy Number")+ 
  ylab("UBE2H expression")+
  theme_classic()
# DepMap.UBE2H.RNA.7qCNboxplot.pdf
# Figure S7M
mean(subset(DepMap_RNA_Aneu, X7q == "0" )$UBE2H..7328..x) - mean(subset(DepMap_RNA_Aneu, X7q == "1" )$UBE2H..7328..x)
x<- t.test(subset(DepMap_RNA_Aneu, X7q == "0")$UBE2H..7328..x, subset(DepMap_RNA_Aneu, X7q == "1")$UBE2H..7328..x)
x$p.value # neutral gain p= 6.047952e-12
mean(subset(DepMap_RNA_Aneu, X7q == "0")$UBE2H..7328..x) - mean(subset(DepMap_RNA_Aneu, X7q == "-1" )$UBE2H..7328..x)
x<- t.test(subset(DepMap_RNA_Aneu, X7q == "0")$UBE2H..7328..x, subset(DepMap_RNA_Aneu, X7q == "-1")$UBE2H..7328..x)
x$p.value # neutral loss 0.02158043
 

# correlation of UBE2H expression and aneuploidy score in 7q neutral cells only
ggplot(subset(DepMap_RNA_Aneu, X7q==0), aes(y=log2(UBE2H..7328..x), x=Aneu.ploidy))+
  geom_point()+
  geom_smooth(method="lm")+
  ylab("UBE2H Expression")+ 
  xlab("DepMap: Aneuploidy Score")+
  theme_classic()
# DepMap.UBE2H.RNA.Aneuploidy.ploidy_Neutral7q.pdf
# Figure S7N
cor.test(log2(subset(DepMap_RNA_Aneu, X7q==0)$UBE2H..7328..x), subset(DepMap_RNA_Aneu, X7q==0)$Aneu.ploidy)
# p=0.01246 
# corr= 0.1352008 




 ##### DepMap Aneuploidy score & Protein expression ####
# Correlate Aneuploidy score with protein abundance in DepMap: 
# are mitochondrial metabolism proteins downregulated in aneuploid cells? 
ColNames.Proteins<- colnames(DepMap_Protein_Aneu)[7:length(DepMap_Protein_Aneu)]

DepMap.Corr.Aneu.Protein<- data.frame(Protein_Id= character(), 
                               Corr= numeric(), 
                               pvalue= numeric())
for (i in ColNames.Proteins){ 
  x<- cor.test(DepMap_Protein_Aneu$Aneu.ploidy, DepMap_Protein_Aneu[,i])
  DepMap.Corr.Aneu.Protein<- rbind(DepMap.Corr.Aneu.Protein, data.frame(Protein_Id= i, 
                                                          Corr= x$estimate, 
                                                          pvalue= x$p.value ))
}
#get gene name from protein ID
DepMap.Corr.Aneu.Protein$Protein_Id2<- gsub("\\|",".", DepMap.Corr.Aneu.Protein$Protein_Id)
DepMap.Corr.Aneu.Protein$Protein_Id2[1:10]
DepMap.Corr.Aneu.Protein<- merge(Protein_ProID, DepMap.Corr.Aneu.Protein, by.x="Protein_Id", by.y="Protein_Id2")


# Plot protein abundance correlation with  
ggplot(DepMap.Corr.Aneu.Protein, aes(x=Corr))+
  geom_density(data = subset(DepMap.Corr.Aneu.Protein, Gene_Symbol %in% RNA_Processing_rRNA ), color="deepskyblue4", size=2) + #00BFC4
  geom_density(data = subset(DepMap.Corr.Aneu.Protein, Gene_Symbol %in% Ribosomal_Genes_noMT$Approved.symbol ), color="#C77CFF", size=2) + #C77CFF
  geom_density(data = subset(DepMap.Corr.Aneu.Protein, Gene_Symbol %in% Proteosome_Genes$Approved.symbol ), color="chartreuse1", size=2) + #7CAE00
  geom_density(data = subset(DepMap.Corr.Aneu.Protein, Gene_Symbol %in% RNA_Processing_Spliceosome ), color="deepskyblue1", size=2) + #00BFC4
  geom_density(data = subset(DepMap.Corr.Aneu.Protein, Gene_Symbol %in% MitoTranslationTranscription ), color="gold2", size=2) +
  geom_density(data = subset(DepMap.Corr.Aneu.Protein, Gene_Symbol %in% MitoElectronChainnAssembly), color = "#F8766D", size=2) + 
  geom_density(data = subset(DepMap.Corr.Aneu.Protein, !Gene_Symbol %in% 
                               c(MitoElectronChainnAssembly, MitoTranslationTranscription, RNA_Processing_rRNA, Ribosomal_Genes_noMT$Approved.symbol, Proteosome_Genes$Approved.symbol, RNA_Processing_Spliceosome)), 
               size=2, color = "black")+ # all other genes in black
  xlab("Correlation (DepMap):\n Aneuploidy score & protein abundance")+ 
  xlim(-0.5, 0.5)+  
  geom_hline(yintercept=0)+
  geom_vline(xintercept=0)+
  theme_classic()
# DepMap.proteinAbundance.Corr.AneuScore.pdf
# Supplementary Figure 3E

# t-test to see if DepMap protein abundance significantly correlates with aneuploidy score: 
t.test(subset(DepMap.Corr.Aneu.Protein, Gene_Symbol %in% RNA_Processing_rRNA)$Corr,  subset(DepMap.Corr.Aneu.Protein, ! Gene_Symbol %in% RNA_Processing_rRNA)$Corr)
# rRNA processing: p = 2.101e-09 downreg
t.test(subset(DepMap.Corr.Aneu.Protein, Gene_Symbol %in% Ribosomal_Genes_noMT$Approved.symbol)$Corr,  subset(DepMap.Corr.Aneu.Protein, ! Gene_Symbol %in% Ribosomal_Genes_noMT$Approved.symbol)$Corr)$p.value
# Ribosome: p = 5.682144e-24 downreg
t.test(subset(DepMap.Corr.Aneu.Protein, Gene_Symbol %in% Proteosome_Genes$Approved.symbol)$Corr,  subset(DepMap.Corr.Aneu.Protein, ! Gene_Symbol %in% Proteosome_Genes$Approved.symbol)$Corr)
# Proteosome: p = 1.359e-10 downreg
t.test(subset(DepMap.Corr.Aneu.Protein, Gene_Symbol %in% RNA_Processing_Spliceosome)$Corr,  subset(DepMap.Corr.Aneu.Protein, ! Gene_Symbol %in% RNA_Processing_Spliceosome)$Corr)
# Spliceosome: p = 0.05299 NS
t.test(subset(DepMap.Corr.Aneu.Protein, Gene_Symbol %in% MitoTranslationTranscription)$Corr,  subset(DepMap.Corr.Aneu.Protein, ! Gene_Symbol %in% MitoTranslationTranscription)$Corr)
# MitoTrans: p = 0.1353 NS
t.test(subset(DepMap.Corr.Aneu.Protein, Gene_Symbol %in% MitoElectronChainnAssembly)$Corr,  subset(DepMap.Corr.Aneu.Protein, ! Gene_Symbol %in% MitoElectronChainnAssembly)$Corr)
# MitoECA: p = 0.005027 downregulated in aneuploid cancer cell lines. 


 ##### DepMap Aneuploidy Score with RNA expression, gene groups ####
##  Find gene (groups) that are up or down regulated in aneuploid cells: 
setwd(Dependency)
DepMap_transcription<- read.csv("OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv", header=TRUE)

MitochiondrialGenome_Genes<- read.csv("group-1972.csv", skip = 1 , header=TRUE)#note: these genes are not in the DepMap or paired CRISPR data


# Merge aneuploidy score and DepMap CRISPR data: 
DepMap_expression_AneuScore<- merge(Aneuploidy_Scores, DepMap_transcription, by.x="DepMap_ID", by.y="ModelID")


#Plot a significant correlation gene, UBE2H: 
ggplot(DepMap_expression_AneuScore, aes(y= UBE2H..7328., x=Score.divPloidy))+
  geom_point()+
  theme_classic()+
  geom_smooth(method="lm")+
  xlab("Aneuploidy Score")+
  ylab("Expression")

cor.test(DepMap_expression_AneuScore$UBE2H..7328., DepMap_expression_AneuScore$Score.divPloidy, method = "pearson")
# cor =  0.1681761
# p-value= 1.395e-07




### Now calculate and plot gene beta scores with aneuploidy score

DepMap_RNA.An<-data.frame(GeneID= character(), 
                          Gene = character(), 
                          Mean.Expression = numeric(),
                          SampleNum = numeric(),
                          
                          Aneu.RNA.Corr.P= numeric(),
                          Aneu.RNA.Corr.Coef= numeric(), 
                          Aneu.divPloidy.RNA.Corr.P= numeric(),
                          Aneu.divPloidy.RNA.Corr.Coef= numeric(), 
                          RNA.Ploidy.P=numeric(),
                          RNA.Ploidy.Coef=numeric()
)

for (i in 11:length(DepMap_expression_AneuScore)){
  x<- cor.test(DepMap_expression_AneuScore[,i], DepMap_expression_AneuScore$Score.divPloidy, method = "pearson")
  y<- cor.test(DepMap_expression_AneuScore[,i], DepMap_expression_AneuScore$Aneuploidy.score, method = "pearson")
  z<- cor.test(DepMap_expression_AneuScore[,i], DepMap_expression_AneuScore$Ploidy, method = "pearson")
  
  DepMap_RNA.An<- rbind(DepMap_RNA.An, data.frame(GeneID= colnames(DepMap_expression_AneuScore)[i],
                                                  Gene = sub("[..].*", "", colnames(DepMap_expression_AneuScore)[i] ),
                                                  Mean.Expression = mean(na.omit(DepMap_expression_AneuScore[,i])),
                                                  SampleNum= length(DepMap_expression_AneuScore[,i][!is.na(DepMap_expression_AneuScore[,i])]), 
                                                  
                                                  Aneu.RNA.Corr.P= y$p.value, 
                                                  Aneu.RNA.Corr.Coef= y$estimate,
                                                  
                                                  AneudivPloidy.RNA.Corr.P= x$p.value, 
                                                  AneudivPloidy.RNA.Corr.Coef = x$estimate, 
                                                  
                                                  RNA.Ploidy.P= z$p.value, 
                                                  RNA.Ploidy.Coef = z$estimate) ) 
} 

# Length of Genes analyzed: 19215 Genes
length(DepMap_RNA.An$GeneID)
DepMap_RNA.An<- DepMap_RNA.An[order(DepMap_RNA.An$AneudivPloidy.RNA.Corr.Coef),]


DepMap_RNA.An$SigDiff<- "FALSE"
BonferroniCorrectedPvalue<- 0.05/length(DepMap_RNA.An$AneudivPloidy.RNA.Corr.P)
for (i in 1:length(DepMap_RNA.An$GeneID)){
  if (DepMap_RNA.An$AneudivPloidy.RNA.Corr.P[i] < BonferroniCorrectedPvalue ){
    DepMap_RNA.An$SigDiff[i]<- "TRUE"
  }
}

DepMap_RNA.An$Gene<- sub("\\.\\..*", "", DepMap_RNA.An$GeneID) #Two periods \\.\\. followed by any character .* replace with nothing
DepMap_RNA.An$Gene<- sub("[.]", "-", DepMap_RNA.An$Gene) # replace single period with -. 
head(DepMap_RNA.An)

length(DepMap_RNA.An$Gene)# 19215


# plot RNA correlation with aneuploidy score per group: 
DepMap_RNA.An$Group<- "Other"
for (i in 1:length(DepMap_RNA.An$Group)){
  if (DepMap_RNA.An$Gene[i] %in% MitochiondrialGenome_Genes$Approved.symbol ){
    DepMap_RNA.An$Group[i]<- "Mitochondrial genome"
  } else if (DepMap_RNA.An$Gene[i] %in% MitoElectronChainnAssembly){
    DepMap_RNA.An$Group[i]<- "Mito ETC"
  } else if (DepMap_RNA.An$Gene[i] %in% MitoTranslationTranscription){
    DepMap_RNA.An$Group[i]<- "Mito Transcription Translation"
  }else if (DepMap_RNA.An$Gene[i] %in% RNA_Processing_rRNA){
    DepMap_RNA.An$Group[i]<- "rRNA Processing"
  }else if (DepMap_RNA.An$Gene[i] %in% Ribosomal_Genes_noMT$Approved.symbol){
    DepMap_RNA.An$Group[i]<- "Ribosomal Genes"
  }else if (DepMap_RNA.An$Gene[i] %in% Proteosome_Genes$Approved.symbol){
    DepMap_RNA.An$Group[i]<- "Proteosome"
  }else if (DepMap_RNA.An$Gene[i] %in% RNA_Processing_Spliceosome){
    DepMap_RNA.An$Group[i]<- "Spliceosome"
  }
}
DepMap_RNA.An$Group<- factor(DepMap_RNA.An$Group, levels= c("Other", "Mitochondrial genome", "Mito ETC", "Mito Transcription Translation",  "rRNA Processing","Ribosomal Genes", "Proteosome", "Spliceosome"))

ggplot(DepMap_RNA.An, aes(y=Aneu.RNA.Corr.Coef, 
                          x=Group, fill= Group))+
  geom_hline(yintercept=0)+
  geom_boxplot(outlier.shape = NA)+
  scale_fill_manual(values=c("black", "red3", "#F8766D","gold2", "deepskyblue4","#C77CFF","chartreuse1","deepskyblue1"))+ 
  theme_classic()+
  ylim(-0.5,0.5)+
  xlab("")+
  ylab("Correlation (DepMap):\n Gene expression & Aneuploid Arms")
# 5x4
# Plot.DepMap.RNACorrAneuploidyArms_All.pdf


# Aneuploidy score:
ggplot(DepMap_RNA.An, aes(y=AneudivPloidy.RNA.Corr.Coef, 
                          x=Group, fill= Group))+
  geom_hline(yintercept=0)+
  geom_boxplot(outlier.shape = NA)+
  scale_fill_manual(values=c("black", "red3", "#F8766D","gold2", "deepskyblue4","#C77CFF","chartreuse1","deepskyblue1"))+ 
  theme_classic()+  
  ylim(-0.5,0.5)+
  xlab("")+
  ylab("Correlation (DepMap):\n Gene expression & Aneuploidy Score")
# 5x4
# Plot.DepMap.RNACorrPloidy_All.pdf
# Supplementary Figure S3F

t.test(subset(DepMap_RNA.An, Group == "Other")$AneudivPloidy.RNA.Corr.Coef, 
       subset(DepMap_RNA.An, Group == "Mitochondrial genome")$AneudivPloidy.RNA.Corr.Coef) # 1.379e-08 Higher
t.test(subset(DepMap_RNA.An, Group == "Other")$AneudivPloidy.RNA.Corr.Coef, 
       subset(DepMap_RNA.An, Group == "Mito ETC")$AneudivPloidy.RNA.Corr.Coef) # 0.09557 mean lower NS
t.test(subset(DepMap_RNA.An, Group == "Other")$AneudivPloidy.RNA.Corr.Coef, 
       subset(DepMap_RNA.An, Group == "Mito Transcription Translation")$AneudivPloidy.RNA.Corr.Coef) # 0.3712, mean higher

t.test(subset(DepMap_RNA.An, Group == "Other")$AneudivPloidy.RNA.Corr.Coef, 
       subset(DepMap_RNA.An, Group == "rRNA Processing")$AneudivPloidy.RNA.Corr.Coef) #5.902e-08, mean lower
t.test(subset(DepMap_RNA.An, Group == "Other")$AneudivPloidy.RNA.Corr.Coef, 
       subset(DepMap_RNA.An, Group == "Ribosomal Genes")$AneudivPloidy.RNA.Corr.Coef)$p.value #8.258149e-20
t.test(subset(DepMap_RNA.An, Group == "Other")$AneudivPloidy.RNA.Corr.Coef, 
       subset(DepMap_RNA.An, Group == "Proteosome")$AneudivPloidy.RNA.Corr.Coef) # 0.001938, mean higher
t.test(subset(DepMap_RNA.An, Group == "Other")$AneudivPloidy.RNA.Corr.Coef, 
       subset(DepMap_RNA.An, Group == "Spliceosome")$AneudivPloidy.RNA.Corr.Coef) # 7.159e-08, mean lower

subset(DepMap_RNA.An, Group == "Mitochondrial genome" & SigDiff==TRUE)
subset(DepMap_RNA.An, Group == "Mito ETC" & SigDiff==TRUE)
subset(DepMap_RNA.An, Group == "Mito Translation Transcription" & SigDiff==TRUE)
subset(DepMap_RNA.An, Group == "rRNA Processing" & SigDiff==TRUE)
subset(DepMap_RNA.An, Group == "Ribosomal Genes" & SigDiff==TRUE)
subset(DepMap_RNA.An, Group == "Proteosome" & SigDiff==TRUE)
subset(DepMap_RNA.An, Group == "Spliceosome" & SigDiff==TRUE)

#examples of MT gene correlation: 
subset(DepMap_RNA.An, Group == "Mitochondrial genome")





ggplot(DepMap_expression_AneuScore, aes(y=MT.ND4L..4539., 
                                        x=Score.divPloidy))+
  geom_point()+
  theme_classic()+
  geom_smooth(method="lm")+
  xlab("Aneuploidy Score")+
  ylab("MT-ND4L Expression")
# plot.DepMap.expression.AneuScore.pdf
cor.test(DepMap_expression_AneuScore$MT.ND4L..4539., 
         DepMap_expression_AneuScore$Score.divPloidy, method = "pearson")$p.value
# cor =  0.2917194
# p-value= 1.834657e-20


ggplot(DepMap_expression_AneuScore, aes(y= (MT.ND4L..4539.), 
                                        x=Ploidy))+
  geom_point()+
  theme_classic()+
  geom_smooth(method="lm")+
  xlab("Ploidy")+
  ylab("MT-ND4L Expression")
# plot.DepMap.expression.Ploidy.pdf
cor.test(DepMap_expression_AneuScore$MT.ND4L..4539., 
         DepMap_expression_AneuScore$Ploidy, method = "pearson")
# cor =  0.08301323
# p-value= 0.009731





# plot UBE2H RNA correlation with mitochondrial ETC groups: 
DepMap_RNA.An$Group2<- "Other"
for (i in 1:length(DepMap_RNA.An$Group)){
  if (DepMap_RNA.An$Gene[i] %in% c("NDUFS1","NDUFS2","NDUFS3","NDUFS7","NDUFS8","NDUFV1","NDUFV2", "NDUFAB1","NDUFA1",
                                   "NDUFA2","NDUFA3","NDUFA5","NDUFA6","NDUFA7","NDUFA8","NDUFA9","NDUFA10","NDUFA11","NDUFA12",
                                   "NDUFA13","NDUFB1","NDUFB2","NDUFB3","NDUFB4","NDUFB5","NDUFB6","NDUFB7","NDUFB8","NDUFB9","NDUFB10",
                                   "NDUFB11","NDUFC1","NDUFC2","NDUFS4","NDUFS5","NDUFS6","NDUFV3") ){
    DepMap_RNA.An$Group2[i]<- "ETC I"
  } else if (DepMap_RNA.An$Gene[i] %in% c("SDHA", "SDHB", "SDHC", "SDHD")){
    DepMap_RNA.An$Group2[i]<- "ETC II"
  } else if (DepMap_RNA.An$Gene[i] %in% c("UQCRB","UQCRQ","UQCRC1","UQCRC2","MT-CYB","CYC1","UQCRFS1","UQCRH","UQCR10","UQCR11")){
    DepMap_RNA.An$Group2[i]<- "ETC III"
  }else if (DepMap_RNA.An$Gene[i] %in% c("COX4I1","COX4I2","COX5A","COX5B","COX6A1", "COX6A2","COX6B1","COX6B2","COX6C", "COX7A1","COX7A2","COX7B","COX7B2","COX7C",
                                         "COX8A","COX8C","MT-CO1","MT-CO2","MT-CO3")){
    DepMap_RNA.An$Group2[i]<- "ETC IV"
  }else if (DepMap_RNA.An$Gene[i] %in% c("ATP5F1A", "ATP5F1B", "ATP5F1C",  "ATP5F1D", "ATP5F1E", "ATP5MC1", "ATP5MC2", "ATP5MC3", "ATP5ME", "ATP5MF", "ATP5MG", "ATP5MJ", "ATP5MK", 
                                         "MT-ATP6", "MT-ATP8", "ATP5PB", "ATP5PD", "ATP5PF", "ATP5PO", "ATP5IF1")){
    DepMap_RNA.An$Group2[i]<- "ETC V"
  }else if (DepMap_RNA.An$Gene[i] %in% c("ATPAF1","ATPAF2","BCS1L","CHCHD7","CMC1","CMC2", "COA1",  "COA3", "COA4","COA5", "COA6",  "COA7", "COA8",
                                         "COX7A2L","COX10", "COX11", "COX14","COX15", "COX16",  "COX17", "COX18", "COX19", "COX20",
                                         "DMAC1", "DMAC2", "FDXR", "FDX2", "FMC1", "FOXRED1","HCCS","HIGD1A","HIGD2A","LRPPRC", "LYRM7",
                                         "NDUFAF1", "NDUFAF2", "NDUFAF3", "NDUFAF4", "NDUFAF5", "NDUFAF6","NDUFAF7",  "NDUFAF8","NUBPL",
                                         "OCIAD2", "OXA1L", "PET100","PET117","PNKD", "SCO1", "SCO2",  "SDHAF1","SDHAF2","SDHAF3", "SDHAF4", "SFXN4",
                                         "SMIM20", "SURF1","TACO1", "TIMMDC1","TMEM70","TMEM126A", "TMEM177","TMEM186","TMEM223","TMEM242","TTC19", 
                                         "UQCC1","UQCC2","UQCC3", "UQCC4","UQCC5","UQCC6", 
                                         "ACAD9","COA1","ECSIT","NDUFAF1","TMEM126B","TMEM186")){
    DepMap_RNA.An$Group2[i]<- "Assembly"
  }else if (DepMap_RNA.An$Gene[i] %in% c("BOLA3", "FDXR", "FXN", "GLRX5", "IBA57", "ISCA1", "ISCA2", "ISCU", "LYRM4", "NFS1", "NFU1", "NUBPL")){
    DepMap_RNA.An$Group2[i]<- "Iron-Sulfide Assembly" 
  } else if  (DepMap_RNA.An$Gene[i] %in% c("MRPL1", "MRPL2", "MRPL3", "MRPL4", "MRPL9", "MRPL10", "MRPL11", "MRPL12", "MRPL13", "MRPL14", "MRPL15", "MRPL16", "MRPL17", 
                                           "MRPL18", "MRPL19", "MRPL20", "MRPL21", "MRPL22", "MRPL23", "MRPL24", "MRPL27", "MRPL28", "MRPL30", "MRPL32", "MRPL33", "MRPL34", "MRPL35", "MRPL36", 
                                           "MRPL37", "MRPL38", "MRPL39", "MRPL40", "MRPL41", "MRPL42", "MRPL43", "MRPL44", "MRPL45", "MRPL46", "MRPL47", "MRPL48", "MRPL49", "MRPL50", 
                                           "MRPL51", "MRPL52", "MRPL53", "MRPL54", "MRPL55", "MRPL58", "MRPL57", "GADD45GIP1", "MRPS30", "MRPS18A", "MRPS2", "MRPS5", "MRPS6", "MRPS7", "MRPS9", "MRPS10", 
                                           "MRPS11", "MRPS12", "MRPS14", "MRPS15", "MRPS16", "MRPS17", "MRPS18C", "MRPS21", "MRPS22", "MRPS23", "MRPS24", "MRPS25", "MRPS26", "MRPS27", 
                                           "MRPS28", "DAP3", "MRPS31", "MRPS33", "MRPS34", "MRPS35", "CHCHD1", "AURKAIP1", "PTCD3", "MRPS18B")){
    DepMap_RNA.An$Group2[i]<- "Mitochondrial ribosome"
  } else if  (DepMap_RNA.An$Gene[i] %in% c("FASTK", "FASTKD1", "FASTKD2", "FASTKD3", "TBRG4", "FASTKD5", "MTRFR", "MTRF1", "MTRF1L", "POLRMT", "TFAM", "TFB2M","TRMT10C","HSD17B10","PRORP")){
    DepMap_RNA.An$Group2[i]<- "Mitochondrial transcription"
  }
}


DepMap_RNA.An$Group2<- factor(DepMap_RNA.An$Group2, levels= c("Other", "ETC I", "ETC II", "ETC III",  "ETC IV","ETC V", "Assembly", "Iron-Sulfide Assembly", "Mitochondrial ribosome", "Mitochondrial transcription"))

ggplot(DepMap_RNA.An, aes(y=AneudivPloidy.RNA.Corr.Coef, 
                          x=Group2, fill= Group2))+
  geom_hline(yintercept= median(subset(DepMap_RNA.An, Group2=="Other")$AneudivPloidy.RNA.Corr.Coef, na.rm=TRUE) )+
  geom_boxplot()+
  scale_fill_manual(values=c("black", "#F8766D","red1","red2","red3","red4", "firebrick2", "firebrick3", "gold2","gold4"))+ 
  theme_classic()+
  xlab("")+
  ylab("Correlation (DepMap):\n Gene expression & Aneuploidy Score")
# 5x4
# Plot.DepMap.RNACorAn_MitoType.pdf

t.test(subset(DepMap_RNA.An, Group2 == "Other")$AneudivPloidy.RNA.Corr.Coef, 
       subset(DepMap_RNA.An, Group2 == "ETC I")$AneudivPloidy.RNA.Corr.Coef) #  NS
t.test(subset(DepMap_RNA.An, Group2 == "Other")$AneudivPloidy.RNA.Corr.Coef, 
       subset(DepMap_RNA.An, Group2 == "ETC II")$AneudivPloidy.RNA.Corr.Coef) #  NS
t.test(subset(DepMap_RNA.An, Group2 == "Other")$AneudivPloidy.RNA.Corr.Coef, 
       subset(DepMap_RNA.An, Group2 == "ETC III")$AneudivPloidy.RNA.Corr.Coef) # NS
t.test(subset(DepMap_RNA.An, Group2 == "Other")$AneudivPloidy.RNA.Corr.Coef, 
       subset(DepMap_RNA.An, Group2 == "ETC IV")$AneudivPloidy.RNA.Corr.Coef) # NS
t.test(subset(DepMap_RNA.An, Group2 == "Other")$AneudivPloidy.RNA.Corr.Coef, 
       subset(DepMap_RNA.An, Group2 == "ETC V")$AneudivPloidy.RNA.Corr.Coef) # NS
t.test(subset(DepMap_RNA.An, Group2 == "Other")$AneudivPloidy.RNA.Corr.Coef, 
       subset(DepMap_RNA.An, Group2 == "Assembly")$AneudivPloidy.RNA.Corr.Coef) #NS
t.test(subset(DepMap_RNA.An, Group2 == "Other")$AneudivPloidy.RNA.Corr.Coef, 
       subset(DepMap_RNA.An, Group2 == "Iron-Sulfide Assembly")$AneudivPloidy.RNA.Corr.Coef)$p.value #  NS
t.test(subset(DepMap_RNA.An, Group2 == "Other")$AneudivPloidy.RNA.Corr.Coef, 
       subset(DepMap_RNA.An, Group2 == "Mitochondrial ribosome")$AneudivPloidy.RNA.Corr.Coef)$p.value #  NS
t.test(subset(DepMap_RNA.An, Group2 == "Other")$AneudivPloidy.RNA.Corr.Coef, 
       subset(DepMap_RNA.An, Group2 == "Mitochondrial transcription")$AneudivPloidy.RNA.Corr.Coef)$p.value #  NS



 ##### DepMap Aneuploidy Score with Protein expression, gene groups ####
##  Find genes that are up or down regulated in aneuploid cells: 


# Merge aneuploidy score and DepMap CRISPR data: 
DepMap_Prot_Aneu<- merge(Aneuploidy_Scores, ID_Protein_effect_4, by.x="DepMap_ID", by.y="Broad_ID")


#Plot a significant correlation gene, UBE2H: 
ggplot(DepMap_Prot_Aneu, aes(x= Aneuploidy.score, y=`sp|P03886|NU1M_HUMAN`))+
  geom_point()+
  theme_classic()+
  geom_smooth(method="lm")+
  xlab("Aneuploidy Score")+
  ylab("Protein ")

cor.test(DepMap_Prot_Aneu$Score.divPloidy, 
         DepMap_Prot_Aneu$`sp|P03886|NU1M_HUMAN`, method = "pearson")
# cor =  -
# p-value= 




### Now calculate and plot gene beta scores with aneuploidy score

DepMap_Prot.Cor.Aneu<-data.frame(ProteinID= character(), 
                                  Protein = character(), 
                                  Mean.Expression = numeric(),
                                  SampleNum = numeric(),
                                  
                                  Prot.Aneu.P=numeric(),
                                  Prot.Aneu.Coef=numeric()
)

for (i in 7:length(DepMap_Prot_Aneu)){
  x<- cor.test(DepMap_Prot_Aneu[,i], DepMap_Prot_Aneu$Score.divPloidy, method = "pearson")
  y<- cor.test(DepMap_Prot_Aneu[,i], DepMap_Prot_Aneu$Aneuploidy.score, method = "pearson")
  z<- cor.test(DepMap_Prot_Aneu[,i], DepMap_Prot_Aneu$Ploidy, method = "pearson")
  
  DepMap_Prot.Cor.Aneu<- rbind(DepMap_Prot.Cor.Aneu, data.frame(ProteinID= colnames(DepMap_Prot_Aneu)[i],
                                                                  Protein = sub(".*[|]", "", colnames(DepMap_Prot_Aneu)[i]) ,
                                                                  Mean.Expression = mean(na.omit(DepMap_Prot_Aneu[,i])),
                                                                  SampleNum= length(DepMap_Prot_Aneu[,i][!is.na(DepMap_Prot_Aneu[,i])]), 
                                                                  
                                                                  Prot.AneudivPloidy.P= x$p.value, 
                                                                  Prot.AneudivPloidy.Coef= x$estimate, 
                                                                
                                                                Prot.AneuArm.P= y$p.value, 
                                                                Prot.AneuArm.Coef= y$estimate,
                                                                
                                                                Prot.Ploidy.P= z$p.value, 
                                                                Prot.Ploidy.Coef= z$estimate) ) 
} 

# Length of Genes analyzed: 19215 Genes
length(DepMap_Prot.Cor.Aneu$ProteinID)
DepMap_Prot.Cor.Aneu<- DepMap_Prot.Cor.Aneu[order(DepMap_Prot.Cor.Aneu$Prot.AneudivPloidy.Coef),]


DepMap_Prot.Cor.Aneu$ProteinID<- str_replace_all(DepMap_Prot.Cor.Aneu$ProteinID, "[|]", ".")
DepMap_Prot.Cor.Aneu$ProteinID<- str_replace_all(DepMap_Prot.Cor.Aneu$ProteinID, "[-]", ".")

DepMap_Prot.Cor.Aneu<- merge(Protein_ProID, DepMap_Prot.Cor.Aneu, by.x="Protein_Id", by.y="ProteinID")

DepMap_Prot.Cor.Aneu$Protein<- sub("_.*", "", DepMap_Prot.Cor.Aneu$Protein) # replace _HUMAN (_.*) with nothing
DepMap_Prot.Cor.Aneu<- subset(DepMap_Prot.Cor.Aneu, ! Gene_Symbol== "UBE2H") # remove UBE2H


DepMap_Prot.Cor.Aneu$SigDiff<- "FALSE"
for (i in 1:length(DepMap_Prot.Cor.Aneu$Protein)){
  if (DepMap_Prot.Cor.Aneu$Prot.AneudivPloidy.P[i]< (0.05/length(DepMap_Prot.Cor.Aneu$Prot.AneudivPloidy.P)) ){
    DepMap_Prot.Cor.Aneu$SigDiff[i]<- "TRUE"
  }
}

head(DepMap_Prot.Cor.Aneu)

length(DepMap_Prot.Cor.Aneu$Gene_Symbol)# 12723


# plot RNA correlation with aneuploidy score per group: 
DepMap_Prot.Cor.Aneu$Group<- "Other"
for (i in 1:length(DepMap_Prot.Cor.Aneu$Group)){
  if (DepMap_Prot.Cor.Aneu$Gene_Symbol[i] %in% MitochiondrialGenome_Genes$Approved.symbol ){
    DepMap_Prot.Cor.Aneu$Group[i]<- "Mitochondrial genome"
  } else if (DepMap_Prot.Cor.Aneu$Gene_Symbol[i] %in% MitoElectronChainnAssembly){
    DepMap_Prot.Cor.Aneu$Group[i]<- "Mito ETC"
  } else if (DepMap_Prot.Cor.Aneu$Gene_Symbol[i] %in% MitoTranslationTranscription){
    DepMap_Prot.Cor.Aneu$Group[i]<- "Mito Transcription Translation"
  }else if (DepMap_Prot.Cor.Aneu$Gene_Symbol[i] %in% RNA_Processing_rRNA){
    DepMap_Prot.Cor.Aneu$Group[i]<- "rRNA Processing"
  }else if (DepMap_Prot.Cor.Aneu$Gene_Symbol[i] %in% Ribosomal_Genes_noMT$Approved.symbol){
    DepMap_Prot.Cor.Aneu$Group[i]<- "Ribosomal Genes"
  }else if (DepMap_Prot.Cor.Aneu$Gene_Symbol[i] %in% Proteosome_Genes$Approved.symbol){
    DepMap_Prot.Cor.Aneu$Group[i]<- "Proteosome"
  }else if (DepMap_Prot.Cor.Aneu$Gene_Symbol[i] %in% RNA_Processing_Spliceosome){
    DepMap_Prot.Cor.Aneu$Group[i]<- "Spliceosome"
  }
}
DepMap_Prot.Cor.Aneu$Group<- factor(DepMap_Prot.Cor.Aneu$Group, 
                                     levels= c("Other", "Mitochondrial genome", "Mito ETC", 
                                               "Mito Transcription Translation",  "rRNA Processing",
                                               "Ribosomal Genes", "Proteosome", "Spliceosome"))

ggplot(DepMap_Prot.Cor.Aneu, aes(y=Prot.AneudivPloidy.Coef, 
                                  x=Group, fill= Group))+
  geom_hline(yintercept= 0 )+
  geom_boxplot(outlier.shape = NA)+
  scale_fill_manual(values=c("black", "red3", "#F8766D","gold2", "deepskyblue4","#C77CFF","chartreuse1","deepskyblue1"))+ 
  theme_classic()+
  ylim(-0.5, 0.5)+
  xlab("")+
  ylab("Correlation (DepMap):\n Protein & Aneuploid Chromosomes")
# 5x4
# Plot.DepMap.ProteinCorAneuScore_All.pdf
# Supplementary Figure S3G


t.test(subset(DepMap_Prot.Cor.Aneu, Group == "Other")$Prot.AneudivPloidy.Coef, 
       subset(DepMap_Prot.Cor.Aneu, Group == "Mitochondrial genome")$Prot.AneudivPloidy.Coef) # 0.04177 * 7 = NS
t.test(subset(DepMap_Prot.Cor.Aneu, Group == "Other")$Prot.AneudivPloidy.Coef, 
       subset(DepMap_Prot.Cor.Aneu, Group == "Mito ETC")$Prot.AneudivPloidy.Coef)$p.value # 2.060426e-05 * 7 = 0.0001442298
t.test(subset(DepMap_Prot.Cor.Aneu, Group == "Other")$Prot.AneudivPloidy.Coef, 
       subset(DepMap_Prot.Cor.Aneu, Group == "Mito Transcription Translation")$Prot.AneudivPloidy.Coef)$p.value # 0.07138248 NS * 7 = NS

t.test(subset(DepMap_Prot.Cor.Aneu, Group == "Other")$Prot.AneudivPloidy.Coef, 
       subset(DepMap_Prot.Cor.Aneu, Group == "rRNA Processing")$Prot.AneudivPloidy.Coef)$p.value #1.791241e-15 * 7 = 1.253869e-14
t.test(subset(DepMap_Prot.Cor.Aneu, Group == "Other")$Prot.AneudivPloidy.Coef, 
       subset(DepMap_Prot.Cor.Aneu, Group == "Ribosomal Genes")$Prot.AneudivPloidy.Coef)$p.value # 7.391463e-19 * 7 = 5.174024e-18
t.test(subset(DepMap_Prot.Cor.Aneu, Group == "Other")$Prot.AneudivPloidy.Coef, 
       subset(DepMap_Prot.Cor.Aneu, Group == "Proteosome")$Prot.AneudivPloidy.Coef)$p.value # 3.388955e-10 * 7 = 2.372269e-09
t.test(subset(DepMap_Prot.Cor.Aneu, Group == "Other")$Prot.AneudivPloidy.Coef, 
       subset(DepMap_Prot.Cor.Aneu, Group == "Spliceosome")$Prot.AneudivPloidy.Coef)$p.value # 0.002145613 * 7 = 0.01501929


subset(DepMap_Prot.Cor.Aneu, Group == "Mitochondrial genome" & SigDiff==TRUE)
subset(DepMap_Prot.Cor.Aneu, Group == "Mito ETC" & SigDiff==TRUE) 
subset(DepMap_Prot.Cor.Aneu, Group == "Mito Transcription Translation" & SigDiff==TRUE)
subset(DepMap_Prot.Cor.Aneu, Group == "rRNA Processing" & SigDiff==TRUE)
subset(DepMap_Prot.Cor.Aneu, Group == "Ribosomal Genes" & SigDiff==TRUE) #RPL22L1
subset(DepMap_Prot.Cor.Aneu, Group == "Proteosome" & SigDiff==TRUE)
subset(DepMap_Prot.Cor.Aneu, Group == "Spliceosome" & SigDiff==TRUE)


ggplot(DepMap_Prot_Aneu, aes(y=DepMap_Prot_Aneu$Score.divPloidy, 
                              x=DepMap_Prot_Aneu$`sp|Q00059|TFAM_HUMAN` ))+
  geom_point()+
  theme_classic()+
  geom_smooth(method="lm")+
  xlab("Aneuploidy Score")+
  ylab("TFAM Protein")
# plot.DepMap.Protein.Aneu.TFAM.pdf
cor.test(DepMap_Prot_Aneu$`sp|Q00059|TFAM_HUMAN`, 
         DepMap_Prot_Aneu$Score.divPloidy, method = "pearson")
# cor =  
# p-value= 



#### Bokenkamp: Gene expression in adapted aneuploid clones (passage 0 versus passage 50) ####
## Author: Klaske Schukken
# Analyze data from Bokenkamp et al. 2025
# 2025.07.11

# for Suplementary Figure S3 

#### INTRO AND LIBRARIES ####

## Folders 
## !! Update these locations to pathway where data was downloaded: 
# Folder with CRISPR screening and proteomics datasets: 
SchukkenData<- "/Volumes/Schukken_SSD/CRISPR SCREEN PAPER/R Code Data Files/Data"
# Folder with dependency datasets: 
Dependency<- "/Volumes/Schukken_SSD/CRISPR SCREEN PAPER/R Code Data Files/Dependency files"
# Folder with results: 
ResultsFile<- "/Volumes/Schukken_SSD/CRISPR SCREEN PAPER/R Code Data Files/Results"


### libraries
library(ggplot2) 
library(reshape2)
library(tidyverse)
library(readxl)
library(ggpubr)
library("cowplot")
library(plyr)
library('gprofiler2')
library(xlsx)
library("viridis")  # color packet
library(readr)

#### Get data ####
# Analyzing Bokenkamp et al 2025 for mitochondrial metabolism gene expression changes in p50 aneuploid clones
setwd(Dependency)

Bokenkamp_protein<- read_excel("44318_2025_372_moesm5_esm.xlsx", 
                               sheet = "cell_line_comparisons")


Proteosome_Genes<- read.csv("group-690.csv", header=TRUE)

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
Ribosomal_Genes<- read.csv("group-1054.csv", skip = 1 , header=TRUE) #HGNC Ribosomal Proteins
Ribosomal_Genes_noMT<- subset(Ribosomal_Genes, Group %in% c("L ribosomal proteins", "S ribosomal proteins"))
RNA_Processing<- read_tsv("GO_term_summary_20250720_134600.txt")
RNA_Processing1<-unique(toupper(RNA_Processing$Symbol))
RNA_Processing_Spliceosome<- unique(toupper(subset(RNA_Processing, `Annotated Term` %in% 
                                                     c("mRNA splicing, via spliceosome", 
                                                       "regulation of mRNA splicing, via spliceosome") )$Symbol ) )
RNA_Processing_rRNA<- unique(toupper(subset(RNA_Processing, `Annotated Term` %in% 
                                              c("rRNA processing", "regulation of rRNA processing"))$Symbol))




### Plot Protein dosage changes P0 compared to WT: #### 
## P0 aneuploid vs WT

Bokenkamp_protein$meanFCp0WT<- (Bokenkamp_protein$`FC_Hte5-HCT116` +
                                  Bokenkamp_protein$`FC_Htr5-HCT116` +
                                  Bokenkamp_protein$`FC_Rtr5-RPE1` +
                                  Bokenkamp_protein$`FC_Rtr21-RPE1`)/4 
Bokenkamp_protein$meanpvalue_p0WT<- (Bokenkamp_protein$`p_value_Htr5-HCT116` +
                                       Bokenkamp_protein$`p_value_Hte5-HCT116` +
                                       Bokenkamp_protein$`p_value_Rtr5-RPE1` +
                                       Bokenkamp_protein$`p_value_Rtr21-RPE1`)/4

ggplot(Bokenkamp_protein, aes(y=- log2(meanpvalue_p0WT), 
                              x= (meanFCp0WT)))+
  geom_point(color= "Black")+ 
  geom_density2d(color="white")+
  geom_point(data = subset(Bokenkamp_protein, hgnc_symbol %in% MitoElectronChainnAssembly), color="#F8766D")+ #F8766D
  geom_point(data = subset(Bokenkamp_protein, hgnc_symbol %in% MitoTranslationTranscription), color="gold2")+   #7CAE00
  theme_classic()+
  theme(legend.position="none")+
  xlab(paste("log2 protein abundance fold change Aneu to WT"))+
  ylab(paste("-log2(p-value)"))+
  ggtitle(paste("Foldchange protein dosage Aneu to WT \n Bokenkamp et al. (2025) "))
# Boekenkamp.Protein.p0WT.MitoMet.pdf
# Figure S3G


ggplot(Bokenkamp_protein, aes(y=- log2(meanpvalue_p0WT), 
                              x= (meanFCp0WT)))+
  geom_point(color= "Black")+ 
  geom_density2d(color="white")+
  geom_point(data = subset(Bokenkamp_protein, hgnc_symbol %in% RNA_Processing_rRNA ), color="deepskyblue4") + #00BFC4
  geom_point(data = subset(Bokenkamp_protein, hgnc_symbol %in% Ribosomal_Genes_noMT$Approved.symbol ), color="#C77CFF") + #C77CFF
  #geom_point(data = subset(Bokenkamp_protein, hgnc_symbol %in% Proteosome_Genes$Approved.symbol ), color="chartreuse1") + #7CAE00
  #geom_point(data = subset(Bokenkamp_protein, hgnc_symbol %in% RNA_Processing_Spliceosome ), color="deepskyblue1") + #00BFC4
  theme_classic()+
  theme(legend.position="none")+
  xlab(paste("log2 protein abundance fold change Aneu to WT"))+
  ylab(paste("-log2(p-value)"))+
  ggtitle(paste("Foldchange protein dosage Aneu to WT \n Bokenkamp et al. (2025) "))
# Boekenkamp.Protein.p0WT.RiborRNA.pdf
# Figure S3D

t.test(subset(Bokenkamp_protein, ! hgnc_symbol %in% Proteosome_Genes$Approved.symbol)$meanFCp0WT, 
       subset(Bokenkamp_protein, hgnc_symbol %in% Proteosome_Genes$Approved.symbol)$meanFCp0WT ) # NS
t.test(subset(Bokenkamp_protein, ! hgnc_symbol %in% RNA_Processing_Spliceosome)$meanFCp0WT, 
       subset(Bokenkamp_protein, hgnc_symbol %in% RNA_Processing_Spliceosome)$meanFCp0WT ) # 2.955e-13 (downregulation)

t.test(subset(Bokenkamp_protein, ! hgnc_symbol %in% MitoElectronChainnAssembly)$meanFCp0WT, 
       subset(Bokenkamp_protein, hgnc_symbol %in% MitoElectronChainnAssembly)$meanFCp0WT ) # 1.222e-08 (downregulation)
t.test(subset(Bokenkamp_protein, ! hgnc_symbol %in% MitoTranslationTranscription)$meanFCp0WT, 
       subset(Bokenkamp_protein, hgnc_symbol %in% MitoTranslationTranscription)$meanFCp0WT ) #4.441e-05 (downregulation)
x<-t.test(subset(Bokenkamp_protein, ! hgnc_symbol %in% RNA_Processing_rRNA)$meanFCp0WT, 
          subset(Bokenkamp_protein, hgnc_symbol %in% RNA_Processing_rRNA)$meanFCp0WT ) # 9.496386e-23 (downregulation)
x
x$p.value
x<-t.test(subset(Bokenkamp_protein, ! hgnc_symbol %in% Ribosomal_Genes_noMT$Approved.symbol)$meanFCp0WT, 
          subset(Bokenkamp_protein, hgnc_symbol %in% Ribosomal_Genes_noMT$Approved.symbol)$meanFCp0WT ) #1.208278e-21 (downregulation)
x
x$p.value


# newly aneuploid clones downregulate Spliceosome protein, Ribosomal proteins, rRNA processing, and mitochondrial metabolism proteins, relative to WT
# VERY Strong downregulation of Ribosome and rRNA



### Plot Protein dosage changes P50 compared to p0: #### 


# P50 relative to P0. Adaptation to aneuploidy at P50? Protein expression changes? 

Bokenkamp_protein$meanFCp50p0<- (Bokenkamp_protein$`FC_Htr5p50-Htr5` +
                                        Bokenkamp_protein$`FC_Hte5p50-Hte5` +
                                        Bokenkamp_protein$`FC_Rtr5p50-Rtr5` +
                                        Bokenkamp_protein$`FC_Rtr21p50-Rtr21`)/4 
Bokenkamp_protein$meanpvalue_p50p0<- (Bokenkamp_protein$`p_value_Htr5p50-Htr5` +
                                        Bokenkamp_protein$`p_value_Hte5p50-Hte5` +
                                        Bokenkamp_protein$`p_value_Rtr5p50-Rtr5`+
                                        Bokenkamp_protein$`p_value_Rtr21p50-Rtr21`)/4

ggplot(Bokenkamp_protein, aes(y=- log2(relevance_p_value), 
                                 x= (meanFCp50p0)))+
  geom_point(color= "Black")+ 
  geom_density_2d(color="white", bins=10)+
  geom_point(data = subset(Bokenkamp_protein, hgnc_symbol %in% MitoElectronChainnAssembly), color="#F8766D")+ #F8766D
  geom_point(data = subset(Bokenkamp_protein, hgnc_symbol %in% MitoTranslationTranscription), color="gold2")+   #7CAE00
  theme_classic()+
  theme(legend.position="none")+
  xlab(paste("log2 protein abundance fold change P50 to P0"))+
  ylab(paste("-log2(p-value)"))+
  ggtitle(paste("Foldchange protein dosage p50 to p0 \n Bokenkamp et al. (2025) "))
# Boekenkamp.Protein.p50p0.MitoMetab.pdf
# Figure S3H

ggplot(Bokenkamp_protein, aes(y=- log2(relevance_p_value), 
                                 x= (meanFCp50p0)))+
  geom_point(color= "Black")+ 
  geom_density2d(color="white")+
  geom_point(data = subset(Bokenkamp_protein, hgnc_symbol %in% RNA_Processing_rRNA ), color="deepskyblue4") + #00BFC4
  geom_point(data = subset(Bokenkamp_protein, hgnc_symbol %in% Ribosomal_Genes_noMT$Approved.symbol ), color="#C77CFF") + #C77CFF
  #geom_point(data = subset(Bokenkamp_protein, hgnc_symbol %in% Proteosome_Genes$Approved.symbol ), color="chartreuse1") + #7CAE00
  #geom_point(data = subset(Bokenkamp_protein, hgnc_symbol %in% RNA_Processing_Spliceosome ), color="deepskyblue1") + #00BFC4
  theme_classic()+
  theme(legend.position="none")+
  xlab(paste("log2 protein abundance fold change P50 to P0"))+
  ylab(paste("-log2(p-value)"))+
  ggtitle(paste("Foldchange protein dosage p50 to p0 \n Bokenkamp et al. (2025) "))
# Boekenkamp.Protein.p50p0.RiborRNA.pdf
# Figure S3E

#t-test for pathways: 
t.test(subset(Bokenkamp_protein, ! hgnc_symbol %in% Proteosome_Genes$Approved.symbol)$meanFCp50p0, 
       subset(Bokenkamp_protein, hgnc_symbol %in% Proteosome_Genes$Approved.symbol)$meanFCp50p0 ) #NS
t.test(subset(Bokenkamp_protein, ! hgnc_symbol %in% RNA_Processing_Spliceosome)$meanFCp50p0, 
       subset(Bokenkamp_protein, hgnc_symbol %in% RNA_Processing_Spliceosome)$meanFCp50p0 ) # 0.001218 (positive)

t.test(subset(Bokenkamp_protein, ! hgnc_symbol %in% MitoElectronChainnAssembly)$meanFCp50p0, 
          subset(Bokenkamp_protein, hgnc_symbol %in% MitoElectronChainnAssembly)$meanFCp50p0 ) #0.0002962 (negative) 
t.test(subset(Bokenkamp_protein, ! hgnc_symbol %in% MitoTranslationTranscription)$meanFCp50p0, 
           subset(Bokenkamp_protein, hgnc_symbol %in% MitoTranslationTranscription)$meanFCp50p0 ) # 0.01696 (negative)
t.test(subset(Bokenkamp_protein, ! hgnc_symbol %in% RNA_Processing_rRNA)$meanFCp50p0, 
           subset(Bokenkamp_protein, hgnc_symbol %in% RNA_Processing_rRNA)$meanFCp50p0 ) # 3.243e-10 (positive)
t.test(subset(Bokenkamp_protein, ! hgnc_symbol %in% Ribosomal_Genes_noMT$Approved.symbol)$meanFCp50p0, 
           subset(Bokenkamp_protein, hgnc_symbol %in% Ribosomal_Genes_noMT$Approved.symbol)$meanFCp50p0 ) # 5.694e-06 (Positive) 

# Proteosome NS
# Spliceosome 0.001218 P
# Mito electron transport chain 0.0002962 N
# Mito transcription translation 0.01696 N
# rRNA processing: 3.243e-10 P
# Ribosome: 5.694e-06 P


# Check other metabolic pathways: 
ggplot(Bokenkamp_protein, aes(y=- log2(relevance_p_value), 
                              x= (meanFCp50p0)))+
  geom_point(color= "Black")+ 
  geom_density2d(color="white")+  
  #geom_point(data = subset(Bokenkamp_protein, hgnc_symbol %in% MitoElectronChainnAssembly), color="#F8766D")+ #F8766D
  theme_classic()+
  theme(legend.position="none")+
  xlab(paste("log2 protein abundance fold change P50 to P0"))+
  ylab(paste("-log2(p-value)"))+
  ggtitle(paste("Foldchange protein dosage p50 to p0 \n Bokenkamp et al. (2025) "))
# Boekenkamp.Protein.p50p0.Metabolism.pdf



# Ribosomes, rRNA processing show a strong protein upregulation at p50. adaptation to aneuploidy
# Spliceosome genes also show a slight upregulated at p50 compared to p0. adaptation to aneuploidy
# Mitochondrial metabolism genes show a slight but significant downregulation at p50 relative to p0.  (more at p0) 
# general (slight) downregulation on metabolism proteins




### Plot Protein dosage changes P50 compared to WT: #### 
Bokenkamp_protein$meanFCp50WT<- (Bokenkamp_protein$`FC_Hte5p50-HCT116` +
                                     Bokenkamp_protein$`FC_Htr5p50-HCT116` +
                                     Bokenkamp_protein$`FC_Rtr5p50-RPE1`+
                                     Bokenkamp_protein$`FC_Rtr21p50-RPE1`)/4 
Bokenkamp_protein$meanpvalue_p50WT<- (Bokenkamp_protein$`p_value_Htr5p50-HCT116`+
                                          Bokenkamp_protein$`p_value_Hte5p50-HCT116` +
                                          Bokenkamp_protein$`p_value_Rtr5p50-RPE1` +
                                          Bokenkamp_protein$`p_value_Rtr21p50-RPE1`)/4

ggplot(Bokenkamp_protein, aes(y=- log2(meanpvalue_p50WT), 
                                 x= (meanFCp50WT)))+
  geom_point(color= "Black")+ 
  geom_density2d(color="white")+
  geom_point(data = subset(Bokenkamp_protein, hgnc_symbol %in% MitoElectronChainnAssembly), color="#F8766D")+ #F8766D
  geom_point(data = subset(Bokenkamp_protein, hgnc_symbol %in% MitoTranslationTranscription), color="gold2")+   #7CAE00
  theme_classic()+
  theme(legend.position="none")+
  xlab(paste("log2 protein abundance fold change Aneup50 to WT"))+
  ylab(paste("-log2(p-value)"))+
  ggtitle(paste("Foldchange protein dosage Aneup50 to WT \n Bokenkamp et al. (2025) "))
# Boekenkamp.Protein.p50WT.MitoMet.pdf
# Figure S3I

ggplot(Bokenkamp_protein, aes(y=- log2(meanpvalue_p50WT), 
                                 x= (meanFCp50WT)))+
  geom_point(color= "Black")+ 
  geom_density2d(color="white")+
  geom_point(data = subset(Bokenkamp_protein, hgnc_symbol %in% RNA_Processing_rRNA ), color="deepskyblue4") + #00BFC4
  geom_point(data = subset(Bokenkamp_protein, hgnc_symbol %in% Ribosomal_Genes_noMT$Approved.symbol ), color="#C77CFF") + #C77CFF
  #geom_point(data = subset(Bokenkamp_protein, hgnc_symbol %in% Proteosome_Genes$Approved.symbol ), color="chartreuse1") + #7CAE00
  #geom_point(data = subset(Bokenkamp_protein, hgnc_symbol %in% RNA_Processing_Spliceosome ), color="deepskyblue1") + #00BFC4
  theme_classic()+
  theme(legend.position="none")+
  xlab(paste("log2 protein abundance fold change Aneup50 to WT"))+
  ylab(paste("-log2(p-value)"))+
  ggtitle(paste("Foldchange protein dosage Aneup50 to WT \n Bokenkamp et al. (2025) "))
# Boekenkamp.Protein.p50WT.RiborRNA.pdf
# Figure S3F

t.test(subset(Bokenkamp_protein, ! hgnc_symbol %in% Proteosome_Genes$Approved.symbol)$meanFCp50WT, 
       subset(Bokenkamp_protein, hgnc_symbol %in% Proteosome_Genes$Approved.symbol)$meanFCp50WT ) #0.01145 (upreg?) (not major hit)
t.test(subset(Bokenkamp_protein, ! hgnc_symbol %in% RNA_Processing_Spliceosome)$meanFCp50WT, 
       subset(Bokenkamp_protein, hgnc_symbol %in% RNA_Processing_Spliceosome)$meanFCp50WT ) # 1.648e-05 (negative)

t.test(subset(Bokenkamp_protein, ! hgnc_symbol %in% MitoElectronChainnAssembly)$meanFCp50WT, 
          subset(Bokenkamp_protein, hgnc_symbol %in% MitoElectronChainnAssembly)$meanFCp50WT ) # 1.055e-10 (negative)
t.test(subset(Bokenkamp_protein, ! hgnc_symbol %in% MitoTranslationTranscription)$meanFCp50WT, 
           subset(Bokenkamp_protein, hgnc_symbol %in% MitoTranslationTranscription)$meanFCp50WT ) #2.357e-07 (negative)
t.test(subset(Bokenkamp_protein, ! hgnc_symbol %in% RNA_Processing_rRNA)$meanFCp50WT, 
           subset(Bokenkamp_protein, hgnc_symbol %in% RNA_Processing_rRNA)$meanFCp50WT ) # 8.338e-11 (negative)
x<- t.test(subset(Bokenkamp_protein, ! hgnc_symbol %in% Ribosomal_Genes_noMT$Approved.symbol)$meanFCp50WT, 
           subset(Bokenkamp_protein, hgnc_symbol %in% Ribosomal_Genes_noMT$Approved.symbol)$meanFCp50WT ) #1.338554e-16 (negative)
x$p.value
x

# Strong downregulation of Spliceosome, Ribosome, rRNA processing, and mitochondrial metabolism genes at p50 aneuploid relative to WT 
# Proteosome NS
### UBE2H ####

Bokenkamp_protein$meanH.p50p0<- rowMeans(Bokenkamp_protein[,c(12,17)])
Bokenkamp_protein$meanR.p50p0<- rowMeans(Bokenkamp_protein[,c(22,27)])

Bokenkamp_protein$meanH.p0WT<- rowMeans(Bokenkamp_protein[,c(32,37)])
Bokenkamp_protein$meanR.p0WT<- rowMeans(Bokenkamp_protein[,c(42,47)])

Bokenkamp_protein$meanH.p50WT<- rowMeans(Bokenkamp_protein[,c(52,57)])
Bokenkamp_protein$meanR.p50WT<- rowMeans(Bokenkamp_protein[,c(62,67)])

ggplot(Bokenkamp_protein, aes(y= meanH.p50p0, 
                              x= meanR.p50p0 ))+
  geom_point(color= "Black")+ 
  geom_density2d(color="white")+
  geom_point(data = subset(Bokenkamp_protein, hgnc_symbol %in% c("UBE2H")), color="red")+   #7CAE00
  theme_classic()+
  geom_hline(yintercept = 0)+ #add lines at y=0
  geom_vline(xintercept = 0)+ #add lines at x=0
  theme(legend.position="none")+
  xlab(paste("RPE1"))+
  ylab(paste("HCT116"))+
  ggtitle(paste("Foldchange protein p50 to p0 "))
# Boekenkamp.Protein.p50p0.MitoMet.pdf


ggplot(Bokenkamp_protein, aes(y= meanH.p0WT, 
                              x= meanR.p0WT ))+
  geom_point(color= "Black")+ 
  geom_density2d(color="white")+
  geom_point(data = subset(Bokenkamp_protein, hgnc_symbol %in% c("UBE2H")), color="red")+   #7CAE00
  theme_classic()+
  geom_hline(yintercept = 0)+ #add lines at y=0
  geom_vline(xintercept = 0)+ #add lines at x=0
  theme(legend.position="none")+
  xlab(paste("RPE1"))+
  ylab(paste("HCT116"))+
  ggtitle(paste("Foldchange protein p0 to WT "))
# Boekenkamp.Protein.p0WT.MitoMet.pdf


ggplot(Bokenkamp_protein, aes(y= meanH.p50WT, 
                              x= meanR.p50WT ))+
  geom_point(color= "Black")+ 
  geom_density2d(color="white")+
  geom_point(data = subset(Bokenkamp_protein, hgnc_symbol %in% c("UBE2H")), color="red")+   #7CAE00
  theme_classic()+
  geom_hline(yintercept = 0)+ #add lines at y=0
  geom_vline(xintercept = 0)+ #add lines at x=0
  theme(legend.position="none")+
  xlab(paste("RPE1"))+
  ylab(paste("HCT116"))+
  ggtitle(paste("Foldchange protein p50 to WT "))
# Boekenkamp.Protein.p50WT.MitoMet.pdf



ggplot(Bokenkamp_protein, aes(x= meanH.p0WT, 
                              y= meanH.p50WT ))+
  geom_point(color= "Black")+ 
  geom_density2d(color="white")+
  geom_point(data = subset(Bokenkamp_protein, hgnc_symbol %in% c("UBE2H")), color="red")+   #7CAE00
  theme_classic()+
  geom_hline(yintercept = 0)+ #add lines at y=0
  geom_vline(xintercept = 0)+ #add lines at x=0
  theme(legend.position="none")+
  xlab(paste("HCT116 P0"))+
  ylab(paste("HCT116 P50"))+
  ggtitle(paste("Foldchange protein Aneuploid to WT "))
# Boekenkamp.Protein.p50WT.MitoMet.pdf




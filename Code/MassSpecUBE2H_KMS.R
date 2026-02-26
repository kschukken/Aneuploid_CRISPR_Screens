##### Mass Spec analysis: UBE2HKO and AZ3146 ####

### 2026-01-13
# Author: Klaske M. Schukken 

# Mass Spec data (raw) was analyzed by Yale Proteomics Core and imported into *** software. 
# Analysis on protein expression per sample was performed and a weight-based p-value was assigned to 
# calculate differences between gene expression per group. 
# Protein expression per sample, average differences per group and corresponding p-values were exported from ***. 


# AZ3146 is an MPS1 inhibitor drug shown to induce chromosome missegregation and lead to aneuploidy. 

# Three protein samples of each group was analyzed: 
## HCT116 WT
## HCT116 UBE2H-KO clone 2
## HCT116 WT + AZ3146 (96 hour treatment, 1uM) 
## HCT116 UBE2H-KO + AZ3146 (96 hour treatment, 1uM) 


##### Intro and libraries ##### 
# Folders
## !! Update these locations to pathway where data was downloaded: 
# Folder with CRISPR screening and proteomics datasets: 
SchukkenData<- "/Volumes/Schukken_SSD/CRISPR SCREEN PAPER/R Code Data Files/Data"
# Folder with dependency datasets: 
Dependency<- "/Volumes/Schukken_SSD/CRISPR SCREEN PAPER/R Code Data Files/Dependency files"
# Folder with results: 
ResultsFile<- "/Volumes/Schukken_SSD/CRISPR SCREEN PAPER/R Code Data Files/Results"




# libraries 
library(ggplot2) 
library(reshape2)
library(tidyverse)
library(readxl)
library(ggpubr)
library(plyr)
library('gprofiler2')
library(xlsx)
library("viridis")  # color packet
library(readr)




##### Get Data #####

setwd(SchukkenData)
MS<- read_excel("Data_FC_pvalue_HCT116.UBE2H.AZ3146.xlsx")


## Get pathways ## 

# MitoElectronChainnAssembly             color="#F8766D"
# MitoTranslationTranscription           color="gold2"
# RNA_Processing_rRNA                    color="deepskyblue4"
# Ribosomal_Genes_noMT$Approved.symbol   color="#C77CFF"
# Proteosome_Genes$Approved.symbol       color="chartreuse1"
# RNA_Processing_Spliceosome             color="deepskyblue1"

setwd(Dependency)



Proteosome_Genes<- read.csv("group-690.csv", header=TRUE)
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



##### Format data #####
head(MS)
colnames(MS) #191 columns
# [1] "Gene Symbol...1" 
# [7] "Accession" 
# [8]  "Description"  
# [77] "GO Accessions"                                             
# [78]"Entrez Gene ID"                                           
# [79] "Ensembl Gene ID"  

# [94] "Abundance Ratio (log2): (UBE2HKO_AZ) / (AZ)"              
#[95] "Abundance Ratio (log2): (UBE2HKO_AZ) / (UBE2HKO)"          
#[96] "Abundance Ratio (log2): (AZ) / (WT)"                      
#[97] "Abundance Ratio (log2): (UBE2HKO) / (WT)" 

#[102] "Abundance Ratio Adj. P-Value: (UBE2HKO_AZ) / (AZ)"        
#[103] "Abundance Ratio Adj. P-Value: (UBE2HKO_AZ) / (UBE2HKO)"    
#[104] "Abundance Ratio Adj. P-Value: (AZ) / (WT)"         
#[105] "Abundance Ratio Adj. P-Value: (UBE2HKO) / (WT)"     

# [114]"Abundances (Grouped): AZ"                                 [115] "Abundances (Grouped): UBE2HKO"
# [116] "Abundances (Grouped): UBE2HKO_AZ"                        [117] "Abundances (Grouped): WT"  

# [126] "Abundance: F4: Sample, AZ"                                
#[127] "Abundance: F6: Sample, AZ"                                 "Abundance: F9: Sample, AZ"                                
#[129] "Abundance: F3: Sample, UBE2HKO"                            "Abundance: F7: Sample, UBE2HKO"                           
#[131] "Abundance: F12: Sample, UBE2HKO"                           "Abundance: F2: Sample, UBE2HKO_AZ"                        
#[133] "Abundance: F8: Sample, UBE2HKO_AZ"                         "Abundance: F10: Sample, UBE2HKO_AZ"                       
#[135] "Abundance: F1: Sample, WT"                                 "Abundance: F5: Sample, WT"                                
#[137] "Abundance: F11: Sample, WT" 

MS2<-MS[,c(1,7,8,78,79,94,95,96,97,102,103,104,105,114,115,116,117,126:137) ]
MS2<- data.frame(MS2)



### Plot UBE2H abundance ##### 
## check to see if UBE2H is knocked out in UBE2H-KO clones 
# it is much lower or undetected in my UBE2H-KO clones. 
# UBE2H-KO 2/3 undetected, UBE2H-KO +AZ3146 1/3 undetected. 

UBE2HMS<- subset(MS2, Gene.Symbol...1=="UBE2H")
UBE2HMS<- t(UBE2HMS)
UBE2HMS[is.na(UBE2HMS)] <- 0
UBE2HMS<- UBE2HMS[18:29]

PlotUBE2H<- data.frame(Sample =  factor(c("AZ3146", "AZ3146", "AZ3146", "UBE2H-KO", "UBE2H-KO", "UBE2H-KO", 
                                    "UBE2HKO + AZ3145", "UBE2HKO + AZ3145", "UBE2HKO + AZ3145", "WT", "WT", "WT"), level=c("WT", "AZ3146", "UBE2H-KO", "UBE2HKO + AZ3145")),
  Abundance =  as.numeric(UBE2HMS), 
  Detected = c(TRUE,TRUE,TRUE, FALSE,FALSE,TRUE, TRUE,TRUE,FALSE, TRUE,TRUE,TRUE)) 
PlotUBE2H

ggplot(PlotUBE2H, aes(x=Sample, 
                y= Abundance ))+
  geom_boxplot()+
  geom_point()+
  geom_point(data = subset(PlotUBE2H, Detected == FALSE), color = "Red")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+ 
  xlab(paste(""))+
  ylab(paste("UBE2H Abundance"))+
  ggtitle(paste("UBE2H Abundance per sample"))
# plot.MS.UBE2H_Abundance.pdf
# 4x5 portrait
# Supplementary Figure S6F




### Plot pathways ####
colnames(MS2)

#AZ3146
ggplot(MS2, aes(x=Abundance.Ratio..log2....AZ.....WT., 
                  y= -log2(Abundance.Ratio.Adj..P.Value...AZ.....WT.)))+ 
  geom_point()+ 
  theme_classic()+
  geom_point(data = subset(MS2, Gene.Symbol...1 %in% MitoElectronChainnAssembly), color = "#F8766D") + 
  geom_point(data = subset(MS2, Gene.Symbol...1 %in% MitoTranslationTranscription ), color="gold2") + 
  #geom_point(data = subset(MS2, Gene.Symbol...1 %in% RNA_Processing_rRNA ), color="deepskyblue4") + #00BFC4
  #geom_point(data = subset(MS2, Gene.Symbol...1 %in% Ribosomal_Genes_noMT$Approved.symbol ), color="#C77CFF") + #C77CFF
  #geom_point(data = subset(MS2, Gene.Symbol...1 %in% Proteosome_Genes$Approved.symbol ), color="chartreuse1") + #7CAE00
  #geom_point(data = subset(MS2, Gene.Symbol...1 %in% RNA_Processing_Spliceosome ), color="deepskyblue1") + #00BFC4
  geom_density2d(color="white")+
  xlab(paste("Abundance Ratio (log2): AZ/WT"))+
  ylab(paste("-log2(p-value)"))+
  ggtitle(paste("Protein abundance changes AZ3146"))
# plot.MS.AZ.WT_Mito.pdf
# Supplementary Figure S6G

# 4x4 

# Significance 0.01/24 = 0.0004166667 (< 4E-4) 
t.test(subset(MS2, Gene.Symbol...1 %in% MitoElectronChainnAssembly)$Abundance.Ratio..log2....AZ.....WT., 
       subset(MS2, ! Gene.Symbol...1 %in% MitoElectronChainnAssembly)$Abundance.Ratio..log2....AZ.....WT. ) # NS diff =-0.0132556
t.test(subset(MS2, Gene.Symbol...1 %in% MitoTranslationTranscription)$Abundance.Ratio..log2....AZ.....WT., 
       subset(MS2, ! Gene.Symbol...1 %in% MitoTranslationTranscription)$Abundance.Ratio..log2....AZ.....WT. ) # slight Downreg 2.688e-05, diff = -0.1325763
t.test(subset(MS2, Gene.Symbol...1 %in% RNA_Processing_rRNA)$Abundance.Ratio..log2....AZ.....WT., 
       subset(MS2, ! Gene.Symbol...1 %in% RNA_Processing_rRNA)$Abundance.Ratio..log2....AZ.....WT. ) # NS
t.test(subset(MS2, Gene.Symbol...1 %in% Ribosomal_Genes_noMT$Approved.symbol)$Abundance.Ratio..log2....AZ.....WT., 
       subset(MS2, ! Gene.Symbol...1 %in% Ribosomal_Genes_noMT$Approved.symbol)$Abundance.Ratio..log2....AZ.....WT. ) # STRONG Downreg 6.167589e-18, diff = -0.1793359
x$p.value
x
t.test(subset(MS2, Gene.Symbol...1 %in% Proteosome_Genes$Approved.symbol)$Abundance.Ratio..log2....AZ.....WT., 
       subset(MS2, ! Gene.Symbol...1 %in% Proteosome_Genes$Approved.symbol)$Abundance.Ratio..log2....AZ.....WT. ) # downreg? 0.000705 NS
t.test(subset(MS2, Gene.Symbol...1 %in% RNA_Processing_Spliceosome)$Abundance.Ratio..log2....AZ.....WT., 
       subset(MS2, ! Gene.Symbol...1 %in% RNA_Processing_Spliceosome)$Abundance.Ratio..log2....AZ.....WT. ) #  downreg? 0.0309 NS
x$p.value


# UBE2H-KO 
ggplot(MS2, aes(x=Abundance.Ratio..log2....UBE2HKO.....WT., 
                y= -log2(Abundance.Ratio.Adj..P.Value...UBE2HKO.....WT.)))+
  geom_point()+ 
  theme_classic()+
  geom_point(data = subset(MS2, Gene.Symbol...1 %in% MitoElectronChainnAssembly), color = "#F8766D") + 
  geom_point(data = subset(MS2, Gene.Symbol...1 %in% MitoTranslationTranscription ), color="gold2") + 
  #geom_point(data = subset(MS2, Gene.Symbol...1 %in% RNA_Processing_rRNA ), color="deepskyblue4") + #00BFC4
  #geom_point(data = subset(MS2, Gene.Symbol...1 %in% Ribosomal_Genes_noMT$Approved.symbol ), color="#C77CFF") + #C77CFF
  #geom_point(data = subset(MS2, Gene.Symbol...1 %in% Proteosome_Genes$Approved.symbol ), color="chartreuse1") + #7CAE00
  #geom_point(data = subset(MS2, Gene.Symbol...1 %in% RNA_Processing_Spliceosome ), color="deepskyblue1") + #00BFC4
  geom_density2d(color="white")+
  xlab(paste("Abundance Ratio (log2): UBE2H-KO/WT"))+
  ylab(paste("-log2(p-value)"))+
  ggtitle(paste("Protein abundance changes UBE2H-KO"))
# plot.MS.UBE2HKO.WT_Mito.pdf
# Supplementary Figure S6H
# 4x4
# Mitochondrial down-regulation
# spliceosome upregulation

# Significance 0.01/24 = 0.0004166667 (< 4E-4) 
t.test(subset(MS2, Gene.Symbol...1 %in% MitoElectronChainnAssembly)$Abundance.Ratio..log2....UBE2HKO.....WT., 
       subset(MS2, ! Gene.Symbol...1 %in% MitoElectronChainnAssembly)$Abundance.Ratio..log2....UBE2HKO.....WT. ) # Downreg 1.315e-05, diff = -0.4736259
t.test(subset(MS2, Gene.Symbol...1 %in% MitoTranslationTranscription)$Abundance.Ratio..log2....UBE2HKO.....WT., 
       subset(MS2, ! Gene.Symbol...1 %in% MitoTranslationTranscription)$Abundance.Ratio..log2....UBE2HKO.....WT. ) # STRONG downreg 9.874215e-20, diff = -0.3778242
t.test(subset(MS2, Gene.Symbol...1 %in% RNA_Processing_rRNA)$Abundance.Ratio..log2....UBE2HKO.....WT., 
       subset(MS2, ! Gene.Symbol...1 %in% RNA_Processing_rRNA)$Abundance.Ratio..log2....UBE2HKO.....WT. ) # 0.006255 Upreg (NS)
t.test(subset(MS2, Gene.Symbol...1 %in% Ribosomal_Genes_noMT$Approved.symbol)$Abundance.Ratio..log2....UBE2HKO.....WT., 
       subset(MS2, ! Gene.Symbol...1 %in% Ribosomal_Genes_noMT$Approved.symbol)$Abundance.Ratio..log2....UBE2HKO.....WT. ) # STRONG Downreg 8.627754e-17, diff = -0.2683223
t.test(subset(MS2, Gene.Symbol...1 %in% Proteosome_Genes$Approved.symbol)$Abundance.Ratio..log2....UBE2HKO.....WT., 
       subset(MS2, ! Gene.Symbol...1 %in% Proteosome_Genes$Approved.symbol)$Abundance.Ratio..log2....UBE2HKO.....WT. ) # Upreg 3.659e-05
x<-t.test(subset(MS2, Gene.Symbol...1 %in% RNA_Processing_Spliceosome)$Abundance.Ratio..log2....UBE2HKO.....WT., 
       subset(MS2, ! Gene.Symbol...1 %in% RNA_Processing_Spliceosome)$Abundance.Ratio..log2....UBE2HKO.....WT. ) # STRONG upreg 6.613802e-20
x$p.value


# UBE2H-KO + AZ3146 relative to AZ
ggplot(MS2, aes(x=Abundance.Ratio..log2....UBE2HKO_AZ.....AZ., 
                y= -log2(Abundance.Ratio.Adj..P.Value...UBE2HKO_AZ.....AZ.)))+
  geom_point()+ 
  theme_classic()+
  geom_point(data = subset(MS2, Gene.Symbol...1 %in% MitoElectronChainnAssembly), color = "#F8766D") + 
  geom_point(data = subset(MS2, Gene.Symbol...1 %in% MitoTranslationTranscription ), color="gold2") + 
  #geom_point(data = subset(MS2, Gene.Symbol...1 %in% RNA_Processing_rRNA ), color="deepskyblue4") + #00BFC4
  #geom_point(data = subset(MS2, Gene.Symbol...1 %in% Ribosomal_Genes_noMT$Approved.symbol ), color="#C77CFF") + #C77CFF
  #geom_point(data = subset(MS2, Gene.Symbol...1 %in% Proteosome_Genes$Approved.symbol ), color="chartreuse1") + #7CAE00
  #geom_point(data = subset(MS2, Gene.Symbol...1 %in% RNA_Processing_Spliceosome ), color="deepskyblue1") + #00BFC4
  geom_density2d(color="white")+
  xlab(paste("Abundance Ratio (log2): UBE2H-KO+AZ/AZ"))+
  ylab(paste("-log2(p-value)"))+
  ggtitle(paste("Protein abundance changes UBE2H-KO + AZ relative to AZ"))
# plot.MS.UBE2HKOAZ.AZ_Pathways.pdf
# plot.MS.UBE2HKOAZ.AZ_Mito.pdf
# Supplementary Figure S6I
# 4x4
# Mitochondrial down-regulation
# spliceosome upregulation

# Significance 0.01/24 = 0.0004166667 (< 4E-4) 
t.test(subset(MS2, Gene.Symbol...1 %in% MitoElectronChainnAssembly)$Abundance.Ratio..log2....UBE2HKO_AZ.....AZ., 
       subset(MS2, ! Gene.Symbol...1 %in% MitoElectronChainnAssembly)$Abundance.Ratio..log2....UBE2HKO_AZ.....AZ. ) # Downreg 8.754e-06.  diff = -0.3680787
t.test(subset(MS2, Gene.Symbol...1 %in% MitoTranslationTranscription)$Abundance.Ratio..log2....UBE2HKO_AZ.....AZ., 
       subset(MS2, ! Gene.Symbol...1 %in% MitoTranslationTranscription)$Abundance.Ratio..log2....UBE2HKO_AZ.....AZ. ) # Downreg 1.123e-10   diff= -0.2628351
t.test(subset(MS2, Gene.Symbol...1 %in% RNA_Processing_rRNA)$Abundance.Ratio..log2....UBE2HKO_AZ.....AZ., 
       subset(MS2, ! Gene.Symbol...1 %in% RNA_Processing_rRNA)$Abundance.Ratio..log2....UBE2HKO_AZ.....AZ. ) # 0.002081 Upreg (NS)
t.test(subset(MS2, Gene.Symbol...1 %in% Ribosomal_Genes_noMT$Approved.symbol)$Abundance.Ratio..log2....UBE2HKO_AZ.....AZ., 
       subset(MS2, ! Gene.Symbol...1 %in% Ribosomal_Genes_noMT$Approved.symbol)$Abundance.Ratio..log2....UBE2HKO_AZ.....AZ. ) # STRONG Downreg 7.459e-16, diff = -0.2945533
t.test(subset(MS2, Gene.Symbol...1 %in% Proteosome_Genes$Approved.symbol)$Abundance.Ratio..log2....UBE2HKO_AZ.....AZ., 
       subset(MS2, ! Gene.Symbol...1 %in% Proteosome_Genes$Approved.symbol)$Abundance.Ratio..log2....UBE2HKO_AZ.....AZ. ) # Upreg 8.782e-06
t.test(subset(MS2, Gene.Symbol...1 %in% RNA_Processing_Spliceosome)$Abundance.Ratio..log2....UBE2HKO_AZ.....AZ., 
       subset(MS2, ! Gene.Symbol...1 %in% RNA_Processing_Spliceosome)$Abundance.Ratio..log2....UBE2HKO_AZ.....AZ. ) # upreg 1.504e-10




## Effect of AZ in UBE2HKO: UBE2HKO +AZ/ UBE2HKO
ggplot(MS2, aes(x=Abundance.Ratio..log2....UBE2HKO_AZ.....UBE2HKO., 
                y= -log2(Abundance.Ratio.Adj..P.Value...UBE2HKO_AZ.....UBE2HKO.)))+
  geom_point()+ 
  theme_classic()+
  geom_point(data = subset(MS2, Gene.Symbol...1 %in% MitoElectronChainnAssembly), color = "#F8766D") + 
  geom_point(data = subset(MS2, Gene.Symbol...1 %in% MitoTranslationTranscription ), color="gold2") + 
  #geom_point(data = subset(MS2, Gene.Symbol...1 %in% RNA_Processing_rRNA ), color="deepskyblue4") + #00BFC4
  #geom_point(data = subset(MS2, Gene.Symbol...1 %in% Ribosomal_Genes_noMT$Approved.symbol ), color="#C77CFF") + #C77CFF
  #geom_point(data = subset(MS2, Gene.Symbol...1 %in% Proteosome_Genes$Approved.symbol ), color="chartreuse1") + #7CAE00
  #geom_point(data = subset(MS2, Gene.Symbol...1 %in% RNA_Processing_Spliceosome ), color="deepskyblue1") + #00BFC4
  geom_density2d(color="white")+
  xlab(paste("Abundance Ratio (log2): UBE2H-KO+AZ/UBE2HKO"))+
  ylab(paste("-log2(p-value)"))+
  ggtitle(paste("Protein abundance changes UBE2H-KO + AZ relative to UBE2HKO"))
# plot.MS.UBE2HKOAZ.UBE2H_Pathways.pdf
# plot.MS.UBE2HKOAZ.UBE2H_Mito.pdf
# 4x4
# downregulation Ribosomes 
# upreg mitochondria

# Significance 0.01/24 = 0.0004166667 (< 4E-4) 
t.test(subset(MS2, Gene.Symbol...1 %in% MitoElectronChainnAssembly)$Abundance.Ratio..log2....UBE2HKO_AZ.....UBE2HKO., 
       subset(MS2, ! Gene.Symbol...1 %in% MitoElectronChainnAssembly)$Abundance.Ratio..log2....UBE2HKO_AZ.....UBE2HKO. ) # upreg 2.472e-05, difference = 0.114745
t.test(subset(MS2, Gene.Symbol...1 %in% MitoTranslationTranscription)$Abundance.Ratio..log2....UBE2HKO_AZ.....UBE2HKO., 
       subset(MS2, ! Gene.Symbol...1 %in% MitoTranslationTranscription)$Abundance.Ratio..log2....UBE2HKO_AZ.....UBE2HKO. ) # NS
t.test(subset(MS2, Gene.Symbol...1 %in% RNA_Processing_rRNA)$Abundance.Ratio..log2....UBE2HKO_AZ.....UBE2HKO., 
       subset(MS2, ! Gene.Symbol...1 %in% RNA_Processing_rRNA)$Abundance.Ratio..log2....UBE2HKO_AZ.....UBE2HKO. ) # NA
t.test(subset(MS2, Gene.Symbol...1 %in% Ribosomal_Genes_noMT$Approved.symbol)$Abundance.Ratio..log2....UBE2HKO_AZ.....UBE2HKO., 
       subset(MS2, ! Gene.Symbol...1 %in% Ribosomal_Genes_noMT$Approved.symbol)$Abundance.Ratio..log2....UBE2HKO_AZ.....UBE2HKO. ) # STRONG Downreg 4.311e-15
t.test(subset(MS2, Gene.Symbol...1 %in% Proteosome_Genes$Approved.symbol)$Abundance.Ratio..log2....UBE2HKO_AZ.....UBE2HKO., 
       subset(MS2, ! Gene.Symbol...1 %in% Proteosome_Genes$Approved.symbol)$Abundance.Ratio..log2....UBE2HKO_AZ.....UBE2HKO. ) #NS
t.test(subset(MS2, Gene.Symbol...1 %in% RNA_Processing_Spliceosome)$Abundance.Ratio..log2....UBE2HKO_AZ.....UBE2HKO., 
       subset(MS2, ! Gene.Symbol...1 %in% RNA_Processing_Spliceosome)$Abundance.Ratio..log2....UBE2HKO_AZ.....UBE2HKO. ) # 0.01229 lower



### Plot density pathways #####

ggplot(MS2, aes(x=Abundance.Ratio..log2....AZ.....WT.))+
  geom_density(data = subset(MS2, !Gene.Symbol...1 %in% c(MitoElectronChainnAssembly, MitoTranslationTranscription, Ribosomal_Genes_noMT$Approved.symbol)), linewidth =1)+ 
  theme_classic()+
  geom_density(data = subset(MS2, Gene.Symbol...1 %in% MitoElectronChainnAssembly), color = "#F8766D", linewidth =1) + 
  geom_density(data = subset(MS2, Gene.Symbol...1 %in% MitoTranslationTranscription ), color="gold2", linewidth =1) + 
  #geom_density(data = subset(MS2, Gene.Symbol...1 %in% Ribosomal_Genes_noMT$Approved.symbol ), color="#C77CFF", linewidth =1) + #C77CFF
  xlab(paste("Abundance Ratio (log2): AZ/WT"))+
  ylab(paste("Density"))+
  xlim(-2,2)+
  ggtitle(paste("Protein abundance changes AZ"))
# plot.MS.AZ_Pathways_Density.pdf
# 4x4
# Figure S6G

ggplot(MS2, aes(x=Abundance.Ratio..log2....UBE2HKO.....WT.))+
  geom_density(data = subset(MS2, !Gene.Symbol...1 %in% c(MitoElectronChainnAssembly, MitoTranslationTranscription, Ribosomal_Genes_noMT$Approved.symbol)), linewidth =1)+ 
  theme_classic()+
  geom_density(data = subset(MS2, Gene.Symbol...1 %in% MitoElectronChainnAssembly), color = "#F8766D",linewidth =1) + 
  geom_density(data = subset(MS2, Gene.Symbol...1 %in% MitoTranslationTranscription ), color="gold2",linewidth =1) + 
  #geom_density(data = subset(MS2, Gene.Symbol...1 %in% Ribosomal_Genes_noMT$Approved.symbol ), color="#C77CFF",linewidth =1) + #C77CFF
  xlab(paste("Abundance Ratio (log2): UBE2H-KO/WT"))+
  ylab(paste("Density"))+
  xlim(-2,2)+
  ggtitle(paste("Protein abundance changes UBE2H-KO"))
# plot.MS.UBE2HKO_Pathways_Density.pdf
# Figure S6H

ggplot(MS2, aes(x=Abundance.Ratio..log2....UBE2HKO_AZ.....AZ.))+
  geom_density(data = subset(MS2, !Gene.Symbol...1 %in% c(MitoElectronChainnAssembly, MitoTranslationTranscription, Ribosomal_Genes_noMT$Approved.symbol)), linewidth =1)+ 
  theme_classic()+
  geom_density(data = subset(MS2, Gene.Symbol...1 %in% MitoElectronChainnAssembly), color = "#F8766D", linewidth =1) + 
  geom_density(data = subset(MS2, Gene.Symbol...1 %in% MitoTranslationTranscription ), color="gold2", linewidth =1) + 
  #geom_density(data = subset(MS2, Gene.Symbol...1 %in% Ribosomal_Genes_noMT$Approved.symbol ), color="#C77CFF", linewidth =1) + #C77CFF
  xlab(paste("Abundance Ratio (log2): UBE2H-KO+AZ/AZ"))+
  ylab(paste("Density"))+
  xlim(-2,2)+
  ggtitle(paste("Protein abundance changes UBE2H-KO + AZ relative to AZ"))
# plot.MS.UBE2HKOAZ.AZ_Pathways_Density.pdf
# Figure S6I


ggplot(MS2, aes(x=Abundance.Ratio..log2....UBE2HKO_AZ.....UBE2HKO.))+
  geom_density(data = subset(MS2, !Gene.Symbol...1 %in% c(MitoElectronChainnAssembly, MitoTranslationTranscription, Ribosomal_Genes_noMT$Approved.symbol)), linewidth =1)+ 
  theme_classic()+
  geom_density(data = subset(MS2, Gene.Symbol...1 %in% MitoElectronChainnAssembly), color = "#F8766D",linewidth =1) + 
  geom_density(data = subset(MS2, Gene.Symbol...1 %in% MitoTranslationTranscription ), color="gold2",linewidth =1) + 
  #geom_density(data = subset(MS2, Gene.Symbol...1 %in% Ribosomal_Genes_noMT$Approved.symbol ), color="#C77CFF",linewidth =1) + #C77CFF
  xlab(paste("Abundance Ratio (log2): UBE2H-KO+AZ/UBE2HKO"))+
  ylab(paste("Density"))+
  xlim(-1,1)+
  ggtitle(paste("Protein abundance changes UBE2H-KO + AZ relative to UBE2HKO"))
# plot.MS.UBE2HKOAZ.UBE2HKO_Pathways_Density.pdf


### G-profiler ####

# Find top and bottom 15% of proteins for each comparison
MS2_UBE_Down <- quantile(MS2$Abundance.Ratio..log2....UBE2HKO.....WT., 0.20, na.rm = TRUE) #bottom 15 of protein <= this value
MS2_UBE_Up <- quantile(MS2$Abundance.Ratio..log2....UBE2HKO.....WT., 0.80, na.rm = TRUE) #Top 15 of protein >= this value

MS2_AZ_Down <- quantile(MS2$Abundance.Ratio..log2....AZ.....WT., 0.20, na.rm = TRUE) #bottom 15 of protein <= this value
MS2_AZ_Up <- quantile(MS2$Abundance.Ratio..log2....AZ.....WT., 0.80, na.rm = TRUE) #Top 15 of protein >= this value

MS2_UBEAZ.AZ_Down <- quantile(MS2$Abundance.Ratio..log2....UBE2HKO_AZ.....AZ., 0.20, na.rm = TRUE) #bottom 15 of protein <= this value
MS2_UBEAZ.AZ_Up <- quantile(MS2$Abundance.Ratio..log2....UBE2HKO_AZ.....AZ., 0.80, na.rm = TRUE) #Top 15 of protein >= this value

MS2_UBEAZ.UBE_Down <- quantile(MS2$Abundance.Ratio..log2....UBE2HKO_AZ.....UBE2HKO., 0.20, na.rm = TRUE) #bottom 15 of protein <= this value
MS2_UBEAZ.UBE_Up <- quantile(MS2$Abundance.Ratio..log2....UBE2HKO_AZ.....UBE2HKO., 0.80, na.rm = TRUE) #Top 15 of protein >= this value



# AZ3146 DOWNREGULATED PROTEINS 
# SIGNIFICANT: Ribosome, Mitochondrial Metabolism
# G-Profiler: 
g.MS.AZ.DOWNReg<- gost(
  subset(MS2, MS2$Abundance.Ratio..log2....AZ.....WT. < MS2_AZ_Down)$Gene,
  organism = "hsapiens",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = c("g_SCS"),
  domain_scope = "custom",
  custom_bg = as.character(MS2$Gene.Symbol...1),
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE
)

g.MS.AZ.DOWNReg<-g.MS.AZ.DOWNReg$result
g.MS.AZ.DOWNReg<-g.MS.AZ.DOWNReg[order(g.MS.AZ.DOWNReg$p_value),]

##Upregulated protein in UBE2HKO cells 
MS.AZ.D<-g.MS.AZ.DOWNReg[order(g.MS.AZ.DOWNReg$p_value),]
MS.AZ.D$term_name
MS.AZ.D<-MS.AZ.D[c(1,2,20, 3,4,7, 5,9),]
MS.AZ.D$term_name<- factor(MS.AZ.D$term_name, levels= MS.AZ.D$term_name)
MS.AZ.D$Termtype<-c( 
  "Ribosomal", "Ribosomal", "Ribosomal",
  "Mitochondrial Metabolism","Mitochondrial Metabolism","Mitochondrial Metabolism",
  "Mitochondrial Ribosome", "Mitochondrial Ribosome"
)

ggplot(MS.AZ.D, 
       aes(x=term_name, y=-log2(p_value), fill=Termtype))+
  geom_bar(stat="identity")+
  ylab("-log2 (P-Value)") +
  xlab("")+
  ggtitle ("Proteins enriched: Downregulated upon AZ3146 treatment")+
  theme_classic()+
  coord_flip()
## 
# plot.gProfiler.MassSpec.AZ3146.Downreg.pdf
# 10x5
# Figure 6H



# AZ3146 UPREGULATED 
# SIGNIFICANT: transcription, CIN, P53
# G-Profiler: 
g.MS.AZ.UPReg<- gost(
  subset(MS2, MS2$Abundance.Ratio..log2....AZ.....WT. > MS2_AZ_Up)$Gene,
  organism = "hsapiens",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = c("g_SCS"),
  domain_scope = "custom",
  custom_bg = as.character(MS2$Gene.Symbol...1),
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE
)

g.MS.AZ.UPReg<-g.MS.AZ.UPReg$result
g.MS.AZ.UPReg<-g.MS.AZ.UPReg[order(g.MS.AZ.UPReg$p_value),]

##Upregulated protein in AZ cells 
MS.AZ.U<-g.MS.AZ.UPReg[order(g.MS.AZ.UPReg$p_value),]
MS.AZ.U$term_name
MS.AZ.U<-MS.AZ.U[c(2,3,4, 11, 13, 14),]
MS.AZ.U$term_name<- factor(MS.AZ.U$term_name, levels= MS.AZ.U$term_name)
MS.AZ.U$Termtype<-c( 
  "Transcription", "Transcription", "Transcription","Cell Cycle",  "P53", "Chromosome Instability"
)

ggplot(MS.AZ.U, 
       aes(x=term_name, y=-log2(p_value), fill=Termtype))+
  geom_bar(stat="identity")+
  ylab("-log2 (P-Value)") +
  xlab("")+
  ggtitle ("Proteins enriched: Upregulated upon AZ3146 treatment")+
  theme_classic()+
  coord_flip()
## 
# plot.gProfiler.MassSpec.AZ3146.Upreg.pdf
# 10x5 
# Figure S6J




# UBE2H-KO DOWNREGULATED 
#  SIGNIFICANT: Mitochondria, RNA processing
# G-Profiler: 
g.MS.UBE2HKO.DownReg<- gost(
  subset(MS2, MS2$Abundance.Ratio..log2....UBE2HKO.....WT. < MS2_UBE_Down)$Gene,
  organism = "hsapiens",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = c("g_SCS"),
  domain_scope = "custom",
  custom_bg = as.character(MS2$Gene.Symbol...1),
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE
)

g.MS.UBE2HKO.DownReg2<-g.MS.UBE2HKO.DownReg$result
g.MS.UBE2HKO.DownReg2<-g.MS.UBE2HKO.DownReg2[order(g.MS.UBE2HKO.DownReg2$p_value),]

##Downregulated protein in UBE2HKO cells 
MS.UBE.D<-g.MS.UBE2HKO.DownReg2[order(g.MS.UBE2HKO.DownReg2$p_value),]
MS.UBE.D<- subset(MS.UBE.D, !source %in% c("HPA", "TF"))
MS.UBE.D$term_name[1:150]
MS.UBE.D<-MS.UBE.D[c(1,2,3, 60,76,83, 96,103),]
MS.UBE.D$term_name<- factor(MS.UBE.D$term_name, levels= MS.UBE.D$term_name)
MS.UBE.D$Termtype<-c("Mitochondrial Metabolism","Mitochondrial Metabolism","Mitochondrial Metabolism",
                     "RNA processing", "RNA processing", "RNA processing", 
                     "Mitochondrial Ribosome",  "Mitochondrial Ribosome"
)

ggplot(MS.UBE.D, 
       aes(x=term_name, y=-log2(p_value), fill=Termtype))+
  geom_bar(stat="identity")+
  ylab("-log2 (P-Value)") +
  xlab("")+
  ggtitle ("Proteins enriched: Downregulated upon UBE2H-KO")+
  theme_classic()+
  coord_flip()
## 
# plot.gProfiler.MassSpec.UBE2HKO.Downreg.pdf
# 10x5
# Figure 6I


# UBE2H-KO UPREGULATED 
# SIGNIFICANT: Cell cycle
# UBE2H-KO cells grew slower and were still in growth phase, while WT were near full. 
# G-Profiler: 

g.MS.UBE2HKO.UPReg<- gost(
  subset(MS2, MS2$Abundance.Ratio..log2....UBE2HKO.....WT. > MS2_UBE_Up)$Gene,
  organism = "hsapiens",
  ordered_query = TRUE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = c("g_SCS"),
  domain_scope = "custom",
  custom_bg = as.character(MS2$Gene.Symbol...1),
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE
)

g.MS.UBE2HKO.UPReg<-g.MS.UBE2HKO.UPReg$result
g.MS.UBE2HKO.UPReg<-g.MS.UBE2HKO.UPReg[order(g.MS.UBE2HKO.UPReg$p_value),]

##Upregulated protein in UBE2HKO cells 
MS.UBE.U<-g.MS.UBE2HKO.UPReg[order(g.MS.UBE2HKO.UPReg$p_value),]
MS.UBE.U<- subset(MS.UBE.U, !source %in% c("HPA", "TF"))
MS.UBE.U$term_name
MS.UBE.U<-MS.UBE.U[c(1,5,12,20, 23),]
MS.UBE.U$term_name<- factor(MS.UBE.U$term_name, levels= MS.UBE.U$term_name)
MS.UBE.U$Termtype<-c( 
  "Cell Cycle", "Cell Cycle", "Cell Cycle", "Cell Cycle", "Cell Cycle"
)

ggplot(MS.UBE.U, 
       aes(x=term_name, y=-log2(p_value), fill=Termtype))+
  geom_bar(stat="identity")+
  ylab("-log2 (P-Value)") +
  xlab("")+
  ggtitle ("Proteins enriched: Upregulated upon UBE2H-KO")+
  theme_classic()+
  coord_flip()
## 
# plot.gProfiler.MassSpec.UBE2HKO.Upreg.pdf
# 10x5
# Figure S6K





# UBE2H-KO in AZ3146 background DOWNREGULATED PROTEINS 
# SIGNIFICANT: Mitochondrial Metabolism, Ribosome
# G-Profiler: 
g.MS.UBE2HKOAZ.AZ.DOWNReg<- gost(
  subset(MS2, MS2$Abundance.Ratio..log2....UBE2HKO_AZ.....AZ. < MS2_UBEAZ.AZ_Down)$Gene,
  organism = "hsapiens",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = c("g_SCS"),
  domain_scope = "custom",
  custom_bg = as.character(MS2$Gene.Symbol...1),
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE
)

g.MS.UBE2HKOAZ.AZ.DOWNReg<-g.MS.UBE2HKOAZ.AZ.DOWNReg$result
g.MS.UBE2HKOAZ.AZ.DOWNReg<-g.MS.UBE2HKOAZ.AZ.DOWNReg[order(g.MS.UBE2HKOAZ.AZ.DOWNReg$p_value),]

##Upregulated protein in UBE2HKO cells 
MS.UBEAZ.AZ.D<-g.MS.UBE2HKOAZ.AZ.DOWNReg[order(g.MS.UBE2HKOAZ.AZ.DOWNReg$p_value),]
MS.UBEAZ.AZ.D<- subset(MS.UBEAZ.AZ.D, !source %in% c("HPA", "TF"))
MS.UBEAZ.AZ.D$term_name
MS.UBEAZ.AZ.D<-MS.UBEAZ.AZ.D[c(1,2,16, 85,144, 90,148),]
MS.UBEAZ.AZ.D$term_name<- factor(MS.UBEAZ.AZ.D$term_name, levels= MS.UBEAZ.AZ.D$term_name)
MS.UBEAZ.AZ.D$Termtype<-c( 
  "Mitochondrial Metabolism", "Mitochondrial Metabolism", "Mitochondrial Metabolism",
  "Mitochondrial Ribosome", "Mitochondrial Ribosome", "Ribosome","Ribosome"
)

ggplot(MS.UBEAZ.AZ.D, 
       aes(x=term_name, y=-log2(p_value), fill=Termtype))+
  geom_bar(stat="identity")+
  ylab("-log2 (P-Value)") +
  xlab("")+
  ggtitle ("Proteins enriched: Downregulated upon UBE2H-KO in AZ3146 background")+
  theme_classic()+
  coord_flip()
## 
# plot.gProfiler.MassSpec.UBE2HKOAZ.AZ.Downreg.pdf
# 10x5
# Figure 6J



# UBE2H-KO in AZ3146 background UPREGULATED PROTEINS 
# SIGNIFICANT: cell cycle, nucleoplasm, spliceosome, 
# G-Profiler: 
g.MS.UBE2HKOAZ.AZ.UPReg<- gost(
  subset(MS2, MS2$Abundance.Ratio..log2....UBE2HKO_AZ.....AZ. > MS2_UBEAZ.AZ_Up)$Gene,
  organism = "hsapiens",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = c("g_SCS"),
  domain_scope = "custom",
  custom_bg = as.character(MS2$Gene.Symbol...1),
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE
)

g.MS.UBE2HKOAZ.AZ.UPReg<-g.MS.UBE2HKOAZ.AZ.UPReg$result
g.MS.UBE2HKOAZ.AZ.UPReg<-g.MS.UBE2HKOAZ.AZ.UPReg[order(g.MS.UBE2HKOAZ.AZ.UPReg$p_value),]

##Upregulated protein in UBE2HKO cells 
MS.UBEAZ.AZ.U<-g.MS.UBE2HKOAZ.AZ.UPReg[order(g.MS.UBE2HKOAZ.AZ.UPReg$p_value),]
MS.UBEAZ.AZ.U<- subset(MS.UBEAZ.AZ.U, !source %in% c("HPA", "TF"))
MS.UBEAZ.AZ.U$term_name
MS.UBEAZ.AZ.U<-MS.UBEAZ.AZ.U[c(2,5,7, 35,39),]
MS.UBEAZ.AZ.U$term_name<- factor(MS.UBEAZ.AZ.U$term_name, levels= MS.UBEAZ.AZ.U$term_name)
MS.UBEAZ.AZ.U$Termtype<-c( 
  "Cell Cycle", "Cell Cycle", "Cell Cycle", 
  "spliceosome", "spliceosome"
)

ggplot(MS.UBEAZ.AZ.U, 
       aes(x=term_name, y=-log2(p_value), fill=Termtype))+
  geom_bar(stat="identity")+
  ylab("-log2 (P-Value)") +
  xlab("")+
  ggtitle ("Proteins enriched: Upregulated upon UBE2H-KO in AZ3146 background")+
  theme_classic()+
  coord_flip()
## 
# plot.gProfiler.MassSpec.UBE2HKOAZ.AZ.Upreg.pdf
# 10x5
# Figure S6L




# AZ3145 in UBE2H-KO cells UPREGULATED PROTEINS 
# SIGNIFICANT: mitochondrial metabolism
# G-Profiler: 
g.MS.UBE2HKOAZ.UBE2H.UPReg<- gost(
  subset(MS2, MS2$Abundance.Ratio..log2....UBE2HKO_AZ.....UBE2HKO. > MS2_UBEAZ.UBE_Up)$Gene,
  organism = "hsapiens",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = c("g_SCS"),
  domain_scope = "custom",
  custom_bg = as.character(MS2$Gene.Symbol...1),
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE
)

g.MS.UBE2HKOAZ.UBE2H.UPReg<-g.MS.UBE2HKOAZ.UBE2H.UPReg$result
g.MS.UBE2HKOAZ.UBE2H.UPReg<-g.MS.UBE2HKOAZ.UBE2H.UPReg[order(g.MS.UBE2HKOAZ.UBE2H.UPReg$p_value),]

## Upregulated protein upon AZ3146 in UBE2HKO cells 
MS.UBEAZ.UBE.U<-g.MS.UBE2HKOAZ.UBE2H.UPReg[order(g.MS.UBE2HKOAZ.UBE2H.UPReg$p_value),]
MS.UBEAZ.UBE.U<- subset(MS.UBEAZ.UBE.U, !source %in% c("HPA", "TF"))
MS.UBEAZ.UBE.U$term_name
MS.UBEAZ.UBE.U<-MS.UBEAZ.UBE.U[c(2,4,5,10,22),]
MS.UBEAZ.UBE.U$term_name<- factor(MS.UBEAZ.UBE.U$term_name, levels= MS.UBEAZ.UBE.U$term_name)
MS.UBEAZ.UBE.U$Termtype<-c( 
  "Mitochondria", "Mitochondria", "Mitochondria", "Mitochondria",  "Mitochondria"
)

ggplot(MS.UBEAZ.UBE.U, 
       aes(x=term_name, y=-log2(p_value), fill=Termtype))+
  geom_bar(stat="identity")+
  ylab("-log2 (P-Value)") +
  xlab("")+
  ggtitle ("Proteins enriched: Upregulated upon AZ3146 in UBE2H-KO background")+
  theme_classic()+
  coord_flip()
## 
# plot.gProfiler.MassSpec.UBE2HKOAZ.UBE.Upreg.pdf
# 10x5




# AZ3145 in UBE2H-KO cells DOWNREGULATED PROTEINS 
# SIGNIFICANT: Ribosomes, rRNA pocessing
# G-Profiler: 
g.MS.UBE2HKOAZ.UBE2H.DownReg<- gost(
  subset(MS2, MS2$Abundance.Ratio..log2....UBE2HKO_AZ.....UBE2HKO. < MS2_UBEAZ.UBE_Down)$Gene,
  organism = "hsapiens",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = c("g_SCS"),
  domain_scope = "custom",
  custom_bg = as.character(MS2$Gene.Symbol...1),
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE
)

g.MS.UBE2HKOAZ.UBE2H.DownReg<-g.MS.UBE2HKOAZ.UBE2H.DownReg$result
g.MS.UBE2HKOAZ.UBE2H.DownReg<-g.MS.UBE2HKOAZ.UBE2H.DownReg[order(g.MS.UBE2HKOAZ.UBE2H.DownReg$p_value),]

##Downregulated protein upon AZ3146 in UBE2HKO cells 
MS.UBEAZ.UBE.D<-g.MS.UBE2HKOAZ.UBE2H.DownReg[order(g.MS.UBE2HKOAZ.UBE2H.DownReg$p_value),]
MS.UBEAZ.UBE.D<- subset(MS.UBEAZ.UBE.D, !source %in% c("HPA", "TF"))
MS.UBEAZ.UBE.D$term_name
MS.UBEAZ.UBE.D<-MS.UBEAZ.UBE.D[c(1,2,3,5, 12,22,24,27),]
MS.UBEAZ.UBE.D$term_name<- factor(MS.UBEAZ.UBE.D$term_name, levels= MS.UBEAZ.UBE.D$term_name)
MS.UBEAZ.UBE.D$Termtype<-c( 
  "Ribosomes","Ribosomes","Ribosomes","Ribosomes", 
  "rRNA processing","rRNA processing","rRNA processing","rRNA processing"
)

ggplot(MS.UBEAZ.UBE.D, 
       aes(x=term_name, y=-log2(p_value), fill=Termtype))+
  geom_bar(stat="identity")+
  ylab("-log2 (P-Value)") +
  xlab("")+
  ggtitle ("Proteins enriched: Downregulated upon AZ3146 in UBE2H-KO background")+
  theme_classic()+
  coord_flip()
## 
# plot.gProfiler.MassSpec.UBE2HKOAZ.UBE.Downreg.pdf
# 10x5


setwd(ResultsFile)
#write.xlsx(g.MS.UBE2HKO.DownReg2[,1:13], "gProfiler.MassSPec.UBE2HKO.DownRegulated.xlsx")
#write.xlsx(g.MS.UBE2HKO.UPReg[,1:13], "gProfiler.MassSPec.UBE2HKO.UpRegulated.xlsx")

#write.xlsx(g.MS.AZ.DOWNReg[,1:13], "gProfiler.MassSPec.AZ.DownRegulated.xlsx")
#write.xlsx(g.MS.AZ.UPReg[,1:13], "gProfiler.MassSPec.AZ.UPRegulated.xlsx")

#write.xlsx(g.MS.UBE2HKOAZ.AZ.DOWNReg[,1:13], "gProfiler.MassSPec.UBE2HKOAZ.AZ.DownRegulated.xlsx")
#write.xlsx(g.MS.UBE2HKOAZ.AZ.UPReg[,1:13], "gProfiler.MassSPec.UBE2HKOAZ_AZ.UPRegulated.xlsx")

#write.xlsx(g.MS.UBE2HKOAZ.UBE2H.DownReg[,1:13], "gProfiler.MassSPec.UBE2HKOAZ_UBE2H.DownRegulated.xlsx")
#write.xlsx(g.MS.UBE2HKOAZ.UBE2H.UPReg[,1:13], "gProfiler.MassSPec.UBE2HKOAZ_UBE2H.UPRegulated.xlsx")

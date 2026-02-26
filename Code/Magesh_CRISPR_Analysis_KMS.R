#### CRISPR Screen data analysis from Magesh et al. (Emma Watson lab) ####
##   Mageshetal_CRISPR_Analysis_CC.R

## data from Magesh R.Y. et al. 2025
##  From Emma Watson lab
## University of Massachusetts Chem medical school 
# Human Mammary epithelial cells. 
# Paired aneuploidy screens



### INTRO AND LIBRARIES ####
## November 2024
### Watson lab screen data 

# Data from: 
# Magesh RY et al. (2024) An aneuploidy epistasis map reveals metabolic vulnerabilities 
# associated with supernumerary chromosomes in cancer. BioRXIV preprint
# https://www.biorxiv.org/content/10.1101/2024.09.30.615609v1


#Publication: 
#Rayna Y Magesh,  Arshia N Kaur ,  Faith N Keller ,  Abdulrazak Frederick , Tenzin Tseyang ,  John A Haley , Alejandra M Rivera-Nieves, Anthony C Liang ,  David A Guertin , Jessica B Spinelli , Stephen J Elledge , Emma V Watson
# Aneuploidy generates enhanced nucleotide dependency and sensitivity to metabolic perturbation
# Genes Dev (2025)  Jun 2;39(11-12):770-786. doi: 10.1101/gad.352512.124.
#PMID: 40324880 PMCID: PMC12128873 DOI: 10.1101/gad.352512.124

# Human Mammary Epithelial cells

## Bioarchive paper goes into metabolism. increased metabolism in aneuploid and 4N aneuploid cells, especially. 


# Aneuploid cell lines: 
# 1: Trisomy 15 and 20 and Monosomy 10p 
# 2: Trisomy 9, 12 and 20 and Tet8q
# 3: Trisomy 10, 12, 13, 20 and X
# 4: Trisomy 5 and 17 and Tetrasomy 20

# All have gain 20. 
# Gain of 12 and 17 occur twice. 
# ILK 20q gain
# 9p loss



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
library("Hmisc")
library(corrplot)
library("PerformanceAnalytics")



### GET DATA #####
# Human Mammary epithelial cells


setwd(Dependency)

Dropouts.Watson <- read.csv("Magesh_etal_HMEC_CRISPR_aneupPD6.vs.dipPD6.csv")
# 19,005 genes 
# Whole genome library

colnames(Dropouts.Watson)



#### Rename all Excel "Date" gene names
# Make function to rename "Sept-1" genes and "SEPT1" to their updated names "SEPTIN1" etc. for all 
#   genes who's names are misinterpretted as dates in Excel. 
RenameExcelGenes2<- function(GeneMLEFile){
  GeneMLEFile$gene<- as.character(GeneMLEFile$gene)
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '1-Dec'] <- "DELEC1"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == 'DEC1'] <- "DELEC1"
  
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '15-Sep'] <- "SELENOF"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '1-Sep'] <- "SEPTIN1"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '2-Sep'] <- "SEPTIN2"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '3-Sep'] <- "SEPTIN3"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '4-Sep'] <- "SEPTIN4"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '5-Sep'] <- "SEPTIN5"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '6-Sep'] <- "SEPTIN6"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '7-Sep'] <- "SEPTIN7"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '8-Sep'] <- "SEPTIN8"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '9-Sep'] <- "SEPTIN9"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '10-Sep'] <- "SEPTIN10"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '11-Sep'] <- "SEPTIN11"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '12-Sep'] <- "SEPTIN12"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '14-Sep'] <- "SEPTIN14"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == 'SEP15'] <- "SELENOF"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == 'SEPT1'] <- "SEPTIN1"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == 'SEPT2'] <- "SEPTIN2"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == 'SEPT3'] <- "SEPTIN3"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == 'SEPT4'] <- "SEPTIN4"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == 'SEPT5'] <- "SEPTIN5"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == 'SEPT6'] <- "SEPTIN6"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == 'SEPT7'] <- "SEPTIN7"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == 'SEPT8'] <- "SEPTIN8"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == 'SEPT9'] <- "SEPTIN9"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == 'SEPT10'] <- "SEPTIN10"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == 'SEPT11'] <- "SEPTIN11"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == 'SEPT12'] <- "SEPTIN12"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == 'SEPT14'] <- "SEPTIN14"
  
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '1-Mar'] <- "MARCHF1"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '2-Mar'] <- "MARCHF2"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '3-Mar'] <- "MARCHF3"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '4-Mar'] <- "MARCHF4"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '5-Mar'] <- "MARCHF5"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '6-Mar'] <- "MARCHF6"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '7-Mar'] <- "MARCHF7"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '8-Mar'] <- "MARCHF8"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '9-Mar'] <- "MARCHF9"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '10-Mar'] <- "MARCHF10"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '11-Mar'] <- "MARCHF11"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == 'MARCH1'] <- "MARCHF1"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == 'MARCH2'] <- "MARCHF2"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == 'MARCH3'] <- "MARCHF3"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == 'MARCH4'] <- "MARCHF4"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == 'MARCH5'] <- "MARCHF5"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == 'MARCH6'] <- "MARCHF6"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == 'MARCH7'] <- "MARCHF7"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == 'MARCH8'] <- "MARCHF8"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == 'MARCH9'] <- "MARCHF9"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == 'MARCH10'] <- "MARCHF10"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == 'MARCH11'] <- "MARCHF11"
  
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '43525'] <- 'MARCHF1'
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '43534'] <- 'MARCHF10'
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '43535'] <- 'MARCHF11'
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '43526'] <- 'MARCHF2'
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '43527'] <- 'MARCHF3'
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '43528'] <- 'MARCHF4'
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '43529'] <- 'MARCHF5'
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '43530'] <- 'MARCHF6'
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '43531'] <- 'MARCHF7'
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '43532'] <- 'MARCHF8'
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '43533'] <- 'MARCHF9'
  
  return(GeneMLEFile)
} 

Dropouts.Watson<- RenameExcelGenes2(Dropouts.Watson) 




setwd(Dependency)
Proteosome_Genes<- read.csv("group-690.csv", header=TRUE)
Ribosomal_Genes<- read.csv("group-1054.csv", skip = 1 , header=TRUE) #HGNC Ribosomal Proteins
Ribosomal_Genes_noMT<- subset(Ribosomal_Genes, Group %in% c("L ribosomal proteins", "S ribosomal proteins"))
Olfactory_Genes<- read.csv("group-141.csv", skip = 1 , header=TRUE)

RNA_Processing<- read_tsv("GO_term_summary_20250720_134600.txt")
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



### Plot Magesh difference p.value #####
length(Dropouts.Watson$gene) # 19,005 genes

# they already calculated log2 fold change and p-value for all 4 aneuploids, as well as for each aneuploid clone



## Highlight hits

HighlightTheseGenes<- c("UBE2H")
HighlightTheseGenes<- c("HSF1")
HighlightTheseGenes<- c("FOSL1")
HighlightTheseGenes<- c("FOSL1", "JUNB", "JUND", "JUN", "FOSL2", "FOS", "FOSB")

HighlightTheseGenes<- subset(Dropouts.Watson, log2FC_meta< -1 & Pvalue_meta< 0.05)$gene

# plot significant genes: 
ggplot(Dropouts.Watson, aes(y=-log(Pvalue_meta,2), x=log2FC_meta, 
                                      color=gene %in% HighlightTheseGenes))+
  geom_point()+
  geom_point(data = subset(Dropouts.Watson, gene %in% HighlightTheseGenes))+
  geom_density_2d(color = "white")+
  scale_color_manual(values=c("Black", "red"))+ 
  theme_classic()+
  theme(legend.position="none")+
  xlab("Beta delta (aneuploid - euploid)")+
  ylab("-log2(p-value)")
# 5x4
# plot.Magesh_HMEF_Diff.pvalue.pdf


# plot pathways: 
ggplot(Dropouts.Watson, aes(y=-log(Pvalue_meta,2), x=log2FC_meta))+
  geom_point()+
  geom_point(data = subset(Dropouts.Watson, gene %in% MitoElectronChainnAssembly), color="#F8766D")+ #F8766D
  geom_point(data = subset(Dropouts.Watson, gene %in% MitoTranslationTranscription), color="gold2")+   #7CAE00
  geom_point(data = subset(Dropouts.Watson, gene %in% RNA_Processing_rRNA ), color="deepskyblue4") + #00BFC4
  geom_point(data = subset(Dropouts.Watson, gene %in% Ribosomal_Genes_noMT$Approved.symbol ), color="#C77CFF") + #C77CFF
  geom_point(data = subset(Dropouts.Watson, gene %in% RNA_Processing_Spliceosome ), color="deepskyblue1") + #00BFC4
  geom_point(data = subset(Dropouts.Watson, gene %in% Proteosome_Genes$Approved.symbol ), color="chartreuse1") + #7CAE00
  geom_density_2d(color = "white")+
  theme(legend.position="none")+
  theme_classic()+
  xlab("log2 FC")+
  ylab("-log2(p-value)")
# 5x4
# Figure 4D
# plot.Magesh.Watson_HMEF.pairedP.Diff_AllPathway.pdf


# plot pathways as density: 
ggplot(Dropouts.Watson, aes( x=log2FC_meta))+
  geom_density(size=2, data = subset(Dropouts.Watson, gene %in% MitoElectronChainnAssembly), color="#F8766D")+ #F8766D
  geom_density(size=2, data = subset(Dropouts.Watson, gene %in% MitoTranslationTranscription), color="gold2")+   #7CAE00
  geom_density(size=2, data = subset(Dropouts.Watson, gene %in% RNA_Processing_rRNA ), color="deepskyblue4") + #00BFC4
  geom_density(size=2, data = subset(Dropouts.Watson, gene %in% Ribosomal_Genes_noMT$Approved.symbol ), color="#C77CFF") + #C77CFF
  geom_density(size=2, data = subset(Dropouts.Watson, gene %in% Proteosome_Genes$Approved.symbol ), color="chartreuse1") + #7CAE00
  geom_density(size=2, data = subset(Dropouts.Watson, gene %in% RNA_Processing_Spliceosome ), color="deepskyblue1") + #00BFC4
  geom_density(size=2)+
  theme(legend.position="none")+
  theme_classic()+
  xlab("log2 FC")+
  ylab("-log2(p-value)")
# 5x4
# plot.Magesh.Watson_HMEF.DENSITY_AllPathway.pdf


# t-test for pathways
t.test(subset(Dropouts.Watson, ! gene %in% Proteosome_Genes$Approved.symbol)$log2FC_meta, 
       subset(Dropouts.Watson, gene %in% Proteosome_Genes$Approved.symbol)$log2FC_meta ) #0.09832
t.test(subset(Dropouts.Watson, ! gene %in% RNA_Processing_Spliceosome)$log2FC_meta, 
       subset(Dropouts.Watson, gene %in% RNA_Processing_Spliceosome)$log2FC_meta ) # 0.0005969
t.test(subset(Dropouts.Watson, ! gene %in% MitoElectronChainnAssembly)$log2FC_meta, 
       subset(Dropouts.Watson, gene %in% MitoElectronChainnAssembly)$log2FC_meta ) #9.932e-15
t.test(subset(Dropouts.Watson, ! gene %in% MitoTranslationTranscription)$log2FC_meta, 
       subset(Dropouts.Watson, gene %in% MitoTranslationTranscription)$log2FC_meta ) #4.489e-12
t.test(subset(Dropouts.Watson, ! gene %in% Ribosomal_Genes_noMT$Approved.symbol)$log2FC_meta, 
       subset(Dropouts.Watson, gene %in% Ribosomal_Genes_noMT$Approved.symbol)$log2FC_meta ) #0.003039
t.test(subset(Dropouts.Watson, ! gene %in% RNA_Processing_rRNA)$log2FC_meta, 
       subset(Dropouts.Watson, gene %in% RNA_Processing_rRNA)$log2FC_meta ) #3.633e-06

# proteosome       0.098
# Spliceosome      5.9E-4
# MitoETC.         9.9E-15
# MitoTrans.       4.5E-12
# Ribosome.        0.003
# RNA processing   3.6E-6



### plot Magesh RANK ####
Dropouts.Watson$AneuploidRank<- rank(Dropouts.Watson$log2FC_meta)


ggplot(Dropouts.Watson, aes(x=AneuploidRank, y=log2FC_meta, color= AneuploidRank < 1901))+
  geom_point()+
  geom_point(data = subset(Dropouts.Watson, AneuploidRank < 1901))+
  scale_color_manual(values=c("Black", "Red"))+ 
  theme_classic()+
  theme(legend.position="none")+
  ylab("log2 FC:\n hMEC")+
  xlab("Magesh et al. (2024): Rank") 
# 5x3  
# Figure 4B
# plot.Magesh.Watson.Rank_Top10.pdf 

subset(Dropouts.Watson,  AneuploidRank < 1901)$gene

mean(subset(Dropouts.Watson, gene %in% MitoTranslationTranscription)$log2FC_meta ) # -0.2511305
mean(subset(Dropouts.Watson, gene %in% MitoElectronChainnAssembly)$log2FC_meta ) # -0.302472

# Write top ten
setwd(ResultsFile)
write.csv(subset(Dropouts.Watson, AneuploidRank < 1901)$gene, "Magesh.Watson.RANK.Top10.csv")
write.csv(subset(Dropouts.Watson), "Magesh.Watson.RANK.csv")


# Top Top hits: 
subset(Dropouts.Watson,log2FC_meta< -1.25 & Pvalue_meta< 0.01)$gene
# [1] "DHODH"         "YRDC"          "CAD"           "UMPS"          "SMEK1"         "MAT2A"         "KEAP1"         "DTYMK"         "MYBL2"        
#[10] "PGD"           "RPE"           "PITPNB"        "CKS1B"         "GMPPB"         "TKT"           "UBE2H"         "OGT"           "UQCRFS1"      
#[19] "POGZ"          "PGS1"          "MRE11A"        "POLR2C"        "CCS"           "ZPR1"          "COPS5"         "POP5"          "LSS"          
#[28] "LRR1"          "NAPA"          "RTFDC1"        "RUVBL2"        "RAD1"          "PI4KA"         "CYP51A1"       "SEPHS1"        "TAF15"        
#[37] "NPIPB3"        "SOD1"          "HMGCS1"        "GAPDH"         "ALG1L"         "SNRNP200"      "RP11-812E19.9"

# UBE2H in top hits. 

### G Profiler: Magesh #####

Dropouts.Watson
### terms enriched in toxic to aneuploid cells
g.Antoxic.Watson<- gost(
  subset(Dropouts.Watson, log2FC_meta< -0.2)$gene,
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
  custom_bg = as.character(Dropouts.Watson$gene),
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE
)
g.Antoxic.Watson_result<-g.Antoxic.Watson$result


# terms enriched in toxic to aneuploid cells
P.An.Wa<-g.Antoxic.Watson_result[order(g.Antoxic.Watson_result$p_value),]
P.An.Wa<- subset(P.An.Wa, ! source %in% ("HPA"))
P.An.Wa$term_name[1:300]
P.An.Wa<-P.An.Wa[c(3,5,7,93,97,136,  8,10, 27,37,40,  36,42,51, 151),]
P.An.Wa$term_name<- factor(P.An.Wa$term_name, levels= P.An.Wa$term_name)
P.An.Wa$Termtype<-c("Metabolism", "Metabolism", "Metabolism", "Metabolism", "Metabolism", "Metabolism", 
                    
                    "Protein Complex", "Protein Complex", 
                    "Purine Ribonucleotide binding", "Purine Ribonucleotide binding", "Purine Ribonucleotide binding",
                    "RNA processing", "RNA processing","RNA processing",
                    "Ribosomal") 

ggplot(P.An.Wa, 
       aes(x=term_name, y=-log2(p_value), fill=Termtype))+
  geom_bar(stat="identity")+
  ylab("-log2 (P-Value)") +
  xlab("Key RNA terms")+
  ggtitle ("Top biological terms enriched: toxic to Aneuploid, Watson")+
  theme_classic()+
  coord_flip()
# plot.GProfiler.Watson.ANvEU.pdf
# 8x5


setwd(ResultsFile)
write.csv(g.Antoxic.Watson_result[1:13], "gprofiler.Aneuploidtoxic.Magesh.csv")


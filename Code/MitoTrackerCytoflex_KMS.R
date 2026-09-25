### Analyze Mitotracker Cytoflex data ######
# date: 2026-May 26
#Klaske Schukken

## The purpose of this code is to take Cytoflex data from cells treated with Hoechst and MitoTracker and measure the ammount of Mitotracker, cell size, Hoechst in different cells. 
## Main question: Does UBE2H-KO alter the ammount of mitochondria? Does UBE2H-Ko alter the ammount of mitochondria when treated with AZ3146? 
## How many cells/well? What was the growth rate? Assuming all cells had 100k cells/well plates. 

setwd("/Users/christinescaduto/Documents/Klaske/Klaske/Digital lab notebook/Hit Gene follow up/UBE2H/MitoTracker/Mitotracker cytoflex")



### INTRO AND LIBRARIES ####
# Folders
## !! Update these locations to pathway where data was downloaded: 
# Folder with CRISPR screening datasets and data generated for this paper: 
SchukkenData<- "/Volumes/Schukken_SSD/CRISPR SCREEN PAPER/R Code Data Files/Data"
# Folder with dependency datasets: 
Dependency<- "/Volumes/Schukken_SSD/CRISPR SCREEN PAPER/R Code Data Files/Dependency files"
# Folder with results: 
ResultsFile<- "/Volumes/Schukken_SSD/CRISPR SCREEN PAPER/R Code Data Files/Results"


library(ggplot2)
library(stringr)


### Get Data ####
setwd(SchukkenData)

# HCT116 WT, KO clones and KO + bulk overexpression cells
HCT116.518<- read.csv("260518_HCT116_Mitotracker_data_Mean.csv", skip = 2) # start on row 3
HCT116.525<- read.csv("260525_HCT116_Mitotracker.csv", skip = 2) #May 25 2026 mitotracker cytoflex data 
HCT116.531<- read.csv("260531_HCT116_MitoTracker3.csv", skip =2)
HCT116.609<- read.csv("2026.06.09_HCT116_Mito.csv", skip =2)

# HCT116 WT, KO clones and KO + single cell overexpression clones: (Cleaner data)
HCT116.817<- read.csv("2026.08.17_Cytoflex_KOOEc.csv", skip =2) # Aug 17, HCT116 KO + OE clones
HCT116.825<- read.csv("2026.08.25_Cytoflex_KOOEc.csv", skip =2) # Aug 17, HCT116 KO + OE clones

# Cal51
Cal51.603<- read.csv("260603_Cal51_MitoTracker.csv", skip =2)

### Data format #####

# Data type:
# Cytoflex data for cells. looking at 30ul of 200ul total 
# FSC: Forward scatter. How big is cell? 
# SSC: Side Scatter. how dense is cell?
# Violet610: Hoechst stain. DNA 
# PE: Red fluorescence. MitoTracker. 

# Note: Hoechst staining is fairly consistent WITHIN experiemnts, but across experiments there is a lot of variability.
# Potentially due to slight differences in length of time with hoechst staining between experiments. 
# Use cell size as a control instead. 

#[1] "Tube.Name."                      "Sample.ID."                      "All.Events"                     
#[4] "Cell.size.Events"                "Hoechst.Events"                  "medium.Mito.Events"             
#[7] "low.Mito.Events"                 "High.Mito.Events"                "All.Events....Total"            
#[10] "Cell.size...Total"               "Hoechst...Total"                 "medium.Mito...Total"            
#[13] "low.Mito...Total"                "High.Mito...Total"               "All.Events....Parent"           
#[16] "Cell.size...Parent"              "Hoechst...Parent"                "medium.Mito...Parent"           
#[19] "low.Mito...Parent"               "High.Mito...Parent"              "All.Events..Median.PE.A"        
#[22] "Cell.size.Median.PE.A"           "Hoechst.Median.PE.A"             "medium.Mito.Median.PE.A"        
#[25] "low.Mito.Median.PE.A"            "High.Mito.Median.PE.A"           "All.Events..Mean.PE.A"          
#[28] "Cell.size.Mean.PE.A"             "Hoechst.Mean.PE.A"               "medium.Mito.Mean.PE.A"          
#[31] "low.Mito.Mean.PE.A"              "High.Mito.Mean.PE.A"             "All.Events..rSD.PE.A"           
#[34] "Cell.size.rSD.PE.A"              "Hoechst.rSD.PE.A"                "medium.Mito.rSD.PE.A"           
#[37] "low.Mito.rSD.PE.A"               "High.Mito.rSD.PE.A"              "All.Events..CV.PE.A"            
#[40] "Cell.size.CV.PE.A"               "Hoechst.CV.PE.A"                 "medium.Mito.CV.PE.A"            
#[43] "low.Mito.CV.PE.A"                "High.Mito.CV.PE.A"               "All.Events..Mean.FSC.A"         
#[46] "Cell.size.Mean.FSC.A"            "Hoechst.Mean.FSC.A"              "medium.Mito.Mean.FSC.A"         
#[49] "low.Mito.Mean.FSC.A"             "High.Mito.Mean.FSC.A"            "All.Events..Median.FSC.A"       
#[52] "Cell.size.Median.FSC.A"          "Hoechst.Median.FSC.A"            "medium.Mito.Median.FSC.A"       
#[55] "low.Mito.Median.FSC.A"           "High.Mito.Median.FSC.A"          "All.Events..rSD.FSC.A"          
#[58] "Cell.size.rSD.FSC.A"             "Hoechst.rSD.FSC.A"               "medium.Mito.rSD.FSC.A"          
#[61] "low.Mito.rSD.FSC.A"              "High.Mito.rSD.FSC.A"             "All.Events..CV.FSC.A"           
#[64] "Cell.size.CV.FSC.A"              "Hoechst.CV.FSC.A"                "medium.Mito.CV.FSC.A"           
#[67] "low.Mito.CV.FSC.A"               "High.Mito.CV.FSC.A"              "All.Events..SD.FSC.A"           
#[70] "Cell.size.SD.FSC.A"              "Hoechst.SD.FSC.A"                "medium.Mito.SD.FSC.A"           
#[73] "low.Mito.SD.FSC.A"               "High.Mito.SD.FSC.A"              "All.Events..rCV.FSC.A"          
#[76] "Cell.size.rCV.FSC.A"             "Hoechst.rCV.FSC.A"               "medium.Mito.rCV.FSC.A"          
#[79] "low.Mito.rCV.FSC.A"              "High.Mito.rCV.FSC.A"             "All.Events..GeoMean.FSC.A"      
#[82] "Cell.size.GeoMean.FSC.A"         "Hoechst.GeoMean.FSC.A"           "medium.Mito.GeoMean.FSC.A"      
#[85] "low.Mito.GeoMean.FSC.A"          "High.Mito.GeoMean.FSC.A"         "All.Events..Mean.SSC.A"         
#[88] "Cell.size.Mean.SSC.A"            "Hoechst.Mean.SSC.A"              "medium.Mito.Mean.SSC.A"         
#[91] "low.Mito.Mean.SSC.A"             "High.Mito.Mean.SSC.A"            "All.Events..GeoMean.SSC.A"      
#[94] "Cell.size.GeoMean.SSC.A"         "Hoechst.GeoMean.SSC.A"           "medium.Mito.GeoMean.SSC.A"      
#[97] "low.Mito.GeoMean.SSC.A"          "High.Mito.GeoMean.SSC.A"         "All.Events..Median.SSC.A"       
#[100] "Cell.size.Median.SSC.A"          "Hoechst.Median.SSC.A"            "medium.Mito.Median.SSC.A"       
#[103] "low.Mito.Median.SSC.A"           "High.Mito.Median.SSC.A"          "All.Events..rCV.SSC.A"          
#[106] "Cell.size.rCV.SSC.A"             "Hoechst.rCV.SSC.A"               "medium.Mito.rCV.SSC.A"          
#[109] "low.Mito.rCV.SSC.A"              "High.Mito.rCV.SSC.A"             "All.Events..rSD.SSC.A"          
#[112] "Cell.size.rSD.SSC.A"             "Hoechst.rSD.SSC.A"               "medium.Mito.rSD.SSC.A"          
#[115] "low.Mito.rSD.SSC.A"              "High.Mito.rSD.SSC.A"             "All.Events..CV.SSC.A"           
#[118] "Cell.size.CV.SSC.A"              "Hoechst.CV.SSC.A"                "medium.Mito.CV.SSC.A"           
#[121] "low.Mito.CV.SSC.A"               "High.Mito.CV.SSC.A"              "All.Events..SD.SSC.A"           
#[124] "Cell.size.SD.SSC.A"              "Hoechst.SD.SSC.A"                "medium.Mito.SD.SSC.A"           
#[127] "low.Mito.SD.SSC.A"               "High.Mito.SD.SSC.A"              "All.Events..Mean.Violet610.A"   
#[130] "Cell.size.Mean.Violet610.A"      "Hoechst.Mean.Violet610.A"        "medium.Mito.Mean.Violet610.A"   
#[133] "low.Mito.Mean.Violet610.A"       "High.Mito.Mean.Violet610.A"      "All.Events..GeoMean.Violet610.A"
#[136] "Cell.size.GeoMean.Violet610.A"   "Hoechst.GeoMean.Violet610.A"     "medium.Mito.GeoMean.Violet610.A"
#[139] "low.Mito.GeoMean.Violet610.A"    "High.Mito.GeoMean.Violet610.A"   "All.Events..Median.Violet610.A" 
#[142] "Cell.size.Median.Violet610.A"    "Hoechst.Median.Violet610.A"      "medium.Mito.Median.Violet610.A" 
#[145] "low.Mito.Median.Violet610.A"     "High.Mito.Median.Violet610.A"    "All.Events..rCV.Violet610.A"    
#[148] "Cell.size.rCV.Violet610.A"       "Hoechst.rCV.Violet610.A"         "medium.Mito.rCV.Violet610.A"    
#[151] "low.Mito.rCV.Violet610.A"        "High.Mito.rCV.Violet610.A"       "All.Events..rSD.Violet610.A"    
#[154] "Cell.size.rSD.Violet610.A"       "Hoechst.rSD.Violet610.A"         "medium.Mito.rSD.Violet610.A"    
#[157] "low.Mito.rSD.Violet610.A"        "High.Mito.rSD.Violet610.A"       "All.Events..CV.Violet610.A"     
#[160] "Cell.size.CV.Violet610.A"        "Hoechst.CV.Violet610.A"          "medium.Mito.CV.Violet610.A"     
#[163] "low.Mito.CV.Violet610.A"         "High.Mito.CV.Violet610.A"        "All.Events..SD.Violet610.A"     
#[166] "Cell.size.SD.Violet610.A"        "Hoechst.SD.Violet610.A"          "medium.Mito.SD.Violet610.A"     
#[169] "low.Mito.SD.Violet610.A"         "High.Mito.SD.Violet610.A"        "All.Events..GeoMean.PE.A"       
#[172] "Cell.size.GeoMean.PE.A"          "Hoechst.GeoMean.PE.A"            "medium.Mito.GeoMean.PE.A"       
#[175] "low.Mito.GeoMean.PE.A"           "High.Mito.GeoMean.PE.A"          "All.Events..rCV.PE.A"           
#[178] "Cell.size.rCV.PE.A"              "Hoechst.rCV.PE.A"                "medium.Mito.rCV.PE.A"           
#[181] "low.Mito.rCV.PE.A"               "High.Mito.rCV.PE.A"              "All.Events..SD.PE.A"            
#[184] "Cell.size.SD.PE.A"               "Hoechst.SD.PE.A"                 "medium.Mito.SD.PE.A"            
#[187] "low.Mito.SD.PE.A"                "High.Mito.SD.PE.A"               


### HCT116 MERGE all 4: KO and bulk OE cells #### 
# HCT116 cells treated with 1uM AZ3146 for 96 hours
HCT116.518<- read.csv("260518_HCT116_Mitotracker_data_Mean.csv", skip = 2) # start on row 3
HCT116.525<- read.csv("HCT116_Mitotracker_26.05.25.csv", skip = 2) #May 25 2026 mitotracker cytoflex data 
HCT116.531<- read.csv("260531_HCT116_MitoTracker3.csv", skip =2)
HCT116.609<- read.csv("2026.06.09_HCT116_Mito.csv", skip =2)

## Assign Sample ID names: 
HCT116.518$Sample.ID.<- str_sub(HCT116.518$Tube.Name., start=4, end = -4)
HCT116.525$Sample.ID.<- str_sub(HCT116.525$Tube.Name., start=4, end = -4)
HCT116.531$Sample.ID.<- str_sub(HCT116.531$Tube.Name., start=4, end = -4)
HCT116.609$Sample.ID.<- str_sub(HCT116.609$Tube.Name., start=4, end = -4)

# make relative mitoTracker: Mitotracker (PE) fluorescence . mean of WT Mitotracker
HCT116.518$RelativeMito<- HCT116.518$Hoechst.Mean.PE.A/ mean(subset(HCT116.518, Sample.ID. == "WT")$Hoechst.Mean.PE.A)
HCT116.525$RelativeMito<- HCT116.525$Hoechst.Mean.PE.A/ mean(subset(HCT116.525, Sample.ID. == "WT")$Hoechst.Mean.PE.A)
HCT116.531$RelativeMito<- HCT116.531$Hoechst.Mean.PE.A/ mean(subset(HCT116.531, Sample.ID. == "WT")$Hoechst.Mean.PE.A)
HCT116.609$RelativeMito<- HCT116.609$Hoechst.Mean.PE.A/ mean(subset(HCT116.609, Sample.ID. == "WT")$Hoechst.Mean.PE.A)


HCT116.Merge<- rbind(HCT116.525, HCT116.518[c(1:9,13:33,37:48),]) # I am removing KO3 from test #1 
# because I do not trust that I added the same ammount of mitotracker to samples. make new media batch. less Mitotracker? 
HCT116.Merge<- rbind(HCT116.Merge, HCT116.531)
HCT116.Merge<- rbind(HCT116.Merge, HCT116.609)



# Assign Sample ID names: 
# Remove the last threee characters from "Tube.Name." (..._01, ..._02, ..._03, etc.) and 
# and remove the first three characters from "Tube.Name." (01_...) 
HCT116.Merge$Sample.ID.[HCT116.Merge$Sample.ID. == "H.OE.AZ"] <- "OE_AZ"
HCT116.Merge$Sample.ID.[HCT116.Merge$Sample.ID. == "KO2.AZ"] <- "KO2_AZ"
HCT116.Merge$Sample.ID.[HCT116.Merge$Sample.ID. == "KO3.OE"] <- "KO3_OE"
HCT116.Merge$Sample.ID.[HCT116.Merge$Sample.ID. == "KO3.OE.AZ"] <- "KO3_OE_AZ"
HCT116.Merge$Sample.ID.<- factor(HCT116.Merge$Sample.ID., levels= c("WT", "WT_AZ", "OE", "OE_AZ", 
                                                                    "KO1", "KO1_AZ",  "KO1_OE","KO1_OE_AZ",
                                                                    "KO2","KO2_AZ", "KO2_OE", "KO2_OE_AZ",
                                                                    "KO3", "KO3_AZ", "KO3_OE", "KO3_OE_AZ"
))
HCT116.Merge$Drug<- grepl( "AZ", HCT116.Merge$Sample.ID., fixed = TRUE)

# rename that one wierdly names collumn HCT116 + UBE2H-OE + AZ3146



# now merge KO clones as a single "type" 
# first name the May 25, then may 18 samples: 
HCT116.Merge$Type<- c("WT","WT","WT", 
                      "KO","KO","KO", "KO","KO","KO", "KO","KO","KO",
                      "OE","OE","OE",
                      "KO+OE","KO+OE","KO+OE", "KO+OE","KO+OE","KO+OE", "KO+OE","KO+OE","KO+OE", 
                      "WT+AZ","WT+AZ","WT+AZ", 
                      "KO+AZ","KO+AZ","KO+AZ","KO+AZ","KO+AZ","KO+AZ","KO+AZ","KO+AZ","KO+AZ",
                      "KO+OE+AZ","KO+OE+AZ","KO+OE+AZ","KO+OE+AZ","KO+OE+AZ","KO+OE+AZ","KO+OE+AZ","KO+OE+AZ","KO+OE+AZ", 
                      "OE+AZ","OE+AZ","OE+AZ", 
                      
                      # Now for May 18 samples: 
                      "WT","WT","WT", 
                      "KO","KO","KO", "KO","KO","KO",                      #    "KO","KO","KO",
                      "OE","OE","OE",
                      "KO+OE","KO+OE","KO+OE", "KO+OE","KO+OE","KO+OE", "KO+OE","KO+OE","KO+OE", 
                      "WT+AZ","WT+AZ","WT+AZ", 
                      "KO+AZ","KO+AZ","KO+AZ","KO+AZ","KO+AZ","KO+AZ",      #    "KO+AZ","KO+AZ","KO+AZ",
                      "OE+AZ","OE+AZ","OE+AZ",
                      "KO+OE+AZ","KO+OE+AZ","KO+OE+AZ","KO+OE+AZ","KO+OE+AZ","KO+OE+AZ","KO+OE+AZ","KO+OE+AZ","KO+OE+AZ",
                      
                      # Now for May 31 samples: 
                      "WT","WT","WT", 
                      "KO","KO","KO", "KO","KO","KO", "KO","KO","KO",
                      "OE","OE","OE",
                      "KO+OE","KO+OE","KO+OE", "KO+OE","KO+OE","KO+OE", "KO+OE","KO+OE","KO+OE", 
                      "WT+AZ","WT+AZ","WT+AZ", 
                      "KO+AZ","KO+AZ","KO+AZ","KO+AZ","KO+AZ","KO+AZ", "KO+AZ","KO+AZ","KO+AZ",
                      "OE+AZ","OE+AZ","OE+AZ",
                      "KO+OE+AZ","KO+OE+AZ","KO+OE+AZ","KO+OE+AZ","KO+OE+AZ","KO+OE+AZ","KO+OE+AZ","KO+OE+AZ","KO+OE+AZ", 
                      
                      # Now for June 9 samples: 
                      "WT","WT","WT", 
                      "KO","KO","KO", "KO","KO","KO", "KO","KO","KO",
                      "OE","OE","OE",
                      "KO+OE","KO+OE","KO+OE", "KO+OE","KO+OE","KO+OE", "KO+OE","KO+OE","KO+OE", 
                      "WT+AZ","WT+AZ","WT+AZ", 
                      "KO+AZ","KO+AZ","KO+AZ","KO+AZ","KO+AZ","KO+AZ", "KO+AZ","KO+AZ","KO+AZ",
                      "OE+AZ","OE+AZ","OE+AZ",
                      "KO+OE+AZ","KO+OE+AZ","KO+OE+AZ","KO+OE+AZ","KO+OE+AZ","KO+OE+AZ","KO+OE+AZ","KO+OE+AZ","KO+OE+AZ")


HCT116.Merge$Type<- factor(HCT116.Merge$Type, levels= c("WT", "WT+AZ","OE", "OE+AZ", 
                                                        "KO", "KO+AZ","KO+OE", "KO+OE+AZ"
))



# Cell growth (96 hours)
## Note that 30ul from a total of 200ul is measured. so multiuple # cells by 200/30
## Count number of cells: cell sized & Hoechst positive. 
MeanWTCount<- mean(subset(HCT116.Merge, Sample.ID. == "WT")$Hoechst.Events)
HCT116.Merge$RelativeCount<- HCT116.Merge$Hoechst.Events/MeanWTCount
HCT116.Merge$Doubling<- HCT116.Merge$Hoechst.Events*200/(30*100000*4) # doubling: final/initial 
# Final is number of Hoechst + cells *200/30 (only analyzed 30ul from 200ul total), started with 100k cells, grew for 4 days so divide by 4 to get doubling per day


ggplot(HCT116.Merge, aes(x=HCT116.Merge$Sample.ID., y=Doubling, color=Type ))+
  geom_point()+
  scale_color_manual(values=c("cyan3","cyan3", "cyan4", "cyan4", "purple1","purple1", "purple4",  "purple4"), guide="none")+
  stat_summary(fun.y=mean, geom="point", shape="-", size=10, color="red", fill="red") +
  ylab("Doubling per day")+
  xlab("")+  
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
# plot.HCT116.Merge.Doubling.pdf


####
## plot mean MitoTracker (PE) 
ggplot(HCT116.Merge, aes(x=Sample.ID., y=Hoechst.Mean.PE.A , fill = Type))+
  geom_boxplot()+
  geom_point()+
  scale_fill_manual(values=c("cyan3","cyan3", "cyan4", "cyan4", "purple1","purple1", "purple4",  "purple4"), guide="none")+
  geom_hline(yintercept = mean(subset(HCT116.Merge, Sample.ID. == "WT")$Hoechst.Mean.PE.A) ) + # average WT mitotracker expression
  #stat_summary(fun.y=mean, geom="point", shape="-", size=10, color="red", fill="red") +
  ylab("Mean Mitotracker")+
  xlab("")+  
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
# plot.HCT116.Merge.Mito.pdf


# Relative MitoTracker (Mitotracker divided by mean Mitotracker of WT per batch): 
ggplot(HCT116.Merge, aes(x=Sample.ID., y=RelativeMito , color = Type))+
  geom_point()+
  scale_color_manual(values=c("cyan3","cyan3", "cyan4", "cyan4", "purple1","purple1", "purple4",  "purple4"), guide="none")+
  geom_hline(yintercept = mean(subset(HCT116.Merge, Sample.ID. == "WT")$RelativeMito) ) + # average WT mitotracker expression
  stat_summary(fun.y=mean, geom="point", shape="-", size=10, color="red", fill="red") +
  ylab("Relative Mitotracker")+
  xlab("")+  
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
# plot.HCT116.Merge.RelativeMito.pdf


## plot mean  cell size
HCT116.Merge$RelativeSize<- HCT116.Merge$Hoechst.Mean.FSC.A/ mean(subset(HCT116.Merge, Sample.ID. == "WT")$Hoechst.Mean.FSC.A)
ggplot(HCT116.Merge, aes(x=Sample.ID., y=RelativeSize, fill= Type))+
  geom_boxplot()+
  geom_point()+
  #scale_color_manual(values=c("cyan3", "cyan4", "purple1", "purple4", "cyan3", "cyan4", "purple1", "purple4"))+
  scale_fill_manual(values=c("cyan3","cyan3", "cyan4", "cyan4", "purple1","purple1", "purple4",  "purple4"), guide="none")+
  geom_hline(yintercept = mean(subset(HCT116.Merge, Sample.ID. == "WT")$RelativeSize) ) + # average WT mitotracker expression
  #stat_summary(fun.y=mean, geom="point", shape="-", size=10, color="red", fill="red") +
  ylab("Relative Size")+
  xlab("")+  
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
# UBE2H-KO reduces cell size (KO+OE rescues cell size), OE increases cell size. 
# plot.HCT116.Merge.Size.pdf




### pvalues PAIRED T-TEST!!!! (Due to variability in technical replicates)
## MitoTracker:  AZ increase significant? 
t.test(subset(HCT116.Merge, Sample.ID. == "WT")$Hoechst.Mean.PE.A, 
       subset(HCT116.Merge, Sample.ID. == "WT_AZ")$Hoechst.Mean.PE.A, paired=TRUE)
#WT vs AZ: 0.006296
t.test(subset(HCT116.Merge, Sample.ID. == "OE")$Hoechst.Mean.PE.A, 
       subset(HCT116.Merge, Sample.ID. == "OE_AZ")$Hoechst.Mean.PE.A, paired=TRUE)
#OE  vs AZ: 0.01437
t.test(subset(HCT116.Merge, Sample.ID. == "KO1")$Hoechst.Mean.PE.A, 
       subset(HCT116.Merge, Sample.ID. == "KO1_AZ")$Hoechst.Mean.PE.A, paired=TRUE)
#KO1  vs AZ: 0.0003994
t.test(subset(HCT116.Merge, Sample.ID. == "KO2")$Hoechst.Mean.PE.A, 
       subset(HCT116.Merge, Sample.ID. == "KO2_AZ")$Hoechst.Mean.PE.A, paired=TRUE)
#KO2 vs AZ: 0.001109 
t.test(subset(HCT116.Merge, Sample.ID. == "KO3")$Hoechst.Mean.PE.A, 
       subset(HCT116.Merge, Sample.ID. == "KO3_AZ")$Hoechst.Mean.PE.A, paired=TRUE)
#KO3 vs AZ: 0.0001564
t.test(subset(HCT116.Merge, Sample.ID. == "KO1_OE")$Hoechst.Mean.PE.A, 
       subset(HCT116.Merge, Sample.ID. == "KO1_OE_AZ")$Hoechst.Mean.PE.A, paired=TRUE)
#KO1OE vs AZ: 9.825e-05
t.test(subset(HCT116.Merge, Sample.ID. == "KO2_OE")$Hoechst.Mean.PE.A, 
       subset(HCT116.Merge, Sample.ID. == "KO2_OE_AZ")$Hoechst.Mean.PE.A, paired=TRUE)
#KO2OE vs AZ: 6.406e-05
t.test(subset(HCT116.Merge, Sample.ID. == "KO3_OE")$Hoechst.Mean.PE.A, 
       subset(HCT116.Merge, Sample.ID. == "KO3_OE_AZ")$Hoechst.Mean.PE.A, paired=TRUE)$p.value
#KO3OE vs AZ: 1.314679e-06



## MitoTracker:  AZ mito treated greater than WT+AZ? 
t.test(subset(HCT116.Merge, Sample.ID. == "OE_AZ")$Hoechst.Mean.PE.A, 
       subset(HCT116.Merge, Sample.ID. == "WT_AZ")$Hoechst.Mean.PE.A, paired=TRUE)
#OE AZ  vs WTAZ: 0.01952 LESS
t.test(subset(HCT116.Merge, Sample.ID. == "KO1_AZ")$Hoechst.Mean.PE.A, 
       subset(HCT116.Merge, Sample.ID. == "WT_AZ")$Hoechst.Mean.PE.A, paired=TRUE)
#KO1AZ  vs WTAZ: 9.75e-06
t.test(subset(HCT116.Merge, Sample.ID. == "KO2_AZ")$Hoechst.Mean.PE.A, 
       subset(HCT116.Merge, Sample.ID. == "WT_AZ")$Hoechst.Mean.PE.A, paired=TRUE)
#KO2AZ vs WTAZ: 2.805e-05 
t.test(subset(HCT116.Merge, Sample.ID. == "KO3_AZ")$Hoechst.Mean.PE.A, 
       subset(HCT116.Merge, Sample.ID. == "WT_AZ")$Hoechst.Mean.PE.A[4:12], paired=TRUE)# t-test between paired samples of technical replicates
#KO3AZ vs WTAZ: 0.6627 NS
t.test(subset(HCT116.Merge, Sample.ID. == "KO1_OE_AZ")$Hoechst.Mean.PE.A, 
       subset(HCT116.Merge, Sample.ID. == "WT_AZ")$Hoechst.Mean.PE.A, paired=TRUE)
#KO1OEAZ vs WTAZ: 0.03635 
t.test(subset(HCT116.Merge, Sample.ID. == "KO2_OE_AZ")$Hoechst.Mean.PE.A, 
       subset(HCT116.Merge, Sample.ID. == "WT_AZ")$Hoechst.Mean.PE.A, paired=TRUE)
#KO2OEAZ vs WTAZ: 0.0546 NS
t.test(subset(HCT116.Merge, Sample.ID. == "KO3_OE_AZ")$Hoechst.Mean.PE.A, 
       subset(HCT116.Merge, Sample.ID. == "WT_AZ")$Hoechst.Mean.PE.A, paired=TRUE)$p.value
#KO3OEAZ vs WTAZ: 0.8292314 NS



## MitoTracker:  KO+AZ va KO+OE+AZ? 
t.test(subset(HCT116.Merge, Sample.ID. == "KO1_AZ")$Hoechst.Mean.PE.A, 
       subset(HCT116.Merge, Sample.ID. == "KO1_OE_AZ")$Hoechst.Mean.PE.A, paired=TRUE)
#KO1 AZ  vs KO1.oe.AZ: 0.2303 NS lower
t.test(subset(HCT116.Merge, Sample.ID. == "KO2_AZ")$Hoechst.Mean.PE.A, 
       subset(HCT116.Merge, Sample.ID. == "KO2_OE_AZ")$Hoechst.Mean.PE.A, paired=TRUE)
#KO2 AZ  vs KO2.oe.AZ: 0.04487 
t.test(subset(HCT116.Merge, Sample.ID. == "KO3_AZ")$Hoechst.Mean.PE.A, 
       subset(HCT116.Merge, Sample.ID. == "KO3_OE_AZ")$Hoechst.Mean.PE.A[4:12], paired=TRUE) # t-test between paired samples of technical replicates
#KO3 AZ  vs KO3.oe.AZ: 0.1799  lower



# change in mitotracker upon AZ3146
HCT116.Merge$MitoTracker.AZ.EtOH<- c(NA,NA,NA, NA,NA,NA, NA,NA,NA, NA,NA,NA, #HCT116 WT KO
                                     NA,NA,NA, NA,NA,NA, NA,NA,NA, NA,NA,NA, 
                                     HCT116.Merge$Hoechst.Mean.PE.A[25:36]/HCT116.Merge$Hoechst.Mean.PE.A[1:12], #HCT116 WT & KO: AZ treated divided by EtOH treated
                                     HCT116.Merge$Hoechst.Mean.PE.A[46:48]/HCT116.Merge$Hoechst.Mean.PE.A[13:15],# HCT116 OE
                                     HCT116.Merge$Hoechst.Mean.PE.A[37:45]/HCT116.Merge$Hoechst.Mean.PE.A[16:24], # HCT116 OE+KO
                                     
                                     NA,NA,NA, NA,NA,NA, NA,NA,NA, 
                                     NA,NA,NA, NA,NA,NA, NA,NA,NA, NA,NA,NA, 
                                     HCT116.Merge$Hoechst.Mean.PE.A[70:90]/HCT116.Merge$Hoechst.Mean.PE.A[49:69], #HCT116: AZ treated divided by EtOH treated
                                     
                                     NA,NA,NA, NA,NA,NA, NA,NA,NA, NA,NA,NA, #HCT116 WT KO
                                     NA,NA,NA, NA,NA,NA, NA,NA,NA, NA,NA,NA, 
                                     HCT116.Merge$Hoechst.Mean.PE.A[115:138]/HCT116.Merge$Hoechst.Mean.PE.A[91:114], 
                                     
                                     NA,NA,NA, NA,NA,NA, NA,NA,NA, NA,NA,NA, #HCT116 WT KO
                                     NA,NA,NA, NA,NA,NA, NA,NA,NA, NA,NA,NA, 
                                     HCT116.Merge$Hoechst.Mean.PE.A[163:186]/HCT116.Merge$Hoechst.Mean.PE.A[139:162]) #HCT116 WT & KO: AZ treated divided by EtOH treated)
ggplot(subset(HCT116.Merge, Sample.ID. %in% c("WT_AZ", "OE_AZ","KO1_AZ", "KO2_AZ", "KO3_AZ", "KO1_OE_AZ", "KO2_OE_AZ", "KO3_OE_AZ")) , aes(x=Sample.ID., y=MitoTracker.AZ.EtOH, color= Type ))+
  geom_point()+
  scale_color_manual(values=c("cyan3", "cyan4", "purple1", "purple4"))+
  stat_summary(fun.y=mean, geom="point", shape="-", size=10, color="red", fill="red") +
  geom_hline(yintercept = mean(subset(HCT116.Merge, Sample.ID. == "WT")$MitoTracker.AZ.EtOH) ) + # average WT mitotracker expression
  ylab("Change in MitoTracker upon AZ3146")+
  xlab("")+  
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
# UBE2H-KO Increases Mitotracker. OE reduces it. AZ3146 increases it. 



### HCT116 Aug 17 2026 :KO + OE clones #### 
# HCT116 cells treated with 1uM AZ3146 for 96 hours
# HCT116 KO plus UBE2H overexpression single cell clones

colnames(HCT116.817)

# Assign groups: 
# Remove the last threee characters from "Tube.Name." (..._01, ..._02, ..._03, etc.) and 
# and remove the first three characters from "Tube.Name." (01_...) 
HCT116.817$Sample.ID.<- str_sub(HCT116.817$Tube.Name., start=4, end = -4)
HCT116.817$Sample.ID.<- c( "WT" ,          "WT"    ,       "WT"    ,       "WT.AZ" ,       "WT.AZ"  ,      "WT.AZ"   ,     "OE"      ,     "OE"    ,      
   "OE"  ,         "OE.AZ"  ,      "OE.AZ" ,       "OE.AZ"  ,      "KO1"     ,     "KO1"     ,     "KO1"     ,     "KO1.AZ"   ,   
   "KO1.AZ",       "KO1.AZ" ,      "KO2"  ,        "KO2"   ,       "KO2"    ,      "KO2.AZ"  ,     "KO2.AZ"  ,     "KO2.AZ"   ,   
   "KO3" ,         "KO3"   ,       "KO3"    ,      "KO3.AZ"  ,     "KO3.AZ"  ,     "KO3.AZ"  ,     "KO1.OEc1"  ,   "KO1.OEc1"   , 
   "KO1.OEc1",     "KO1.OEc1.AZ",  "KO1.OEc1.AZ" , "KO1.OEc1.AZ",  "KO2.OEc1" ,    "KO2.OEc1" ,    "KO2.OEc1"   ,  "KO2.OEc1.AZ" ,
   "KO2.OEc1.AZ",  "KO2.OEc1.AZ",  "KO2.OEc2" ,    "KO2.OEc2"   ,  "KO2.OEc2" ,    "KO2.OEc2.AZ" , "KO2.OEc2.AZ" , "KO2.OEc2.AZ" ,
   "KO3.OEc2",     "KO3.OEc2"  ,   "KO3.OEc2" ,    "KO3.OEc2.AZ" , "KO3.OEc2.AZ" , "KO3.OEc2.AZ" , "KO3.OEc5"  ,   "KO3.OEc5"  , 
   "KO3.OEc5",     "KO3.OEc5.AZ",  "KO3.OEc5.AZ",  "KO3.OEc5.AZ")
HCT116.817$Sample.ID.<- factor(HCT116.817$Sample.ID., levels= c("WT", "WT.AZ", "OE", "OE.AZ", 
                                                                "KO1", "KO1.AZ", "KO1.OEc1", "KO1.OEc1.AZ",
                                                                "KO2", "KO2.AZ", "KO2.OEc1",  "KO2.OEc1.AZ","KO2.OEc2",  "KO2.OEc2.AZ",
                                                                "KO3",  "KO3.AZ",  "KO3.OEc2",  "KO3.OEc2.AZ",  "KO3.OEc5",  "KO3.OEc5.AZ"
))


# now merge KO clones as a single "type" 
HCT116.817$Type<- c("WT","WT","WT", "WT+AZ","WT+AZ","WT+AZ", 
                    "OE","OE","OE", "OE+AZ","OE+AZ","OE+AZ",
                    "KO","KO","KO", "KO+AZ","KO+AZ","KO+AZ",
                    "KO","KO","KO", "KO+AZ","KO+AZ","KO+AZ",
                    "KO","KO","KO", "KO+AZ","KO+AZ","KO+AZ",
                    "KO+OE","KO+OE","KO+OE",  "KO+OE+AZ","KO+OE+AZ","KO+OE+AZ",
                    "KO+OE","KO+OE","KO+OE", "KO+OE+AZ","KO+OE+AZ","KO+OE+AZ",
                    "KO+OE","KO+OE","KO+OE", "KO+OE+AZ","KO+OE+AZ","KO+OE+AZ",
                    "KO+OE","KO+OE","KO+OE", "KO+OE+AZ","KO+OE+AZ","KO+OE+AZ",
                    "KO+OE","KO+OE","KO+OE", "KO+OE+AZ","KO+OE+AZ","KO+OE+AZ"
)
HCT116.817$Type<- factor(HCT116.817$Type, levels= c("WT", "WT+AZ", "OE", "OE+AZ", "KO", "KO+AZ",  "KO+OE", "KO+OE+AZ"))



# Cell growth (96 hours)
## Note that 30ul from a total of 200ul is measured. so multiuple # cells by 200/30
## Count number of cells: cell sized & Hoechst positive. 
MeanWTCount<- mean(subset(HCT116.817, Sample.ID. == "WT")$Hoechst.Events)
HCT116.817$RelativeCount<- HCT116.817$Hoechst.Events/MeanWTCount
HCT116.817$Doubling<- HCT116.817$Hoechst.Events*200/(30*100000*4) # doubling: final/initial 
# Final is number of Hoechst + cells *200/30 (only analyzed 30ul from 200ul total), started with 100k cells, grew for 4 days so divide by 4 to get doubling per day



ggplot(HCT116.817, aes(x=HCT116.817$Sample.ID., y=Doubling, color=Type ))+
  geom_point()+
  scale_color_manual(values=c("cyan3","cyan3", "cyan4", "cyan4", "purple1","purple1", "purple4",  "purple4"), guide="none")+
  stat_summary(fun.y=mean, geom="point", shape="-", size=10, color="red", fill="red") +
  ylab("Doubling per day")+
  xlab("")+  
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
# plot.HCT116.OEc.Doubling.pdf


####
## plot mean mitoTracker (PE) 
HCT116.817$RelativeMito<- HCT116.817$Hoechst.Mean.PE.A/ mean(subset(HCT116.817, Sample.ID. == "WT")$Hoechst.Mean.PE.A)

ggplot(HCT116.817, aes(x=HCT116.817$Sample.ID., y=HCT116.817$Hoechst.Mean.PE.A, color=Type ))+
  geom_point()+
  scale_color_manual(values=c("cyan3","cyan3", "cyan4", "cyan4", "purple1","purple1", "purple4",  "purple4"), guide="none")+
  ylab("Mean Mitotracker")+
  geom_hline(yintercept = mean(subset(HCT116.817, Sample.ID. == "WT")$Hoechst.Mean.PE.A) ) + # average WT mitotracker expression
  #geom_hline(yintercept = mean(subset(HCT116.817, Sample.ID. == "WT_AZ")$Hoechst.Mean.PE.A) ) + # average WT mitotracker expression
  stat_summary(fun.y=mean, geom="point", shape="-", size=10, color="red", fill="red") +
  xlab("")+  
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
#plot.HCT116.OEc.Mito.pdf

HCT116.817$AZMito<- c( HCT116.817$Hoechst.Mean.PE.A[1:6]/mean(HCT116.817$Hoechst.Mean.PE.A[1:3]), 
                       HCT116.817$Hoechst.Mean.PE.A[7:12]/mean(HCT116.817$Hoechst.Mean.PE.A[7:9]), 
                       
                       HCT116.817$Hoechst.Mean.PE.A[13:18]/mean(HCT116.817$Hoechst.Mean.PE.A[13:15]), 
                       HCT116.817$Hoechst.Mean.PE.A[19:24]/mean(HCT116.817$Hoechst.Mean.PE.A[19:21]), 
                       HCT116.817$Hoechst.Mean.PE.A[25:30]/mean(HCT116.817$Hoechst.Mean.PE.A[25:27]), 
                       
                       HCT116.817$Hoechst.Mean.PE.A[31:36]/mean(HCT116.817$Hoechst.Mean.PE.A[31:33]), 
                       HCT116.817$Hoechst.Mean.PE.A[37:42]/mean(HCT116.817$Hoechst.Mean.PE.A[37:39]), 
                       HCT116.817$Hoechst.Mean.PE.A[43:48]/mean(HCT116.817$Hoechst.Mean.PE.A[43:45]), 
                       HCT116.817$Hoechst.Mean.PE.A[49:54]/mean(HCT116.817$Hoechst.Mean.PE.A[49:51]), 
                       HCT116.817$Hoechst.Mean.PE.A[55:60]/mean(HCT116.817$Hoechst.Mean.PE.A[55:57])
                       )
ggplot(subset(HCT116.817, Sample.ID. %in% c("WT.AZ", "OE.AZ", "KO1.AZ", "KO2.AZ", "KO3.AZ", "KO1.OEc1.AZ", "KO2.OEc1.AZ","KO2.OEc2.AZ","KO3.OEc2.AZ","KO3.OEc5.AZ")), 
       aes(x=Sample.ID., y=AZMito, color=Type ))+
  geom_point()+
  scale_color_manual(values=c("cyan3", "cyan4", "purple1",  "purple4"), guide="none")+
  ylab("Relative Mitotracker upon AZ3146")+
  geom_hline(yintercept = mean(subset(HCT116.817, Sample.ID. == "WT")$AZMito) ) + # average WT mitotracker expression
  stat_summary(fun.y=mean, geom="point", shape="-", size=10, color="red", fill="red") +
  xlab("")+  
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
#plot.HCT116.OEc.MitoAZ.pdf


## plot mean  cell size
ggplot(HCT116.817, aes(x=HCT116.817$Sample.ID., y=HCT116.817$Hoechst.Mean.FSC.A, color=Type ))+
  geom_point()+
  scale_color_manual(values=c("cyan3","cyan3", "cyan4", "cyan4", "purple1","purple1", "purple4",  "purple4"), guide="none")+
  stat_summary(fun.y=mean, geom="point", shape="-", size=10, color="red", fill="red") +
  geom_hline(yintercept = mean(subset(HCT116.817, Sample.ID. == "WT")$Hoechst.Mean.FSC.A) ) + # average WT mitotracker expression
  ylab("Mean Size")+
  xlab("")+  
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
# UBE2H-KO reduces cell size (KO+OE rescues cell size), OE increases cell size. 
#plot.HCT116.OEc.Size.pdf




## plot mean  Hoechst
meanDNA<- mean(subset(HCT116.817, Sample.ID. == "WT")$Hoechst.Mean.Violet610.A)
HCT116.817$RelativeHoechst <- HCT116.817$Hoechst.Mean.Violet610.A/meanDNA
ggplot(HCT116.817, aes(x=HCT116.817$Sample.ID., y=HCT116.817$RelativeHoechst , color=Type))+
  geom_point()+
  scale_color_manual(values=c("cyan3","cyan3", "cyan4", "cyan4", "purple1","purple1", "purple4",  "purple4"), guide="none")+
  geom_hline(yintercept = mean(subset(HCT116.817, Sample.ID. == "WT")$RelativeHoechst) ) + # average WT mitotracker expression
  stat_summary(fun.y=mean, geom="point", shape="-", size=10, color="red", fill="red") +
  ylab("Mean Hoechst Fluoescence")+
  xlab("")+  
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
# HCT116.Hoechst.KOOEc_relative.pdf
# UBE2H-KO + AZ significantly increase DNA/Hoechst content


### MitoTracker per DNA Hoechst 
HCT116.817$Mito.div.Hoechst<- HCT116.817$Hoechst.Mean.PE.A/ HCT116.817$Hoechst.Mean.Violet610.A
ggplot(HCT116.817, aes(x=HCT116.817$Sample.ID., y=HCT116.817$Mito.div.Hoechst, color= Type ))+
  geom_point()+
  stat_summary(fun.y=mean, geom="point", shape="-", size=10, color="red", fill="red") +
  scale_color_manual(values=c("cyan3","cyan3", "cyan4", "cyan4", "purple1","purple1", "purple4",  "purple4"), guide="none")+
  geom_hline(yintercept = mean(subset(HCT116.817, Sample.ID. == "WT")$Mito.div.Hoechst) ) + # average WT mitotracker expression
  ylab("MitoTracker/DNA (Hoechst)")+
  xlab("")+  
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
# UBE2H-KO Increases Mitotracker/DNA. OE reduces it. AZ3146 increases it. 


### MitoTracker per cell size 
HCT116.817$Mito.div.FSC<- HCT116.817$Hoechst.Mean.PE.A/ HCT116.817$Hoechst.Mean.FSC.A
ggplot(HCT116.817, aes(x=HCT116.817$Sample.ID., y=HCT116.817$Mito.div.FSC, color= Type ))+
  geom_point()+
  stat_summary(fun.y=mean, geom="point", shape="-", size=10, color="red", fill="red") +
  scale_color_manual(values=c("cyan3","cyan3", "cyan4", "cyan4", "purple1","purple1", "purple4",  "purple4"), guide="none")+
  geom_hline(yintercept = mean(subset(HCT116.817, Sample.ID. == "WT")$Mito.div.FSC) ) + # average WT mitotracker expression
  ylab("MitoTracker/Cell Size")+
  xlab("")+  
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
# UBE2H-KO Increases Mitotracker/Size. OE reduces it. AZ3146 increases it. 
#plot.HCT116.OEc.MitoSize.pdf





### pvalues
## MitoTracker:  AZ increase significant? 
t.test(subset(HCT116.817, Sample.ID. == "WT")$Hoechst.Mean.PE.A, subset(HCT116.817, Sample.ID. == "WT.AZ")$Hoechst.Mean.PE.A)
#WT vs AZ: 0.04671
t.test(subset(HCT116.817, Sample.ID. == "OE")$Hoechst.Mean.PE.A, subset(HCT116.817, Sample.ID. == "OE.AZ")$Hoechst.Mean.PE.A)
#OE  vs AZ: 0.1442
t.test(subset(HCT116.817, Sample.ID. == "KO1")$Hoechst.Mean.PE.A, subset(HCT116.817, Sample.ID. == "KO1.AZ")$Hoechst.Mean.PE.A)
#KO1  vs AZ: 0.0005526
t.test(subset(HCT116.817, Sample.ID. == "KO2")$Hoechst.Mean.PE.A, subset(HCT116.817, Sample.ID. == "KO2.AZ")$Hoechst.Mean.PE.A)
#KO2 vs AZ:  0.00754
t.test(subset(HCT116.817, Sample.ID. == "KO3")$Hoechst.Mean.PE.A, subset(HCT116.817, Sample.ID. == "KO3.AZ")$Hoechst.Mean.PE.A)
#KO3 vs AZ: 0.00424
t.test(subset(HCT116.817, Sample.ID. == "KO1.OEc1")$Hoechst.Mean.PE.A, subset(HCT116.817, Sample.ID. == "KO1.OEc1.AZ")$Hoechst.Mean.PE.A)
#KO1OE vs AZ: 0.8959
t.test(subset(HCT116.817, Sample.ID. == "KO2.OEc1")$Hoechst.Mean.PE.A, subset(HCT116.817, Sample.ID. == "KO2.OEc1.AZ")$Hoechst.Mean.PE.A)
#KO2OE vs AZ: 0.0001035 (Very low UBE2H)
t.test(subset(HCT116.817, Sample.ID. == "KO2.OEc2")$Hoechst.Mean.PE.A, subset(HCT116.817, Sample.ID. == "KO2.OEc2.AZ")$Hoechst.Mean.PE.A)
#KO2OE vs AZ: 0.1521
t.test(subset(HCT116.817, Sample.ID. == "KO3.OEc2")$Hoechst.Mean.PE.A, subset(HCT116.817, Sample.ID. == "KO3.OEc2.AZ")$Hoechst.Mean.PE.A)
#KO3OE vs AZ: 0.5759836
t.test(subset(HCT116.817, Sample.ID. == "KO3.OEc5")$Hoechst.Mean.PE.A, subset(HCT116.817, Sample.ID. == "KO3.OEc5.AZ")$Hoechst.Mean.PE.A)
#KO3OE vs AZ: 0.9372619



## MitoTracker:  AZ mito treated greater than WT+AZ? 
t.test(subset(HCT116.817, Sample.ID. == "OE.AZ")$Hoechst.Mean.PE.A, subset(HCT116.817, Sample.ID. == "WT.AZ")$Hoechst.Mean.PE.A)
#OE AZ  vs WTAZ: 0.7907 NS
t.test(subset(HCT116.817, Sample.ID. == "KO1.AZ")$Hoechst.Mean.PE.A, subset(HCT116.817, Sample.ID. == "WT.AZ")$Hoechst.Mean.PE.A)
#KO1AZ  vs WTAZ: 0.001845
t.test(subset(HCT116.817, Sample.ID. == "KO2.AZ")$Hoechst.Mean.PE.A, subset(HCT116.817, Sample.ID. == "WT.AZ")$Hoechst.Mean.PE.A)
#KO2AZ vs WTAZ:  0.004474
t.test(subset(HCT116.817, Sample.ID. == "KO3.AZ")$Hoechst.Mean.PE.A, subset(HCT116.817, Sample.ID. == "WT.AZ")$Hoechst.Mean.PE.A)
#KO3AZ vs WTAZ: 0.2946 NS
t.test(subset(HCT116.817, Sample.ID. == "KO1.OEc1.AZ")$Hoechst.Mean.PE.A, subset(HCT116.817, Sample.ID. == "WT.AZ")$Hoechst.Mean.PE.A)
#KO1OEAZ vs WTAZ: 4.567e-05
t.test(subset(HCT116.817, Sample.ID. == "KO2.OEc1.AZ")$Hoechst.Mean.PE.A, subset(HCT116.817, Sample.ID. == "WT.AZ")$Hoechst.Mean.PE.A)
#KO2OEAZ vs WTAZ: 0.0001444
t.test(subset(HCT116.817, Sample.ID. == "KO2.OEc2.AZ")$Hoechst.Mean.PE.A, subset(HCT116.817, Sample.ID. == "WT.AZ")$Hoechst.Mean.PE.A)
#KO2OEAZ vs WTAZ: 0.3941 NS
t.test(subset(HCT116.817, Sample.ID. == "KO3.OEc2.AZ")$Hoechst.Mean.PE.A, subset(HCT116.817, Sample.ID. == "WT.AZ")$Hoechst.Mean.PE.A)
#KO3OEAZ vs WTAZ: 0.01439694
t.test(subset(HCT116.817, Sample.ID. == "KO3.OEc5.AZ")$Hoechst.Mean.PE.A, subset(HCT116.817, Sample.ID. == "WT.AZ")$Hoechst.Mean.PE.A)
#KO3OEAZ vs WTAZ: 0.2778 NS


## MitoTracker:  KO+AZ va KO+OE+AZ? 
t.test(subset(HCT116.817, Sample.ID. == "KO1.AZ")$Hoechst.Mean.PE.A, subset(HCT116.817, Sample.ID. == "KO1.OEc1.AZ")$Hoechst.Mean.PE.A)
#KO1 AZ  vs KO1.oe.AZ: 0.1311 NS
t.test(subset(HCT116.817, Sample.ID. == "KO2.AZ")$Hoechst.Mean.PE.A, subset(HCT116.817, Sample.ID. == "KO2.OEc1.AZ")$Hoechst.Mean.PE.A)
#KO2 AZ  vs KO2.oe.AZ: 0.005946
t.test(subset(HCT116.817, Sample.ID. == "KO2.AZ")$Hoechst.Mean.PE.A, subset(HCT116.817, Sample.ID. == "KO2.OEc2.AZ")$Hoechst.Mean.PE.A)
#KO2 AZ  vs KO2.oe.AZ: 0.001923
t.test(subset(HCT116.817, Sample.ID. == "KO3.AZ")$Hoechst.Mean.PE.A, subset(HCT116.817, Sample.ID. == "KO3.OEc2.AZ")$Hoechst.Mean.PE.A)
#KO3 AZ  vs KO3.oe.AZ: 0.01033
t.test(subset(HCT116.817, Sample.ID. == "KO3.AZ")$Hoechst.Mean.PE.A, subset(HCT116.817, Sample.ID. == "KO3.OEc5.AZ")$Hoechst.Mean.PE.A)
#KO3 AZ  vs KO3.oe.AZ: 0.1733 NS





### HCT116 Aug 25 2026 :KO + OE clones #### 
# HCT116 cells treated with 1uM AZ3146 for 96 hours
# HCT116 KO plus UBE2H overexpression single cell clones

colnames(HCT116.825)

# Assign groups: 
# Remove the last threee characters from "Tube.Name." (..._01, ..._02, ..._03, etc.) and 
# and remove the first three characters from "Tube.Name." (01_...) 
HCT116.825$Sample.ID.<- str_sub(HCT116.825$Tube.Name., start=4, end = -4)
HCT116.825$Sample.ID.<- c( "WT" ,          "WT"    ,       "WT"    ,      "KO1"     ,     "KO1"     ,     "KO1"     ,  
                           "KO2"  ,        "KO2"   ,       "KO2"    ,     "KO3" ,         "KO3"   ,       "KO3"    ,    
                           "WT.AZ" ,       "WT.AZ"  ,      "WT.AZ"   ,       "KO1.AZ"   ,  "KO1.AZ",       "KO1.AZ" ,     
                           "KO2.AZ"  ,     "KO2.AZ"  ,     "KO2.AZ"   ,    "KO3.AZ"  ,     "KO3.AZ"  ,     "KO3.AZ"  ,   
                           "KO1.OEc1"  ,   "KO1.OEc1"   , "KO1.OEc1",     "KO1.OEc1.AZ",  "KO1.OEc1.AZ" , "KO1.OEc1.AZ",  
                           "KO2.OEc1" ,    "KO2.OEc1" ,    "KO2.OEc1"   ,  "KO2.OEc1.AZ" ,"KO2.OEc1.AZ",  "KO2.OEc1.AZ", 
                           "KO2.OEc2" ,    "KO2.OEc2"   ,  "KO2.OEc2" ,    "KO2.OEc2.AZ" , "KO2.OEc2.AZ" , "KO2.OEc2.AZ" ,
                           "KO3.OEc2",     "KO3.OEc2"  ,   "KO3.OEc2" ,    "KO3.OEc2.AZ" , "KO3.OEc2.AZ" , "KO3.OEc2.AZ" , 
                           "KO3.OEc5"  ,   "KO3.OEc5"  ,  "KO3.OEc5",     "KO3.OEc5.AZ",  "KO3.OEc5.AZ",  "KO3.OEc5.AZ",   
                           "OE"      ,     "OE"    ,        "OE"  ,         "OE.AZ"  ,      "OE.AZ" ,       "OE.AZ"  )
HCT116.825$Sample.ID.<- factor(HCT116.825$Sample.ID., levels= c("WT", "WT.AZ", "OE", "OE.AZ", 
                                                                "KO1", "KO1.AZ", "KO1.OEc1", "KO1.OEc1.AZ",
                                                                "KO2", "KO2.AZ", "KO2.OEc1",  "KO2.OEc1.AZ","KO2.OEc2",  "KO2.OEc2.AZ",
                                                                "KO3",  "KO3.AZ",  "KO3.OEc2",  "KO3.OEc2.AZ",  "KO3.OEc5",  "KO3.OEc5.AZ"
))


# now merge KO clones as a single "type" 
HCT116.825$Type<- c("WT","WT","WT", "KO","KO","KO", 
                    "KO","KO","KO", "KO","KO","KO", 
                    "WT+AZ","WT+AZ","WT+AZ",  "KO+AZ","KO+AZ","KO+AZ",
                    "KO+AZ","KO+AZ","KO+AZ", "KO+AZ","KO+AZ","KO+AZ",
                    
                    "KO+OE","KO+OE","KO+OE", "KO+OE+AZ","KO+OE+AZ","KO+OE+AZ",
                    "KO+OE","KO+OE","KO+OE", "KO+OE+AZ","KO+OE+AZ","KO+OE+AZ",
                    "KO+OE","KO+OE","KO+OE", "KO+OE+AZ","KO+OE+AZ","KO+OE+AZ",
                    "KO+OE","KO+OE","KO+OE", "KO+OE+AZ","KO+OE+AZ","KO+OE+AZ",
                    "KO+OE","KO+OE","KO+OE", "KO+OE+AZ","KO+OE+AZ","KO+OE+AZ", 
                    "OE","OE","OE", "OE+AZ","OE+AZ","OE+AZ"
)
HCT116.825$Type<- factor(HCT116.825$Type, levels= c("WT", "WT+AZ", "OE", "OE+AZ", "KO", "KO+AZ",  "KO+OE", "KO+OE+AZ"))



# Cell growth (96 hours)
## Note that 30ul from a total of 200ul is measured. so multiuple # cells by 200/30
## Count number of cells: cell sized & Hoechst positive. 
MeanWTCount<- mean(subset(HCT116.825, Sample.ID. == "WT")$Hoechst.Events)
HCT116.825$RelativeCount<- HCT116.825$Hoechst.Events/MeanWTCount
HCT116.825$Doubling<- HCT116.825$Hoechst.Events*200/(30*100000*4) # doubling: final/initial 
# Final is number of Hoechst + cells *200/30 (only analyzed 30ul from 200ul total), started with 100k cells, grew for 4 days so divide by 4 to get doubling per day



ggplot(HCT116.825, aes(x=HCT116.825$Sample.ID., y=Doubling, color=Type ))+
  geom_point()+
  scale_color_manual(values=c("cyan3","cyan3", "cyan4", "cyan4", "purple1","purple1", "purple4",  "purple4"), guide="none")+
  stat_summary(fun.y=mean, geom="point", shape="-", size=10, color="red", fill="red") +
  ylab("Doubling per day")+
  xlab("")+  
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
# plot.HCT116.OEc.Doubling_2.pdf


####
## plot mean mitoTracker (PE) 
HCT116.825$RelativeMito<- HCT116.825$Hoechst.Mean.PE.A/ mean(subset(HCT116.825, Sample.ID. == "WT")$Hoechst.Mean.PE.A)

ggplot(HCT116.825, aes(x=HCT116.825$Sample.ID., y=HCT116.825$Hoechst.Mean.PE.A, color=Type ))+
  geom_point()+
  scale_color_manual(values=c("cyan3","cyan3", "cyan4", "cyan4", "purple1","purple1", "purple4",  "purple4"), guide="none")+
  ylab("Mean Mitotracker")+
  geom_hline(yintercept = mean(subset(HCT116.825, Sample.ID. == "WT")$Hoechst.Mean.PE.A) ) + # average WT mitotracker expression
  #geom_hline(yintercept = mean(subset(HCT116.825, Sample.ID. == "WT_AZ")$Hoechst.Mean.PE.A) ) + # average WT mitotracker expression
  stat_summary(fun.y=mean, geom="point", shape="-", size=10, color="red", fill="red") +
  xlab("")+  
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
#plot.HCT116.OEc.Mito_2.pdf

HCT116.825$AZMito<- c( HCT116.825$Hoechst.Mean.PE.A[c(1:3)]/mean(HCT116.825$Hoechst.Mean.PE.A[1:3]), 
                       HCT116.825$Hoechst.Mean.PE.A[c(4:6)]/mean(HCT116.825$Hoechst.Mean.PE.A[4:6]), 
                       HCT116.825$Hoechst.Mean.PE.A[c(7:9)]/mean(HCT116.825$Hoechst.Mean.PE.A[7:9]), 
                       HCT116.825$Hoechst.Mean.PE.A[c(10:12)]/mean(HCT116.825$Hoechst.Mean.PE.A[10:12]), 
                       HCT116.825$Hoechst.Mean.PE.A[c(13:15)]/mean(HCT116.825$Hoechst.Mean.PE.A[1:3]), 
                       HCT116.825$Hoechst.Mean.PE.A[c(16:18)]/mean(HCT116.825$Hoechst.Mean.PE.A[4:6]), 
                       HCT116.825$Hoechst.Mean.PE.A[c(19:21)]/mean(HCT116.825$Hoechst.Mean.PE.A[7:9]), 
                       HCT116.825$Hoechst.Mean.PE.A[c(22:24)]/mean(HCT116.825$Hoechst.Mean.PE.A[10:12]), 
                       
                       HCT116.825$Hoechst.Mean.PE.A[25:30]/mean(HCT116.825$Hoechst.Mean.PE.A[25:27]), 
                       HCT116.825$Hoechst.Mean.PE.A[31:36]/mean(HCT116.825$Hoechst.Mean.PE.A[31:33]), 
                       HCT116.825$Hoechst.Mean.PE.A[37:42]/mean(HCT116.825$Hoechst.Mean.PE.A[37:39]), 
                       HCT116.825$Hoechst.Mean.PE.A[43:48]/mean(HCT116.825$Hoechst.Mean.PE.A[43:45]), 
                       HCT116.825$Hoechst.Mean.PE.A[49:54]/mean(HCT116.825$Hoechst.Mean.PE.A[49:51]), 
                       
                       HCT116.825$Hoechst.Mean.PE.A[55:60]/mean(HCT116.825$Hoechst.Mean.PE.A[55:57])
)
ggplot(subset(HCT116.825, Sample.ID. %in% c("WT.AZ", "OE.AZ", "KO1.AZ", "KO2.AZ", "KO3.AZ", "KO1.OEc1.AZ", "KO2.OEc1.AZ","KO2.OEc2.AZ","KO3.OEc2.AZ","KO3.OEc5.AZ")), 
       aes(x=Sample.ID., y=AZMito, color=Type ))+
  geom_point()+
  scale_color_manual(values=c("cyan3", "cyan4", "purple1",  "purple4"), guide="none")+
  ylab("Relative Mitotracker upon AZ3146")+
  geom_hline(yintercept = mean(subset(HCT116.825, Sample.ID. == "WT")$AZMito) ) + # average WT mitotracker expression
  stat_summary(fun.y=mean, geom="point", shape="-", size=10, color="red", fill="red") +
  xlab("")+  
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
#plot.HCT116.OEc.MitoAZ_2.pdf


## plot mean  cell size
ggplot(HCT116.825, aes(x=HCT116.825$Sample.ID., y=HCT116.825$Hoechst.Mean.FSC.A, color=Type ))+
  geom_point()+
  scale_color_manual(values=c("cyan3","cyan3", "cyan4", "cyan4", "purple1","purple1", "purple4",  "purple4"), guide="none")+
  stat_summary(fun.y=mean, geom="point", shape="-", size=10, color="red", fill="red") +
  geom_hline(yintercept = mean(subset(HCT116.825, Sample.ID. == "WT")$Hoechst.Mean.FSC.A) ) + # average WT mitotracker expression
  ylab("Mean Size")+
  xlab("")+  
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
# UBE2H-KO reduces cell size (KO+OE rescues cell size), OE increases cell size. 
#plot.HCT116.OEc.Size_2.pdf




## plot mean  Hoechst
meanDNA<- mean(subset(HCT116.825, Sample.ID. == "WT")$Hoechst.Mean.Violet610.A)
HCT116.825$RelativeHoechst <- HCT116.825$Hoechst.Mean.Violet610.A/meanDNA
ggplot(HCT116.825, aes(x=HCT116.825$Sample.ID., y=HCT116.825$RelativeHoechst , color=Type))+
  geom_point()+
  scale_color_manual(values=c("cyan3","cyan3", "cyan4", "cyan4", "purple1","purple1", "purple4",  "purple4"), guide="none")+
  geom_hline(yintercept = mean(subset(HCT116.825, Sample.ID. == "WT")$RelativeHoechst) ) + # average WT mitotracker expression
  stat_summary(fun.y=mean, geom="point", shape="-", size=10, color="red", fill="red") +
  ylab("Mean Hoechst Fluoescence")+
  xlab("")+  
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
# HCT116.Hoechst.KOOEc_relative_2.pdf
# UBE2H-KO + AZ significantly increase DNA/Hoechst content


### MitoTracker per DNA Hoechst 
HCT116.825$Mito.div.Hoechst<- HCT116.825$Hoechst.Mean.PE.A/ HCT116.825$Hoechst.Mean.Violet610.A
ggplot(HCT116.825, aes(x=HCT116.825$Sample.ID., y=HCT116.825$Mito.div.Hoechst, color= Type ))+
  geom_point()+
  stat_summary(fun.y=mean, geom="point", shape="-", size=10, color="red", fill="red") +
  scale_color_manual(values=c("cyan3","cyan3", "cyan4", "cyan4", "purple1","purple1", "purple4",  "purple4"), guide="none")+
  geom_hline(yintercept = mean(subset(HCT116.825, Sample.ID. == "WT")$Mito.div.Hoechst) ) + # average WT mitotracker expression
  ylab("MitoTracker/DNA (Hoechst)")+
  xlab("")+  
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
# UBE2H-KO Increases Mitotracker/DNA. OE reduces it. AZ3146 increases it. 


### MitoTracker per cell size 
HCT116.825$Mito.div.FSC<- HCT116.825$Hoechst.Mean.PE.A/ HCT116.825$Hoechst.Mean.FSC.A
ggplot(HCT116.825, aes(x=HCT116.825$Sample.ID., y=HCT116.825$Mito.div.FSC, color= Type ))+
  geom_point()+
  stat_summary(fun.y=mean, geom="point", shape="-", size=10, color="red", fill="red") +
  scale_color_manual(values=c("cyan3","cyan3", "cyan4", "cyan4", "purple1","purple1", "purple4",  "purple4"), guide="none")+
  geom_hline(yintercept = mean(subset(HCT116.825, Sample.ID. == "WT")$Mito.div.FSC) ) + # average WT mitotracker expression
  ylab("MitoTracker/Cell Size")+
  xlab("")+  
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
# UBE2H-KO Increases Mitotracker/Size. OE reduces it. AZ3146 increases it. 
#plot.HCT116.OEc.MitoSize_2.pdf





### p-values
## MitoTracker:  AZ increase significant? 
t.test(subset(HCT116.825, Sample.ID. == "WT")$Hoechst.Mean.PE.A, subset(HCT116.825, Sample.ID. == "WT.AZ")$Hoechst.Mean.PE.A)
#WT vs AZ: 0.000938 ***
t.test(subset(HCT116.825, Sample.ID. == "OE")$Hoechst.Mean.PE.A, subset(HCT116.825, Sample.ID. == "OE.AZ")$Hoechst.Mean.PE.A)
#OE  vs AZ: 0.04328 *
t.test(subset(HCT116.825, Sample.ID. == "KO1")$Hoechst.Mean.PE.A, subset(HCT116.825, Sample.ID. == "KO1.AZ")$Hoechst.Mean.PE.A)
#KO1  vs AZ: 0.005597 *
t.test(subset(HCT116.825, Sample.ID. == "KO2")$Hoechst.Mean.PE.A, subset(HCT116.825, Sample.ID. == "KO2.AZ")$Hoechst.Mean.PE.A)
#KO2 vs AZ:  0.004761 **
t.test(subset(HCT116.825, Sample.ID. == "KO3")$Hoechst.Mean.PE.A, subset(HCT116.825, Sample.ID. == "KO3.AZ")$Hoechst.Mean.PE.A)
#KO3 vs AZ: 0.001234 **
t.test(subset(HCT116.825, Sample.ID. == "KO1.OEc1")$Hoechst.Mean.PE.A, subset(HCT116.825, Sample.ID. == "KO1.OEc1.AZ")$Hoechst.Mean.PE.A)
#KO1OE vs AZ: 0.09341 NS
t.test(subset(HCT116.825, Sample.ID. == "KO2.OEc1")$Hoechst.Mean.PE.A, subset(HCT116.825, Sample.ID. == "KO2.OEc1.AZ")$Hoechst.Mean.PE.A)
#KO2OE vs AZ: 0.06018 NS
t.test(subset(HCT116.825, Sample.ID. == "KO2.OEc2")$Hoechst.Mean.PE.A, subset(HCT116.825, Sample.ID. == "KO2.OEc2.AZ")$Hoechst.Mean.PE.A)
#KO2OE vs AZ: 0.7271 NS
t.test(subset(HCT116.825, Sample.ID. == "KO3.OEc2")$Hoechst.Mean.PE.A, subset(HCT116.825, Sample.ID. == "KO3.OEc2.AZ")$Hoechst.Mean.PE.A)
#KO3OE vs AZ: 0.09509 NS
t.test(subset(HCT116.825, Sample.ID. == "KO3.OEc5")$Hoechst.Mean.PE.A, subset(HCT116.825, Sample.ID. == "KO3.OEc5.AZ")$Hoechst.Mean.PE.A)
#KO3OE vs AZ: 0.0008758 ***



## MitoTracker:  AZ mito treated greater than WT+AZ? 
t.test(subset(HCT116.825, Sample.ID. == "OE.AZ")$Hoechst.Mean.PE.A, subset(HCT116.825, Sample.ID. == "WT.AZ")$Hoechst.Mean.PE.A)
#OE AZ  vs WTAZ: 0.0004524 ***
t.test(subset(HCT116.825, Sample.ID. == "KO1.AZ")$Hoechst.Mean.PE.A, subset(HCT116.825, Sample.ID. == "WT.AZ")$Hoechst.Mean.PE.A)
#KO1AZ  vs WTAZ: 0.007219 *
t.test(subset(HCT116.825, Sample.ID. == "KO2.AZ")$Hoechst.Mean.PE.A, subset(HCT116.825, Sample.ID. == "WT.AZ")$Hoechst.Mean.PE.A)
#KO2AZ vs WTAZ:  0.002769 **
t.test(subset(HCT116.825, Sample.ID. == "KO3.AZ")$Hoechst.Mean.PE.A, subset(HCT116.825, Sample.ID. == "WT.AZ")$Hoechst.Mean.PE.A)
#KO3AZ vs WTAZ: 0.007076 *
t.test(subset(HCT116.825, Sample.ID. == "KO1.OEc1.AZ")$Hoechst.Mean.PE.A, subset(HCT116.825, Sample.ID. == "WT.AZ")$Hoechst.Mean.PE.A)
#KO1OEAZ vs WTAZ: 0.00417 **
t.test(subset(HCT116.825, Sample.ID. == "KO2.OEc1.AZ")$Hoechst.Mean.PE.A, subset(HCT116.825, Sample.ID. == "WT.AZ")$Hoechst.Mean.PE.A)
#KO2OEAZ vs WTAZ: 0.05188 NS
t.test(subset(HCT116.825, Sample.ID. == "KO2.OEc2.AZ")$Hoechst.Mean.PE.A, subset(HCT116.825, Sample.ID. == "WT.AZ")$Hoechst.Mean.PE.A)
#KO2OEAZ vs WTAZ: 0.008035 *
t.test(subset(HCT116.825, Sample.ID. == "KO3.OEc2.AZ")$Hoechst.Mean.PE.A, subset(HCT116.825, Sample.ID. == "WT.AZ")$Hoechst.Mean.PE.A)
#KO3OEAZ vs WTAZ: 0.0003061 ***
t.test(subset(HCT116.825, Sample.ID. == "KO3.OEc5.AZ")$Hoechst.Mean.PE.A, subset(HCT116.825, Sample.ID. == "WT.AZ")$Hoechst.Mean.PE.A)
#KO3OEAZ vs WTAZ: 0.0002291  ***


## MitoTracker:  KO+AZ va KO+OE+AZ? 
t.test(subset(HCT116.825, Sample.ID. == "KO1.AZ")$Hoechst.Mean.PE.A, subset(HCT116.825, Sample.ID. == "KO1.OEc1.AZ")$Hoechst.Mean.PE.A)
#KO1 AZ  vs KO1.oe.AZ: 0.006047 *
t.test(subset(HCT116.825, Sample.ID. == "KO2.AZ")$Hoechst.Mean.PE.A, subset(HCT116.825, Sample.ID. == "KO2.OEc1.AZ")$Hoechst.Mean.PE.A)
#KO2 AZ  vs KO2.oe.AZ: 0.00014 ***
t.test(subset(HCT116.825, Sample.ID. == "KO2.AZ")$Hoechst.Mean.PE.A, subset(HCT116.825, Sample.ID. == "KO2.OEc2.AZ")$Hoechst.Mean.PE.A)
#KO2 AZ  vs KO2.oe.AZ: 0.0004326 ***
t.test(subset(HCT116.825, Sample.ID. == "KO3.AZ")$Hoechst.Mean.PE.A, subset(HCT116.825, Sample.ID. == "KO3.OEc2.AZ")$Hoechst.Mean.PE.A)
#KO3 AZ  vs KO3.oe.AZ: 2.849e-05 ***
t.test(subset(HCT116.825, Sample.ID. == "KO3.AZ")$Hoechst.Mean.PE.A, subset(HCT116.825, Sample.ID. == "KO3.OEc5.AZ")$Hoechst.Mean.PE.A)
#KO3 AZ  vs KO3.oe.AZ: 0.001175 **



### Change in MitoTracker upon AZ3146
## MitoTracker:  AZ mito treated greater than WT+AZ? 
t.test(subset(HCT116.825, Sample.ID. == "OE.AZ")$AZMito, 
       subset(HCT116.825, Sample.ID. == "WT.AZ")$AZMito)
#OE AZ  vs WTAZ: 0.02639
t.test(subset(HCT116.825, Sample.ID. == "KO1.AZ")$AZMito, 
       subset(HCT116.825, Sample.ID. == "WT.AZ")$AZMito)
#KO1AZ  vs WTAZ: 0.01056
t.test(subset(HCT116.825, Sample.ID. == "KO2.AZ")$AZMito, 
       subset(HCT116.825, Sample.ID. == "WT.AZ")$AZMito)
#KO2AZ vs WTAZ:  0.03329
t.test(subset(HCT116.825, Sample.ID. == "KO3.AZ")$AZMito, 
       subset(HCT116.825, Sample.ID. == "WT.AZ")$AZMito)
#KO3AZ vs WTAZ: 0.0009187
t.test(subset(HCT116.825, Sample.ID. == "KO1.OEc1.AZ")$AZMito, 
       subset(HCT116.825, Sample.ID. == "WT.AZ")$AZMito)
#KO1OEAZ vs WTAZ: 3.232e-05
t.test(subset(HCT116.825, Sample.ID. == "KO2.OEc1.AZ")$AZMito, 
       subset(HCT116.825, Sample.ID. == "WT.AZ")$AZMito)
#KO2OEAZ vs WTAZ: 0.62
t.test(subset(HCT116.825, Sample.ID. == "KO2.OEc2.AZ")$AZMito, 
       subset(HCT116.825, Sample.ID. == "WT.AZ")$AZMito)
#KO2OEAZ vs WTAZ:0.007682
t.test(subset(HCT116.825, Sample.ID. == "KO3.OEc2.AZ")$AZMito, 
       subset(HCT116.825, Sample.ID. == "WT.AZ")$AZMito)
#KO3OEAZ vs WTAZ: 0.02934
t.test(subset(HCT116.825, Sample.ID. == "KO3.OEc5.AZ")$AZMito, 
       subset(HCT116.825, Sample.ID. == "WT.AZ")$AZMito)
#KO3OEAZ vs WTAZ: 0.1222


## MitoTracker:  KO+AZ vs KO+OE+AZ? 
t.test(subset(HCT116.825, Sample.ID. == "KO1.AZ")$AZMito, 
       subset(HCT116.825, Sample.ID. == "KO1.OEc1.AZ")$AZMito)
#KO1 AZ  vs KO1.oe.AZ: 0.005496
t.test(subset(HCT116.825, Sample.ID. == "KO2.AZ")$AZMito, 
       subset(HCT116.825, Sample.ID. == "KO2.OEc1.AZ")$AZMito)
#KO2 AZ  vs KO2.oe.AZ: 0.03959
t.test(subset(HCT116.825, Sample.ID. == "KO2.AZ")$AZMito, 
       subset(HCT116.825, Sample.ID. == "KO2.OEc2.AZ")$AZMito)
#KO2 AZ  vs KO2.oe.AZ: 0.001089
t.test(subset(HCT116.825, Sample.ID. == "KO3.AZ")$AZMito, 
       subset(HCT116.825, Sample.ID. == "KO3.OEc2.AZ")$AZMito)
#KO3 AZ  vs KO3.oe.AZ: 0.0003137
t.test(subset(HCT116.825, Sample.ID. == "KO3.AZ")$AZMito, 
       subset(HCT116.825, Sample.ID. == "KO3.OEc5.AZ")$AZMito)
#KO3 AZ  vs KO3.oe.AZ: 0.002511



### HCT116 MERGE: KO OE clones  #### 
# HCT116 cells treated with 1uM AZ3146 for 96 hours
# HCT116 KO plus UBE2H overexpression single cell clones

# Note: HCT116 KO1 _ OE clone 1 is a tetraploid cell. Much larger and twice as much DNA as the others. not odd that it also has more mitotracker. 

colnames(HCT116.817)[1:20]
colnames(HCT116.825)[1:20]

columnnames<- c("Tube.Name.", "Sample.ID.", "All.Events" , "Cell.size.Events", "Hoechst.Events", 
             "Type", "RelativeCount", "Doubling", "RelativeMito", "Hoechst.Mean.PE.A", "AZMito", 
             "Hoechst.Mean.FSC.A", "Hoechst.Mean.Violet610.A", "Mito.div.Hoechst", "Mito.div.FSC", "RelativeHoechst")

HCT116.KOOEc<- rbind(HCT116.817[, columnnames], HCT116.825[, columnnames])


ggplot(HCT116.KOOEc, aes(x=Sample.ID., y=Doubling, color=Type ))+
  geom_point()+
  scale_color_manual(values=c("cyan3","cyan3", "cyan4", "cyan4", "purple1","purple1", "purple4",  "purple4"), guide="none")+
  stat_summary(fun.y=mean, geom="point", shape="-", size=10, color="red", fill="red") +
  ylab("Doubling per day")+
  xlab("")+  
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
# plot.HCT116.OEc.Doubling_Merge.pdf


## plot mean mitoTracker (PE) 
ggplot(HCT116.KOOEc, aes(x=Sample.ID., y=RelativeMito, color=Type ))+
  geom_point()+
  scale_color_manual(values=c("cyan3","cyan3", "cyan4", "cyan4", "purple1","purple1", "purple4",  "purple4"), guide="none")+
  ylab("Mean Mitotracker")+
  geom_hline(yintercept = mean(subset(HCT116.KOOEc, Sample.ID. == "WT")$RelativeMito) ) + # average WT mitotracker expression
  #geom_hline(yintercept = mean(subset(HCT116.KOOEc, Sample.ID. == "WT_AZ")$Hoechst.Mean.PE.A) ) + # average WT mitotracker expression
  stat_summary(fun.y=mean, geom="point", shape="-", size=10, color="red", fill="red") +
  xlab("")+  
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
#plot.HCT116.OEc.Mito_Merge.pdf


## No AZ
ggplot(subset(HCT116.KOOEc, Sample.ID. %in% c("WT", "OE", "KO1", "KO2", "KO3", "KO1.OEc1", "KO2.OEc1", "KO2.OEc2", "KO3.OEc2", "KO3.OEc5")), 
       aes(x=Sample.ID., y=RelativeMito, color=Type ))+
  geom_point()+
  scale_color_manual(values=c("cyan3", "cyan4",  "purple1", "purple4"), guide="none")+
  ylab("Mean Mitotracker")+
  geom_hline(yintercept = mean(subset(HCT116.KOOEc, Sample.ID. == "WT")$RelativeMito) ) + # average WT mitotracker expression
  #geom_hline(yintercept = mean(subset(HCT116.KOOEc, Sample.ID. == "WT_AZ")$Hoechst.Mean.PE.A) ) + # average WT mitotracker expression
  stat_summary(fun.y=mean, geom="point", shape="-", size=10, color="red", fill="red") +
  xlab("")+  
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
#plot.HCT116.OEc.Mito_Merge.pdf




# change in mitotracker upon AZ3146
ggplot(subset(HCT116.KOOEc, Sample.ID. %in% c("WT.AZ", "OE.AZ", "KO1.AZ", "KO2.AZ", "KO3.AZ", "KO1.OEc1.AZ", "KO2.OEc1.AZ","KO2.OEc2.AZ","KO3.OEc2.AZ","KO3.OEc5.AZ")), 
       aes(x=Sample.ID., y=AZMito, color=Type ))+
  geom_point()+
  scale_color_manual(values=c("cyan3", "cyan4", "purple1",  "purple4"), guide="none")+
  ylab("Change in Mitotracker\nupon AZ3146")+
  geom_hline(yintercept = mean(subset(HCT116.KOOEc, Sample.ID. == "WT")$AZMito) ) + # average WT mitotracker expression
  stat_summary(fun.y=mean, geom="point", shape="-", size=10, color="red", fill="red") +
  xlab("")+  
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
# plot.HCT116.OEc.MitoAZ_Merge.pdf
## Figure 6D


## plot mean  cell size
ggplot(HCT116.KOOEc, aes(x=Sample.ID., y=Hoechst.Mean.FSC.A, color=Type ))+
  geom_point()+
  scale_color_manual(values=c("cyan3","cyan3", "cyan4", "cyan4", "purple1","purple1", "purple4",  "purple4"), guide="none")+
  stat_summary(fun.y=mean, geom="point", shape="-", size=10, color="red", fill="red") +
  geom_hline(yintercept = mean(subset(HCT116.KOOEc, Sample.ID. == "WT")$Hoechst.Mean.FSC.A) ) + # average WT mitotracker expression
  ylab("Mean Size")+
  xlab("")+  
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
# UBE2H-KO reduces cell size (KO+OE rescues cell size), OE increases cell size. 
#plot.HCT116.OEc.Size_Merge.pdf


## plot relative  Hoechst
ggplot(HCT116.KOOEc, aes(x=Sample.ID., y=RelativeHoechst , color=Type))+
  geom_point()+
  scale_color_manual(values=c("cyan3","cyan3", "cyan4", "cyan4", "purple1","purple1", "purple4",  "purple4"), guide="none")+
  geom_hline(yintercept = mean(subset(HCT116.KOOEc, Sample.ID. == "WT")$RelativeHoechst) ) + # average WT mitotracker expression
  stat_summary(fun.y=mean, geom="point", shape="-", size=10, color="red", fill="red") +
  ylab("Relative Hoechst Fluoescence")+
  xlab("")+ 
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
# HCT116.Hoechst.KOOEc_relative_Merge.pdf
# UBE2H-KO + AZ significantly increase DNA/Hoechst content





### pvalues
## MitoTracker:  AZ increase significant? 
t.test(subset(HCT116.KOOEc, Sample.ID. == "WT")$Hoechst.Mean.PE.A, subset(HCT116.KOOEc, Sample.ID. == "WT.AZ")$Hoechst.Mean.PE.A)
#WT vs AZ: 0.0002344 ***
t.test(subset(HCT116.KOOEc, Sample.ID. == "OE")$Hoechst.Mean.PE.A, subset(HCT116.KOOEc, Sample.ID. == "OE.AZ")$Hoechst.Mean.PE.A)
#OE  vs AZ: 0.599 NS
t.test(subset(HCT116.KOOEc, Sample.ID. == "KO1")$Hoechst.Mean.PE.A, subset(HCT116.KOOEc, Sample.ID. == "KO1.AZ")$Hoechst.Mean.PE.A)
#KO1  vs AZ: 1.594e-07 ***
t.test(subset(HCT116.KOOEc, Sample.ID. == "KO2")$Hoechst.Mean.PE.A, subset(HCT116.KOOEc, Sample.ID. == "KO2.AZ")$Hoechst.Mean.PE.A)
#KO2 vs AZ:  0.0001338 ***
t.test(subset(HCT116.KOOEc, Sample.ID. == "KO3")$Hoechst.Mean.PE.A, subset(HCT116.KOOEc, Sample.ID. == "KO3.AZ")$Hoechst.Mean.PE.A)
#KO3 vs AZ: 0.00213 **
t.test(subset(HCT116.KOOEc, Sample.ID. == "KO1.OEc1")$Hoechst.Mean.PE.A, subset(HCT116.KOOEc, Sample.ID. == "KO1.OEc1.AZ")$Hoechst.Mean.PE.A)
#KO1OE vs AZ: 0.9123 NS
t.test(subset(HCT116.KOOEc, Sample.ID. == "KO2.OEc1")$Hoechst.Mean.PE.A, subset(HCT116.KOOEc, Sample.ID. == "KO2.OEc1.AZ")$Hoechst.Mean.PE.A)
#KO2OE vs AZ: 0.1503 NS
t.test(subset(HCT116.KOOEc, Sample.ID. == "KO2.OEc2")$Hoechst.Mean.PE.A, subset(HCT116.KOOEc, Sample.ID. == "KO2.OEc2.AZ")$Hoechst.Mean.PE.A)
#KO2OE vs AZ: 0.7269 NS
t.test(subset(HCT116.KOOEc, Sample.ID. == "KO3.OEc2")$Hoechst.Mean.PE.A, subset(HCT116.KOOEc, Sample.ID. == "KO3.OEc2.AZ")$Hoechst.Mean.PE.A)
#KO3OE vs AZ: 0.779 NS
t.test(subset(HCT116.KOOEc, Sample.ID. == "KO3.OEc5")$Hoechst.Mean.PE.A, subset(HCT116.KOOEc, Sample.ID. == "KO3.OEc5.AZ")$Hoechst.Mean.PE.A)
#KO3OE vs AZ: 0.5043 NS



## MitoTracker:  AZ mito treated greater than WT+AZ? 
t.test(subset(HCT116.KOOEc, Sample.ID. == "OE.AZ")$Hoechst.Mean.PE.A, subset(HCT116.KOOEc, Sample.ID. == "WT.AZ")$Hoechst.Mean.PE.A)
#OE AZ  vs WTAZ: 0.109 NS
t.test(subset(HCT116.KOOEc, Sample.ID. == "KO1.AZ")$Hoechst.Mean.PE.A, subset(HCT116.KOOEc, Sample.ID. == "WT.AZ")$Hoechst.Mean.PE.A)
#KO1AZ  vs WTAZ: 3.431e-06 ***
t.test(subset(HCT116.KOOEc, Sample.ID. == "KO2.AZ")$Hoechst.Mean.PE.A, subset(HCT116.KOOEc, Sample.ID. == "WT.AZ")$Hoechst.Mean.PE.A)
#KO2AZ vs WTAZ:  3.057e-05 ***
t.test(subset(HCT116.KOOEc, Sample.ID. == "KO3.AZ")$Hoechst.Mean.PE.A, subset(HCT116.KOOEc, Sample.ID. == "WT.AZ")$Hoechst.Mean.PE.A)
#KO3AZ vs WTAZ: 0.007182. *
t.test(subset(HCT116.KOOEc, Sample.ID. == "KO1.OEc1.AZ")$Hoechst.Mean.PE.A, subset(HCT116.KOOEc, Sample.ID. == "WT.AZ")$Hoechst.Mean.PE.A)
#KO1OEAZ vs WTAZ: 0.1451 NS
t.test(subset(HCT116.KOOEc, Sample.ID. == "KO2.OEc1.AZ")$Hoechst.Mean.PE.A, subset(HCT116.KOOEc, Sample.ID. == "WT.AZ")$Hoechst.Mean.PE.A)
#KO2OEAZ vs WTAZ: 0.958 NS
t.test(subset(HCT116.KOOEc, Sample.ID. == "KO2.OEc2.AZ")$Hoechst.Mean.PE.A, subset(HCT116.KOOEc, Sample.ID. == "WT.AZ")$Hoechst.Mean.PE.A)
#KO2OEAZ vs WTAZ: 0.2702 NS
t.test(subset(HCT116.KOOEc, Sample.ID. == "KO3.OEc2.AZ")$Hoechst.Mean.PE.A, subset(HCT116.KOOEc, Sample.ID. == "WT.AZ")$Hoechst.Mean.PE.A)
#KO3OEAZ vs WTAZ: 0.01651. *
t.test(subset(HCT116.KOOEc, Sample.ID. == "KO3.OEc5.AZ")$Hoechst.Mean.PE.A, subset(HCT116.KOOEc, Sample.ID. == "WT.AZ")$Hoechst.Mean.PE.A)
#KO3OEAZ vs WTAZ: 0.0718 NS


## MitoTracker:  KO+AZ vs KO+OE+AZ? 
t.test(subset(HCT116.KOOEc, Sample.ID. == "KO1.AZ")$Hoechst.Mean.PE.A, subset(HCT116.KOOEc, Sample.ID. == "KO1.OEc1.AZ")$Hoechst.Mean.PE.A)
#KO1 AZ  vs KO1.oe.AZ: 0.08978 NS
t.test(subset(HCT116.KOOEc, Sample.ID. == "KO2.AZ")$Hoechst.Mean.PE.A, subset(HCT116.KOOEc, Sample.ID. == "KO2.OEc1.AZ")$Hoechst.Mean.PE.A)
#KO2 AZ  vs KO2.oe.AZ: 5.06e-06
t.test(subset(HCT116.KOOEc, Sample.ID. == "KO2.AZ")$Hoechst.Mean.PE.A, subset(HCT116.KOOEc, Sample.ID. == "KO2.OEc2.AZ")$Hoechst.Mean.PE.A)
#KO2 AZ  vs KO2.oe.AZ: 1.024e-06
t.test(subset(HCT116.KOOEc, Sample.ID. == "KO3.AZ")$Hoechst.Mean.PE.A, subset(HCT116.KOOEc, Sample.ID. == "KO3.OEc2.AZ")$Hoechst.Mean.PE.A)
#KO3 AZ  vs KO3.oe.AZ: 0.006899
t.test(subset(HCT116.KOOEc, Sample.ID. == "KO3.AZ")$Hoechst.Mean.PE.A, subset(HCT116.KOOEc, Sample.ID. == "KO3.OEc5.AZ")$Hoechst.Mean.PE.A)
#KO3 AZ  vs KO3.oe.AZ: 0.01996



### Change in MitoTracker upon AZ3146
## MitoTracker:  AZ mito treated greater than WT+AZ? 
t.test(subset(HCT116.KOOEc, Sample.ID. == "OE.AZ")$AZMito, 
       subset(HCT116.KOOEc, Sample.ID. == "WT.AZ")$AZMito)
#OE AZ  vs WTAZ: 0.0003949
t.test(subset(HCT116.KOOEc, Sample.ID. == "KO1.AZ")$AZMito, 
       subset(HCT116.KOOEc, Sample.ID. == "WT.AZ")$AZMito)
#KO1AZ  vs WTAZ: 0.0009943
t.test(subset(HCT116.KOOEc, Sample.ID. == "KO2.AZ")$AZMito, 
       subset(HCT116.KOOEc, Sample.ID. == "WT.AZ")$AZMito)
#KO2AZ vs WTAZ:  0.001144
t.test(subset(HCT116.KOOEc, Sample.ID. == "KO3.AZ")$AZMito, 
       subset(HCT116.KOOEc, Sample.ID. == "WT.AZ")$AZMito)
#KO3AZ vs WTAZ: 0.1689

t.test(subset(HCT116.KOOEc, Sample.ID. == "KO1.OEc1.AZ")$AZMito, 
       subset(HCT116.KOOEc, Sample.ID. == "WT.AZ")$AZMito)
#KO1OEAZ vs WTAZ: 4.279e-08
t.test(subset(HCT116.KOOEc, Sample.ID. == "KO2.OEc1.AZ")$AZMito, 
       subset(HCT116.KOOEc, Sample.ID. == "WT.AZ")$AZMito)
#KO2OEAZ vs WTAZ: 0.3947
t.test(subset(HCT116.KOOEc, Sample.ID. == "KO2.OEc2.AZ")$AZMito, 
       subset(HCT116.KOOEc, Sample.ID. == "WT.AZ")$AZMito)
#KO2OEAZ vs WTAZ: 5.384e-06
t.test(subset(HCT116.KOOEc, Sample.ID. == "KO3.OEc2.AZ")$AZMito, 
       subset(HCT116.KOOEc, Sample.ID. == "WT.AZ")$AZMito)
#KO3OEAZ vs WTAZ: 0.0003307
t.test(subset(HCT116.KOOEc, Sample.ID. == "KO3.OEc5.AZ")$AZMito, 
       subset(HCT116.KOOEc, Sample.ID. == "WT.AZ")$AZMito)
#KO3OEAZ vs WTAZ: 0.1975


## MitoTracker:  KO+AZ vs KO+OE+AZ? 
t.test(subset(HCT116.KOOEc, Sample.ID. == "KO1.AZ")$AZMito, 
       subset(HCT116.KOOEc, Sample.ID. == "KO1.OEc1.AZ")$AZMito)
#KO1 AZ  vs KO1.oe.AZ: 0.0001692
t.test(subset(HCT116.KOOEc, Sample.ID. == "KO2.AZ")$AZMito, 
       subset(HCT116.KOOEc, Sample.ID. == "KO2.OEc1.AZ")$AZMito)
#KO2 AZ  vs KO2.oe.AZ: 0.0004316
t.test(subset(HCT116.KOOEc, Sample.ID. == "KO2.AZ")$AZMito, 
       subset(HCT116.KOOEc, Sample.ID. == "KO2.OEc2.AZ")$AZMito)
#KO2 AZ  vs KO2.oe.AZ: 8.009e-07
t.test(subset(HCT116.KOOEc, Sample.ID. == "KO3.AZ")$AZMito, 
       subset(HCT116.KOOEc, Sample.ID. == "KO3.OEc2.AZ")$AZMito)
#KO3 AZ  vs KO3.oe.AZ: 0.01977
t.test(subset(HCT116.KOOEc, Sample.ID. == "KO3.AZ")$AZMito, 
       subset(HCT116.KOOEc, Sample.ID. == "KO3.OEc5.AZ")$AZMito)
#KO3 AZ  vs KO3.oe.AZ: 0.06234



### Cal51  June 3:  #### 
### Cal51 cells treated with 0.25uM AZ3146 for 96 hours 


colnames(Cal51.603)


# Assign groups: 
# Remove the last threee characters from "Tube.Name." (..._01, ..._02, ..._03, etc.) and 
# and remove the first three characters from "Tube.Name." (01_...) 
Cal51.603$Sample.ID.<- str_sub(Cal51.603$Tube.Name., start=4, end = -4)
Cal51.603$Sample.ID.<- factor(Cal51.603$Sample.ID., levels= c("WT", "WT_AZ", 
                                                              "KO2", "KO2_AZ",  "KO7",  "KO7_AZ"
))

Cal51.603$Type<- c("WT", "WT", "WT", "WT", "WT",  "WT+AZ", "WT+AZ", "WT+AZ",  "WT+AZ", "WT+AZ", "WT+AZ",
                   "KO", "KO", "KO", "KO+AZ", "KO+AZ", "KO+AZ",
                   "KO", "KO", "KO", "KO+AZ", "KO+AZ", "KO+AZ", "WT")
Cal51.603$Type<- factor(Cal51.603$Type, levels= c("WT", "WT+AZ", "KO", "KO+AZ"))


# Cell growth (96 hours)
## Note that 30ul from a total of 200ul is measured. so multiuple # cells by 200/30
## Count number of cells: cell sized & Hoechst positive. 
MeanWTCount<- mean(subset(Cal51.603, Sample.ID. == "WT")$Hoechst.Events)
Cal51.603$RelativeCount<- Cal51.603$Hoechst.Events/MeanWTCount
Cal51.603$Doubling<- Cal51.603$Hoechst.Events*200/(30*100000*4) # doubling: final/initial 
# Final is number of Hoechst + cells *200/30 (only analyzed 30ul from 200ul total), started with 100k cells, grew for 4 days so divide by 4 to get doubling per day

ggplot(Cal51.603, aes(x=Cal51.603$Sample.ID., y=RelativeCount, color= Type ))+
  geom_point()+
  scale_color_manual(values=c("cyan3", "cyan3", "purple1",  "purple1"), guide="none")+
  #scale_fill_manual(values=c("cyan3", "cyan3",  "purple1", "purple1"))+
  stat_summary(fun.y=mean, geom="point", shape="-", size=10, color="red", fill="red") +
  ylab("Relative Cell count")+
  xlab("")+  
  theme_classic()


ggplot(Cal51.603, aes(x=Cal51.603$Sample.ID., y=Doubling, color= Type ))+
  geom_point()+
  scale_color_manual(values=c("cyan3", "cyan3", "purple1",  "purple1"), guide="none")+
  #scale_fill_manual(values=c("cyan3", "cyan3",  "purple1", "purple1"))+
  stat_summary(fun.y=mean, geom="point", shape="-", size=10, color="red", fill="red") +
  ylab("Doubling per day")+
  xlab("")+  
  theme_classic()
#plot.Cal51.AZ.Doubling.250_3.pdf


####
## plot mean mitoTracker (PE) 

ggplot(Cal51.603, aes(x=Cal51.603$Sample.ID., y=Cal51.603$Hoechst.Mean.PE.A, color=Type ))+
  geom_point()+
  scale_color_manual(values=c("cyan3", "cyan3", "purple1",  "purple1"), guide="none")+
  #scale_fill_manual(values=c("cyan3", "cyan3",  "purple1", "purple1"))+
  geom_hline(yintercept = mean(subset(Cal51.603, Sample.ID. == "WT")$Hoechst.Mean.PE.A) ) + # average WT mitotracker expression
  #geom_hline(yintercept = mean(subset(Cal51.603, Sample.ID. == "WT_AZ")$Hoechst.Mean.PE.A) ) + # average WT mitotracker expression
  stat_summary(fun.y=mean, geom="point", shape="-", size=10, color="red", fill="red") +
  ylab("Mean Mitotracker")+
  xlab("")+  
  theme_classic()
# KO and AZ increase mitotracker. KO increases it more than AZ
# plot.Cal51.AZ.Mitotracker.250_3.pdf

Cal51.603$Change.Mito<- c(Cal51.603$Hoechst.Mean.PE.A[c(1:5)]/mean(Cal51.603$Hoechst.Mean.PE.A[c(1:5,24)]), 
                          Cal51.603$Hoechst.Mean.PE.A[6:11]/mean(Cal51.603$Hoechst.Mean.PE.A[c(1:5,24)]), 
                          Cal51.603$Hoechst.Mean.PE.A[12:14]/mean(Cal51.603$Hoechst.Mean.PE.A[12:14]), 
                          Cal51.603$Hoechst.Mean.PE.A[15:17]/mean(Cal51.603$Hoechst.Mean.PE.A[12:14]), 
                          Cal51.603$Hoechst.Mean.PE.A[18:20]/mean(Cal51.603$Hoechst.Mean.PE.A[18:20]), 
                          Cal51.603$Hoechst.Mean.PE.A[21:23]/mean(Cal51.603$Hoechst.Mean.PE.A[18:20]), 
                          Cal51.603$Hoechst.Mean.PE.A[c(24)]/mean(Cal51.603$Hoechst.Mean.PE.A[c(1:5,24)])) 

ggplot(subset(Cal51.603, Type %in% c("WT+AZ", "KO+AZ")), 
       aes(x=Sample.ID., y=Change.Mito, color=Type ))+
  geom_point()+
  scale_color_manual(values=c("cyan3", "purple1",  "purple1"), guide="none")+
  #scale_fill_manual(values=c("cyan3", "purple1", "purple1"))+
  geom_hline(yintercept = mean(subset(Cal51.603, Sample.ID. == "WT")$Change.Mito) ) + # average WT mitotracker expression
  #geom_hline(yintercept = mean(subset(Cal51.603, Sample.ID. == "WT_AZ")$Hoechst.Mean.PE.A) ) + # average WT mitotracker expression
  stat_summary(fun.y=mean, geom="point", shape="-", size=10, color="red", fill="red") +
  ylab("Change in Mitotracker") +
  xlab("") +  
  theme_classic()
# plot.Cal51.AZ.Mitotracker.250_ChangeinMito.pdf
# Figure 6E

## plot mean cell size
ggplot(Cal51.603, aes(x=Cal51.603$Sample.ID., y=Cal51.603$Hoechst.Mean.FSC.A, fill=Type ))+
  geom_boxplot()+
  geom_point()+
  scale_fill_manual(values=c("cyan3", "cyan3",  "purple1", "purple1"))+
  geom_hline(yintercept = mean(subset(Cal51.603, Sample.ID. == "WT")$Hoechst.Mean.FSC.A) ) + # average WT mitotracker expression
  #stat_summary(fun.y=mean, geom="point", shape="-", size=10, color="red", fill="red") +
  ylab("Mean Size")+
  xlab("")+  
  theme_classic()
# UBE2H-KO increases cell size. AZ inceases it? 
# plot.Cal51.AZ.Size.250_2.pdf


## plot mean  Hoechst
Cal51.603$relative.Hoechst<- Cal51.603$Hoechst.Mean.Violet610.A/mean(subset(Cal51.603, Sample.ID. == "WT")$Hoechst.Mean.Violet610.A)
ggplot(Cal51.603, aes(x=Cal51.603$Sample.ID., y=Cal51.603$relative.Hoechst, fill=Type ))+
  geom_boxplot()+
  geom_point()+
  scale_fill_manual(values=c("cyan3", "cyan3",  "purple1", "purple1"))+
  geom_hline(yintercept = mean(subset(Cal51.603, Sample.ID. == "WT")$relative.Hoechst) ) + # average WT mitotracker expression
  #stat_summary(fun.y=mean, geom="point", shape="-", size=10, color="red", fill="red") +
  ylab("Relative  Hoechst  Fluoescence")+
  xlab("")+  
  theme_classic()
# UBE2H-KO + AZ significantly increase DNA/Hoechst content
# plot.Cal51.AZ.Hoechst.250.pdf

### MitoTracker per DNA Hoechst 
Cal51.603$Mito.div.Hoechst<- Cal51.603$Hoechst.Mean.PE.A/ Cal51.603$Hoechst.Mean.Violet610.A
ggplot(Cal51.603, aes(x=Cal51.603$Sample.ID., y=Cal51.603$Mito.div.Hoechst, fill=Type ))+
  geom_boxplot()+
  geom_point()+
  scale_fill_manual(values=c("cyan3", "cyan3",  "purple1", "purple1"))+
  geom_hline(yintercept = mean(subset(Cal51.603, Sample.ID. == "WT")$Mito.div.Hoechst) ) + # average WT mitotracker expression
  #stat_summary(fun.y=mean, geom="point", shape="-", size=10, color="red", fill="red") +
  ylab("MitoTracker/DNA (Hoechst)")+
  xlab("")+  
  theme_classic()
# UBE2H-KO Increases Mitotracker/DNA.  AZ3146 decreases it? 


### MitoTracker per cell size 
Cal51.603$Mito.div.FSC<- Cal51.603$Hoechst.Mean.PE.A/ Cal51.603$Hoechst.Mean.FSC.A
ggplot(Cal51.603, aes(x=Sample.ID., y=Mito.div.FSC, fill=Type ))+
  geom_boxplot()+
  geom_point()+
  scale_fill_manual(values=c("cyan3", "cyan3",  "purple1", "purple1"))+
  geom_hline(yintercept = mean(subset(Cal51.603, Sample.ID. == "WT")$Mito.div.FSC) ) + # average WT mitotracker expression
  #stat_summary(fun.y=mean, geom="point", shape="-", size=10, color="red", fill="red") +
  ylab("MitoTracker/Cell Size")+
  xlab("")+  
  theme_classic()
# UBE2H-KO Increases Mitotracker/Size.  AZ3146 increases it. 
# plot.Cal51.AZ.MitoCellSize.250.pdf



### pvalues
## MitoTracker:  AZ increase significant? 
t.test(subset(Cal51.603, Sample.ID. == "WT")$Hoechst.Mean.PE.A, 
       subset(Cal51.603, Sample.ID. == "WT_AZ")$Hoechst.Mean.PE.A)
#WT vs AZ: 0.01114
t.test(subset(Cal51.603, Sample.ID. == "KO2")$Hoechst.Mean.PE.A, 
       subset(Cal51.603, Sample.ID. == "KO2_AZ")$Hoechst.Mean.PE.A)
#KO2  vs KO2AZ: 0.000558
t.test(subset(Cal51.603, Sample.ID. == "KO7")$Hoechst.Mean.PE.A, 
       subset(Cal51.603, Sample.ID. == "KO7_AZ")$Hoechst.Mean.PE.A)
#KO7  vs KO7AZ: 0.0003054

## MitoTracker:  AZ mito treated greater than WT+AZ? 
t.test(subset(Cal51.603, Sample.ID. == "WT_AZ")$Hoechst.Mean.PE.A, 
       subset(Cal51.603, Sample.ID. == "KO2_AZ")$Hoechst.Mean.PE.A)
#WTAZ vs KO2AZ: 5.904e-05
t.test(subset(Cal51.603, Sample.ID. == "WT_AZ")$Hoechst.Mean.PE.A, 
       subset(Cal51.603, Sample.ID. == "KO7_AZ")$Hoechst.Mean.PE.A)
#WTAZ vs KO7AZ: 0.002922



## CHANGE in MITOTRACKER: 
t.test(subset(Cal51.603, Sample.ID. == "WT_AZ")$Change.Mito, 
       subset(Cal51.603, Sample.ID. == "KO2_AZ")$Change.Mito)
#WTAZ vs KO2AZ: 0.4553
t.test(subset(Cal51.603, Sample.ID. == "WT_AZ")$Change.Mito, 
       subset(Cal51.603, Sample.ID. == "KO7_AZ")$Change.Mito)
#WTAZ vs KO7AZ: 0.005739


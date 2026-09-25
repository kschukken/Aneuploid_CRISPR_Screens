### Plot colony formation, cell death, Rotenone drug curve, and Seahorse ####
## Author: Klaske M. Schukken 
## Date: September 11, 2025

# This file to to take colony formation numbers and plot them. 
# to see if UBE2H-KO clones are more sensitive to drugs than WT clones 


### INTRO AND LIBRARIES ####
# Folders
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
library(plyr)
library('gprofiler2')
library(xlsx)
library("viridis")  # color packet
library(drc) # Drug curve library


### Rotenone Mitochondrial drug analysis ####
setwd(SchukkenData)

RotentoneIC50 <- read_excel("Rotenone_IC50.xlsx", sheet="Sheet4")

RotentoneIC50<- subset(RotentoneIC50, ! `Concentration (nM)`==0)

# five runs: 
# run 1: concentrations too high, all cells dead
# run 2: Good
# run 3: ok. only 3 samples, not 5
# run 4: Good. clean replicates
# run 5: ok. but SNU1 C114 had outlier points. High cell number at high concentration?? not seen in other replicates. 



# Load required libraries



# fit_multi <- drm(response ~ dose, curveid = cell_line, data = df, fct = LL.4())

RotentoneRun4<- subset(RotentoneIC50, Run ==4)[,1:3]
colnames(RotentoneRun4)<- c("dose", "Response", "CellLine")
RotentoneRun4$CellLine<- factor(RotentoneRun4$CellLine, level = c("SNU1 C13", "SNU1 C12", "SNU1 C24", "SNU1 C111", "SNU1 C114"))
RotentoneRun4$condition<-paste0(RotentoneRun4$CellLine,"_", RotentoneRun4$dose)
RotentoneRun4<- data.frame(RotentoneRun4)

# 1. Calculate the Mean and SEM for each unique Dose and Cell Line condition
summary_df<- data.frame()
for (i in unique(RotentoneRun4$condition) ){
  a<- subset(RotentoneRun4, condition == i)
  average = mean(a$Response)
  SEM  = sd(a$Response) # SEM = SD/srt(n)
  
  summary_df<- rbind(summary_df, data.frame(
    CellLine= a$CellLine, 
    dose= a$dose, 
    mean_percent= average, 
    sem_percent=SEM, 
    lower_sem= average -SEM, 
    upper_sem= average + SEM
  ))
}

# 2. Fit the multi-group log-logistic model
fit_multi <- drm(
  Response ~ dose, 
  curveid = CellLine, 
  data = RotentoneRun4, 
  fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50"))
)


# 3. Create the smooth layout grid for predictions
smooth_doses <- exp(seq(log(min(RotentoneRun4$dose)), log(max(RotentoneRun4$dose)), length.out = 200))
grid_df <- expand.grid(dose = smooth_doses, CellLine = levels(RotentoneRun4$CellLine) )

# 4. Predict smooth curve lines (we no longer need 'interval = "confidence"')
predictions <- predict(fit_multi, newdata = grid_df)

# Combine the grid layout with the predicted curve points
predicted_df <- cbind(grid_df, Response = predictions)


# 5. Build the ggplot
ggplot() +
  # A. Overlay the smooth regression line calculated by drc
  geom_line(data = predicted_df, 
            aes(x = dose, y = Response, color = CellLine), 
            size = 1.1) +
  
  # B. Add SD Error Bars around those means
  # width = 0.05 adds small horizontal caps to the top and bottom of the error bars
  geom_errorbar(data = summary_df, 
                aes(x = dose, ymin = lower_sem, ymax = upper_sem, color = CellLine), 
                width = 0.05, size = 0.8, color="grey30") +
  
  # C. Plot Mean Data Points (instead of raw jittered replicates)
  geom_point(data = summary_df, 
             aes(x = dose, y = mean_percent, color = CellLine), 
             size = 3) +
  #ylim(0,1.2)+
  scale_x_log10() + 
  labs(
    x = "Rotenone [nM]",
    y = "Survival (%)",
    color = "Cell Lines"
  ) +
  theme_classic()
# plot.Rotenone.SNU1_Run4.pdf
# Supplementary Figure 3A


## To find the significance between cell line curves: ##
# Compare your 4 aneuploid clones back to the 1 control clones
# The 'percVec = 50' means we are comparing the EC50 values
# The 'adjustment = "dunnett"' corrects the P-values for making multiple comparisons

# Dunnet adjustment for multiple tests: 
EDcomp(fit_multi, percVec = c(50, 50), adjustment = "dunnett")

# Run 4 p-values: 
#                             Estimate  Std. Error     t-value     p-value
#SNU1 C111/SNU1 C13:50/50   2.4808e-01  1.4722e-01 -5.1074e+00  1.9742e-06
#SNU1 C114/SNU1 C13:50/50   5.6128e-01  1.9276e-01 -2.2760e+00  2.5358e-02
#SNU1 C12/SNU1 C13:50/50    6.0781e-01  2.4241e-01 -1.6179e+00  1.0939e-01
#SNU1 C13/SNU1 C24:50/50    3.3451e+00  1.4604e+00  1.6058e+00  1.1202e-01



#IC50: 
summary(fit_multi)


# Summary IC50: 
# Run#: 2,   3,   4,   5
# C13: 323, 355, 379, 315 (ALL SIMILAR, ~350nM) 
# C12: 185, 67,  230, 281 (All except #3 ~200nM)
# C24: 3,    2,  113, 0.1 (All except #4 ~2nM)
#C111: 164,   ,   94, 218 (All except #4 ~200nM)
#C114: 257,   ,  212, 654 (All except #5 ~200nM)

# Run 2 and 4 are the most consistent. 


#                        Estimate Std. Error  t-value   p-value    
# run2: 
#ED50:SNU1 C12         185.882100  32.651071  5.6930 1.742e-07 ***
#ED50:SNU1 C13         323.592553 347.859469  0.9302 0.3548816    
#ED50:SNU1 C24           2.654970  26.410636  0.1005 0.9201630    
#ED50:SNU1 C111        164.708299  32.482364  5.0707 2.290e-06 ***
#ED50:SNU1 C114        257.393458  34.759818  7.4049 8.785e-11 ***

# Run 3
#ED50:SNU1 C12         67.802480 110.451633  0.6139 0.5420326    
#ED50:SNU1 C13        355.604591  99.672554  3.5677 0.0007943 ***
#ED50:SNU1 C24          2.050068  13.575089  0.1510 0.8805581    

# Run 4
#ED50:SNU1 C12         230.536904  90.866653  2.5371 0.0130039 *  
#ED50:SNU1 C13         379.289521  23.079561 16.4340 < 2.2e-16 ***
#ED50:SNU1 C24         113.386209  49.017466  2.3132 0.0231282 *  
#ED50:SNU1 C111         94.094908  55.544767  1.6940 0.0939196 .  
#ED50:SNU1 C114        212.887058  71.953875  2.9587 0.0040014 ** 

# run 5
# ED50:SNU1 C12         281.242403  64.790465  4.3408 3.897e-05 ***
# ED50:SNU1 C13         315.399369  58.758977  5.3677 6.806e-07 ***
# ED50:SNU1 C24           0.109989   0.202957  0.5419   0.58928    
# ED50:SNU1 C111        218.219993 158.022745  1.3809   0.17092    
# ED50:SNU1 C114        652.872755 106.411942  6.1353 2.600e-08 *** 


# Dose-response curves were analyzed using nonlinear regression with a four-parameter log-logistic model 
# using the drc package in R. Pairwise comparisons of the IC50 values were performed using the EDcomp function, and 
# p-values were corrected for multiple tests using Dunnett’s test.  



### Get Colony formation Data ####
# Scaned plates for scientific reproducibility, but manually counted number of colonies per well. Using Fijito visualize and count. 
# 6 well plates

Cal51.UBE2H.AZ<- data.frame( CellLine = c("WT", "WT","WT","WT","WT","WT",
                                          "UBE2H-KO_c2", "UBE2H-KO_c2",  "UBE2H-KO_c2",  "UBE2H-KO_c2",  "UBE2H-KO_c2",  "UBE2H-KO_c2", 
                                          "UBE2H-KO_c7", "UBE2H-KO_c7", "UBE2H-KO_c7", "UBE2H-KO_c7", "UBE2H-KO_c7", "UBE2H-KO_c7"), 
                             AZ3146_250nM = c("-","-","-","+","+","+", 
                                              "-","-","-","+","+","+", 
                                              "-","-","-","+","+","+"), 
                             Colony = c(61,75,63,34,26,32, 
                                        31,41,43,19,17,10, 
                                        67,40,53,21,9,9))

Cal51.UBE2H.AZ2<- data.frame( CellLine = c("WT", "WT","WT","WT","WT","WT",
                                           "UBE2H-KO_c2", "UBE2H-KO_c2",  "UBE2H-KO_c2",  "UBE2H-KO_c2",  "UBE2H-KO_c2",  "UBE2H-KO_c2", 
                                           "UBE2H-KO_c7", "UBE2H-KO_c7", "UBE2H-KO_c7", "UBE2H-KO_c7", "UBE2H-KO_c7", "UBE2H-KO_c7"), 
                              AZ3146_250nM = c("-","-","-","+","+","+", 
                                               "-","-","-","+","+","+", 
                                               "-","-","-","+","+","+"), 
                              Colony = c(86,89,75,61,50,56,
                                         46,34,30,5,20,12,
                                         38,42,50,24,27,19))

Cal51.UBE2H.AZ3<- data.frame( CellLine = c("WT", "WT","WT","WT","WT","WT",
                                           "UBE2H-KO_c2", "UBE2H-KO_c2",  "UBE2H-KO_c2",  "UBE2H-KO_c2",  "UBE2H-KO_c2",  "UBE2H-KO_c2", 
                                           "UBE2H-KO_c7", "UBE2H-KO_c7", "UBE2H-KO_c7", "UBE2H-KO_c7", "UBE2H-KO_c7", "UBE2H-KO_c7"), 
                              AZ3146_250nM = c("-","-","-","+","+","+", 
                                               "-","-","-","+","+","+", 
                                               "-","-","-","+","+","+"), 
                              Colony = c(145,146,156, 129,139,157,
                                         15,13,21,7,11,7,
                                         93,108,116, 48,67,64))

HCT116.UBE2H.AZ<- data.frame( CellLine = c("WT","WT","WT",  "WT","WT","WT",
                                           "UBE2H-KO_c1", "UBE2H-KO_c1", "UBE2H-KO_c1", "UBE2H-KO_c1", "UBE2H-KO_c1", "UBE2H-KO_c1",
                                           "UBE2H-KO_c2", "UBE2H-KO_c2", "UBE2H-KO_c2", "UBE2H-KO_c2", "UBE2H-KO_c2", "UBE2H-KO_c2",
                                           "UBE2H-KO_c3", "UBE2H-KO_c3", "UBE2H-KO_c3", "UBE2H-KO_c3", "UBE2H-KO_c3", "UBE2H-KO_c3"), 
                              AZ3146_250nM = c("-","-","-","+","+","+", 
                                               "-","-","-","+","+","+",
                                               "-","-","-","+","+","+",
                                               "-","-","-","+","+","+"), 
                              Colony = c(209,229,208,195,205,221, 42,60,50,16,12,14, 75,79,89,57,52,57,89,73,96,37,40,43))



HCT116.UBE2H.KOOE<- data.frame( CellLine = c("WT","WT","WT",  "WT","WT","WT",
                                             "UBE2H-KO_c1", "UBE2H-KO_c1", "UBE2H-KO_c1", "UBE2H-KO_c1", "UBE2H-KO_c1", "UBE2H-KO_c1",
                                             "UBE2H-KO_c2", "UBE2H-KO_c2", "UBE2H-KO_c2", "UBE2H-KO_c2", "UBE2H-KO_c2", "UBE2H-KO_c2",
                                             "UBE2H-KO_c3", "UBE2H-KO_c3", "UBE2H-KO_c3", "UBE2H-KO_c3", "UBE2H-KO_c3", "UBE2H-KO_c3", 
                                             "UBE2H-OE","UBE2H-OE","UBE2H-OE",  "UBE2H-OE","UBE2H-OE","UBE2H-OE",
                                             "UBE2H-KO_c1+OE", "UBE2H-KO_c1+OE", "UBE2H-KO_c1+OE", "UBE2H-KO_c1+OE", "UBE2H-KO_c1+OE", "UBE2H-KO_c1+OE",
                                             "UBE2H-KO_c2+OE", "UBE2H-KO_c2+OE", "UBE2H-KO_c2+OE", "UBE2H-KO_c2+OE", "UBE2H-KO_c2+OE", "UBE2H-KO_c2+OE",
                                             "UBE2H-KO_c3+OE", "UBE2H-KO_c3+OE", "UBE2H-KO_c3+OE", "UBE2H-KO_c3+OE", "UBE2H-KO_c3+OE", "UBE2H-KO_c3+OE"), 
                                AZ3146_250nM = c("-","-","-","+","+","+", 
                                                 "-","-","-","+","+","+",
                                                 "-","-","-","+","+","+",
                                                 "-","-","-","+","+","+", 
                                                 "-","-","-","+","+","+", 
                                                 "-","-","-","+","+","+",
                                                 "-","-","-","+","+","+",
                                                 "-","-","-","+","+","+"), 
                                Colony = c(155,129,149,128,114,116, 
                                           23,22,36,9,8,9, 
                                           73,41,49,36,30,32,
                                           89,82,98,51,47,50, 
                                           198,173,167,163,171,160,
                                           63,50,47,33,22,26,
                                           75,84,64,61,58,63,
                                           71,95,66,67,56,57))


HCT116.UBE2H.KOOE2<- data.frame( CellLine = c("WT","WT","WT",  "WT","WT","WT",
                                              "UBE2H-KO_c1", "UBE2H-KO_c1", "UBE2H-KO_c1", "UBE2H-KO_c1", "UBE2H-KO_c1", "UBE2H-KO_c1",
                                              "UBE2H-KO_c2", "UBE2H-KO_c2", "UBE2H-KO_c2", "UBE2H-KO_c2", "UBE2H-KO_c2", "UBE2H-KO_c2",
                                              "UBE2H-KO_c3", "UBE2H-KO_c3", "UBE2H-KO_c3", "UBE2H-KO_c3", "UBE2H-KO_c3", "UBE2H-KO_c3", 
                                              "UBE2H-OE","UBE2H-OE","UBE2H-OE",  "UBE2H-OE","UBE2H-OE","UBE2H-OE",
                                              "UBE2H-KO_c1+OE", "UBE2H-KO_c1+OE", "UBE2H-KO_c1+OE", "UBE2H-KO_c1+OE", "UBE2H-KO_c1+OE", "UBE2H-KO_c1+OE",
                                              "UBE2H-KO_c2+OE", "UBE2H-KO_c2+OE", "UBE2H-KO_c2+OE", "UBE2H-KO_c2+OE", "UBE2H-KO_c2+OE", "UBE2H-KO_c2+OE",
                                              "UBE2H-KO_c3+OE", "UBE2H-KO_c3+OE", "UBE2H-KO_c3+OE", "UBE2H-KO_c3+OE", "UBE2H-KO_c3+OE", "UBE2H-KO_c3+OE"), 
                                 AZ3146_250nM = c("-","-","-","+","+","+", 
                                                  "-","-","-","+","+","+",
                                                  "-","-","-","+","+","+",
                                                  "-","-","-","+","+","+", 
                                                  "-","-","-","+","+","+", 
                                                  "-","-","-","+","+","+",
                                                  "-","-","-","+","+","+",
                                                  "-","-","-","+","+","+"), 
                                 Colony = c(207,231,251,209,181,195,
                                            97,116,125,54,41,71,
                                            77,97,76,42,49,40,
                                            149,182,157,121,90,96, 
                                            206,207,178,210,202,179,
                                            101,117,112,77,102,87,
                                            78,69,69,55,63,58,
                                            104,113,83,76,79,88))

#### Reversine
HCT116.Rev<- data.frame( CellLine = c("WT","WT","WT",  "WT","WT","WT",
                                      "UBE2H-KO_c1", "UBE2H-KO_c1", "UBE2H-KO_c1", "UBE2H-KO_c1", "UBE2H-KO_c1", "UBE2H-KO_c1",
                                      "UBE2H-KO_c2", "UBE2H-KO_c2", "UBE2H-KO_c2", "UBE2H-KO_c2", "UBE2H-KO_c2", "UBE2H-KO_c2",
                                      "UBE2H-KO_c3", "UBE2H-KO_c3", "UBE2H-KO_c3", "UBE2H-KO_c3", "UBE2H-KO_c3", "UBE2H-KO_c3", 
                                      "UBE2H-OE","UBE2H-OE","UBE2H-OE",  "UBE2H-OE","UBE2H-OE","UBE2H-OE",
                                      "UBE2H-KO_c1+OE", "UBE2H-KO_c1+OE", "UBE2H-KO_c1+OE", "UBE2H-KO_c1+OE", "UBE2H-KO_c1+OE", "UBE2H-KO_c1+OE",
                                      "UBE2H-KO_c2+OE", "UBE2H-KO_c2+OE", "UBE2H-KO_c2+OE", "UBE2H-KO_c2+OE", "UBE2H-KO_c2+OE", "UBE2H-KO_c2+OE",
                                      "UBE2H-KO_c3+OE", "UBE2H-KO_c3+OE", "UBE2H-KO_c3+OE", "UBE2H-KO_c3+OE", "UBE2H-KO_c3+OE", "UBE2H-KO_c3+OE"), 
                         Reversine_100nM = c("-","-","-","+","+","+", 
                                             "-","-","-","+","+","+",
                                             "-","-","-","+","+","+",
                                             "-","-","-","+","+","+", 
                                             "-","-","-","+","+","+", 
                                             "-","-","-","+","+","+",
                                             "-","-","-","+","+","+",
                                             "-","-","-","+","+","+"), 
                         Colony = c(179,173,160, 69,73,81, 
                                    54,56,41,2,6,6, 
                                    86,63,70,5,1,5, 
                                    0,0,0,0,0,0,
                                    205,188,177,125,104,107, 
                                    75,52,54, 3,6,6, 
                                    96,84,79, 25,26,30, 
                                    137,88,110, 26,15,17))
HCT116.Rev2<- data.frame( CellLine = c("WT","WT","WT",  "WT","WT","WT",
                                       "UBE2H-KO_c1", "UBE2H-KO_c1", "UBE2H-KO_c1", "UBE2H-KO_c1", "UBE2H-KO_c1", "UBE2H-KO_c1",
                                       "UBE2H-KO_c2", "UBE2H-KO_c2", "UBE2H-KO_c2", "UBE2H-KO_c2", "UBE2H-KO_c2", "UBE2H-KO_c2",
                                       "UBE2H-KO_c3", "UBE2H-KO_c3", "UBE2H-KO_c3", "UBE2H-KO_c3", "UBE2H-KO_c3", "UBE2H-KO_c3", 
                                       "UBE2H-OE","UBE2H-OE","UBE2H-OE",  "UBE2H-OE","UBE2H-OE","UBE2H-OE",
                                       "UBE2H-KO_c1+OE", "UBE2H-KO_c1+OE", "UBE2H-KO_c1+OE", "UBE2H-KO_c1+OE", "UBE2H-KO_c1+OE", "UBE2H-KO_c1+OE",
                                       "UBE2H-KO_c2+OE", "UBE2H-KO_c2+OE", "UBE2H-KO_c2+OE", "UBE2H-KO_c2+OE", "UBE2H-KO_c2+OE", "UBE2H-KO_c2+OE",
                                       "UBE2H-KO_c3+OE", "UBE2H-KO_c3+OE", "UBE2H-KO_c3+OE", "UBE2H-KO_c3+OE", "UBE2H-KO_c3+OE", "UBE2H-KO_c3+OE"), 
                          Reversine_100nM = c("-","-","-","+","+","+", 
                                              "-","-","-","+","+","+",
                                              "-","-","-","+","+","+",
                                              "-","-","-","+","+","+", 
                                              "-","-","-","+","+","+", 
                                              "-","-","-","+","+","+",
                                              "-","-","-","+","+","+",
                                              "-","-","-","+","+","+"), 
                          Colony = c(183, 196, 170, 192, 194, 187, 
                                     163,132,129,102,112,109,
                                     84,63,72,  49,32,44,
                                     249, 254, 238,  189, 130, 165, 
                                     233, 177, 173,  205, 176, 189, 
                                     113, 128, 111,  110, 112, 130, 
                                     104, 111, 92,   105, 99, 98, 
                                     205, 177, 174,   174, 161, 158))


HCT116.Rev3<- data.frame( CellLine = c("WT","WT","WT",  "WT","WT","WT",
                                       "UBE2H-KO_c1", "UBE2H-KO_c1", "UBE2H-KO_c1", "UBE2H-KO_c1", "UBE2H-KO_c1", "UBE2H-KO_c1",
                                       "UBE2H-KO_c2", "UBE2H-KO_c2", "UBE2H-KO_c2", "UBE2H-KO_c2", "UBE2H-KO_c2", "UBE2H-KO_c2",
                                       "UBE2H-KO_c3", "UBE2H-KO_c3", "UBE2H-KO_c3", "UBE2H-KO_c3", "UBE2H-KO_c3", "UBE2H-KO_c3", 
                                       "UBE2H-OE","UBE2H-OE","UBE2H-OE",  "UBE2H-OE","UBE2H-OE","UBE2H-OE",
                                       "UBE2H-KO_c1+OE", "UBE2H-KO_c1+OE", "UBE2H-KO_c1+OE", "UBE2H-KO_c1+OE", "UBE2H-KO_c1+OE", "UBE2H-KO_c1+OE",
                                       "UBE2H-KO_c2+OE", "UBE2H-KO_c2+OE", "UBE2H-KO_c2+OE", "UBE2H-KO_c2+OE", "UBE2H-KO_c2+OE", "UBE2H-KO_c2+OE",
                                       "UBE2H-KO_c3+OE", "UBE2H-KO_c3+OE", "UBE2H-KO_c3+OE", "UBE2H-KO_c3+OE", "UBE2H-KO_c3+OE", "UBE2H-KO_c3+OE"), 
                          Reversine_100nM = c("-","-","-","+","+","+", 
                                              "-","-","-","+","+","+",
                                              "-","-","-","+","+","+",
                                              "-","-","-","+","+","+", 
                                              "-","-","-","+","+","+", 
                                              "-","-","-","+","+","+",
                                              "-","-","-","+","+","+",
                                              "-","-","-","+","+","+"), 
                          Colony = c(365,415,385,  365, 341, 351, # WT 1 top half # 379,330,322, 309,357,310, #307, 349, 333,   329, 306, 340,
                                     126, 119, 137,   74,  63,  81,  
                                     68,  71,  79,    42,  32,  46, 
                                     194, 180, 204,   112, 132, 110,
                                     256, 299, 306,   226, 239, 241, 
                                     121, 115, 118,   87,  97,  80, 
                                     89,  94,  75,    72,  54,  63, 
                                     134, 116, 125,   96,  88,  94))


Cal51.Rev<- data.frame( CellLine = c("WT", "WT","WT","WT","WT","WT",
                                     "UBE2H-KO_c2", "UBE2H-KO_c2",  "UBE2H-KO_c2",  "UBE2H-KO_c2",  "UBE2H-KO_c2",  "UBE2H-KO_c2", 
                                     "UBE2H-KO_c7", "UBE2H-KO_c7", "UBE2H-KO_c7", "UBE2H-KO_c7", "UBE2H-KO_c7", "UBE2H-KO_c7"), 
                        Reversine_100nM = c("-","-","-","+","+","+", 
                                            "-","-","-","+","+","+", 
                                            "-","-","-","+","+","+"), 
                        Colony = c(170,181,185, 74,49,58, 
                                   66,61,63, 0,0,0, 
                                   130,112,118,4,8,12))

Cal51.Rev2<- data.frame( CellLine = c("WT", "WT","WT","WT","WT","WT",
                                      "WT", "WT","WT","WT","WT","WT",
                                      "UBE2H-KO_c2", "UBE2H-KO_c2",  "UBE2H-KO_c2",  "UBE2H-KO_c2",  "UBE2H-KO_c2",  "UBE2H-KO_c2", 
                                      "UBE2H-KO_c7", "UBE2H-KO_c7", "UBE2H-KO_c7", "UBE2H-KO_c7", "UBE2H-KO_c7", "UBE2H-KO_c7"), 
                         Reversine_100nM = c("-","-","-","+","+","+", 
                                             "-","-","-","+","+","+", 
                                             "-","-","-","+","+","+", 
                                             "-","-","-","+","+","+"), 
                         Colony = c(336, 370, 344,  304, 340, 296,
                                    356,313,372,  262,281,279, 
                                    327,353,294,   224,205,222,
                                    373,263,330,  182,214,181))

Cal51.Rev3<- data.frame( CellLine = c("WT", "WT","WT", "WT","WT","WT",
                                      "WT", "WT","WT", "WT","WT","WT",
                                      "UBE2H-KO_c2", "UBE2H-KO_c2",  "UBE2H-KO_c2",  "UBE2H-KO_c2",  "UBE2H-KO_c2",  "UBE2H-KO_c2", 
                                      "UBE2H-KO_c2", "UBE2H-KO_c2",  "UBE2H-KO_c2",  "UBE2H-KO_c2",  "UBE2H-KO_c2",  "UBE2H-KO_c2", 
                                      "UBE2H-KO_c7", "UBE2H-KO_c7", "UBE2H-KO_c7", "UBE2H-KO_c7", "UBE2H-KO_c7", "UBE2H-KO_c7", 
                                      "UBE2H-KO_c7", "UBE2H-KO_c7", "UBE2H-KO_c7", "UBE2H-KO_c7", "UBE2H-KO_c7", "UBE2H-KO_c7"), 
                         Reversine_100nM = c("-","-","-","+","+","+", 
                                             "-","-","-","+","+","+", 
                                             "-","-","-","+","+","+", 
                                             "-","-","-","+","+","+", 
                                             "-","-","-","+","+","+", 
                                             "-","-","-","+","+","+"), 
                         Colony = c(353,390,405, 226,216,224, 
                                    335,363,352,  216,231,229, 
                                    348,339,342,  106,156,129, 
                                    317,358,324,  144,145,140, 
                                    313,274,313,  80, 89, 76, 
                                    318,229,267,  93,72,79))

# BAY1217389
HCT116.Bay<- data.frame( CellLine = c("WT","WT","WT",  "WT","WT","WT",
                                      "UBE2H-KO_c1", "UBE2H-KO_c1", "UBE2H-KO_c1", "UBE2H-KO_c1", "UBE2H-KO_c1", "UBE2H-KO_c1",
                                      "UBE2H-KO_c2", "UBE2H-KO_c2", "UBE2H-KO_c2", "UBE2H-KO_c2", "UBE2H-KO_c2", "UBE2H-KO_c2",
                                      "UBE2H-KO_c3", "UBE2H-KO_c3", "UBE2H-KO_c3", "UBE2H-KO_c3", "UBE2H-KO_c3", "UBE2H-KO_c3", 
                                      "UBE2H-OE","UBE2H-OE","UBE2H-OE",  "UBE2H-OE","UBE2H-OE","UBE2H-OE",
                                      "UBE2H-KO_c1+OE", "UBE2H-KO_c1+OE", "UBE2H-KO_c1+OE", "UBE2H-KO_c1+OE", "UBE2H-KO_c1+OE", "UBE2H-KO_c1+OE",
                                      "UBE2H-KO_c2+OE", "UBE2H-KO_c2+OE", "UBE2H-KO_c2+OE", "UBE2H-KO_c2+OE", "UBE2H-KO_c2+OE", "UBE2H-KO_c2+OE",
                                      "UBE2H-KO_c3+OE", "UBE2H-KO_c3+OE", "UBE2H-KO_c3+OE", "UBE2H-KO_c3+OE", "UBE2H-KO_c3+OE", "UBE2H-KO_c3+OE"), 
                         BAY1218389_10nM = c("-","-","-","+","+","+", 
                                             "-","-","-","+","+","+",
                                             "-","-","-","+","+","+",
                                             "-","-","-","+","+","+", 
                                             "-","-","-","+","+","+", 
                                             "-","-","-","+","+","+",
                                             "-","-","-","+","+","+",
                                             "-","-","-","+","+","+"), 
                         Colony = c(215,241,174,169,129,137, 
                                    213,204,223, 4,4,7, 
                                    71,73,62, 1,2,1, 
                                    322,299,299, 65,66,27, 
                                    322,302,320, 106,189,201,
                                    101,109,94, 56,50,32, 
                                    61,66,54, 27,24,26, 
                                    109,102,109, 48,38,25))

### Plot colonies AZ3146 ####
  #### Cal51 AZ3146 colony formation ####
    ##### August 2025 ####
## Plot Cal51 +/- UBE2H-KO +/- AZ3146 
# Calculate relative colony number: 

Cal51.UBE2H.AZ$Cell.Drug  <-  paste0(Cal51.UBE2H.AZ$CellLine, Cal51.UBE2H.AZ$AZ3146_250nM)
WTMean<- mean(Cal51.UBE2H.AZ$Colony[1:3])
WTc2<- mean(Cal51.UBE2H.AZ$Colony[7:9])
WTc7<- mean(Cal51.UBE2H.AZ$Colony[13:15])

Cal51.UBE2H.AZ$Rel.colony <-  c( Cal51.UBE2H.AZ$Colony[1:6]/WTMean, 
                                 Cal51.UBE2H.AZ$Colony[7:12]/WTc2, 
                                 Cal51.UBE2H.AZ$Colony[13:18]/WTc7)
Cal51.UBE2H.AZ$Cell.Drug <- factor(Cal51.UBE2H.AZ$Cell.Drug, level= unique(Cal51.UBE2H.AZ$Cell.Drug))

# p-value
# one sided t-test
# is AZ3146 more toxic in UBE2H-KO clones than wt? 
# Cal51 UBE2H-KO c2: 
t.test(Cal51.UBE2H.AZ$Rel.colony[4:6], # WT + AZ3146 relative colony number
       Cal51.UBE2H.AZ$Rel.colony[10:12],
       alternative = c("greater") ) # UBE2H-KO clone 2
# 0.2463

# Cal51 UBE2H-KO c7: 
t.test(Cal51.UBE2H.AZ$Rel.colony[4:6], # WT + AZ3146 relative colony number
       Cal51.UBE2H.AZ$Rel.colony[16:18],
       alternative = c("greater") ) # UBE2H-KO clone 7
# 0.04099



# plot
ggplot(Cal51.UBE2H.AZ, aes(x= Cell.Drug, y= Rel.colony))+ 
  geom_boxplot()+
  geom_point()+  
  theme_classic()+
  ylim(0,1.3)+
  ylab("Relative colony formation\n upon 250nM AZ3146")+
  xlab("Cell line")
# Cal51.AZ3146.UBE2HKO.RelativeColony.pdf

Cal51.UBE2H.AZ$color <- c("WT", "WT","WT", "WT", "WT","WT", 
                           "KO", "KO", "KO", "KO", "KO", "KO", 
                           "KO", "KO", "KO", "KO", "KO", "KO")

ggplot(subset(Cal51.UBE2H.AZ, AZ3146_250nM=="+") , aes(x= Cell.Drug, y= Rel.colony, fill = color))+ 
  geom_boxplot()+
  geom_point()+  
  theme_classic()+
  ylim(0, 0.6)+
  ylab("Relative colony formation\n upon 250nM AZ3146")+
  xlab("Cell line")
# Cal51.AZ3146.UBE2HKO.RelativeColony_AZ_Aug2025.pdf





    ##### October 2025 ####
## Plot Cal51 +/- UBE2H-KO +/- AZ3146 
# Calculate relative colony number: 

Cal51.UBE2H.AZ2$Cell.Drug  <-  paste0(Cal51.UBE2H.AZ2$CellLine, Cal51.UBE2H.AZ2$AZ3146_250nM)
WTMean<- mean(Cal51.UBE2H.AZ2$Colony[1:3])
WTc2<- mean(Cal51.UBE2H.AZ2$Colony[7:9])
WTc7<- mean(Cal51.UBE2H.AZ2$Colony[13:15])

Cal51.UBE2H.AZ2$Rel.colony <-  c( Cal51.UBE2H.AZ2$Colony[1:6]/WTMean, 
                                  Cal51.UBE2H.AZ2$Colony[7:12]/WTc2, 
                                  Cal51.UBE2H.AZ2$Colony[13:18]/WTc7)
Cal51.UBE2H.AZ2$Cell.Drug <- factor(Cal51.UBE2H.AZ2$Cell.Drug, level= unique(Cal51.UBE2H.AZ2$Cell.Drug))

# p-value
# one sided t-test
# is AZ3146 more toxic in UBE2H-KO clones than wt? 
# Cal51 UBE2H-KO c2: 
t.test(Cal51.UBE2H.AZ2$Rel.colony[4:6], # WT + AZ3146 relative colony number
       Cal51.UBE2H.AZ2$Rel.colony[10:12],
       alternative = c("greater") ) # UBE2H-KO clone 2
# 0.04772

# Cal51 UBE2H-KO c7: 
t.test(Cal51.UBE2H.AZ2$Rel.colony[4:6], # WT + AZ3146 relative colony number
       Cal51.UBE2H.AZ2$Rel.colony[16:18],
       alternative = c("greater") ) # UBE2H-KO clone 7
# 0.06445

# Cal51 UBE2H-KO both: 
t.test(Cal51.UBE2H.AZ2$Rel.colony[4:6], # WT + AZ3146 relative colony number
       Cal51.UBE2H.AZ2$Rel.colony[c(10:12,16:18)],
       alternative = c("greater") ) # UBE2H-KO clone 7
# 0.01398

Cal51.UBE2H.AZ2$color <- c("WT", "WT","WT", "WT", "WT","WT", 
                           "KO", "KO", "KO", "KO", "KO", "KO", 
                           "KO", "KO", "KO", "KO", "KO", "KO")

# plot
ggplot(Cal51.UBE2H.AZ2, aes(x= Cell.Drug, y= Rel.colony))+ 
  geom_boxplot()+
  geom_point()+  
  theme_classic()+
  ylim(0,1.3)+
  ylab("Relative colony formation\n upon 250nM AZ3146")+
  xlab("Cell line")
# Cal51.AZ3146.UBE2HKO.RelativeColony.pdf


ggplot(subset(Cal51.UBE2H.AZ2, AZ3146_250nM=="+") , aes(x= Cell.Drug, y= Rel.colony, fill =color ))+ 
  geom_boxplot()+
  geom_point()+  
  theme_classic()+
  scale_fill_manual(values= c("mediumorchid2", "turquoise2"))+
  ylim(0, 1)+
  ylab("Relative colony formation\n upon 250nM AZ3146")+
  xlab("Cal51")
# Cal51.AZ3146.UBE2HKO.RelativeColony_AZ.pdf





    ##### *February 2026 ####
## Plot Cal51 +/- UBE2H-KO +/- AZ3146 
# Calculate relative colony number: 
# 250nM AZ3146 

Cal51.UBE2H.AZ3$Cell.Drug  <-  paste0(Cal51.UBE2H.AZ3$CellLine, Cal51.UBE2H.AZ3$AZ3146_250nM)
WTMean<- mean(Cal51.UBE2H.AZ3$Colony[1:3])
WTc2<- mean(Cal51.UBE2H.AZ3$Colony[7:9])
WTc7<- mean(Cal51.UBE2H.AZ3$Colony[13:15])

Cal51.UBE2H.AZ3$Rel.colony <-  c( Cal51.UBE2H.AZ3$Colony[1:6]/WTMean, 
                                  Cal51.UBE2H.AZ3$Colony[7:12]/WTc2, 
                                  Cal51.UBE2H.AZ3$Colony[13:18]/WTc7)
Cal51.UBE2H.AZ3$Cell.Drug <- factor(Cal51.UBE2H.AZ3$Cell.Drug, level= unique(Cal51.UBE2H.AZ3$Cell.Drug))

# p-value
# one sided t-test
# is AZ3146 more toxic in UBE2H-KO clones than wt? 
# Cal51 UBE2H-KO c2: 
t.test(Cal51.UBE2H.AZ3$Rel.colony[4:6], # WT + AZ3146 relative colony number
       Cal51.UBE2H.AZ3$Rel.colony[10:12]) # UBE2H-KO clone 2
# 0.01484

# Cal51 UBE2H-KO c7: 
t.test(Cal51.UBE2H.AZ3$Rel.colony[4:6], # WT + AZ3146 relative colony number
       Cal51.UBE2H.AZ3$Rel.colony[16:18]) # UBE2H-KO clone 7
# 0.007886


Cal51.UBE2H.AZ3$color <- c("WT", "WT","WT", "WT", "WT","WT", 
                           "KO", "KO", "KO", "KO", "KO", "KO", 
                           "KO", "KO", "KO", "KO", "KO", "KO")

# plot
ggplot(Cal51.UBE2H.AZ3, aes(x= Cell.Drug, y= Rel.colony))+ 
  geom_boxplot()+
  geom_point()+  
  theme_classic()+
  ylim(0,1.3)+
  ylab("Relative colony formation\n upon 250nM AZ3146")+
  xlab("Cell line")
# Cal51.AZ3146.UBE2HKO.RelativeColony_Feb2026.pdf


ggplot(subset(Cal51.UBE2H.AZ3, AZ3146_250nM=="+") , aes(x= Cell.Drug, y= Rel.colony, fill =color ))+ 
  geom_boxplot()+
  geom_point()+  
  theme_classic()+
  scale_fill_manual(values= c("mediumorchid2", "turquoise2"))+
  ylim(0, 1)+
  ylab("Relative colony formation\n upon 250nM AZ3146")+
  xlab("Cal51")
# Cal51.AZ3146.UBE2HKO.RelativeColony_AZ_Feb2026.pdf
# Figure 5H





# Merge 2 and 3: 
Cal51.Merge<- rbind(Cal51.UBE2H.AZ2, Cal51.UBE2H.AZ3)
# Cal51 UBE2H-KO c2: 
t.test(Cal51.Merge$Rel.colony[c(4:6, 22:24)], # WT + AZ3146 relative colony number
       Cal51.Merge$Rel.colony[c(10:12, 28:30)]) # UBE2H-KO clone 2
# 0.00374

# Cal51 UBE2H-KO c7: 
t.test(Cal51.Merge$Rel.colony[c(4:6, 22:24)], # WT + AZ3146 relative colony number
       Cal51.Merge$Rel.colony[c(16:18, 34:36)]) # UBE2H-KO clone 7
# 0.0123

ggplot(subset(Cal51.Merge, AZ3146_250nM=="+") , aes(x= Cell.Drug, y= Rel.colony, fill =color ))+ 
  geom_boxplot()+
  geom_point()+  
  theme_classic()+
  scale_fill_manual(values= c("mediumorchid2", "turquoise2"))+
  ylim(0, 1)+
  ylab("Relative colony formation\n upon 250nM AZ3146")+
  xlab("Cal51")
# Cal51.AZ3146.UBE2HKO.RelativeColony_AZ.pdf



  #### HCT116 AZ3146 colony formation ####
    ##### November 19 2025: 250nM  KO + OE ####
# 250nM AZ3146

HCT116.UBE2H.KOOE$Cell.Drug  <-  paste0(HCT116.UBE2H.KOOE$CellLine, HCT116.UBE2H.KOOE$AZ3146_250nM)
HCTMean<- mean(HCT116.UBE2H.KOOE$Colony[1:3])
HCTc1<- mean(HCT116.UBE2H.KOOE$Colony[7:9])
HCTc2<- mean(HCT116.UBE2H.KOOE$Colony[13:15])
HCTc3<- mean(HCT116.UBE2H.KOOE$Colony[19:21])
HCTOE<- mean(HCT116.UBE2H.KOOE$Colony[25:27])
HCTc1OE<- mean(HCT116.UBE2H.KOOE$Colony[31:33])
HCTc2OE<- mean(HCT116.UBE2H.KOOE$Colony[37:39])
HCTc3OE<- mean(HCT116.UBE2H.KOOE$Colony[43:45])

HCT116.UBE2H.KOOE$Rel.colony <-  c( HCT116.UBE2H.KOOE$Colony[1:6]/HCTMean, 
                                    HCT116.UBE2H.KOOE$Colony[7:12]/HCTc1,
                                    HCT116.UBE2H.KOOE$Colony[13:18]/HCTc2 ,
                                    HCT116.UBE2H.KOOE$Colony[19:24]/HCTc3 , 
                                    HCT116.UBE2H.KOOE$Colony[25:30]/HCTOE, 
                                    HCT116.UBE2H.KOOE$Colony[31:36]/HCTc1OE,
                                    HCT116.UBE2H.KOOE$Colony[37:42]/HCTc2OE ,
                                    HCT116.UBE2H.KOOE$Colony[43:48]/HCTc3OE 
)

HCT116.UBE2H.KOOE$Cell.Drug <- factor(HCT116.UBE2H.KOOE$Cell.Drug, 
                                      level= c("WT-", "WT+", "UBE2H-OE-", "UBE2H-OE+", 
                                               "UBE2H-KO_c1-", "UBE2H-KO_c1+","UBE2H-KO_c1+OE-", "UBE2H-KO_c1+OE+",
                                               "UBE2H-KO_c2-", "UBE2H-KO_c2+","UBE2H-KO_c2+OE-", "UBE2H-KO_c2+OE+",
                                               "UBE2H-KO_c3-", "UBE2H-KO_c3+","UBE2H-KO_c3+OE-", "UBE2H-KO_c3+OE+"))

HCT116.UBE2H.KOOE$color<- c("WT", "WT","WT", "WT", "WT","WT", 
                            "KO", "KO", "KO", "KO", "KO", "KO", 
                            "KO", "KO", "KO", "KO", "KO", "KO", 
                            "KO", "KO", "KO", "KO", "KO", "KO", 
                            "OE", "OE", "OE", "OE", "OE", "OE", 
                            "KO + OE","KO + OE","KO + OE",  "KO + OE","KO + OE","KO + OE", 
                            "KO + OE","KO + OE","KO + OE",  "KO + OE","KO + OE","KO + OE", 
                            "KO + OE","KO + OE","KO + OE",  "KO + OE","KO + OE","KO + OE")

# p-value
# one sided t-test: 
# is AZ3146 more toxic in UBE2H-KO clones than wt? 

# HCT116 UBE2H-KO c1: 
t.test(HCT116.UBE2H.KOOE$Rel.colony[4:6], # WT + AZ3146 relative colony number
       HCT116.UBE2H.KOOE$Rel.colony[10:12]) # UBE2H-KO clone 1
#  0.00114

# HCT116 UBE2H-KO c2: 
t.test(HCT116.UBE2H.KOOE$Rel.colony[4:6], # WT + AZ3146 relative colony number
       HCT116.UBE2H.KOOE$Rel.colony[16:18]) # UBE2H-KO clone 2
#  0.00717

# HCT116 UBE2H-KO c3: 
t.test(HCT116.UBE2H.KOOE$Rel.colony[4:6], # WT + AZ3146 relative colony number
       HCT116.UBE2H.KOOE$Rel.colony[22:24]) # UBE2H-KO clone 3
#  0.004922


# Overexpression different from WT? 
# HCT116 UBE2H-OE: 
t.test(HCT116.UBE2H.KOOE$Rel.colony[4:6], # WT + AZ3146 relative colony number
       HCT116.UBE2H.KOOE$Rel.colony[28:30]) # UBE2H-OE
#  0.03
# HCT116 KO1 UBE2H-OE: 
t.test(HCT116.UBE2H.KOOE$Rel.colony[4:6], # WT + AZ3146 relative colony number
       HCT116.UBE2H.KOOE$Rel.colony[34:36]) # UBE2H-OE
#  0.01838
# HCT116 KO2 UBE2H-OE: 
t.test(HCT116.UBE2H.KOOE$Rel.colony[4:6], # WT + AZ3146 relative colony number
       HCT116.UBE2H.KOOE$Rel.colony[40:42]) # UBE2H-OE
#  0.7848
# HCT116 KO3 UBE2H-OE: 
t.test(HCT116.UBE2H.KOOE$Rel.colony[4:6], # WT + AZ3146 relative colony number
       HCT116.UBE2H.KOOE$Rel.colony[46:48]) # UBE2H-OE
#  0.4108


# p-value
# one sided t-test: 
# is AZ3146 Less toxic in UBE2H-OE clones than corresponding KO? 

# HCT116 UBE2H-KO c1 + OE: 
t.test(HCT116.UBE2H.KOOE$Rel.colony[10:12], # WT + AZ3146 relative colony number
       HCT116.UBE2H.KOOE$Rel.colony[34:36],
       alternative = c("less")) # UBE2H-KO clone 1
#  0.04299

# HCT116 UBE2H-KO c2 + OE: 
t.test(HCT116.UBE2H.KOOE$Rel.colony[16:18], # WT + AZ3146 relative colony number
       HCT116.UBE2H.KOOE$Rel.colony[40:42],
       alternative = c("less")) # UBE2H-KO clone 2
#  0.00424

# HCT116 UBE2H-KO c3 + OE: 
t.test(HCT116.UBE2H.KOOE$Rel.colony[22:24], # WT + AZ3146 relative colony number
       HCT116.UBE2H.KOOE$Rel.colony[46:48],
       alternative = c("less")) # UBE2H-KO clone 3
#  0.01498



# plot HCT116: 
ggplot(HCT116.UBE2H.KOOE, aes(x= Cell.Drug, y= Rel.colony))+ 
  geom_boxplot()+
  geom_point()+  
  theme_classic()+
  ylim(0,1)+
  ylab("Relative colony formation\n upon 250nM AZ3146")+
  xlab("Cell line")
# HCT116.AZ3146.UBE2HKO.OE.RelativeColony.pdf


ggplot(subset(HCT116.UBE2H.KOOE, AZ3146_250nM=="+") , aes(x= Cell.Drug, y= Rel.colony, fill = color))+ 
  geom_boxplot()+
  scale_fill_manual(values= c("mediumorchid2", "mediumorchid4", "turquoise4", "turquoise2"))+
  geom_point()+  
  theme_classic()+
  ylim(0,1)+
  ylab("Relative colony formation\n upon 250nM AZ3146")+
  xlab("Cell line")
# HCT116.AZ3146.UBE2HKO.OE.RelativeColony.pdf




    ##### *February 7 2026:  250nM   KO +OE ####
# 250nM AZ3146
HCT116.UBE2H.KOOE2$Cell.Drug  <-  paste0(HCT116.UBE2H.KOOE2$CellLine, HCT116.UBE2H.KOOE2$AZ3146_250nM)
HCTMean<- mean(HCT116.UBE2H.KOOE2$Colony[c(1:3)])
HCTc1<- mean(HCT116.UBE2H.KOOE2$Colony[7:9])
HCTc2<- mean(HCT116.UBE2H.KOOE2$Colony[13:15])
HCTc3<- mean(HCT116.UBE2H.KOOE2$Colony[19:21])
HCTOE<- mean(HCT116.UBE2H.KOOE2$Colony[25:27])
HCTc1OE<- mean(HCT116.UBE2H.KOOE2$Colony[31:33])
HCTc2OE<- mean(HCT116.UBE2H.KOOE2$Colony[37:39])
HCTc3OE<- mean(HCT116.UBE2H.KOOE2$Colony[43:45])

HCT116.UBE2H.KOOE2$Rel.colony <-  c( HCT116.UBE2H.KOOE2$Colony[1:6]/HCTMean, 
                                     HCT116.UBE2H.KOOE2$Colony[7:12]/HCTc1 ,
                                     HCT116.UBE2H.KOOE2$Colony[13:18]/HCTc2 ,
                                     HCT116.UBE2H.KOOE2$Colony[19:24]/HCTc3 , 
                                     HCT116.UBE2H.KOOE2$Colony[25:30]/HCTOE, 
                                     HCT116.UBE2H.KOOE2$Colony[31:36]/HCTc1OE,
                                     HCT116.UBE2H.KOOE2$Colony[37:42]/HCTc2OE ,
                                     HCT116.UBE2H.KOOE2$Colony[43:48]/HCTc3OE 
)

HCT116.UBE2H.KOOE2$Cell.Drug <- factor(HCT116.UBE2H.KOOE2$Cell.Drug, 
                                      level= c("WT-", "WT+", "UBE2H-OE-", "UBE2H-OE+", 
                                               "UBE2H-KO_c1-", "UBE2H-KO_c1+","UBE2H-KO_c1+OE-", "UBE2H-KO_c1+OE+",
                                               "UBE2H-KO_c2-", "UBE2H-KO_c2+","UBE2H-KO_c2+OE-", "UBE2H-KO_c2+OE+",
                                               "UBE2H-KO_c3-", "UBE2H-KO_c3+","UBE2H-KO_c3+OE-", "UBE2H-KO_c3+OE+"))

HCT116.UBE2H.KOOE2$color<- c("WT", "WT","WT", "WT", "WT","WT", 
                            "KO", "KO", "KO", "KO", "KO", "KO", 
                            "KO", "KO", "KO", "KO", "KO", "KO", 
                            "KO", "KO", "KO", "KO", "KO", "KO", 
                            "OE", "OE", "OE", "OE", "OE", "OE", 
                            "KO + OE","KO + OE","KO + OE",  "KO + OE","KO + OE","KO + OE", 
                            "KO + OE","KO + OE","KO + OE",  "KO + OE","KO + OE","KO + OE", 
                            "KO + OE","KO + OE","KO + OE",  "KO + OE","KO + OE","KO + OE")

# p-value
# one sided t-test: 
# is AZ3146 more toxic in UBE2H-KO clones than wt? 

# HCT116 UBE2H-KO c1: 
t.test(HCT116.UBE2H.KOOE2$Rel.colony[c(4:6)], # WT + AZ3146 relative colony number
       HCT116.UBE2H.KOOE2$Rel.colony[10:12]) # UBE2H-KO clone 1
# 0.02781

# HCT116 UBE2H-KO c2: 
t.test(HCT116.UBE2H.KOOE2$Rel.colony[c(4:6)], # WT + AZ3146 relative colony number
       HCT116.UBE2H.KOOE2$Rel.colony[16:18]) # UBE2H-KO clone 2
# 0.002541

# HCT116 UBE2H-KO c3: 
t.test(HCT116.UBE2H.KOOE2$Rel.colony[c(4:6)], # WT + AZ3146 relative colony number
       HCT116.UBE2H.KOOE2$Rel.colony[22:24]) # UBE2H-KO clone 3
# 0.04248


# Overexpression different from WT? 
# HCT116 UBE2H-OE: 
t.test(HCT116.UBE2H.KOOE2$Rel.colony[c(4:6)], # WT + AZ3146 relative colony number
       HCT116.UBE2H.KOOE2$Rel.colony[28:30]) # UBE2H-OE
# 0.06726
# HCT116 KO1 UBE2H-OE: 
t.test(HCT116.UBE2H.KOOE2$Rel.colony[c(4:6)], # WT + AZ3146 relative colony number
       HCT116.UBE2H.KOOE2$Rel.colony[34:36]) # UBE2H-OE
# 0.6052
# HCT116 KO2 UBE2H-OE: 
t.test(HCT116.UBE2H.KOOE2$Rel.colony[c(4:6)], # WT + AZ3146 relative colony number
       HCT116.UBE2H.KOOE2$Rel.colony[40:42]) # UBE2H-OE
# 0.514
# HCT116 KO3 UBE2H-OE: 
t.test(HCT116.UBE2H.KOOE2$Rel.colony[c(4:6)], # WT + AZ3146 relative colony number
       HCT116.UBE2H.KOOE2$Rel.colony[46:48]) # UBE2H-OE
# 0.4815


# p-value
# one sided t-test: 
# is AZ3146 Less toxic in UBE2H-OE clones than corresponding KO? 

# HCT116 UBE2H-KO c1 + OE: 
t.test(HCT116.UBE2H.KOOE2$Rel.colony[10:12], # WT + AZ3146 relative colony number
       HCT116.UBE2H.KOOE2$Rel.colony[34:36],
       alternative = c("less")) # UBE2H-KO clone 1
#  0.01865
# 0.01865

# HCT116 UBE2H-KO c2 + OE: 
t.test(HCT116.UBE2H.KOOE2$Rel.colony[16:18], # WT + AZ3146 relative colony number
       HCT116.UBE2H.KOOE2$Rel.colony[40:42],
       alternative = c("less")) # UBE2H-KO clone 2
#  0.001611
# 0.001611

# HCT116 UBE2H-KO c3 + OE: 
t.test(HCT116.UBE2H.KOOE2$Rel.colony[22:24], # WT + AZ3146 relative colony number
       HCT116.UBE2H.KOOE2$Rel.colony[46:48],
       alternative = c("less")) # UBE2H-KO clone 3
#  0.03487



# plot HCT116: 
ggplot(HCT116.UBE2H.KOOE2, aes(x= Cell.Drug, y= Rel.colony))+ 
  geom_boxplot()+
  geom_point()+  
  theme_classic()+
  ylim(0,1.2)+
  ylab("Relative colony formation\n upon 250nM AZ3146")+
  xlab("Cell line")
# HCT116.AZ3146.UBE2HKO.OE.RelativeColony.pdf


ggplot(subset(HCT116.UBE2H.KOOE2, AZ3146_250nM=="+") , aes(x= Cell.Drug, y= Rel.colony, fill = color))+ 
  geom_boxplot()+
  scale_fill_manual(values= c("mediumorchid2", "mediumorchid4", "turquoise4", "turquoise2"))+
  geom_point()+  
  theme_classic()+
  ylim(0,1.2)+
  ylab("Relative colony formation\n upon 250nM AZ3146")+
  xlab("Cell line")
# HCT116.AZ3146.UBE2HKO.OE.RelativeColony_February.pdf
# Figure 5F


# second run: colony formation is reproducible! 


    ##### Merge Nov Feb ####
# 250nM AZ3146
# both runs similar. REPRODUCIBLE!

# Merge Nov Feb: 
HCT116.Merge<- rbind(HCT116.UBE2H.KOOE, HCT116.UBE2H.KOOE2)

ggplot(subset(HCT116.Merge, AZ3146_250nM=="+") , aes(x= Cell.Drug, y= Rel.colony, fill = color))+ 
  geom_boxplot()+
  scale_fill_manual(values= c("mediumorchid2", "mediumorchid4", "turquoise4", "turquoise2"))+
  geom_point()+  
  theme_classic()+
  ylim(0,1.2)+
  ylab("Relative colony formation\n upon 250nM AZ3146")+
  xlab("Cell line")
# HCT116.AZ3146.UBE2HKO.OE.RelativeColony_February.pdf


# p-value
# one sided t-test: 
# is AZ3146 more toxic in UBE2H-KO clones than wt? 
# HCT116 UBE2H-KO c1: 
t.test(subset(HCT116.Merge, Cell.Drug =="WT+" )$Rel.colony, # WT + AZ3146 relative colony number
       subset(HCT116.Merge, Cell.Drug =="UBE2H-KO_c1+" )$Rel.colony) # UBE2H-KO clone 1
# 0.0001464

# HCT116 UBE2H-KO c2: 
t.test(subset(HCT116.Merge, Cell.Drug =="WT+" )$Rel.colony, # WT + AZ3146 relative colony number
       subset(HCT116.Merge, Cell.Drug =="UBE2H-KO_c2+" )$Rel.colony) # UBE2H-KO clone 1
# 1.566e-05

# HCT116 UBE2H-KO c3: 
t.test(subset(HCT116.Merge, Cell.Drug =="WT+" )$Rel.colony, # WT + AZ3146 relative colony number
       subset(HCT116.Merge, Cell.Drug =="UBE2H-KO_c3+" )$Rel.colony) # UBE2H-KO clone 1
# 0.0001379


# Overexpression different from WT? 
# WT vs HCT116 UBE2H-KO c3 + OE : 
t.test(subset(HCT116.Merge, Cell.Drug =="WT+" )$Rel.colony, # WT + AZ3146 relative colony number
       subset(HCT116.Merge, Cell.Drug =="UBE2H-OE+" )$Rel.colony) # OE
# 0.00819

# WT vs HCT116 UBE2H-KO c1 + OE : 
t.test(subset(HCT116.Merge, Cell.Drug =="WT+" )$Rel.colony, # WT + AZ3146 relative colony number
       subset(HCT116.Merge, Cell.Drug =="UBE2H-KO_c1+OE+" )$Rel.colony) # UBE2H-KO clone 1
# 0.06776 NS
 
# WT vs HCT116 UBE2H-KO c2 + OE : 
t.test(subset(HCT116.Merge, Cell.Drug =="WT+" )$Rel.colony, # WT + AZ3146 relative colony number
       subset(HCT116.Merge, Cell.Drug =="UBE2H-KO_c2+OE+" )$Rel.colony) # UBE2H-KO clone 2
# 0.4303 NS

# WT vs HCT116 UBE2H-KO c3 + OE : 
t.test(subset(HCT116.Merge, Cell.Drug =="WT+" )$Rel.colony, # WT + AZ3146 relative colony number
       subset(HCT116.Merge, Cell.Drug =="UBE2H-KO_c3+OE+" )$Rel.colony) # UBE2H-KO clone 3
# 0.2223


# p-value
# one sided t-test: 
# is AZ3146 Less toxic in UBE2H-OE clones than corresponding KO? 

# HCT116 UBE2H-KO c1 + OE: 
t.test(subset(HCT116.Merge, Cell.Drug =="UBE2H-KO_c1+" )$Rel.colony, # WT + AZ3146 relative colony number
       subset(HCT116.Merge, Cell.Drug =="UBE2H-KO_c1+OE+" )$Rel.colony) # UBE2H-KO clone 3
# 0.02635

# HCT116 UBE2H-KO c2 + OE: 
t.test(subset(HCT116.Merge, Cell.Drug =="UBE2H-KO_c2+" )$Rel.colony, # WT + AZ3146 relative colony number
       subset(HCT116.Merge, Cell.Drug =="UBE2H-KO_c2+OE+" )$Rel.colony) # UBE2H-KO clone 3
# 3.388e-05

# HCT116 UBE2H-KO c3 + OE: 
t.test(subset(HCT116.Merge, Cell.Drug =="UBE2H-KO_c3+" )$Rel.colony, # WT + AZ3146 relative colony number
       subset(HCT116.Merge, Cell.Drug =="UBE2H-KO_c3+OE+" )$Rel.colony) # UBE2H-KO clone 3
# 0.0007313



### Colony formation Reversine #####
  #### HCT116 colony formation REVERSINE ####
     ##### June 2026 ####
# 100nM Reversine
HCT116.Rev2$Cell.Drug  <-  paste0(HCT116.Rev2$CellLine, HCT116.Rev2$Reversine_100nM)

HCTMean<- mean(HCT116.Rev2$Colony[1:3])
HCTc1<- mean(HCT116.Rev2$Colony[7:9])
HCTc2<- mean(HCT116.Rev2$Colony[13:15])
HCTc3<- mean(HCT116.Rev2$Colony[19:21])
HCTOE<- mean(HCT116.Rev2$Colony[25:27])
HCTc1OE<- mean(HCT116.Rev2$Colony[31:33])
HCTc2OE<- mean(HCT116.Rev2$Colony[37:39])
HCTc3OE<- mean(HCT116.Rev2$Colony[43:45])

HCT116.Rev2$Rel.colony <-  c( HCT116.Rev2$Colony[1:6]/HCTMean, 
                              HCT116.Rev2$Colony[7:12]/HCTc1,
                              HCT116.Rev2$Colony[13:18]/HCTc2 ,
                              HCT116.Rev2$Colony[19:24]/HCTc3 , 
                              HCT116.Rev2$Colony[25:30]/HCTOE, 
                              HCT116.Rev2$Colony[31:36]/HCTc1OE,
                              HCT116.Rev2$Colony[37:42]/HCTc2OE ,
                              HCT116.Rev2$Colony[43:48]/HCTc3OE 
)

HCT116.Rev2$Cell.Drug <- factor(HCT116.Rev2$Cell.Drug, 
                               level= c("WT-", "WT+", "UBE2H-OE-", "UBE2H-OE+", 
                                        "UBE2H-KO_c1-", "UBE2H-KO_c1+","UBE2H-KO_c1+OE-", "UBE2H-KO_c1+OE+",
                                        "UBE2H-KO_c2-", "UBE2H-KO_c2+","UBE2H-KO_c2+OE-", "UBE2H-KO_c2+OE+",
                                        "UBE2H-KO_c3-", "UBE2H-KO_c3+","UBE2H-KO_c3+OE-", "UBE2H-KO_c3+OE+"))

HCT116.Rev2$color<- c("WT", "WT","WT", "WT", "WT","WT", 
                     "KO", "KO", "KO", "KO", "KO", "KO", 
                     "KO", "KO", "KO", "KO", "KO", "KO", 
                     "KO", "KO", "KO", "KO", "KO", "KO", 
                     "OE", "OE", "OE", "OE", "OE", "OE", 
                     "KO + OE","KO + OE","KO + OE",  "KO + OE","KO + OE","KO + OE", 
                     "KO + OE","KO + OE","KO + OE",  "KO + OE","KO + OE","KO + OE", 
                     "KO + OE","KO + OE","KO + OE",  "KO + OE","KO + OE","KO + OE")

# p-value
# one sided t-test: 
# is AZ3146 more toxic in UBE2H-KO clones than wt? 

# HCT116 UBE2H-KO c1: 
t.test(HCT116.Rev2$Rel.colony[4:6], # WT + AZ3146 relative colony number
       HCT116.Rev2$Rel.colony[10:12]) # UBE2H-KO clone 1
#  0.001141

# HCT116 UBE2H-KO c2: 
t.test(HCT116.Rev2$Rel.colony[4:6], # WT + AZ3146 relative colony number
       HCT116.Rev2$Rel.colony[16:18]) # UBE2H-KO clone 2
#  0.01849

# HCT116 UBE2H-KO c3: 
t.test(HCT116.Rev2$Rel.colony[4:6], # WT + AZ3146 relative colony number
       HCT116.Rev2$Rel.colony[22:24]) # UBE2H-KO clone 3
#  0.02745


# Overexpression different from WT? 
# HCT116 UBE2H-OE: 
t.test(HCT116.Rev2$Rel.colony[4:6], # WT + AZ3146 relative colony number
       HCT116.Rev2$Rel.colony[28:30]) # UBE2H-OE
#  0.2626
# HCT116 KO1 UBE2H-OE: 
t.test(HCT116.Rev2$Rel.colony[4:6], # WT + AZ3146 relative colony number
       HCT116.Rev2$Rel.colony[34:36]) # UBE2H-OE
#  0.5068
# HCT116 KO2 UBE2H-OE: 
t.test(HCT116.Rev2$Rel.colony[4:6], # WT + AZ3146 relative colony number
       HCT116.Rev2$Rel.colony[40:42]) # UBE2H-OE
#  0.08792
# HCT116 KO3 UBE2H-OE: 
t.test(HCT116.Rev2$Rel.colony[4:6], # WT + AZ3146 relative colony number
       HCT116.Rev2$Rel.colony[46:48]) # UBE2H-OE
#  0.01571


# p-value
# one sided t-test: 
# is Rev Less toxic in UBE2H-OE clones than corresponding KO? 

# HCT116 UBE2H-KO c1 + OE: 
t.test(HCT116.Rev2$Rel.colony[10:12], # WT + AZ3146 relative colony number
       HCT116.Rev2$Rel.colony[34:36],
       alternative = c("less")) # UBE2H-KO clone 1
#  0.01737

# HCT116 UBE2H-KO c2 + OE: 
t.test(HCT116.Rev2$Rel.colony[16:18], # WT + AZ3146 relative colony number
       HCT116.Rev2$Rel.colony[40:42],
       alternative = c("less")) # UBE2H-KO clone 2
#  0.009727

# HCT116 UBE2H-KO c3 + OE: 
t.test(HCT116.Rev2$Rel.colony[22:24], # WT + AZ3146 relative colony number
       HCT116.Rev2$Rel.colony[46:48],
       alternative = c("less")) # UBE2H-KO clone 3
#  0.03166



# plot HCT116: 
ggplot(HCT116.Rev2, aes(x= Cell.Drug, y= Rel.colony))+ 
  geom_boxplot()+
  geom_point()+  
  theme_classic()+
  ylim(0,1.3)+
  ylab("Relative colony formation\n upon 100nM Reversine")+
  xlab("Cell line")
# HCT116.Rev.UBE2HKO.OE.RelativeColony2.pdf


ggplot(subset(HCT116.Rev2, Reversine_100nM=="+") , aes(x= Cell.Drug, y= Rel.colony, fill = color))+ 
  geom_boxplot()+
  scale_fill_manual(values= c("mediumorchid2", "mediumorchid4", "turquoise4", "turquoise2"))+
  geom_point()+  
  theme_classic()+
  ylim(0,1.2)+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ylab("Relative colony formation\n upon 100nM Reversine")+
  xlab("Cell line")
# HCT116.Rev.UBE2HKO.OE.RelativeColony2.pdf






     ##### June 30 2026 ####
# 100nM Reversine
HCT116.Rev3$Cell.Drug  <-  paste0(HCT116.Rev3$CellLine, HCT116.Rev3$Reversine_100nM)

HCTMean<- mean(HCT116.Rev3$Colony[1:3])
HCTc1<- mean(HCT116.Rev3$Colony[7:9])
HCTc2<- mean(HCT116.Rev3$Colony[13:15])
HCTc3<- mean(HCT116.Rev3$Colony[19:21])
HCTOE<- mean(HCT116.Rev3$Colony[25:27])
HCTc1OE<- mean(HCT116.Rev3$Colony[31:33])
HCTc2OE<- mean(HCT116.Rev3$Colony[37:39])
HCTc3OE<- mean(HCT116.Rev3$Colony[43:45])

HCT116.Rev3$Rel.colony <-  c( HCT116.Rev3$Colony[1:6]/HCTMean, 
                              HCT116.Rev3$Colony[7:12]/HCTc1,
                              HCT116.Rev3$Colony[13:18]/HCTc2 ,
                              HCT116.Rev3$Colony[19:24]/HCTc3 , 
                              HCT116.Rev3$Colony[25:30]/HCTOE, 
                              HCT116.Rev3$Colony[31:36]/HCTc1OE,
                              HCT116.Rev3$Colony[37:42]/HCTc2OE ,
                              HCT116.Rev3$Colony[43:48]/HCTc3OE 
)

HCT116.Rev3$Cell.Drug <- factor(HCT116.Rev3$Cell.Drug, 
                                level= c("WT-", "WT+", "UBE2H-OE-", "UBE2H-OE+", 
                                         "UBE2H-KO_c1-", "UBE2H-KO_c1+","UBE2H-KO_c1+OE-", "UBE2H-KO_c1+OE+",
                                         "UBE2H-KO_c2-", "UBE2H-KO_c2+","UBE2H-KO_c2+OE-", "UBE2H-KO_c2+OE+",
                                         "UBE2H-KO_c3-", "UBE2H-KO_c3+","UBE2H-KO_c3+OE-", "UBE2H-KO_c3+OE+"))

HCT116.Rev3$color<- c("WT", "WT","WT", "WT", "WT","WT", 
                      "KO", "KO", "KO", "KO", "KO", "KO", 
                      "KO", "KO", "KO", "KO", "KO", "KO", 
                      "KO", "KO", "KO", "KO", "KO", "KO", 
                      "OE", "OE", "OE", "OE", "OE", "OE", 
                      "KO + OE","KO + OE","KO + OE",  "KO + OE","KO + OE","KO + OE", 
                      "KO + OE","KO + OE","KO + OE",  "KO + OE","KO + OE","KO + OE", 
                      "KO + OE","KO + OE","KO + OE",  "KO + OE","KO + OE","KO + OE")

# p-value
# one sided t-test: 
# is AZ3146 more toxic in UBE2H-KO clones than wt? 

# HCT116 UBE2H-KO c1: 
t.test(HCT116.Rev3$Rel.colony[4:6], # WT + AZ3146 relative colony number
       HCT116.Rev3$Rel.colony[10:12]) # UBE2H-KO clone 1
#  0.006723

# HCT116 UBE2H-KO c2: 
t.test(HCT116.Rev3$Rel.colony[4:6], # WT + AZ3146 relative colony number
       HCT116.Rev3$Rel.colony[16:18]) # UBE2H-KO clone 2
#  0.01759

# HCT116 UBE2H-KO c3: 
t.test(HCT116.Rev3$Rel.colony[4:6], # WT + AZ3146 relative colony number
       HCT116.Rev3$Rel.colony[22:24]) # UBE2H-KO clone 3
#  0.00595


# Overexpression different from WT? 
# HCT116 UBE2H-OE: 
t.test(HCT116.Rev3$Rel.colony[4:6], # WT + AZ3146 relative colony number
       HCT116.Rev3$Rel.colony[28:30]) # UBE2H-OE
#  0.02315
# HCT116 KO1 UBE2H-OE: 
t.test(HCT116.Rev3$Rel.colony[4:6], # WT + AZ3146 relative colony number
       HCT116.Rev3$Rel.colony[34:36]) # UBE2H-OE
#  0.04469
# HCT116 KO2 UBE2H-OE: 
t.test(HCT116.Rev3$Rel.colony[4:6], # WT + AZ3146 relative colony number
       HCT116.Rev3$Rel.colony[40:42]) # UBE2H-OE
#  0.09139 
# HCT116 KO3 UBE2H-OE: 
t.test(HCT116.Rev3$Rel.colony[4:6], # WT + AZ3146 relative colony number
       HCT116.Rev3$Rel.colony[46:48]) # UBE2H-OE
#  0.003273


# p-value
# one sided t-test: 
# is Rev Less toxic in UBE2H-OE clones than corresponding KO? 

# HCT116 UBE2H-KO c1 + OE: 
t.test(HCT116.Rev3$Rel.colony[10:12], # WT + AZ3146 relative colony number
       HCT116.Rev3$Rel.colony[34:36],
       alternative = c("less")) # UBE2H-KO clone 1
#  0.02028

# HCT116 UBE2H-KO c2 + OE: 
t.test(HCT116.Rev3$Rel.colony[16:18], # WT + AZ3146 relative colony number
       HCT116.Rev3$Rel.colony[40:42],
       alternative = c("less")) # UBE2H-KO clone 2
#  0.0471

# HCT116 UBE2H-KO c3 + OE: 
t.test(HCT116.Rev3$Rel.colony[22:24], # WT + AZ3146 relative colony number
       HCT116.Rev3$Rel.colony[46:48],
       alternative = c("less")) # UBE2H-KO clone 3
#  0.02571



# plot HCT116: 
ggplot(HCT116.Rev3, aes(x= Cell.Drug, y= Rel.colony))+ 
  geom_boxplot()+
  geom_point()+  
  theme_classic()+
  ylim(0,1.2)+
  ylab("Relative colony formation\n upon 100nM Reversine")+
  xlab("Cell line")
# HCT116.Rev.UBE2HKO.OE.RelativeColony3.pdf


ggplot(subset(HCT116.Rev3, Reversine_100nM=="+") , aes(x= Cell.Drug, y= Rel.colony, fill = color))+ 
  geom_boxplot()+
  scale_fill_manual(values= c("mediumorchid2", "mediumorchid4", "turquoise4", "turquoise2"))+
  geom_point()+  
  theme_classic()+
  ylim(0,1)+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ylab("Relative colony formation\n upon 100nM Reversine")+
  xlab("Cell line")
# HCT116.Rev.UBE2HKO.OE.RelativeColony3.pdf






     ##### *Merge Rev HCT116 ####
# 100nM Reversine

HCT116.Rev_merge<- rbind(HCT116.Rev2, HCT116.Rev3)

# p-value
# one sided t-test: 
# is AZ3146 more toxic in UBE2H-KO clones than wt? 

# HCT116 UBE2H-KO c1: 
t.test(subset(HCT116.Rev_merge, Cell.Drug == "WT+")$Rel.colony, # WT + AZ3146 relative colony number
       subset(HCT116.Rev_merge, Cell.Drug == "UBE2H-KO_c1+")$Rel.colony) # UBE2H-KO clone 1
#  0.0004716

# HCT116 UBE2H-KO c2: 
t.test(subset(HCT116.Rev_merge, Cell.Drug == "WT+")$Rel.colony, # WT + AZ3146 relative colony number
       subset(HCT116.Rev_merge, Cell.Drug == "UBE2H-KO_c2+")$Rel.colony)$p.value # UBE2H-KO clone 2
#  1.519568e-05

# HCT116 UBE2H-KO c3: 
t.test(subset(HCT116.Rev_merge, Cell.Drug == "WT+")$Rel.colony, # WT + AZ3146 relative colony number
       subset(HCT116.Rev_merge, Cell.Drug == "UBE2H-KO_c3+")$Rel.colony)$p.value # UBE2H-KO clone 3
#  3.580208e-05


# Over-expression different from WT? 
# HCT116 UBE2H-OE: 
t.test(subset(HCT116.Rev_merge, Cell.Drug == "WT+")$Rel.colony, # WT + AZ3146 relative colony number
       subset(HCT116.Rev_merge, Cell.Drug == "UBE2H-OE+")$Rel.colony)$p.value # UBE2H-KO clone 3
#  0.1719431

# HCT116 KO1 OE: 
t.test(subset(HCT116.Rev_merge, Cell.Drug == "WT+")$Rel.colony, # WT + AZ3146 relative colony number
       subset(HCT116.Rev_merge, Cell.Drug == "UBE2H-KO_c1+OE+")$Rel.colony)$p.value # UBE2H-KO clone 1 + OE
#  0.1955405
# HCT116 KO2 OE: 
t.test(subset(HCT116.Rev_merge, Cell.Drug == "WT+")$Rel.colony, # WT + AZ3146 relative colony number
       subset(HCT116.Rev_merge, Cell.Drug == "UBE2H-KO_c2+OE+")$Rel.colony)$p.value # UBE2H-KO clone 2 + OE
#  0.1383959
# HCT116 KO3 OE: 
t.test(subset(HCT116.Rev_merge, Cell.Drug == "WT+")$Rel.colony, # WT + AZ3146 relative colony number
       subset(HCT116.Rev_merge, Cell.Drug == "UBE2H-KO_c3+OE+")$Rel.colony)$p.value # UBE2H-KO clone 3 + OE
#  0.007195056



# p-value
# one sided t-test: 
# is Rev Less toxic in UBE2H-OE clones than corresponding KO? 

# HCT116 KO1 UBE2H-OE: 
t.test(subset(HCT116.Rev_merge, Cell.Drug == "UBE2H-KO_c1+")$Rel.colony, # KO1
       subset(HCT116.Rev_merge, Cell.Drug == "UBE2H-KO_c1+OE+")$Rel.colony,
       alternative = c("less"))$p.value # KO1 + OE
#  0.01465065
# HCT116 KO2 UBE2H-OE: 
t.test(subset(HCT116.Rev_merge, Cell.Drug == "UBE2H-KO_c2+")$Rel.colony, # KO2
       subset(HCT116.Rev_merge, Cell.Drug == "UBE2H-KO_c2+OE+")$Rel.colony,
       alternative = c("less"))$p.value # KO2 + OE
#  0.00180977
# HCT116 KO3 UBE2H-OE: 
t.test(subset(HCT116.Rev_merge, Cell.Drug == "UBE2H-KO_c3+")$Rel.colony, # KO2
       subset(HCT116.Rev_merge, Cell.Drug == "UBE2H-KO_c3+OE+")$Rel.colony,
       alternative = c("less"))$p.value # KO2 + OE
#  0.002564414





# plot HCT116: 
ggplot(HCT116.Rev_merge, aes(x= Cell.Drug, y= Rel.colony))+ 
  geom_boxplot()+
  geom_point()+  
  theme_classic()+
  #ylim(0,1.2)+
  ylab("Relative colony formation\n upon 100nM Reversine")+
  xlab("Cell line")
# HCT116.Rev.UBE2HKO.OE.RelativeColony3.pdf



ggplot(subset(HCT116.Rev_merge, Reversine_100nM=="+") , aes(x= Cell.Drug, y= Rel.colony, fill = color))+ 
  geom_boxplot()+
  scale_fill_manual(values= c("mediumorchid2", "mediumorchid4", "turquoise4", "turquoise2"))+
  geom_point()+  
  theme_classic()+
  ylim(0,1.2)+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ylab("Relative colony formation\n upon 100nM Reversine")+
  xlab("Cell line")
# HCT116.Rev.UBE2HKO.OE.RelativeColony_Merge.pdf
# Figure 5G





  #### Cal51 colony formation REVERSINE ####
     ##### June 2026 ####
## Plot Cal51 +/- UBE2H-KO +/- AZ3146 
## 100nM Reversine
# Calculate relative colony number: 

Cal51.Rev2$Cell.Drug  <-  paste0(Cal51.Rev2$CellLine, Cal51.Rev2$Reversine_100nM)
WTMean<- mean(Cal51.Rev2$Colony[c(1:3, 7:9)])
WTc2<- mean(Cal51.Rev2$Colony[13:15])
WTc7<- mean(Cal51.Rev2$Colony[19:21])

Cal51.Rev2$Rel.colony <-  c( Cal51.Rev2$Colony[1:12]/WTMean, 
                            Cal51.Rev2$Colony[13:18]/WTc2, 
                            Cal51.Rev2$Colony[19:24]/WTc7)
Cal51.Rev2$Cell.Drug <- factor(Cal51.Rev2$Cell.Drug, level= unique(Cal51.Rev2$Cell.Drug))

# p-value
# one sided t-test
# is AZ3146 more toxic in UBE2H-KO clones than wt? 
# Cal51 UBE2H-KO c2: 
t.test(Cal51.Rev2$Rel.colony[c(4:6, 10:12)], # WT + AZ3146 relative colony number
       Cal51.Rev2$Rel.colony[16:18],
       alternative = c("greater") ) # UBE2H-KO clone 2
# 0.001047

# Cal51 UBE2H-KO c7: 
t.test(Cal51.Rev2$Rel.colony[c(4:6, 10:12)], # WT + AZ3146 relative colony number
       Cal51.Rev2$Rel.colony[22:24],
       alternative = c("greater") ) # UBE2H-KO clone 
# 0.001249



# plot
ggplot(Cal51.Rev2, aes(x= Cell.Drug, y= Rel.colony))+ 
  geom_boxplot()+
  geom_point()+  
  theme_classic()+
  ylim(0,1.3)+
  ylab("Relative colony formation\n upon 100nM Reversine")+
  xlab("Cell line")
# Cal51.Rev.UBE2HKO.RelativeColony_2.pdf

Cal51.Rev2$color <- c("WT", "WT","WT", "WT", "WT","WT", 
                      "WT", "WT","WT", "WT", "WT","WT", 
                     "KO", "KO", "KO", "KO", "KO", "KO", 
                     "KO", "KO", "KO", "KO", "KO", "KO")

ggplot(subset(Cal51.Rev2, Reversine_100nM=="+") , aes(x= Cell.Drug, y= Rel.colony, fill = color))+ 
  geom_boxplot()+
  scale_fill_manual(values= c("mediumorchid2", "turquoise2"))+
  geom_point()+  
  theme_classic()+
  ylim(0, 1)+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ylab("Relative colony formation\n upon 100nM Reversine")+
  xlab("Cell line")
# Cal51.Reversine.UBE2HKO.RelativeColony_2.pdf





     ##### *July 2026 ####
## Plot Cal51 +/- UBE2H-KO +/- AZ3146 
## 100nM Reversine
# Calculate relative colony number: 

Cal51.Rev3$Cell.Drug  <-  paste0(Cal51.Rev3$CellLine, Cal51.Rev3$Reversine_100nM)
WTMean<- mean(Cal51.Rev3$Colony[c(1:3, 7:9)])
WTc2<- mean(Cal51.Rev3$Colony[c(13:15, 19:21)])
WTc7<- mean(Cal51.Rev3$Colony[c(25:27, 31:33)])

Cal51.Rev3$Rel.colony <-  c( Cal51.Rev3$Colony[1:12]/WTMean, 
                             Cal51.Rev3$Colony[13:24]/WTc2, 
                             Cal51.Rev3$Colony[25:36]/WTc7)
Cal51.Rev3$Cell.Drug <- factor(Cal51.Rev3$Cell.Drug, level= unique(Cal51.Rev3$Cell.Drug))

# p-value
# one sided t-test
# is AZ3146 more toxic in UBE2H-KO clones than wt? 
# Cal51 UBE2H-KO c2: 
t.test(Cal51.Rev3$Rel.colony[c(4:6, 10:12)], # WT + AZ3146 relative colony number
       Cal51.Rev3$Rel.colony[c(16:18, 22:24)],
       alternative = c("greater") ) # UBE2H-KO clone 2
# 3.785e-05

# Cal51 UBE2H-KO c7: 
t.test(Cal51.Rev3$Rel.colony[c(4:6, 10:12)], # WT + AZ3146 relative colony number
       Cal51.Rev3$Rel.colony[c(28:30, 34:36)],
       alternative = c("greater") ) # UBE2H-KO clone 
# 2.272e-09



# plot
ggplot(Cal51.Rev3, aes(x= Cell.Drug, y= Rel.colony))+ 
  geom_boxplot()+
  geom_point()+  
  theme_classic()+
  ylim(0,1.3)+
  ylab("Relative colony formation\n upon 100nM Reversine")+
  xlab("Cell line")
# Cal51.Rev.UBE2HKO.RelativeColony_3.pdf

Cal51.Rev3$color <- c("WT", "WT","WT", "WT", "WT","WT", 
                      "WT", "WT","WT", "WT", "WT","WT", 
                      "KO", "KO", "KO", "KO", "KO", "KO", 
                      "KO", "KO", "KO", "KO", "KO", "KO", 
                      "KO", "KO", "KO", "KO", "KO", "KO", 
                      "KO", "KO", "KO", "KO", "KO", "KO")

ggplot(subset(Cal51.Rev3, Reversine_100nM=="+") , aes(x= Cell.Drug, y= Rel.colony, fill = color))+ 
  geom_boxplot()+
  scale_fill_manual(values= c("mediumorchid2", "turquoise2"))+
  geom_point()+  
  theme_classic()+
  ylim(0, 1)+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ylab("Relative colony formation\n upon 100nM Reversine")+
  xlab("Cell line")
# Cal51.Reversine.UBE2HKO.RelativeColony_3.pdf
# Figure 5I

# Similar results in both runs. Reproducible! 


### HCT116 colony formation BAY-1217389 ####
     ##### *July 2027 ####
# 10nM BAY-1217389
HCT116.Bay$Cell.Drug  <-  paste0(HCT116.Bay$CellLine, HCT116.Bay$BAY1218389_10nM)

HCTMean<- mean(HCT116.Bay$Colony[1:3])
HCTc1<- mean(HCT116.Bay$Colony[7:9])
HCTc2<- mean(HCT116.Bay$Colony[13:15])
HCTc3<- mean(HCT116.Bay$Colony[19:21])
HCTOE<- mean(HCT116.Bay$Colony[25:27])
HCTc1OE<- mean(HCT116.Bay$Colony[31:33])
HCTc2OE<- mean(HCT116.Bay$Colony[37:39])
HCTc3OE<- mean(HCT116.Bay$Colony[43:45])

HCT116.Bay$Rel.colony <-  c( HCT116.Bay$Colony[1:6]/HCTMean, 
                             HCT116.Bay$Colony[7:12]/HCTc1,
                             HCT116.Bay$Colony[13:18]/HCTc2 ,
                             HCT116.Bay$Colony[19:24]/HCTc3 , 
                             HCT116.Bay$Colony[25:30]/HCTOE, 
                             HCT116.Bay$Colony[31:36]/HCTc1OE,
                             HCT116.Bay$Colony[37:42]/HCTc2OE ,
                             HCT116.Bay$Colony[43:48]/HCTc3OE 
)

HCT116.Bay$Cell.Drug <- factor(HCT116.Bay$Cell.Drug, 
                               level= c("WT-", "WT+", "UBE2H-OE-", "UBE2H-OE+", 
                                        "UBE2H-KO_c1-", "UBE2H-KO_c1+","UBE2H-KO_c1+OE-", "UBE2H-KO_c1+OE+",
                                        "UBE2H-KO_c2-", "UBE2H-KO_c2+","UBE2H-KO_c2+OE-", "UBE2H-KO_c2+OE+",
                                        "UBE2H-KO_c3-", "UBE2H-KO_c3+","UBE2H-KO_c3+OE-", "UBE2H-KO_c3+OE+"))

HCT116.Bay$color<- c("WT", "WT","WT", "WT", "WT","WT", 
                     "KO", "KO", "KO", "KO", "KO", "KO", 
                     "KO", "KO", "KO", "KO", "KO", "KO", 
                     "KO", "KO", "KO", "KO", "KO", "KO", 
                     "OE", "OE", "OE", "OE", "OE", "OE", 
                     "KO + OE","KO + OE","KO + OE",  "KO + OE","KO + OE","KO + OE", 
                     "KO + OE","KO + OE","KO + OE",  "KO + OE","KO + OE","KO + OE", 
                     "KO + OE","KO + OE","KO + OE",  "KO + OE","KO + OE","KO + OE")

# p-value
# two sided t-test: 
# is AZ3146 more toxic in UBE2H-KO clones than wt? 

# HCT116 UBE2H-KO c1: 
t.test(HCT116.Bay$Rel.colony[4:6], # WT + AZ3146 relative colony number
       HCT116.Bay$Rel.colony[10:12]) # UBE2H-KO clone 1
#  0.007229

# HCT116 UBE2H-KO c2: 
t.test(HCT116.Bay$Rel.colony[4:6], # WT + AZ3146 relative colony number
       HCT116.Bay$Rel.colony[16:18]) # UBE2H-KO clone 2
#  0.007122

# HCT116 UBE2H-KO c3: 
t.test(HCT116.Bay$Rel.colony[4:6], # WT + AZ3146 relative colony number
       HCT116.Bay$Rel.colony[22:24]) # UBE2H-KO clone 3
#  0.002779


# Overexpression different from WT? 
# HCT116 UBE2H-OE: 
t.test(HCT116.Bay$Rel.colony[4:6], # WT + AZ3146 relative colony number
       HCT116.Bay$Rel.colony[28:30]) # UBE2H-OE
#  0.2264 NS
# HCT116 KO1 UBE2H-OE: 
t.test(HCT116.Bay$Rel.colony[4:6], # WT + AZ3146 relative colony number
       HCT116.Bay$Rel.colony[34:36]) # UBE2H-OE
#  0.06417 NS
# HCT116 KO2 UBE2H-OE: 
t.test(HCT116.Bay$Rel.colony[4:6], # WT + AZ3146 relative colony number
       HCT116.Bay$Rel.colony[40:42]) # UBE2H-OE
#  0.03819 Lower?
# HCT116 KO3 UBE2H-OE: 
t.test(HCT116.Bay$Rel.colony[4:6], # WT + AZ3146 relative colony number
       HCT116.Bay$Rel.colony[46:48]) # UBE2H-OE
#  0.01594 Lower


# p-value
# one sided t-test: 
# is BAY Less toxic in UBE2H-OE clones than corresponding KO? 

# HCT116 UBE2H-KO c1 + OE: 
t.test(HCT116.Bay$Rel.colony[10:12], # WT + AZ3146 relative colony number
       HCT116.Bay$Rel.colony[34:36],
       alternative = c("less")) # UBE2H-KO clone 1
#  0.01291 Rescue

# HCT116 UBE2H-KO c2 + OE: 
t.test(HCT116.Bay$Rel.colony[16:18], # WT + AZ3146 relative colony number
       HCT116.Bay$Rel.colony[40:42],
       alternative = c("less")) # UBE2H-KO clone 2
#  0.0002357 Rescue

# HCT116 UBE2H-KO c3 + OE: 
t.test(HCT116.Bay$Rel.colony[22:24], # WT + AZ3146 relative colony number
       HCT116.Bay$Rel.colony[46:48],
       alternative = c("less")) # UBE2H-KO clone 3
#  0.04485 Rescue



# plot HCT116: 
ggplot(HCT116.Bay, aes(x= Cell.Drug, y= Rel.colony))+ 
  geom_boxplot()+
  geom_point()+  
  theme_classic()+
  ylim(0,1)+
  ylab("Relative colony formation\n upon 10nM BAY-1217389")+
  xlab("Cell line")
# HCT116.BAY.UBE2HKO.OE.RelativeColony.pdf


ggplot(subset(HCT116.Bay, BAY1218389_10nM=="+") , aes(x= Cell.Drug, y= Rel.colony, fill = color))+ 
  geom_boxplot()+
  scale_fill_manual(values= c("mediumorchid2", "mediumorchid4", "turquoise4", "turquoise2"))+
  geom_point()+  
  theme_classic()+
  ylim(0,1)+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ylab("Relative colony formation\n upon 10nM BAY-1217389")+
  xlab("Cell line")
# HCT116.BAY.UBE2HKO.OE.RelativeColony.pdf
# Supplementary Figure S6D





### Cell death Data ####

# cell death data from try 1, 25 August, 2025
# 48 hours of AZ3146

# colony formation data from October 

# Plot relative colony count compared to cell death (annexin V and DRAQ7) after 48 hours AZ treatment 

Colony.Death<- data.frame(CellLine = c("Cal51 WT", "Cal51 UBE2H-KO c2", "Cal51 UBE2H-KO c7", 
                                       "HCT116", "HCT116 UBE2H-KO c1", "HCT116 UBE2H-KO c2", "HCT116 UBE2H-KO c3"), 
                          UBE2H = c("WT", "KO", "KO",
                                    "WT", "KO", "KO", "KO"),
                          CellType = c("Cal51", "Cal51", "Cal51",
                                       "HCT116", "HCT116", "HCT116", "HCT116"),
                          DiffCellDeath = c(8.91, 7.16, 12.9, 
                                            6.78, 9.72, 6.7, 7.15), 
                          FCcelldeath = c(1.55, 1.43, 1.94, 
                                          1.45, 1.92, 1.54, 1.98), 
                          FCAnnexinV = c(2.11, 1.92, 3.46, 
                                         1.76, 1.93, 1.57, 2.23), 
                          FCDraq7 = c(0.83, 1.16, 0.98 , 
                                      1.4, 2.35, 1.16, 0.94)
)



  #### Data from  October 2025: single well each####
# 96 hours of AZ3146

# Plot relative colony count compared to cell death (annexin V and DRAQ7) after 48 hours AZ treatment 

Colony.Death2<- data.frame(CellLine = c("Cal51 WT", "Cal51 UBE2H-KO c2", "Cal51 UBE2H-KO c7", 
                                        "HCT116", "HCT116 UBE2H-KO c1", "HCT116 UBE2H-KO c2", "HCT116 UBE2H-KO c3"), 
                           UBE2H = c("WT", "KO", "KO",
                                     "WT", "KO", "KO", "KO"),
                           CellType = c("Cal51", "Cal51", "Cal51",
                                        "HCT116", "HCT116", "HCT116", "HCT116"),
                           DiffCellDeath = c(36.08, 39.22, 51.97, 
                                             2.78, 15.74, 11.61, 3.14), 
                           FCcelldeath = c(6.78, 7.73, 5.72,
                                           1.4, 1.67, 1.44, 1.15), 
                           DiffAnnexinV = c(28.22, 28.46, 41.05,
                                            0.73, 5.84, -3.87, 0.2), 
                           DiffDraq7 = c(27.19, 30.71, 44.47,
                                         2.96, 14.14, 12.12, 2.79)
)



### Cell death stain plot ###
# Data from October 20, 2025
# 96 hour AZ3146
# 1 uM AZ3146

t.test(c(36.08, 36.08, 2.78, 2.78,2.78), 
       c(39.22, 51.97, 15.74, 11.61, 3.14), 
       paired=TRUE) # p= 0.04729 , 96 hour analysis 


# Paired test HCT116 & Cal51 
t.test(c(8.91,8.91,6.78, 6.78, 6.78), 
       c(7.16, 12.9, 9.72, 6.7, 7.15), 
       paired=TRUE) # p= 0.3537 , 48 hour analysis 

t.test(c(8.91,8.91,6.78, 6.78, 6.78, 36.08, 36.08, 2.78, 2.78,2.78), 
       c(7.16, 12.9, 9.72, 6.7, 7.15, 39.22, 51.97, 15.74, 11.61, 3.14), 
       paired=TRUE)
# p=0.035, 48 and 96 hour cell death increase, paired data 


  #### *Add cell death staining October 2025 triplicate ####

# cell death data from 10 October 2025
# 96 hours of 1uM AZ3146
# HCT116 and Cal51 in triplicate



### Paired Cell death stain plot ###
# Data from October 20, 2025
# 96 hour 1 uM AZ3146

Colony.Death<- data.frame(CellLine = c("Cal51 WT", "Cal51 WT", "Cal51 WT", 
                                       "Cal51 UBE2H-KO c2","Cal51 UBE2H-KO c2","Cal51 UBE2H-KO c2", 
                                       "Cal51 UBE2H-KO c7", "Cal51 UBE2H-KO c7", "Cal51 UBE2H-KO c7", 
                                       "HCT116", "HCT116", "HCT116", 
                                       "HCT116 UBE2H-KO c1","HCT116 UBE2H-KO c1","HCT116 UBE2H-KO c1", 
                                       "HCT116 UBE2H-KO c2", "HCT116 UBE2H-KO c2", "HCT116 UBE2H-KO c2", 
                                       "HCT116 UBE2H-KO c3", "HCT116 UBE2H-KO c3", "HCT116 UBE2H-KO c3"), 
                          Clone = c("WT", "WT", "WT", 
                                    "KO2", "KO2", "KO2", 
                                    "KO7", "KO7", "KO7",
                                    "WT", "WT", "WT", 
                                    "KO1", "KO1", "KO1",
                                    "KO2", "KO2", "KO2",
                                    "KO3", "KO3", "KO3"),
                          UBE2H = c("WT", "WT", "WT", 
                                    "KO", "KO", "KO", 
                                    "KO", "KO", "KO",
                                    "WT", "WT", "WT", 
                                    "KO", "KO", "KO",
                                    "KO", "KO", "KO",
                                    "KO", "KO", "KO"),
                          CellType = c("Cal51", "Cal51", "Cal51",
                                       "Cal51", "Cal51", "Cal51",
                                       "Cal51", "Cal51", "Cal51",
                                       "HCT116", "HCT116", "HCT116", 
                                       "HCT116", "HCT116", "HCT116",
                                       "HCT116", "HCT116", "HCT116",
                                       "HCT116", "HCT116", "HCT116"),
                          DiffCellDeath = c(57.60, 51.74, 59.75,
                                            61.16,58.07,57.06,
                                            69.11,71.31,73.72,  
                                            5.41, 8.01, 6.08, 
                                            36.32,36.97,35.34,
                                            33.30,23.57,35.11,
                                            18.10,12.16,25.30)
)
Colony.Death$UBE2H<- factor(Colony.Death$Clone, level= c("WT", "KO1","KO2","KO3","KO7"))
#Colony.Death$UBE2H<- factor(Colony.Death$UBE2H, level= c("WT", "KO"))

ggplot(subset(Colony.Death, CellType=="HCT116"), aes(x=UBE2H, y= DiffCellDeath, fill = UBE2H))+
  geom_boxplot()+
  geom_point()+
  theme_classic()+
  scale_fill_manual(values= c("turquoise2", "mediumorchid2","mediumorchid2","mediumorchid2"))+
  ylim(0,40)+
  xlab("HCT116")+
  ylab("Relative cell death upon AZ3146 \n(Annexin V + DRAQ 7)")
# boxplot.CellDeath.HCT116_Oct10.pdf

ggplot(subset(Colony.Death, CellType=="Cal51"), aes(x=UBE2H, y= DiffCellDeath, fill = UBE2H))+
  geom_boxplot()+  
  geom_point()+
  theme_classic()+
  scale_fill_manual(values= c("turquoise2", "mediumorchid2","mediumorchid2"))+
  ylim(0,75)+
  xlab("Cal51")+
  ylab("Increased cell death upon AZ3146 \n(Annexin V + DRAQ 7)")
# boxplot.CellDeath.Cal51_Oct10.pdf
# Figure 5K

t.test(c(5.41, 8.01, 6.08, 
         5.41, 8.01, 6.08, 
         5.41, 8.01, 6.08, 
         57.60, 51.74, 59.75, 
         57.60, 51.74, 59.75), 
       c(36.32,36.97,35.34,
         33.30,23.57,35.11,
         18.10,12.16,25.30,
         61.16,58.07,57.06,
         69.11, 71.31, 73.72), 
       paired=TRUE) # p= 3.781e-05 , 96 hour analysis 


# Paired test HCT116
t.test(c(5.41, 8.01, 6.08), 
       c(36.32,36.97,35.34,
         33.30,23.57,35.11,
         18.10,12.16,25.30), 
       alternative = c("less")) # p= 3.1e-05 , 96 hour analysis 
t.test(c(5.41, 8.01, 6.08), 
       c(36.32,36.97,35.34), 
       alternative = c("less")) # p= 1.5e-05 , 96 hour analysis 
t.test(c(5.41, 8.01, 6.08), 
       c(33.30,23.57,35.11), 
       alternative = c("less")) # p= 0.00881 , 96 hour analysis 
t.test(c(5.41, 8.01, 6.08), 
       c(18.10,12.16,25.30), 
       alternative = c("less")) # p= 0.04074 , 96 hour analysis 

# Does UBE2H-KO increase cell death upon AZ3146? 
# Paired test Cal51
t.test(c(57.60, 51.74, 59.75), 
       c(61.16,58.07,57.06,
         69.11, 71.31, 73.72), 
       alternative = c("less")) # p= 0.02869 , 96 hour analysis 
# Paired test Cal51
t.test(c(57.60, 51.74, 59.75), 
       c(61.16,58.07,57.06), 
       alternative = c("less")) # p= 0.2193 , 96 hour analysis 
# Paired test Cal51
t.test(c(57.60, 51.74, 59.75), 
       c(69.11, 71.31, 73.72), 
       alternative = c("less")) # p= 0.005337 , 96 hour analysis 




  #### *Cell death stain February 2026 (HCT116) #####
## Paired Cell death stain plot ##
# Data from Febuary 1, 2026
# HCT116 +/- KO   96 hour 1 uM AZ3146

Colony.Death.Feb26<- data.frame(CellLine = c("HCT116", "HCT116", "HCT116", 
                                             "HCT116 UBE2H-KO c1","HCT116 UBE2H-KO c1","HCT116 UBE2H-KO c1", 
                                             "HCT116 UBE2H-KO c2", "HCT116 UBE2H-KO c2", "HCT116 UBE2H-KO c2", 
                                             "HCT116 UBE2H-KO c3", "HCT116 UBE2H-KO c3", "HCT116 UBE2H-KO c3"), 
                                Clone = c("WT", "WT", "WT", 
                                          "KO1", "KO1", "KO1",
                                          "KO2", "KO2", "KO2",
                                          "KO3", "KO3", "KO3"),
                                UBE2H = c("WT", "WT", "WT", 
                                          "KO", "KO", "KO",
                                          "KO", "KO", "KO",
                                          "KO", "KO", "KO"),
                                CellType = c("HCT116", "HCT116", "HCT116", 
                                             "HCT116", "HCT116", "HCT116",
                                             "HCT116", "HCT116", "HCT116",
                                             "HCT116", "HCT116", "HCT116"),
                                DiffCellDeath = c(34.1, 13.7, 15.6, 
                                                  50.9, 53.2, 53.3,
                                                  49.7,47.4,  51.9,
                                                  52.3, 45.8, 51.5) 
)
Colony.Death.Feb26$Clone<- factor(Colony.Death.Feb26$Clone, level= c("WT", "KO1","KO2","KO3"))


ggplot(subset(Colony.Death.Feb26, CellType=="HCT116"), aes(x=Clone, y= DiffCellDeath, fill = Clone))+
  geom_boxplot()+
  geom_point()+
  theme_classic()+
  scale_fill_manual(values= c("turquoise2", "mediumorchid2","mediumorchid2","mediumorchid2"))+
  ylim(0,60)+
  xlab("HCT116")+
  ylab("Relative cell death upon AZ3146 \n(Annexin V + DRAQ 7)")
# boxplot.CellDeath.HCT116_Feb1.pdf
# Figure 5J

t.test(c(34.1, 13.7, 15.6), 
       c(50.9, 53.2, 53.3,
         49.7,47.4,  51.9,
         52.3, 45.8, 51.5)) # p= 3.781e-05 , 96 hour analysis 


# Paired test HCT116
t.test(c(34.1, 13.7, 15.6), 
       c(50.9, 53.2, 53.3,
         49.7,47.4,  51.9,
         52.3, 45.8, 51.5), 
       alternative = c("less")) # p= 0.02156, 96 hour analysis all 3 clones
t.test(c(34.1, 13.7, 15.6), 
       c(50.9, 53.2, 53.3), 
       alternative = c("less")) # p= 0.01941, 96 hour analysis clone 1
t.test(c(34.1, 13.7, 15.6), 
       c(49.7,47.4,  51.9), 
       alternative = c("less")) # p= 0.0218, 96 hour analysis clone 2
t.test(c(34.1, 13.7, 15.6), 
       c(52.3, 45.8, 51.5), 
       alternative = c("less")) # p= 0.01885, 96 hour analysis clone 3




### Counting IF nuclei (Multinuclei) #####
# Counted a minimum of 1,000 cells per condition, and a minimum of 3 images per condition
# Images immunofluorescent fixed cells. looked at Hoechst and MitoTracker to identify multi-nucleated cells. 
# did not look at UBE2H antibody staining as this did not aid in cell border identification. 
# counted manually. as FIJI was not great at identifying multi-nucleated cells. 

      ##### HCT116 KO.OE ####
HCT116_IF<- data.frame(cell_line= c("WT", "WT_AZ", "KO1", "KO1_AZ", "KO2", "KO2_AZ", "KO3", "KO3_AZ", 
                                    "OE", "OE_AZ", "KO1.OE", "KO1.OE_AZ", "KO2.OE", "KO2.OE_AZ", "KO3.OE", "KO3.OE_AZ"), 
                       color= c("WT", "WT", "KO", "KO","KO", "KO","KO", "KO", "OE", "OE", "KO.OE", "KO.OE", "KO.OE", "KO.OE","KO.OE", "KO.OE"), 
                      total_Cell= c(1874, 1374, 1247, 1264, 1325, 1228, 1608, 1022, 1463, 1097, 1413, 1264, 1317, 1012, 1491, 1230), 
                      Micronuclei= c(130, 129, 36, 95, 24, 57, 35, 38,  112, 75, 48, 55, 73, 159, 31, 45), 
                      MultiNuclei= c(16, 11, 101, 212, 191, 293, 87, 108,  20, 8, 51, 107, 50, 100, 38, 74))

HCT116_IF$PercentMicro<- HCT116_IF$Micronuclei/HCT116_IF$total_Cell
HCT116_IF$PercentMultiNuclei<- HCT116_IF$MultiNuclei/HCT116_IF$total_Cell
HCT116_IF$cell_line<- factor(HCT116_IF$cell_line, levels = c("WT","WT_AZ", "OE", "OE_AZ", "KO1", "KO1_AZ",  "KO1.OE", "KO1.OE_AZ",
                                                             "KO2", "KO2_AZ", "KO2.OE", "KO2.OE_AZ",  
                                                             "KO3", "KO3_AZ", "KO3.OE", "KO3.OE_AZ"))

ggplot(HCT116_IF, aes(x= cell_line, y= PercentMultiNuclei*100, fill = color))+ 
  geom_bar(stat="identity")+
  scale_fill_manual(values= c("mediumorchid2", "mediumorchid4", "turquoise4", "turquoise2"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ylab("Cells Multi-Nucleited (%)")+
  xlab("")
# HCT116_IF_Multinuclei.pdf

#Chi.squared tests: 
HCT116_IF2<- HCT116_IF[,c("cell_line", "total_Cell", "MultiNuclei")]
chisq.test(subset(HCT116_IF2, cell_line %in% c("WT", "WT_AZ"))[,2:3]) # NS
chisq.test(subset(HCT116_IF2, cell_line %in% c("OE", "OE_AZ"))[,2:3]) # NS
chisq.test(subset(HCT116_IF2, cell_line %in% c("KO1", "KO1_AZ"))[,2:3]) # 8.966e-09
chisq.test(subset(HCT116_IF2, cell_line %in% c("KO1.OE", "KO1.OE_AZ"))[,2:3]) # 8.23e-07
chisq.test(subset(HCT116_IF2, cell_line %in% c("KO2", "KO2_AZ"))[,2:3]) # 6.775e-07
chisq.test(subset(HCT116_IF2, cell_line %in% c("KO2.OE", "KO2.OE_AZ"))[,2:3]) # 4.86e-08
chisq.test(subset(HCT116_IF2, cell_line %in% c("KO3", "KO3_AZ"))[,2:3]) # 7.851e-06
chisq.test(subset(HCT116_IF2, cell_line %in% c("KO3.OE", "KO3.OE_AZ"))[,2:3]) # 2.18e-05

chisq.test(subset(HCT116_IF2, cell_line %in% c("WT", "OE"))[,2:3]) # NS
chisq.test(subset(HCT116_IF2, cell_line %in% c("WT_AZ", "OE_AZ"))[,2:3]) # NS
chisq.test(subset(HCT116_IF2, cell_line %in% c("WT", "KO1"))[,2:3])$p.value # 4.428432e-23
chisq.test(subset(HCT116_IF2, cell_line %in% c("WT_AZ", "KO1_AZ"))[,2:3])$p.value # 2.698333e-41
chisq.test(subset(HCT116_IF2, cell_line %in% c("WT", "KO2"))[,2:3])$p.value # 9.842668e-46
chisq.test(subset(HCT116_IF2, cell_line %in% c("WT_AZ", "KO2_AZ"))[,2:3])$p.value # 8.892641e-59
chisq.test(subset(HCT116_IF2, cell_line %in% c("WT", "KO3"))[,2:3])$p.value # 3.73458e-14
chisq.test(subset(HCT116_IF2, cell_line %in% c("WT_AZ", "KO3_AZ"))[,2:3])$p.value # 1.943171e-24

chisq.test(subset(HCT116_IF2, cell_line %in% c("KO1", "KO1.OE"))[,2:3])$p.value # 3.964315e-06
chisq.test(subset(HCT116_IF2, cell_line %in% c("KO1_AZ", "KO1.OE_AZ"))[,2:3])$p.value # 4.157063e-08
chisq.test(subset(HCT116_IF2, cell_line %in% c("KO2", "KO2.OE"))[,2:3])$p.value # 8.406814e-18
chisq.test(subset(HCT116_IF2, cell_line %in% c("KO2_AZ", "KO2.OE_AZ"))[,2:3])$p.value # 4.174842e-13
chisq.test(subset(HCT116_IF2, cell_line %in% c("KO3", "KO3.OE"))[,2:3])$p.value # 0.0001465414
chisq.test(subset(HCT116_IF2, cell_line %in% c("KO3_AZ", "KO3.OE_AZ"))[,2:3])$p.value # 0.0003783084



      ##### Cal51 KO AZ ####
### Cal51 
# (Cal51 KO2 done on a different day from WT and KO7, less growth? too full?)

Cal51_IF<- data.frame(cell_line= c("WT", "WT_AZ",  "KO2", "KO2_AZ",   "KO7", "KO7_AZ"), 
                       color= c("WT", "WT", "KO", "KO",  "KO", "KO"), 
                       total_Cell= c(1373, 1609, 1958, 2126, 1442, 1104), 
                       Micronuclei= c(68, 497,   49, 309,    53, 336), 
                       MultiNuclei= c(13, 453,   33, 305,    11, 395))

Cal51_IF$PercentMicro<- Cal51_IF$Micronuclei/Cal51_IF$total_Cell
Cal51_IF$PercentMultiNuclei<- Cal51_IF$MultiNuclei/Cal51_IF$total_Cell
Cal51_IF$cell_line<- factor(Cal51_IF$cell_line, levels = c("WT","WT_AZ", "KO2", "KO2_AZ", "KO7", "KO7_AZ"))

ggplot(Cal51_IF, aes(x= cell_line, y= PercentMultiNuclei*100, fill = color))+ 
  geom_bar(stat="identity")+
  scale_fill_manual(values= c("mediumorchid2", "turquoise2"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ylab("Cells Multi-Nucleited (%)")+
  xlab("")
# Cal51_IF_Multinuclei.pdf

#Chi.squared tests: 
# untreated vs AZ treated
Cal51_IFt<- Cal51_IF[,c("cell_line", "total_Cell", "MultiNuclei")]
chisq.test(subset(Cal51_IFt, cell_line %in% c("WT", "WT_AZ"))[,2:3])$p.value # 8.617979e-70
chisq.test(subset(Cal51_IFt, cell_line %in% c("KO7", "KO7_AZ"))[,2:3])$p.value # 3.882236e-90
chisq.test(subset(Cal51_IFt, cell_line %in% c("KO2", "KO2_AZ"))[,2:3])$p.value # 1.526289e-41

# WT AZ vs KO AZ
chisq.test(subset(Cal51_IFt, cell_line %in% c("WT", "KO7"))[,2:3])$p.value # NS
chisq.test(subset(Cal51_IFt, cell_line %in% c("WT_AZ", "KO7_AZ"))[,2:3])$p.value # 0.002779954 KO greater
chisq.test(subset(Cal51_IFt, cell_line %in% c("WT", "KO2"))[,2:3])$p.value # NS 
chisq.test(subset(Cal51_IFt, cell_line %in% c("WT_AZ", "KO2_AZ"))[,2:3])$p.value # 6.047612e-17. KO less




## Repeat Cal51, this time all cells grown at same time. 100k/well 
# and added a lower concentration of AZ3146, as 1uM was very deadly to all cells inccluding controls
### Cal51 
# Update August 20 with 0.25uM AZ3146 as well 

Cal51_IF2<- data.frame(cell_line= c("WT", "WT_AZ_.25", "WT_AZ_1", "KO2", "KO2_AZ_.25","KO2_AZ_1",  "KO7", "KO7_AZ_.25","KO7_AZ_1"), 
                      color= c("WT", "WT", "WT", "KO", "KO", "KO",  "KO", "KO", "KO"), 
                      total_Cell= c(1484,2872,1249,   1661,2482,1366,   1758,1715,1046), 
                      Micronuclei= c(46,215,462,    60,205,426,      41, 208,272), 
                      MultiNuclei= c(71,229,378,    146,269,423,     75,173,417))

Cal51_IF2$PercentMicro<- Cal51_IF2$Micronuclei/Cal51_IF2$total_Cell
Cal51_IF2$PercentMultiNuclei<- Cal51_IF2$MultiNuclei/Cal51_IF2$total_Cell
Cal51_IF2$cell_line<- factor(Cal51_IF2$cell_line, levels = c("WT","KO2","KO7", "WT_AZ_.25", "KO2_AZ_.25","KO7_AZ_.25", "WT_AZ_1","KO2_AZ_1", "KO7_AZ_1"))

ggplot(Cal51_IF2, aes(x= cell_line, y= PercentMultiNuclei*100, fill = color))+ 
  geom_bar(stat="identity")+
  scale_fill_manual(values= c("mediumorchid2",  "turquoise2"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ylab("Cells Multi-Nucleited (%)")+
  xlab("")
# Cal51_IF_Multinuclei_AZlowHigh.pdf


ggplot(subset(Cal51_IF2,cell_line %in% c("WT","KO2","KO7", "WT_AZ_.25", "KO2_AZ_.25","KO7_AZ_.25") ), 
       aes(x= cell_line, y= PercentMultiNuclei*100, fill = color))+ 
  geom_bar(stat="identity")+
  scale_fill_manual(values= c("mediumorchid2",  "turquoise2"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ylab("Cells Multi-Nucleited (%)")+
  xlab("")
# Cal51_IF_Multinuclei_AZ0.25.pdf


Cal51_IF2$cell_line<- factor(Cal51_IF2$cell_line, levels = c("WT", "WT_AZ_.25", "WT_AZ_1",  "KO2", "KO2_AZ_.25","KO2_AZ_1", "KO7", "KO7_AZ_.25", "KO7_AZ_1" ))
ggplot(subset(Cal51_IF2,cell_line %in% c("WT","KO2","KO7", "WT_AZ_.25", "KO2_AZ_.25","KO7_AZ_.25", "WT_AZ_1","KO2_AZ_1", "KO7_AZ_1") ), 
       aes(x= cell_line, y= PercentMultiNuclei*100, fill = color))+ 
  geom_bar(stat="identity")+
  scale_fill_manual(values= c("mediumorchid2",  "turquoise2"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ylab("Cells Multi-Nucleited (%)")+
  xlab("")
# Cal51_IF_Multinuclei_AZ0.25.pdf
# Cal51_IF_Multinuclei_AZlowHigh_WT.KO2.Ko7.pdf


#Chi.squared tests: 
# Untreated vs 0.25 AZ treated
chisq.test(subset(Cal51_IF2t, cell_line %in% c("WT", "WT_AZ_.25"))[,2:3])$p.value # 0.0002818403
chisq.test(subset(Cal51_IF2t, cell_line %in% c("KO2", "KO2_AZ_.25"))[,2:3])$p.value # 0.05779116 NS
chisq.test(subset(Cal51_IF2t, cell_line %in% c("KO7", "KO7_AZ_.25"))[,2:3])$p.value # 8.496384e-10

# Untreated vs 1uM AZ treated
Cal51_IF2t<- Cal51_IF2[,c("cell_line", "total_Cell", "MultiNuclei")]
chisq.test(subset(Cal51_IF2t, cell_line %in% c("WT", "WT_AZ_1"))[,2:3])$p.value # 2.601916e-51
chisq.test(subset(Cal51_IF2t, cell_line %in% c("KO2", "KO2_AZ_1"))[,2:3])$p.value # 3.476235e-37
chisq.test(subset(Cal51_IF2t, cell_line %in% c("KO7", "KO7_AZ_1"))[,2:3])$p.value # 1.320824e-84


# WT  vs KO 
chisq.test(subset(Cal51_IF2t, cell_line %in% c("WT", "KO2"))[,2:3])$p.value # 4.830624e-05 
chisq.test(subset(Cal51_IF2t, cell_line %in% c("WT", "KO7"))[,2:3])$p.value # NS p= 0.553499

# WT AZ vs KO 0.25 AZ
chisq.test(subset(Cal51_IF2t, cell_line %in% c("WT_AZ_.25", "KO2_AZ_.25"))[,2:3])$p.value # 0.001245466
chisq.test(subset(Cal51_IF2t, cell_line %in% c("WT_AZ_.25", "KO7_AZ_.25"))[,2:3])$p.value # 0.02891357 
chisq.test(subset(Cal51_IF2t, cell_line %in% c("WT_AZ_1", "KO2_AZ_1"))[,2:3])$p.value # 0.8079125 NS
chisq.test(subset(Cal51_IF2t, cell_line %in% c("WT_AZ_1", "KO7_AZ_1"))[,2:3])$p.value # 0.0009500675 KO greater




      ##### HCT116 KO.OE single cell clones ####
HCT116_IFclones<- data.frame(cell_line= c("WT", "WT_AZ", "KO1", "KO1_AZ", "KO2", "KO2_AZ", "KO3", "KO3_AZ", "OE", "OE_AZ", 
                                    "KO1.OEc1", "KO1.OEc1_AZ", "KO2.OEc1", "KO2.OEc1_AZ","KO2.OEc2", "KO2.OEc2_AZ", "KO3.OEc2", "KO3.OEc2_AZ", "KO3.OEc5", "KO3.OEc5_AZ"), 
                       color= c("WT", "WT", "KO", "KO","KO", "KO","KO", "KO", "OE", "OE", 
                                "KO.OE", "KO.OE", "KO.OE", "KO.OE","KO.OE", "KO.OE", "KO.OE", "KO.OE", "KO.OE", "KO.OE"), 
                       total_Cell= c(1874, 1374, 1247, 1264, 1325, 1228, 1608, 1022, 1463, 1097, 
                                     2025, 1500, 1427, 1363, 1302,1120, 1448,1443, 1014, 1035), 
                       Micronuclei= c(130, 129, 36, 95, 24, 57, 35, 38,  112, 75, 
                                      134, 89, 49, 35, 61, 99, 41, 68, 64, 60), 
                       MultiNuclei= c(16, 11, 101, 212, 191, 293, 87, 108,  20, 8, 
                                      89, 69, 26, 20, 32, 52, 31, 33, 39, 43))

HCT116_IFclones$PercentMicro<- HCT116_IFclones$Micronuclei/HCT116_IFclones$total_Cell
HCT116_IFclones$PercentMultiNuclei<- HCT116_IFclones$MultiNuclei/HCT116_IFclones$total_Cell
HCT116_IFclones$cell_line<- factor(HCT116_IFclones$cell_line, levels = c("WT","WT_AZ", "OE", "OE_AZ", 
                                                                         "KO1", "KO1_AZ",  "KO1.OEc1", "KO1.OEc1_AZ",
                                                             "KO2", "KO2_AZ", "KO2.OEc1", "KO2.OEc1_AZ",  "KO2.OEc2", "KO2.OEc2_AZ",  
                                                             "KO3", "KO3_AZ", "KO3.OEc2", "KO3.OEc2_AZ", "KO3.OEc5", "KO3.OEc5_AZ"))

ggplot(HCT116_IFclones, aes(x= cell_line, y= PercentMultiNuclei*100, fill = color))+ 
  geom_bar(stat="identity")+
  scale_fill_manual(values= c("mediumorchid2", "mediumorchid4", "turquoise4", "turquoise2"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ylab("Cells Multi-Nucleited (%)")+
  xlab("")
# HCT116_IF_clones_Multinuclei.pdf


#Chi.squared tests: 
HCT116_IFclones2<- HCT116_IFclones[,c("cell_line", "total_Cell", "MultiNuclei")]
chisq.test(subset(HCT116_IFclones2, cell_line %in% c("WT", "WT_AZ"))[,2:3]) # NS
chisq.test(subset(HCT116_IFclones2, cell_line %in% c("OE", "OE_AZ"))[,2:3]) # NS
chisq.test(subset(HCT116_IFclones2, cell_line %in% c("KO1", "KO1_AZ"))[,2:3]) # 8.966e-09
chisq.test(subset(HCT116_IFclones2, cell_line %in% c("KO1.OEc1", "KO1.OEc1_AZ"))[,2:3]) # NS
chisq.test(subset(HCT116_IFclones2, cell_line %in% c("KO2", "KO2_AZ"))[,2:3]) # 6.775e-07
chisq.test(subset(HCT116_IFclones2, cell_line %in% c("KO2.OEc1", "KO2.OEc1_AZ"))[,2:3]) # NS
chisq.test(subset(HCT116_IFclones2, cell_line %in% c("KO2.OEc2", "KO2.OEc2_AZ"))[,2:3]) # 0.006586
chisq.test(subset(HCT116_IFclones2, cell_line %in% c("KO3", "KO3_AZ"))[,2:3]) # 7.851e-06
chisq.test(subset(HCT116_IFclones2, cell_line %in% c("KO3.OEc2", "KO3.OEc2_AZ"))[,2:3]) # NS
chisq.test(subset(HCT116_IFclones2, cell_line %in% c("KO3.OEc5", "KO3.OEc5_AZ"))[,2:3]) # NS


chisq.test(subset(HCT116_IFclones2, cell_line %in% c("WT", "OE"))[,2:3]) # NS
chisq.test(subset(HCT116_IFclones2, cell_line %in% c("WT_AZ", "OE_AZ"))[,2:3]) # NS
chisq.test(subset(HCT116_IFclones2, cell_line %in% c("WT", "KO1"))[,2:3])$p.value # 4.428432e-23
chisq.test(subset(HCT116_IFclones2, cell_line %in% c("WT_AZ", "KO1_AZ"))[,2:3])$p.value # 2.698333e-41
chisq.test(subset(HCT116_IFclones2, cell_line %in% c("WT", "KO2"))[,2:3])$p.value # 9.842668e-46
chisq.test(subset(HCT116_IFclones2, cell_line %in% c("WT_AZ", "KO2_AZ"))[,2:3])$p.value # 8.892641e-59
chisq.test(subset(HCT116_IFclones2, cell_line %in% c("WT", "KO3"))[,2:3])$p.value # 3.73458e-14
chisq.test(subset(HCT116_IFclones2, cell_line %in% c("WT_AZ", "KO3_AZ"))[,2:3])$p.value # 1.943171e-24


chisq.test(subset(HCT116_IFclones2, cell_line %in% c("KO1", "KO1.OEc1"))[,2:3])$p.value # 4.936559e-05
chisq.test(subset(HCT116_IFclones2, cell_line %in% c("KO1_AZ", "KO1.OEc1_AZ"))[,2:3])$p.value # 3.9868e-21

chisq.test(subset(HCT116_IFclones2, cell_line %in% c("KO2", "KO2.OEc1"))[,2:3])$p.value # 2.549834e-29
chisq.test(subset(HCT116_IFclones2, cell_line %in% c("KO2_AZ", "KO2.OEc1_AZ"))[,2:3])$p.value # 1.54198e-53
chisq.test(subset(HCT116_IFclones2, cell_line %in% c("KO2", "KO2.OEc2"))[,2:3])$p.value # 9.411171e-24
chisq.test(subset(HCT116_IFclones2, cell_line %in% c("KO2_AZ", "KO2.OEc2_AZ"))[,2:3])$p.value # 6.942252e-30

chisq.test(subset(HCT116_IFclones2, cell_line %in% c("KO3", "KO3.OEc2"))[,2:3])$p.value #  1.001525e-05
chisq.test(subset(HCT116_IFclones2, cell_line %in% c("KO3_AZ", "KO3.OEc2_AZ"))[,2:3])$p.value # 5.476886e-16
chisq.test(subset(HCT116_IFclones2, cell_line %in% c("KO3", "KO3.OEc5"))[,2:3])$p.value # 0.09947432
chisq.test(subset(HCT116_IFclones2, cell_line %in% c("KO3_AZ", "KO3.OEc5_AZ"))[,2:3])$p.value # 3.439721e-07




### Seahorse analysis #####
setwd(SchukkenData)
SH.HCT.UBE2H<- read_excel("UBE2H.Result.Overview_7.24.xlsx")

SH.HCT.UBE2H2<- read_excel("UBE2H.Result.Overview_8.28.xlsx")

SH.HCT.UBE2H3<- read_excel("UBE2H.Result.Overview_9.08.xlsx")


    #### HCT116 run 1 ####

SH.HCT.UBE2H$AZ3146<- factor(SH.HCT.UBE2H$AZ3146, levels = c(0,1))
SH.HCT.UBE2H$`Group Name`<- factor(SH.HCT.UBE2H$`Group Name`, levels = unique(SH.HCT.UBE2H$`Group Name`))



# AZ3146 decreases Basal Respiration
ggplot(SH.HCT.UBE2H, aes(x=AZ3146, y=`Basal Respiration`, color= `Group Name`))+  
  geom_point()+  
  scale_color_manual(values= c("turquoise2","turquoise2", "turquoise4",  "turquoise4", 
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4",
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4",
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4"))+
  stat_summary(fun = "mean", size= 0.3, color= "red", geom = "crossbar")+
  ylab("Basal Respiration")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
# plot.Seahorse.AZ3146.BasalResp_H1.pdf
t.test(as.numeric(subset(SH.HCT.UBE2H, AZ3146 == "0")$`Basal Respiration`),
       as.numeric(subset(SH.HCT.UBE2H, AZ3146 == "1")$`Basal Respiration`)  ,
       paired=TRUE) # p= 0.003712, difference = -0.1206501



# AZ3146 decreases ATP production
ggplot(SH.HCT.UBE2H, aes(x=AZ3146, y=`ATP Production`, color= `Group Name`))+  
  geom_point()+  
  scale_color_manual(values= c("turquoise2","turquoise2", "turquoise4",  "turquoise4", 
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4",
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4",
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4"))+
  stat_summary(fun = "mean", size= 0.3, color= "red", geom = "crossbar")+
  ylab("ATP Production")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
# plot.Seahorse.AZ3146.ATP_H1.pdf
# AZ3146 (paired all)
t.test(as.numeric(subset(SH.HCT.UBE2H, AZ3146 == "0")$`ATP Production`),
       as.numeric(subset(SH.HCT.UBE2H, AZ3146 == "1")$`ATP Production`)  ,
       paired=TRUE) # p= 0.000622, difference = -0.131879

t.test(as.numeric(subset(SH.HCT.UBE2H, `Group Name` == "H.WT")$`ATP Production`),
       as.numeric(subset(SH.HCT.UBE2H, `Group Name` == "H.OE")$`ATP Production`)) # p=0.1038
t.test(as.numeric(subset(SH.HCT.UBE2H, `Group Name` == "H.WT")$`ATP Production`),
       as.numeric(subset(SH.HCT.UBE2H, `Group Name` == "H.KO1")$`ATP Production`) ) # p=0.008197
t.test(as.numeric(subset(SH.HCT.UBE2H, `Group Name` == "H.WT")$`ATP Production`),
       as.numeric(subset(SH.HCT.UBE2H, `Group Name` == "H.KO2")$`ATP Production`) ) # p=4.985e-06
t.test(as.numeric(subset(SH.HCT.UBE2H, `Group Name` == "H.WT")$`ATP Production`),
       as.numeric(subset(SH.HCT.UBE2H, `Group Name` == "H.KO3")$`ATP Production`) ) # p=0.0001717
t.test(as.numeric(subset(SH.HCT.UBE2H, `Group Name` == "H.KO1")$`ATP Production`),
       as.numeric(subset(SH.HCT.UBE2H, `Group Name` == "H.KO1.OE")$`ATP Production`) ) # p=0.03479
t.test(as.numeric(subset(SH.HCT.UBE2H, `Group Name` == "H.KO2")$`ATP Production`),
       as.numeric(subset(SH.HCT.UBE2H, `Group Name` == "H.KO2.OE")$`ATP Production`) ) # p=0.05225
t.test(as.numeric(subset(SH.HCT.UBE2H, `Group Name` == "H.KO3")$`ATP Production`),
       as.numeric(subset(SH.HCT.UBE2H, `Group Name` == "H.KO3.OE")$`ATP Production`) ) # p=0.4921




# AZ3146 increases Spare Respiratory Capacity
ggplot(SH.HCT.UBE2H, aes(x=AZ3146, y=`Spare Respiratory Capacity`, color= `Group Name`))+  
  geom_point()+  
  scale_color_manual(values= c("turquoise2","turquoise2", "turquoise4",  "turquoise4", 
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4",
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4",
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4"))+
  stat_summary(fun = "mean", size= 0.3, color= "red", geom = "crossbar")+
  ylab("Spare Respiratory Capacity")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
# plot.Seahorse.AZ3146.SRC_H1.pdf
t.test(as.numeric(subset(SH.HCT.UBE2H, AZ3146 == "0")$`Spare Respiratory Capacity`),
       as.numeric(subset(SH.HCT.UBE2H, AZ3146 == "1")$`Spare Respiratory Capacity`)  ,
       paired=TRUE) # p= 0.001403, difference = 0.3159841



# AZ3146 increases Non-Mitochondrial Oxygen Consumption
ggplot(SH.HCT.UBE2H, aes(x=AZ3146, y=`Non-Mitochondrial Oxygen Consumption`, color= `Group Name`))+  
  geom_point()+  
  scale_color_manual(values= c("turquoise2","turquoise2", "turquoise4",  "turquoise4", 
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4",
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4",
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4"))+
  stat_summary(fun = "mean", size= 0.3, color= "red", geom = "crossbar")+
  ylab("Non-Mitochondrial Oxygen Consumption")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
# plot.Seahorse.AZ3146.NMOC_H1.pdf
t.test(as.numeric(subset(SH.HCT.UBE2H, AZ3146 == "0")$`Non-Mitochondrial Oxygen Consumption`),
       as.numeric(subset(SH.HCT.UBE2H, AZ3146 == "1")$`Non-Mitochondrial Oxygen Consumption`)  ,
       paired=TRUE) # p= 0.01422, difference = -0.0668481 (eh. not that significant) 




# AZ3146  Maximal Respiration NS
ggplot(SH.HCT.UBE2H, aes(x=AZ3146, y=`Maximal Respiration`, color= `Group Name`))+  
  geom_point()+  
  scale_color_manual(values= c("turquoise2","turquoise2", "turquoise4",  "turquoise4", 
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4",
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4",
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4"))+
  stat_summary(fun = "mean", size= 0.3, color= "red", geom = "crossbar")+
  ylab("Maximal Respiration")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
# plot.Seahorse.AZ3146.Max.pdf
t.test(as.numeric(subset(SH.HCT.UBE2H, AZ3146 == "0")$`Maximal Respiration`),
       as.numeric(subset(SH.HCT.UBE2H, AZ3146 == "1")$`Maximal Respiration`)  ,
       paired=TRUE) # p= 0.08816, difference = -0.195334 (NS) 



###
# Basal Respiration (UBE2H decreases)
ggplot(SH.HCT.UBE2H, aes(x=`Group Name`, y=`Basal Respiration`, fill=`Group Name`, shape= AZ3146 ))+
  geom_boxplot()+
  scale_fill_manual(values= c("turquoise2","turquoise2", "turquoise4",  "turquoise4", 
                              "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4",
                              "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4",
                              "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4"))+
  geom_point()+   
  ylab("Basal Respiration")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
#plot.Seahorse.HCT116.BaselResp_H1.pdf

# Just look at AZ3146 (not AZ): Basal Respiration: 
ggplot(subset(SH.HCT.UBE2H, AZ3146=="0"), 
       aes(x=`Group Name`, y=`Basal Respiration`, fill=`Group Name` ))+
  geom_boxplot()+
  scale_fill_manual(values= c("turquoise2", "turquoise4", 
                              "mediumorchid2", "mediumorchid4", "mediumorchid2", "mediumorchid4","mediumorchid2", "mediumorchid4"))+
  geom_point()+   
  ylab("Basal Respiration")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
# 5x4
#plot.Seahorse.HCT116.noAZ.BaselResp_H1.pdf
t.test(as.numeric(subset(SH.HCT.UBE2H, `Group Name` == "H.WT")$`Basal Respiration`),
       as.numeric(subset(SH.HCT.UBE2H, `Group Name` == "H.OE")$`Basal Respiration`)) # p=0.5601
t.test(as.numeric(subset(SH.HCT.UBE2H, `Group Name` == "H.WT")$`Basal Respiration`),
       as.numeric(subset(SH.HCT.UBE2H, `Group Name` == "H.KO1")$`Basal Respiration`) ) # p=0.0001639
t.test(as.numeric(subset(SH.HCT.UBE2H, `Group Name` == "H.WT")$`Basal Respiration`),
       as.numeric(subset(SH.HCT.UBE2H, `Group Name` == "H.KO2")$`Basal Respiration`) ) # p=6.211e-07
t.test(as.numeric(subset(SH.HCT.UBE2H, `Group Name` == "H.WT")$`Basal Respiration`),
       as.numeric(subset(SH.HCT.UBE2H, `Group Name` == "H.KO3")$`Basal Respiration`) ) # p=6.556e-05
t.test(as.numeric(subset(SH.HCT.UBE2H, `Group Name` == "H.KO1")$`Basal Respiration`),
       as.numeric(subset(SH.HCT.UBE2H, `Group Name` == "H.KO1.OE")$`Basal Respiration`) ) # p=0.0008237
t.test(as.numeric(subset(SH.HCT.UBE2H, `Group Name` == "H.KO2")$`Basal Respiration`),
       as.numeric(subset(SH.HCT.UBE2H, `Group Name` == "H.KO2.OE")$`Basal Respiration`) ) # p=0.01436
t.test(as.numeric(subset(SH.HCT.UBE2H, `Group Name` == "H.KO3")$`Basal Respiration`),
       as.numeric(subset(SH.HCT.UBE2H, `Group Name` == "H.KO3.OE")$`Basal Respiration`) ) # p=0.3579


# ATP Production (***) (UBE2H decreases)
ggplot(SH.HCT.UBE2H, aes(x=`Group Name`, y=SH.HCT.UBE2H$`ATP Production`, fill=`Group Name`, shape= AZ3146 ))+
  geom_boxplot()+
  scale_fill_manual(values= c("turquoise2","turquoise2", "turquoise4",  "turquoise4", 
                              "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4",
                              "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4",
                              "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4"))+
  geom_point()+   
  ylab("ATP Production")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
#plot.Seahorse.HCT116.ATP_H1.pdf

# Just look at AZ3146 (not AZ): ATP production: 
ggplot(subset(SH.HCT.UBE2H, AZ3146=="0"), 
       aes(x=`Group Name`, y=`ATP Production`, fill=`Group Name`))+
  geom_boxplot()+
  scale_fill_manual(values= c("turquoise2", "turquoise4", 
                              "mediumorchid2", "mediumorchid4", "mediumorchid2", "mediumorchid4","mediumorchid2", "mediumorchid4"))+
  geom_point()+   
  ylab("ATP Production")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
#plot.Seahorse.HCT116.noAZ.ATP_H1.pdf



# Non-Mitochondrial Oxygen Consumption *** !! (UBE2H decreases,  AZ increases)
ggplot(SH.HCT.UBE2H, aes(x=`Group Name`, y=SH.HCT.UBE2H$`Non-Mitochondrial Oxygen Consumption`, fill=`Group Name`, shape= AZ3146 ))+
  geom_boxplot()+
  scale_fill_manual(values= c("turquoise2","turquoise2", "turquoise4",  "turquoise4", 
                              "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4",
                              "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4",
                              "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4"))+
  geom_point()+   
  ylab("Non-Mitochondrial\nOxygen Consumption")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
#plot.Seahorse.HCT116.NMOC_H1.pdf

# Just look at non-AZ3146 (not AZ): Non-Mitochondrial Oxygen Consumption
ggplot(subset(SH.HCT.UBE2H, AZ3146=="0"), 
       aes(x=`Group Name`, y=`Non-Mitochondrial Oxygen Consumption`, fill=`Group Name`))+
  geom_boxplot()+
  scale_fill_manual(values= c("turquoise2", "turquoise4", 
                              "mediumorchid2", "mediumorchid4", "mediumorchid2", "mediumorchid4","mediumorchid2", "mediumorchid4"))+
  geom_point()+   
  ylab("Non-Mitochondrial\nOxygen Consumption")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
#plot.Seahorse.HCT116.noAZ.NMOC_H1.pdf




# Spare Respiratory Capacity (AZ3146 increases SRC)
ggplot(SH.HCT.UBE2H, aes(x=`Group Name`, y=SH.HCT.UBE2H$`Spare Respiratory Capacity`, fill=`Group Name`, shape= AZ3146 ))+
  geom_boxplot()+
  scale_fill_manual(values= c("turquoise2","turquoise2", "turquoise4",  "turquoise4", 
                              "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4",
                              "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4",
                              "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4"))+
  geom_point()+   
  ylab("Spare Respiratory Capacity")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))

# Spare Respiratory Capacity as a % of basal. (maybe? but I think this is just reflection of basal) 
ggplot(SH.HCT.UBE2H, aes(x=`Group Name`, y=SH.HCT.UBE2H$`Spare Respiratory Capacity as a %`, fill=`Group Name`, shape= AZ3146 ))+
  geom_boxplot()+
  scale_fill_manual(values= c("turquoise2","turquoise2", "turquoise4",  "turquoise4", 
                              "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4",
                              "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4",
                              "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4"))+
  geom_point()+   
  ylab("Spare Respiratory Capacity as a % of basal")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))


# Maximal Respiration (eh.. OE and KO reduce max)
ggplot(SH.HCT.UBE2H, aes(x=`Group Name`, y=SH.HCT.UBE2H$`Maximal Respiration`, fill=`Group Name`, shape= AZ3146 ))+
  geom_boxplot()+
  scale_fill_manual(values= c("turquoise2","turquoise2", "turquoise4",  "turquoise4", 
                              "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4",
                              "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4",
                              "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4"))+
  geom_point()+   
  ylab("Maximal Respiration")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
#plot.Seahorse.HCT116.MAX.pdf
ggplot(subset(SH.HCT.UBE2H, AZ3146=="0"), 
       aes(x=`Group Name`, y=`Maximal Respiration`, fill=`Group Name`))+
  geom_boxplot()+
  scale_fill_manual(values= c("turquoise2", "turquoise4", 
                              "mediumorchid2", "mediumorchid4", "mediumorchid2", "mediumorchid4","mediumorchid2", "mediumorchid4"))+
  geom_point()+   
  ylab("Maximal Respiration")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
#plot.Seahorse.HCT116.noAZ.MAX.pdf

# Proton Leak. NS (OE and KO reduce leak, inconsistent)
ggplot(SH.HCT.UBE2H, aes(x=`Group Name`, y=SH.HCT.UBE2H$`Proton Leak`, fill=`Group Name`, shape= AZ3146 ))+
  geom_boxplot()+
  scale_fill_manual(values= c("turquoise2","turquoise2", "turquoise4",  "turquoise4", 
                              "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4",
                              "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4",
                              "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4"))+
  geom_point()+   
  ylab("Proton Leak")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))




    #### HCT116 run 2 ####

SH.HCT.UBE2H2$AZ3146<- factor(SH.HCT.UBE2H2$AZ3146, levels = c(0,1))
SH.HCT.UBE2H2$`Group Name`<- factor(SH.HCT.UBE2H2$`Group Name`, levels = unique(SH.HCT.UBE2H2$`Group Name`))


# AZ3146 decreases Basal Respiration
ggplot(SH.HCT.UBE2H2, aes(x=AZ3146, y=`Basal Respiration`, color= `Group Name`))+  
  geom_point()+  
  scale_color_manual(values= c("turquoise2","turquoise2", "turquoise4",  "turquoise4", 
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4",
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4","mediumorchid4", "mediumorchid4",
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4", "mediumorchid4", "mediumorchid4"))+
  stat_summary(fun = "mean", size= 0.3, color= "red", geom = "crossbar")+
  ylab("Basal Respiration")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
# plot.Seahorse.AZ3146.BasalResp_H2.pdf
t.test(as.numeric(subset(SH.HCT.UBE2H2, AZ3146 == "0")$`Basal Respiration`),
       as.numeric(subset(SH.HCT.UBE2H2, AZ3146 == "1")$`Basal Respiration`)  ,
       paired=TRUE) # p= 0.001034, difference = -0.2334825



# AZ3146 decreases ATP production
ggplot(SH.HCT.UBE2H2, aes(x=AZ3146, y=`ATP Production`, color= `Group Name`))+  
  geom_point()+  
  scale_color_manual(values= c("turquoise2","turquoise2", "turquoise4",  "turquoise4", 
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4",
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4","mediumorchid4", "mediumorchid4",
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4", "mediumorchid4", "mediumorchid4"))+
  stat_summary(fun = "mean", size= 0.3, color= "red", geom = "crossbar")+
  ylab("ATP Production")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
# plot.Seahorse.AZ3146.ATP_H2.pdf
# AZ3146 (paired all)
t.test(as.numeric(subset(SH.HCT.UBE2H2, AZ3146 == "0")$`ATP Production`),
       as.numeric(subset(SH.HCT.UBE2H2, AZ3146 == "1")$`ATP Production`)  ,
       paired=TRUE) # p=0.002995 , difference = -0.2191349


# AZ3146 increases Spare Respiratory Capacity
ggplot(SH.HCT.UBE2H2, aes(x=AZ3146, y=`Spare Respiratory Capacity`, color= `Group Name`))+  
  geom_point()+  
  scale_color_manual(values= c("turquoise2","turquoise2", "turquoise4",  "turquoise4", 
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4",
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4","mediumorchid4", "mediumorchid4",
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4", "mediumorchid4", "mediumorchid4"))+
  stat_summary(fun = "mean", size= 0.3, color= "red", geom = "crossbar")+
  ylab("Spare Respiratory Capacity")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
# plot.Seahorse.AZ3146.SRC_H2.pdf
t.test(as.numeric(subset(SH.HCT.UBE2H2, AZ3146 == "0")$`Spare Respiratory Capacity`),
       as.numeric(subset(SH.HCT.UBE2H2, AZ3146 == "1")$`Spare Respiratory Capacity`)  ,
       paired=TRUE) # p=0.176 , difference = 0.2778039



# AZ3146  Non-Mitochondrial Oxygen Consumption (NS)
ggplot(SH.HCT.UBE2H2, aes(x=AZ3146, y=`Non-Mitochondrial Oxygen Consumption`, color= `Group Name`))+  
  geom_point()+  
  scale_color_manual(values= c("turquoise2","turquoise2", "turquoise4",  "turquoise4", 
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4",
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4","mediumorchid4", "mediumorchid4",
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4", "mediumorchid4", "mediumorchid4"))+
  stat_summary(fun = "mean", size= 0.3, color= "red", geom = "crossbar")+
  ylab("Non-Mitochondrial Oxygen Consumption")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
# plot.Seahorse.AZ3146.NMOC_H2.pdf
t.test(as.numeric(subset(SH.HCT.UBE2H2, AZ3146 == "0")$`Non-Mitochondrial Oxygen Consumption`),
       as.numeric(subset(SH.HCT.UBE2H2, AZ3146 == "1")$`Non-Mitochondrial Oxygen Consumption`)  ,
       paired=TRUE) # p=0.8586 , difference = 0.0056




# AZ3146  Maximal Respiration NS
ggplot(SH.HCT.UBE2H2, aes(x=AZ3146, y=`Maximal Respiration`, color= `Group Name`))+  
  geom_point()+  
  scale_color_manual(values= c("turquoise2","turquoise2", "turquoise4",  "turquoise4", 
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4",
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4","mediumorchid4", "mediumorchid4",
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4", "mediumorchid4", "mediumorchid4"))+
  stat_summary(fun = "mean", size= 0.3, color= "red", geom = "crossbar")+
  ylab("Maximal Respiration")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
# plot.Seahorse.AZ3146.Max.pdf
t.test(as.numeric(subset(SH.HCT.UBE2H2, AZ3146 == "0")$`Maximal Respiration`),
       as.numeric(subset(SH.HCT.UBE2H2, AZ3146 == "1")$`Maximal Respiration`)  ,
       paired=TRUE) # p= 0.8473, difference = -0.044



###
# Basal Respiration (UBE2H decreases)
ggplot(SH.HCT.UBE2H2, aes(x=`Group Name`, y=`Basal Respiration`, fill=`Group Name`, shape= AZ3146 ))+
  geom_boxplot()+
  scale_fill_manual(values= c("turquoise2","turquoise2", "turquoise4",  "turquoise4", 
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4",
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4","mediumorchid4", "mediumorchid4",
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4", "mediumorchid4", "mediumorchid4"))+
  geom_point()+   
  ylab("Basal Respiration")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
#plot.Seahorse.HCT116.BaselResp_H2.pdf
# Just look at AZ3146 (not AZ): Basal Respiration: 

ggplot(subset(SH.HCT.UBE2H2, AZ3146=="0"), 
       aes(x=`Group Name`, y=`Basal Respiration`, fill=`Group Name` ))+
  geom_boxplot()+
  scale_fill_manual(values= c("turquoise2", "turquoise4", 
                              "mediumorchid2", "mediumorchid4", 
                              "mediumorchid2", "mediumorchid4","mediumorchid4",
                              "mediumorchid2", "mediumorchid4", "mediumorchid4"))+
  geom_point()+   
  ylab("Basal Respiration")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
# 5x4
#plot.Seahorse.HCT116.noAZ.BaselResp_H2.pdf

# Difference in Basal: 
ggplot(subset(SH.HCT.UBE2H2, AZ3146=="1"), 
       aes(x=`Group Name`, y=log2(`AZ.Basal.Diff`), fill=`Group Name` ))+
  geom_boxplot()+
  scale_fill_manual(values= c("turquoise2", "turquoise4", 
                              "mediumorchid2", "mediumorchid4", 
                              "mediumorchid2", "mediumorchid4","mediumorchid4",
                              "mediumorchid2", "mediumorchid4", "mediumorchid4"))+
  geom_point()+   
  ylab("log2 FC (Basal Respiration upon AZ)")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
# 5x4
#plot.Seahorse.HCT116.ATP.Diff.pdf
# UBE2H OE reduced AZ/untreated Basal metabolism 
#

t.test(as.numeric(subset(SH.HCT.UBE2H2, `Group Name` == "WT")$`Basal Respiration`),
       as.numeric(subset(SH.HCT.UBE2H2, `Group Name` == "OE")$`Basal Respiration`)) # p=0.09704 NS
t.test(as.numeric(subset(SH.HCT.UBE2H2, `Group Name` == "WT")$`Basal Respiration`),
       as.numeric(subset(SH.HCT.UBE2H2, `Group Name` == "KO1")$`Basal Respiration`) ) # p=0.01844    decrease
t.test(as.numeric(subset(SH.HCT.UBE2H2, `Group Name` == "WT")$`Basal Respiration`),
       as.numeric(subset(SH.HCT.UBE2H2, `Group Name` == "KO2")$`Basal Respiration`) ) # p=3.292e-05  decrease
t.test(as.numeric(subset(SH.HCT.UBE2H2, `Group Name` == "WT")$`Basal Respiration`),
       as.numeric(subset(SH.HCT.UBE2H2, `Group Name` == "KO3")$`Basal Respiration`) ) # p=0.001545   decrease
t.test(as.numeric(subset(SH.HCT.UBE2H2, `Group Name` == "KO1")$`Basal Respiration`),
       as.numeric(subset(SH.HCT.UBE2H2, `Group Name` == "KO1 OE1")$`Basal Respiration`) ) # p= NS
t.test(as.numeric(subset(SH.HCT.UBE2H2, `Group Name` == "KO2")$`Basal Respiration`),
       as.numeric(subset(SH.HCT.UBE2H2, `Group Name` == "KO2 OE1")$`Basal Respiration`) ) # p= NS
t.test(as.numeric(subset(SH.HCT.UBE2H2, `Group Name` == "KO2")$`Basal Respiration`),
       as.numeric(subset(SH.HCT.UBE2H2, `Group Name` == "KO2 OE2")$`Basal Respiration`) ) # p= NS
t.test(as.numeric(subset(SH.HCT.UBE2H2, `Group Name` == "KO3")$`Basal Respiration`),
       as.numeric(subset(SH.HCT.UBE2H2, `Group Name` == "KO3 OE2")$`Basal Respiration`) ) # p=0.04931 decrease....
t.test(as.numeric(subset(SH.HCT.UBE2H2, `Group Name` == "KO3")$`Basal Respiration`), 
       as.numeric(subset(SH.HCT.UBE2H2, `Group Name` == "KO3 OE5")$`Basal Respiration`) ) # p=0.02875 decrease....


# ATP Production (UBE2H decreases)
ggplot(SH.HCT.UBE2H2, aes(x=`Group Name`, y=`ATP Production`, fill=`Group Name`, shape= AZ3146 ))+
  geom_boxplot()+
  scale_fill_manual(values= c("turquoise2","turquoise2", "turquoise4",  "turquoise4", 
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4",
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4","mediumorchid4", "mediumorchid4",
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4", "mediumorchid4", "mediumorchid4"))+
  geom_point()+   
  ylab("ATP Production")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
#plot.Seahorse.HCT116.ATP_H2.pdf

# Just look at AZ3146 (not AZ): ATP production: 
ggplot(subset(SH.HCT.UBE2H2, AZ3146=="0"), 
       aes(x=`Group Name`, y=`ATP Production`, fill=`Group Name`))+
  geom_boxplot()+
  scale_fill_manual(values= c("turquoise2", "turquoise4", 
                              "mediumorchid2", "mediumorchid4", 
                              "mediumorchid2", "mediumorchid4","mediumorchid4",
                              "mediumorchid2", "mediumorchid4","mediumorchid4"))+
  geom_point()+   
  ylab("ATP Production")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
#plot.Seahorse.HCT116.noAZ.ATP_H2.pdf

# Difference in ATP: 
ggplot(subset(SH.HCT.UBE2H3, AZ3146=="1"), 
       aes(x=`Group Name`, y=log2(`AZ.ATP.Diff`), fill=`Group Name` ))+
  geom_boxplot()+
  scale_fill_manual(values= c("turquoise2", "turquoise4", 
                              "mediumorchid2", "mediumorchid4", 
                              "mediumorchid2", "mediumorchid4","mediumorchid4",
                              "mediumorchid2", "mediumorchid4", "mediumorchid4"))+
  geom_point()+   
  ylab("Basal Respiration")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
# 5x4
#plot.Seahorse.HCT116.ATP.Diff.pdf

t.test(as.numeric(subset(SH.HCT.UBE2H2, `Group Name` == "WT")$`ATP Production`),
       as.numeric(subset(SH.HCT.UBE2H2, `Group Name` == "OE")$`ATP Production`)) # p= NS
t.test(as.numeric(subset(SH.HCT.UBE2H2, `Group Name` == "WT")$`ATP Production`),
       as.numeric(subset(SH.HCT.UBE2H2, `Group Name` == "KO1")$`ATP Production`) ) # p=NS
t.test(as.numeric(subset(SH.HCT.UBE2H2, `Group Name` == "WT")$`ATP Production`),
       as.numeric(subset(SH.HCT.UBE2H2, `Group Name` == "KO2")$`ATP Production`) ) # p=0.0001944
t.test(as.numeric(subset(SH.HCT.UBE2H2, `Group Name` == "WT")$`ATP Production`),
       as.numeric(subset(SH.HCT.UBE2H2, `Group Name` == "KO3")$`ATP Production`) ) # p=0.04056
t.test(as.numeric(subset(SH.HCT.UBE2H2, `Group Name` == "KO1")$`ATP Production`),
       as.numeric(subset(SH.HCT.UBE2H2, `Group Name` == "KO1 OE1")$`ATP Production`) ) # p=
t.test(as.numeric(subset(SH.HCT.UBE2H2, `Group Name` == "KO2")$`ATP Production`),
       as.numeric(subset(SH.HCT.UBE2H2, `Group Name` == "KO2 OE1")$`ATP Production`) ) # p=
t.test(as.numeric(subset(SH.HCT.UBE2H2, `Group Name` == "KO2")$`ATP Production`),
       as.numeric(subset(SH.HCT.UBE2H2, `Group Name` == "KO2 OE2")$`ATP Production`) ) # p=
t.test(as.numeric(subset(SH.HCT.UBE2H2, `Group Name` == "KO3")$`ATP Production`),
       as.numeric(subset(SH.HCT.UBE2H2, `Group Name` == "KO3 OE2")$`ATP Production`) ) # p=
t.test(as.numeric(subset(SH.HCT.UBE2H2, `Group Name` == "KO3")$`ATP Production`),
       as.numeric(subset(SH.HCT.UBE2H2, `Group Name` == "KO3 OE5")$`ATP Production`) ) # p=




# Non-Mitochondrial Oxygen Consumption NS (UBE2H decreases,  AZ increases)
ggplot(SH.HCT.UBE2H3, aes(x=`Group Name`, y=`Non-Mitochondrial Oxygen Consumption`, fill=`Group Name`, shape= AZ3146 ))+
  geom_boxplot()+
  scale_fill_manual(values= c("turquoise2","turquoise2", "turquoise4",  "turquoise4", 
                              "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4",
                              "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4","mediumorchid4", "mediumorchid4",
                              "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4", "mediumorchid4", "mediumorchid4"))+
  geom_point()+   
  ylab("Non-Mitochondrial\nOxygen Consumption")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
#plot.Seahorse.HCT116.NMOC.pdf

# Just look at AZ3146 (not AZ): Non-Mitochondrial Oxygen Consumption
ggplot(subset(SH.HCT.UBE2H2, AZ3146=="1"), 
       aes(x=`Group Name`, y=`Non-Mitochondrial Oxygen Consumption`, fill=`Group Name`))+
  geom_boxplot()+
  scale_fill_manual(values= c("turquoise2", "turquoise4", 
                              "mediumorchid2", "mediumorchid4", 
                              "mediumorchid2", "mediumorchid4", "mediumorchid4",
                              "mediumorchid2",  "mediumorchid4","mediumorchid4"))+
  geom_point()+   
  ylab("Non-Mitochondrial\nOxygen Consumption")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
#plot.Seahorse.HCT116.noAZ.NMOC_H2.pdf




# Spare Respiratory Capacity (AZ3146 increases SRC)
ggplot(SH.HCT.UBE2H2, aes(x=`Group Name`, y=`Spare Respiratory Capacity`, 
                          fill=`Group Name`, shape= AZ3146 ))+
  geom_boxplot()+
  scale_fill_manual(values= c("turquoise2","turquoise2", "turquoise4",  "turquoise4", 
                              "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4",
                              "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4","mediumorchid4", "mediumorchid4",
                              "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4", "mediumorchid4", "mediumorchid4"))+
  geom_point()+   
  ylab("Spare Respiratory Capacity")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
# Spare Respiratory Capacity (AZ3146 increases SRC)
ggplot(subset(SH.HCT.UBE2H2, AZ3146=="1"),  
       aes(x=`Group Name`, y=`Spare Respiratory Capacity`, 
                          fill=`Group Name`, shape= AZ3146 ))+
  geom_boxplot()+
  scale_fill_manual(values= c("turquoise2", "turquoise4",  
                              "mediumorchid2", "mediumorchid4",
                              "mediumorchid2","mediumorchid4", "mediumorchid4",
                              "mediumorchid2", "mediumorchid4", "mediumorchid4"))+
  geom_point()+   
  ylab("Spare Respiratory Capacity")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
#plot.Seahorse.HCT116.noAZ.SRC_H2.pdf

# Maximal Respiration (eh.. OE and KO increase max)
ggplot(SH.HCT.UBE2H2, aes(x=`Group Name`, y=`Maximal Respiration`, fill=`Group Name`, shape= AZ3146 ))+
  geom_boxplot()+
  scale_fill_manual(values= c("turquoise2","turquoise2", "turquoise4",  "turquoise4", 
                              "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4",
                              "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4","mediumorchid4", "mediumorchid4",
                              "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4", "mediumorchid4", "mediumorchid4"))+
  geom_point()+   
  ylab("Maximal Respiration")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
#plot.Seahorse.HCT116.MAX.pdf
ggplot(subset(SH.HCT.UBE2H2, AZ3146=="0"), 
       aes(x=`Group Name`, y=`Maximal Respiration`, fill=`Group Name`))+
  geom_boxplot()+
  scale_fill_manual(values= c("turquoise2", "turquoise4", 
                              "mediumorchid2", "mediumorchid4", 
                              "mediumorchid2", "mediumorchid4","mediumorchid4",
                              "mediumorchid2", "mediumorchid4","mediumorchid4"))+
  geom_point()+   
  ylab("Maximal Respiration")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
#plot.Seahorse.HCT116.noAZ.MAX.pdf



    #### HCT116 run 3 ####

SH.HCT.UBE2H3$AZ3146<- factor(SH.HCT.UBE2H3$AZ3146, levels = c(0,1))
SH.HCT.UBE2H3$`Group Name`<- factor(SH.HCT.UBE2H3$`Group Name`, levels = unique(SH.HCT.UBE2H3$`Group Name`))

# SH.HCT.UBE2H3<- subset(SH.HCT.UBE2H3, CellCoverage == "High")

# AZ3146 decreases Basal Respiration
ggplot(SH.HCT.UBE2H3, aes(x=AZ3146, y=`Basal Respiration`, color= `Group Name`))+  
  geom_point()+  
  scale_color_manual(values= c("turquoise2","turquoise2", "turquoise4",  "turquoise4", 
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4",
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4","mediumorchid4", "mediumorchid4",
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4", "mediumorchid4", "mediumorchid4"))+
  stat_summary(fun = "mean", size= 0.3, color= "red", geom = "crossbar")+
  ylab("Basal Respiration")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
# plot.Seahorse.AZ3146.BasalResp_H3.pdf
t.test(as.numeric(subset(SH.HCT.UBE2H3, AZ3146 == 0)$`Basal Respiration`),
       as.numeric(subset(SH.HCT.UBE2H3, AZ3146 == 1)$`Basal Respiration`)  ,
       paired=TRUE) # p= 0.007474, difference = -0.254



# AZ3146 decreases ATP production
ggplot(SH.HCT.UBE2H3, aes(x=AZ3146, y=`ATP Production`, color= `Group Name`))+  
  geom_point()+  
  scale_color_manual(values= c("turquoise2","turquoise2", "turquoise4",  "turquoise4", 
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4",
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4","mediumorchid4", "mediumorchid4",
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4", "mediumorchid4", "mediumorchid4"))+
  stat_summary(fun = "mean", size= 0.3, color= "red", geom = "crossbar")+
  ylab("ATP Production")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
# plot.Seahorse.AZ3146.ATP_H3.pdf
# AZ3146 (paired all)
t.test(as.numeric(subset(SH.HCT.UBE2H3, AZ3146 == "0")$`ATP Production`),
       as.numeric(subset(SH.HCT.UBE2H3, AZ3146 == "1")$`ATP Production`)  ,
       paired=TRUE) # p=0.0002134 , difference = -0.294


# AZ3146 increases Spare Respiratory Capacity
ggplot(SH.HCT.UBE2H3, aes(x=AZ3146, y=`Spare Respiratory Capacity`, color= `Group Name`))+  
  geom_point()+  
  scale_color_manual(values= c("turquoise2","turquoise2", "turquoise4",  "turquoise4", 
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4",
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4","mediumorchid4", "mediumorchid4",
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4", "mediumorchid4", "mediumorchid4"))+
  stat_summary(fun = "mean", size= 0.3, color= "red", geom = "crossbar")+
  ylab("Spare Respiratory Capacity")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
# plot.Seahorse.AZ3146.SRC_H3.pdf
t.test(as.numeric(subset(SH.HCT.UBE2H3, AZ3146 == "0")$`Spare Respiratory Capacity`),
       as.numeric(subset(SH.HCT.UBE2H3, AZ3146 == "1")$`Spare Respiratory Capacity`)  ,
       paired=TRUE) # p= 0.0757 NS, difference = 0.502



# AZ3146  Non-Mitochondrial Oxygen Consumption (NS)
ggplot(SH.HCT.UBE2H3, aes(x=AZ3146, y=`Non-Mitochondrial Oxygen Consumption`, color= `Group Name`))+  
  geom_point()+  
  scale_color_manual(values= c("turquoise2","turquoise2", "turquoise4",  "turquoise4", 
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4",
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4","mediumorchid4", "mediumorchid4",
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4", "mediumorchid4", "mediumorchid4"))+
  stat_summary(fun = "mean", size= 0.3, color= "red", geom = "crossbar")+
  ylab("Non-Mitochondrial Oxygen Consumption")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
# plot.Seahorse.AZ3146.NMOC_H3.pdf
t.test(as.numeric(subset(SH.HCT.UBE2H3, AZ3146 == "0")$`Non-Mitochondrial Oxygen Consumption`),
       as.numeric(subset(SH.HCT.UBE2H3, AZ3146 == "1")$`Non-Mitochondrial Oxygen Consumption`)  ,
       paired=TRUE) # p= 0.73 NS, difference = -0.014




# AZ3146  Maximal Respiration NS
ggplot(SH.HCT.UBE2H3, aes(x=AZ3146, y=`Maximal Respiration`, color= `Group Name`))+  
  geom_point()+  
  scale_color_manual(values= c("turquoise2","turquoise2", "turquoise4",  "turquoise4", 
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4",
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4","mediumorchid4", "mediumorchid4",
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4", "mediumorchid4", "mediumorchid4"))+
  stat_summary(fun = "mean", size= 0.3, color= "red", geom = "crossbar")+
  ylab("Maximal Respiration")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
# plot.Seahorse.AZ3146.Max.pdf
t.test(as.numeric(subset(SH.HCT.UBE2H3, AZ3146 == "0")$`Maximal Respiration`),
       as.numeric(subset(SH.HCT.UBE2H3, AZ3146 == "1")$`Maximal Respiration`)  ,
       paired=TRUE) # p=0.4 NS , difference = 0.2471751



###
# Basal Respiration (UBE2H-KO decreases, AZ decreases)
ggplot(SH.HCT.UBE2H3, aes(x=`Group Name`, y=`Basal Respiration`, fill=`Group Name`, shape= AZ3146 ))+
  geom_boxplot()+
  scale_fill_manual(values= c("turquoise2","turquoise2", "turquoise4",  "turquoise4", 
                              "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4",
                              "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4","mediumorchid4", "mediumorchid4",
                              "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4", "mediumorchid4", "mediumorchid4"))+
  geom_point()+   
  ylab("Basal Respiration")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
#plot.Seahorse.HCT116.BaselResp.pdf

# Just look at AZ3146 (not AZ): Basal Respiration: 
ggplot(subset(SH.HCT.UBE2H3, AZ3146=="0"), 
       aes(x=`Group Name`, y=`Basal Respiration`, fill=`Group Name` ))+
  geom_boxplot()+
  scale_fill_manual(values= c("turquoise2", "turquoise4", 
                              "mediumorchid2", "mediumorchid4", 
                              "mediumorchid2", "mediumorchid4","mediumorchid4",
                              "mediumorchid2", "mediumorchid4", "mediumorchid4"))+
  geom_point()+   
  ylab("Basal Respiration")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
# 5x4
#plot.Seahorse.HCT116.noAZ.BaselResp_H3.pdf

# Difference in Basal Respiration: 
ggplot(subset(SH.HCT.UBE2H3, AZ3146=="1"), 
       aes(x=`Group Name`, y=log2(`AZ.Basal.Diff`), fill=`Group Name` ))+
  geom_boxplot()+
  scale_fill_manual(values= c("turquoise2", "turquoise4", 
                              "mediumorchid2", "mediumorchid4", 
                              "mediumorchid2", "mediumorchid4","mediumorchid4",
                              "mediumorchid2", "mediumorchid4", "mediumorchid4"))+
  geom_point()+   
  ylab("Basal Respiration")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
# 5x4
#plot.Seahorse.HCT116.BaselResp.Diff.pdf

t.test(as.numeric(subset(SH.HCT.UBE2H3, `Group Name` == "WT")$`Basal Respiration`),
       as.numeric(subset(SH.HCT.UBE2H3, `Group Name` == "OE")$`Basal Respiration`)) # p=
t.test(as.numeric(subset(SH.HCT.UBE2H3, `Group Name` == "WT")$`Basal Respiration`),
       as.numeric(subset(SH.HCT.UBE2H3, `Group Name` == "KO1")$`Basal Respiration`) ) # p=0.004734
t.test(as.numeric(subset(SH.HCT.UBE2H3, `Group Name` == "WT")$`Basal Respiration`),
       as.numeric(subset(SH.HCT.UBE2H3, `Group Name` == "KO2")$`Basal Respiration`) ) # p=0.0002006
t.test(as.numeric(subset(SH.HCT.UBE2H3, `Group Name` == "WT")$`Basal Respiration`),
       as.numeric(subset(SH.HCT.UBE2H3, `Group Name` == "KO3")$`Basal Respiration`) ) # p=0.0004135
t.test(as.numeric(subset(SH.HCT.UBE2H3, `Group Name` == "KO1")$`Basal Respiration`),
       as.numeric(subset(SH.HCT.UBE2H3, `Group Name` == "KO1.OE1")$`Basal Respiration`) ) # p= 0.02335 lower
t.test(as.numeric(subset(SH.HCT.UBE2H3, `Group Name` == "KO2")$`Basal Respiration`),
       as.numeric(subset(SH.HCT.UBE2H3, `Group Name` == "KO2.OE1")$`Basal Respiration`) ) # p= NS higher
t.test(as.numeric(subset(SH.HCT.UBE2H3, `Group Name` == "KO2")$`Basal Respiration`),
       as.numeric(subset(SH.HCT.UBE2H3, `Group Name` == "KO2.OE2")$`Basal Respiration`) ) # p= NS lower
t.test(as.numeric(subset(SH.HCT.UBE2H3, `Group Name` == "KO3")$`Basal Respiration`),
       as.numeric(subset(SH.HCT.UBE2H3, `Group Name` == "KO3.OE2")$`Basal Respiration`) ) # p= NS lower
t.test(as.numeric(subset(SH.HCT.UBE2H3, `Group Name` == "KO3")$`Basal Respiration`), 
       as.numeric(subset(SH.HCT.UBE2H3, `Group Name` == "KO3.OE5")$`Basal Respiration`) ) # p= 0.04807 higher


# ATP Production (***) (UBE2H decreases, AZ decreases)
ggplot(SH.HCT.UBE2H3, aes(x=`Group Name`, y=SH.HCT.UBE2H3$`ATP Production`, fill=`Group Name`, shape= AZ3146 ))+
  geom_boxplot()+
  scale_fill_manual(values= c("turquoise2","turquoise2", "turquoise4",  "turquoise4", 
                              "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4",
                              "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4","mediumorchid4", "mediumorchid4",
                              "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4", "mediumorchid4", "mediumorchid4"))+
  geom_point()+   
  ylab("ATP Production")+  
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
#plot.Seahorse.HCT116.ATP.pdf

# Just look at AZ3146 (not AZ): ATP production: 
ggplot(subset(SH.HCT.UBE2H3, AZ3146=="0"), 
       aes(x=`Group Name`, y=`ATP Production`, fill=`Group Name`))+
  geom_boxplot()+
  scale_fill_manual(values= c("turquoise2", "turquoise4", 
                              "mediumorchid2", "mediumorchid4", 
                              "mediumorchid2", "mediumorchid4","mediumorchid4",
                              "mediumorchid2", "mediumorchid4","mediumorchid4"))+
  geom_point()+   
  ylab("ATP Production")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
#plot.Seahorse.HCT116.noAZ.ATP_H3.pdf

# Difference in ATP: 
ggplot(subset(SH.HCT.UBE2H3, AZ3146=="1"), 
       aes(x=`Group Name`, y=log2(`AZ.ATP.Diff`), fill=`Group Name` ))+
  geom_boxplot()+
  scale_fill_manual(values= c("turquoise2", "turquoise4", 
                              "mediumorchid2", "mediumorchid4", 
                              "mediumorchid2", "mediumorchid4","mediumorchid4",
                              "mediumorchid2", "mediumorchid4", "mediumorchid4"))+
  geom_point()+   
  ylab("log2 FC difference in Basal Respiration")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
# 5x4
#plot.Seahorse.HCT116.ATP.Diff.pdf

t.test(as.numeric(subset(SH.HCT.UBE2H3, `Group Name` == "WT")$`ATP Production`),
       as.numeric(subset(SH.HCT.UBE2H3, `Group Name` == "OE")$`ATP Production`)) # p= 
t.test(as.numeric(subset(SH.HCT.UBE2H3, `Group Name` == "WT")$`ATP Production`),
       as.numeric(subset(SH.HCT.UBE2H3, `Group Name` == "KO1")$`ATP Production`) ) # p=0.0312
t.test(as.numeric(subset(SH.HCT.UBE2H3, `Group Name` == "WT")$`ATP Production`),
       as.numeric(subset(SH.HCT.UBE2H3, `Group Name` == "KO2")$`ATP Production`) ) # p=0.0004411
t.test(as.numeric(subset(SH.HCT.UBE2H3, `Group Name` == "WT")$`ATP Production`),
       as.numeric(subset(SH.HCT.UBE2H3, `Group Name` == "KO3")$`ATP Production`) ) # p=0.002122
t.test(as.numeric(subset(SH.HCT.UBE2H3, `Group Name` == "KO1")$`ATP Production`),
       as.numeric(subset(SH.HCT.UBE2H3, `Group Name` == "KO1.OE1")$`ATP Production`) ) # p=0.002316 lower
t.test(as.numeric(subset(SH.HCT.UBE2H3, `Group Name` == "KO2")$`ATP Production`),
       as.numeric(subset(SH.HCT.UBE2H3, `Group Name` == "KO2.OE1")$`ATP Production`) ) # p= NS lower
t.test(as.numeric(subset(SH.HCT.UBE2H3, `Group Name` == "KO2")$`ATP Production`),
       as.numeric(subset(SH.HCT.UBE2H3, `Group Name` == "KO2.OE2")$`ATP Production`) ) # p= NS lower
t.test(as.numeric(subset(SH.HCT.UBE2H3, `Group Name` == "KO3")$`ATP Production`),
       as.numeric(subset(SH.HCT.UBE2H3, `Group Name` == "KO3.OE2")$`ATP Production`) ) # p=0.02101 lower
t.test(as.numeric(subset(SH.HCT.UBE2H3, `Group Name` == "KO3")$`ATP Production`),
       as.numeric(subset(SH.HCT.UBE2H3, `Group Name` == "KO3.OE5")$`ATP Production`) ) # p=0.06534 NS higher



# Non-Mitochondrial Oxygen Consumption NS (UBE2H decreases,  AZ increases)
ggplot(SH.HCT.UBE2H3, aes(x=`Group Name`, y=`Non-Mitochondrial Oxygen Consumption`, fill=`Group Name`, shape= AZ3146 ))+
  geom_boxplot()+
  scale_fill_manual(values= c("turquoise2","turquoise2", "turquoise4",  "turquoise4", 
                              "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4",
                              "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4","mediumorchid4", "mediumorchid4",
                              "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4", "mediumorchid4", "mediumorchid4"))+
  geom_point()+    
  ylab("Non-Mitochondrial\nOxygen Consumption")+ 
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
#plot.Seahorse.HCT116.NMOC.pdf

# Just look at AZ3146 (not AZ): Non-Mitochondrial Oxygen Consumption
ggplot(subset(SH.HCT.UBE2H3, AZ3146=="0"), 
       aes(x=`Group Name`, y=`Non-Mitochondrial Oxygen Consumption`, fill=`Group Name`))+
  geom_boxplot()+
  scale_fill_manual(values= c("turquoise2", "turquoise4", 
                              "mediumorchid2", "mediumorchid4", 
                              "mediumorchid2", "mediumorchid4", "mediumorchid4",
                              "mediumorchid2",  "mediumorchid4","mediumorchid4"))+
  geom_point()+   
  ylab("Non-Mitochondrial\nOxygen Consumption")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
#plot.Seahorse.HCT116.noAZ.NMOC_H3.pdf




# Spare Respiratory Capacity (AZ3146 increases SRC)
ggplot(subset(SH.HCT.UBE2H3, AZ3146=="0"), 
       aes(x=`Group Name`, y=`Spare Respiratory Capacity`, fill=`Group Name`, shape= AZ3146 ))+
  geom_boxplot()+
  scale_fill_manual(values= c("turquoise2","turquoise2", "turquoise4",  "turquoise4", 
                              "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4",
                              "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4","mediumorchid4", "mediumorchid4",
                              "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4", "mediumorchid4", "mediumorchid4"))+
  geom_point()+   
  ylab("Spare Respiratory Capacity")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))

# Spare Respiratory Capacity as a % of basal. (maybe? but I think this is just reflection of basal) 
ggplot(SH.HCT.UBE2H3, aes(x=`Group Name`, y=SH.HCT.UBE2H3$`Spare Respiratory Capacity as a %`, fill=`Group Name`, shape= AZ3146 ))+
  geom_boxplot()+
  scale_fill_manual(values= c("turquoise2","turquoise2", "turquoise4",  "turquoise4", 
                              "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4",
                              "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4","mediumorchid4", "mediumorchid4",
                              "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4", "mediumorchid4", "mediumorchid4"))+
  geom_point()+   
  ylab("Spare Respiratory Capacity as a % of basal")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))


# Maximal Respiration (eh.. OE and KO decrease max)
ggplot(SH.HCT.UBE2H3, aes(x=`Group Name`, y=SH.HCT.UBE2H3$`Maximal Respiration`, fill=`Group Name`, shape= AZ3146 ))+
  geom_boxplot()+
  scale_fill_manual(values= c("turquoise2","turquoise2", "turquoise4",  "turquoise4", 
                              "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4",
                              "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4","mediumorchid4", "mediumorchid4",
                              "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4", "mediumorchid4", "mediumorchid4"))+
  geom_point()+   
  ylab("Maximal Respiration")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
#plot.Seahorse.HCT116.MAX.pdf

ggplot(subset(SH.HCT.UBE2H3, AZ3146=="0"), 
       aes(x=`Group Name`, y=`Maximal Respiration`, fill=`Group Name`))+
  geom_boxplot()+
  scale_fill_manual(values= c("turquoise2", "turquoise4", 
                              "mediumorchid2", "mediumorchid4", 
                              "mediumorchid2", "mediumorchid4","mediumorchid4",
                              "mediumorchid2", "mediumorchid4","mediumorchid4"))+
  geom_point()+   
  ylab("Maximal Respiration")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
#plot.Seahorse.HCT116.noAZ.MAX.pdf

# Proton Leak. NS (OE and KO reduce leak, inconsistent)
ggplot(SH.HCT.UBE2H3, aes(x=`Group Name`, y=SH.HCT.UBE2H3$`Proton Leak`, fill=`Group Name`, shape= AZ3146 ))+
  geom_boxplot()+
  scale_fill_manual(values= c("turquoise2","turquoise2", "turquoise4",  "turquoise4", 
                              "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4",
                              "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4","mediumorchid4", "mediumorchid4",
                              "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4", "mediumorchid4", "mediumorchid4"))+
  geom_point()+   
  ylab("Proton Leak")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))




    #### HCT116 Merge ####
# Data organization: 
SH.HCT.UBE2H_Merge<- rbind(SH.HCT.UBE2H3, SH.HCT.UBE2H2)
SH.HCT.UBE2H_Merge<- rbind(SH.HCT.UBE2H_Merge, SH.HCT.UBE2H[, 1:15])

SH.HCT.UBE2H_Merge$`Group Name` <- sub("H.", "", SH.HCT.UBE2H_Merge$`Group Name`) 
SH.HCT.UBE2H_Merge$`Group Name` <- sub(" ", ".", SH.HCT.UBE2H_Merge$`Group Name`) 
SH.HCT.UBE2H_Merge$`Group Name` <- sub(" \\+ AZ", ".AZ", SH.HCT.UBE2H_Merge$`Group Name`) 
SH.HCT.UBE2H_Merge$`Group Name` <- sub("\\+ AZ", ".AZ", SH.HCT.UBE2H_Merge$`Group Name`) 
SH.HCT.UBE2H_Merge$`Group Name` <- sub("\\.\\.", ".", SH.HCT.UBE2H_Merge$`Group Name`) 
SH.HCT.UBE2H_Merge$`Group Name` <- factor(SH.HCT.UBE2H_Merge$`Group Name`, 
                                         levels = c("WT", "WT.AZ", "OE", "OE.AZ", "KO1", "KO1.AZ", "KO1.OE", "KO1.OE.AZ", "KO1.OE1", "KO1.OE1.AZ",
                                                    "KO2", "KO2.AZ", "KO2.OE", "KO2.OE.AZ", "KO2.OE1", "KO2.OE1.AZ", "KO2.OE2", "KO2.OE2.AZ",
                                                    "KO3", "KO3.AZ", "KO3.OE", "KO3.OE.AZ", "KO3.OE2", "KO3.OE2.AZ", "KO3.OE5", "KO3.OE5.AZ"))

SH.HCT.UBE2H_Merge$Group2<- as.character(SH.HCT.UBE2H_Merge$`Group Name`)
for (i in 1:length(SH.HCT.UBE2H_Merge$Group2)) {
  if (SH.HCT.UBE2H_Merge$Group2[i] %in% c("KO1", "KO2", "KO3")){
    SH.HCT.UBE2H_Merge$Group2[i]<- "KO"
  } else if (SH.HCT.UBE2H_Merge$Group2[i] %in% c("KO1.OE", "KO2.OE", "KO3.OE", "KO1.OE1", "KO2.OE1", "KO3.OE2", "KO2.OE2", "KO3.OE5")){
    SH.HCT.UBE2H_Merge$Group2[i]<- "KO.OE"
  } else if (SH.HCT.UBE2H_Merge$Group2[i] %in% c("KO1.AZ", "KO2.AZ", "KO3.AZ")){
    SH.HCT.UBE2H_Merge$Group2[i]<- "KO.AZ"
  }else if (SH.HCT.UBE2H_Merge$Group2[i] %in% c("KO1.OE.AZ", "KO2.OE.AZ", "KO3.OE.AZ", "KO1.OE1.AZ", "KO2.OE1.AZ", "KO3.OE2.AZ", "KO2.OE2.AZ", "KO3.OE5.AZ")){
    SH.HCT.UBE2H_Merge$Group2[i]<- "KO.OE.AZ"
  }
}
SH.HCT.UBE2H_Merge$Group2 <- factor(SH.HCT.UBE2H_Merge$Group2, 
                                          levels = c("WT", "WT.AZ", "OE", "OE.AZ", "KO", "KO.AZ", "KO.OE", "KO.OE.AZ"))


SH.HCT.UBE2H_Merge$AZ3146<- factor(SH.HCT.UBE2H_Merge$AZ3146, levels = c(0,1))



# AZ3146 decreases Basal Respiration
ggplot(SH.HCT.UBE2H_Merge, aes(x=AZ3146, y=`Basal Respiration`, color= `Group Name`))+  
  geom_point()+  
  scale_color_manual(values= c("turquoise2","turquoise2", "turquoise4",  "turquoise4", 
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4","mediumorchid4", "mediumorchid4",
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4","mediumorchid4", "mediumorchid4","mediumorchid4", "mediumorchid4",
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4", "mediumorchid4", "mediumorchid4","mediumorchid4", "mediumorchid4"))+
  stat_summary(fun = "mean", size= 0.3, color= "red", geom = "crossbar")+
  ylab("Basal Respiration")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
# plot.Seahorse.AZ3146.BasalResp.pdf
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, AZ3146 == 0)$`Basal Respiration`),
       as.numeric(subset(SH.HCT.UBE2H_Merge, AZ3146 == 1)$`Basal Respiration`)  ,
       paired=TRUE) # p= 1.117e-06, difference = -0.2028318


# AZ3146 decreases ATP production
ggplot(subset(SH.HCT.UBE2H_Merge, `Group Name` %in% c("WT", "WT.AZ")), 
       aes(x=AZ3146, y=`ATP Production`, fill= `Group Name`))+  
  geom_boxplot()+
  scale_fill_manual(values= c("turquoise2", "turquoise2"))+
  geom_point()+   
  #stat_summary(fun = "mean", size= 0.3, color= "red", geom = "crossbar")+
  ylab("Basal Respiration")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
# plot.Seahorse.AZ3146.ATP.WT.pdf
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, AZ3146 == "0" & `Group Name` %in% c("WT", "WT.AZ"))$`ATP Production`),
       as.numeric(subset(SH.HCT.UBE2H_Merge, AZ3146 == "1" & `Group Name` %in% c("WT", "WT.AZ"))$`ATP Production`)) 
# p= 0.7471, difference = 0.05

ggplot(subset(SH.HCT.UBE2H_Merge, `Group Name` %in% c("KO1", "KO1.AZ", "KO2", "KO2.AZ", "KO3", "KO3.AZ")), 
       aes(x=AZ3146, y=`ATP Production`, fill= "mediumorchid2"))+  
  geom_boxplot()+
  scale_fill_manual(values= c("mediumorchid2"))+
  geom_point()+   
  #stat_summary(fun = "mean", size= 0.3, color= "red", geom = "crossbar")+
  ylab("Basal Respiration")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
# plot.Seahorse.AZ3146.ATP.KO.pdf
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, AZ3146 == "0" & `Group Name` %in% c("KO1", "KO1.AZ", "KO2", "KO2.AZ", "KO3", "KO3.AZ"))$`ATP Production`),
       as.numeric(subset(SH.HCT.UBE2H_Merge, AZ3146 == "1" & `Group Name` %in% c("KO1", "KO1.AZ", "KO2", "KO2.AZ", "KO3", "KO3.AZ"))$`ATP Production`)) 
# p= 0.06002, difference = 0.23

ggplot(subset(SH.HCT.UBE2H_Merge, `Group Name` %in% c("KO1.OE1", "KO1.OE1.AZ", 
                                                      "KO2.OE1", "KO2.OE1.AZ", "KO2.OE2", "KO2.OE2.AZ", 
                                                      "KO3.OE2", "KO3.OE2.AZ", "KO3.OE3", "KO3.OE5.AZ" )), 
       aes(x=AZ3146, y=`ATP Production`, fill= "mediumorchid4"))+  
  geom_boxplot()+
  scale_fill_manual(values= c("mediumorchid4","mediumorchid4"))+
  geom_point()+   
  #stat_summary(fun = "mean", size= 0.3, color= "red", geom = "crossbar")+
  ylab("Basal Respiration")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
# plot.Seahorse.AZ3146.ATP.KOOE.pdf
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, AZ3146 == "0" & `Group Name` %in% c("KO1.OE1", "KO1.OE1.AZ",  "KO2.OE1", "KO2.OE1.AZ", "KO2.OE2", "KO2.OE2.AZ", 
                                                                                 "KO3.OE2", "KO3.OE2.AZ", "KO3.OE3", "KO3.OE5.AZ" ))$`ATP Production`),
       as.numeric(subset(SH.HCT.UBE2H_Merge, AZ3146 == "1" & `Group Name` %in% c("KO1.OE1", "KO1.OE1.AZ",  "KO2.OE1", "KO2.OE1.AZ", "KO2.OE2", "KO2.OE2.AZ", 
                                                                                 "KO3.OE2", "KO3.OE2.AZ", "KO3.OE3", "KO3.OE5.AZ" ))$`ATP Production`)) 
# p= 0.02717, difference = 0.25



# AZ3146 decreases ATP production
ggplot(SH.HCT.UBE2H_Merge, aes(x=AZ3146, y=`ATP Production`, color= `Group Name`))+  
  geom_point()+  
  scale_color_manual(values= c("turquoise2","turquoise2", "turquoise4",  "turquoise4", 
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4","mediumorchid4", "mediumorchid4",
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4","mediumorchid4", "mediumorchid4","mediumorchid4", "mediumorchid4",
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4","mediumorchid4", "mediumorchid4", "mediumorchid4", "mediumorchid4"))+
  stat_summary(fun = "mean", size= 0.3, color= "red", geom = "crossbar")+
  ylab("ATP Production")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
# plot.Seahorse.AZ3146.ATP.pdf
# AZ3146 (paired all)
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, AZ3146 == "0")$`ATP Production`),
       as.numeric(subset(SH.HCT.UBE2H_Merge, AZ3146 == "1")$`ATP Production`)  ,
       paired=TRUE) # p= 1.772e-08, difference = -0.2149754


# AZ3146 Non-Mitochondrial Oxygen Consumption (NS)

# AZ3146  Maximal Respiration NS

# AZ3146  Protein Leak NS

#### 
ggplot(subset(SH.HCT.UBE2H_Merge, AZ3146=="0"),  
       aes(x=Group2, y=`Basal Respiration`, fill=Group2))+
  geom_boxplot()+
  scale_fill_manual(values= c("turquoise2", "turquoise4", "mediumorchid2", "mediumorchid4"))+
  geom_point()+   
  ylab("Basal Respiration")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
#plot.Seahorse.HCT116.BaselResp.pdf

ggplot(subset(SH.HCT.UBE2H_Merge, AZ3146=="0" & Group2 %in% c("WT", "KO")),  
       aes(x=`Group Name`, y=`Basal Respiration`, fill=Group2))+
  geom_boxplot()+
  scale_fill_manual(values= c("turquoise2",  "mediumorchid2","mediumorchid2","mediumorchid2"))+
  geom_point()+   
  ylab("Basal Respiration")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
#plot.Seahorse.HCT116.BaselResp_WT.KO.pdf
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, Group2 == "WT")$`Basal Respiration`), 
       as.numeric(subset(SH.HCT.UBE2H_Merge, Group2 == "KO")$`Basal Respiration`) ) # p= 


ggplot(subset(SH.HCT.UBE2H_Merge, Group2 %in% c("WT", "WT.AZ", "KO", "KO.AZ")),  
       aes(x=`Group Name`, y=`Basal Respiration`, fill=Group2))+
  geom_boxplot()+
  scale_fill_manual(values= c("turquoise2", "turquoise2",  "mediumorchid2","mediumorchid2","mediumorchid2", "mediumorchid2","mediumorchid2","mediumorchid2"))+
  geom_point()+   
  ylab("Basal Respiration")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
#plot.Seahorse.HCT116.BaselResp_WT.KO.pdf
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "WT")$`Basal Respiration`), 
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO1")$`Basal Respiration`) ) # p= 
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "WT")$`Basal Respiration`), 
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO2")$`Basal Respiration`) ) # p= 
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "WT")$`Basal Respiration`), 
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO3")$`Basal Respiration`) ) # p= 


ggplot(subset(SH.HCT.UBE2H_Merge, AZ3146=="0" ),  
       aes(x=Group2, y=`ATP Production`, fill=Group2))+
  geom_boxplot()+
  scale_fill_manual(values= c("turquoise2", "turquoise4", "mediumorchid2", "mediumorchid4"))+
  geom_point()+   
  ylab("ATP production")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
#plot.Seahorse.HCT116.ATP.pdf
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, Group2 == "WT")$`ATP Production`), 
       as.numeric(subset(SH.HCT.UBE2H_Merge, Group2 == "KO")$`ATP Production`) ) # p= 3.49e-08
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, Group2 == "WT")$`ATP Production`), 
       as.numeric(subset(SH.HCT.UBE2H_Merge, Group2 == "KO.OE")$`ATP Production`) ) # p= 2.343e-08
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, Group2 == "KO")$`ATP Production`), 
       as.numeric(subset(SH.HCT.UBE2H_Merge, Group2 == "KO.OE")$`ATP Production`) ) # p= NS

ggplot(subset(SH.HCT.UBE2H_Merge, Group2 %in% c("WT", "WT.AZ", "KO", "KO.AZ")),  
       aes(x=`Group Name`, y=`ATP Production`, fill=Group2))+
  geom_boxplot()+
  scale_fill_manual(values= c("turquoise2", "turquoise2",  "mediumorchid2","mediumorchid2","mediumorchid2","mediumorchid2","mediumorchid2","mediumorchid2"))+
  geom_point()+   
  ylab("ATP Production")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
#plot.Seahorse.HCT116.ATP_WT.KO.AZ.pdf
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "WT")$`ATP Production`), 
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO1")$`ATP Production`) ) # p= 0.0008841
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "WT")$`ATP Production`), 
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO2")$`ATP Production`) ) # p= 1.432e-11
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "WT")$`ATP Production`), 
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO3")$`ATP Production`) ) # p= 3.889e-07

t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "WT")$`ATP Production`), 
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "WT.AZ")$`ATP Production`) ) # p= NS
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO1")$`ATP Production`), 
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO1.AZ")$`ATP Production`) ) # p= 0.01898
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO2")$`ATP Production`), 
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO2.AZ")$`ATP Production`) ) # p= 0.1031
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO3")$`ATP Production`), 
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO3.AZ")$`ATP Production`) ) # p= 0.5078



ggplot(subset(SH.HCT.UBE2H_Merge, AZ3146=="0"),  aes(x=Group2, y=`Non-Mitochondrial Oxygen Consumption`, fill=Group2))+
  geom_boxplot()+
  scale_fill_manual(values= c("turquoise2", "turquoise4", "mediumorchid2", "mediumorchid4"))+
  geom_point()+   
  ylab("Non-Mitochondrial Oxygen Consumption")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
#plot.Seahorse.HCT116.NMOC.pdf
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, Group2 == "WT")$`Non-Mitochondrial Oxygen Consumption`), 
       as.numeric(subset(SH.HCT.UBE2H_Merge, Group2 == "KO")$`Non-Mitochondrial Oxygen Consumption`) ) # p= 0.0008607
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, Group2 == "WT")$`Non-Mitochondrial Oxygen Consumption`), 
       as.numeric(subset(SH.HCT.UBE2H_Merge, Group2 == "OE")$`Non-Mitochondrial Oxygen Consumption`) ) # p= NS
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, Group2 == "WT")$`Non-Mitochondrial Oxygen Consumption`), 
       as.numeric(subset(SH.HCT.UBE2H_Merge, Group2 == "KO.OE")$`Non-Mitochondrial Oxygen Consumption`) ) # p= 0.002207
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, Group2 == "KO")$`Non-Mitochondrial Oxygen Consumption`), 
       as.numeric(subset(SH.HCT.UBE2H_Merge, Group2 == "KO.OE")$`Non-Mitochondrial Oxygen Consumption`) ) # p= NS



###
# Basal Respiration (UBE2H-KO decreases, AZ decreases if UBE2H-KO)
ggplot(SH.HCT.UBE2H_Merge, aes(x=`Group Name`, y=`Basal Respiration`, fill=`Group Name`, shape= AZ3146 ))+
  geom_boxplot()+
  scale_fill_manual(values= c("turquoise2","turquoise2", "turquoise4",  "turquoise4", 
                              "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4","mediumorchid4", "mediumorchid4",
                              "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4","mediumorchid4", "mediumorchid4","mediumorchid4", "mediumorchid4",
                              "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4","mediumorchid4", "mediumorchid4", "mediumorchid4", "mediumorchid4"))+
  geom_point()+   
  ylab("Basal Respiration")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
#plot.Seahorse.HCT116.BaselResp.pdf

# Just look at AZ3146 (not AZ): Basal Respiration: 
ggplot(subset(SH.HCT.UBE2H_Merge, AZ3146=="0"), 
       aes(x=`Group Name`, y=`Basal Respiration`, fill=`Group Name` ))+
  geom_boxplot()+
  scale_fill_manual(values= c("turquoise2", "turquoise4", 
                              "mediumorchid2", "mediumorchid4","mediumorchid4", 
                              "mediumorchid2", "mediumorchid4","mediumorchid4","mediumorchid4",
                              "mediumorchid2", "mediumorchid4","mediumorchid4", "mediumorchid4"))+
  geom_point()+   
  ylab("Basal Respiration")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
# 5x4
#plot.Seahorse.HCT116.noAZ.BaselResp.pdf

t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "WT")$`Basal Respiration`),
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "OE")$`Basal Respiration`)) # p= lower
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "WT")$`Basal Respiration`),
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO1")$`Basal Respiration`) ) # p=4.691e-06
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "WT")$`Basal Respiration`),
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO2")$`Basal Respiration`) ) # p=4.005e-12
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "WT")$`Basal Respiration`),
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO3")$`Basal Respiration`) ) # p=1.625e-08
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO1")$`Basal Respiration`), 
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO1.OE")$`Basal Respiration`) ) # p= OE higher
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO1")$`Basal Respiration`),
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO1.OE1")$`Basal Respiration`) ) # p= OE higher
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO2")$`Basal Respiration`), 
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO2.OE")$`Basal Respiration`) ) # p= lower
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO2")$`Basal Respiration`),
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO2.OE1")$`Basal Respiration`) ) # p= OE higher
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO2")$`Basal Respiration`),
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO2.OE2")$`Basal Respiration`) ) # p= OE higher
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO3")$`Basal Respiration`), 
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO3.OE")$`Basal Respiration`) ) # p= lower
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO3")$`Basal Respiration`),
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO3.OE2")$`Basal Respiration`) ) # p= lower
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO3")$`Basal Respiration`), 
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO3.OE5")$`Basal Respiration`) ) # p= OE higher



# ATP Production (***) (UBE2H decreases, AZ decreases)
ggplot(SH.HCT.UBE2H_Merge, aes(x=`Group Name`, y=SH.HCT.UBE2H_Merge$`ATP Production`, fill=`Group Name`, shape= AZ3146 ))+
  geom_boxplot()+
  scale_fill_manual(values= c("turquoise2","turquoise2", "turquoise4",  "turquoise4", 
                              "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4", "mediumorchid4", "mediumorchid4",
                              "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4", "mediumorchid4", "mediumorchid4","mediumorchid4", "mediumorchid4",
                              "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4","mediumorchid4", "mediumorchid4", "mediumorchid4", "mediumorchid4"))+
  geom_point()+   
  ylab("ATP Production")+  
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
#plot.Seahorse.HCT116.ATP_Merge.pdf

# Just look at AZ3146 (not AZ): ATP production: 
ggplot(subset(SH.HCT.UBE2H_Merge, AZ3146=="0"), 
       aes(x=`Group Name`, y=`ATP Production`, fill=`Group Name`))+
  geom_boxplot()+
  scale_fill_manual(values= c("turquoise2", "turquoise4", 
                              "mediumorchid2", "mediumorchid4","mediumorchid4", 
                              "mediumorchid2", "mediumorchid4","mediumorchid4","mediumorchid4",
                              "mediumorchid2", "mediumorchid4","mediumorchid4","mediumorchid4"))+
  geom_point()+   
  ylab("ATP Production")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
#plot.Seahorse.HCT116.noAZ.ATP_Merge.pdf


# Just look at AZ3146 (not AZ): ATP production: 
ggplot(subset(SH.HCT.UBE2H_Merge, AZ3146=="0" & ! `Group Name` %in% c("KO1.OE", "KO2.OE", "KO3.OE")), 
       aes(x=`Group Name`, y=`ATP Production`, fill=`Group Name`))+
  geom_boxplot()+
  scale_fill_manual(values= c("turquoise2", "turquoise4", 
                              "mediumorchid2", "mediumorchid4",
                              "mediumorchid2", "mediumorchid4","mediumorchid4",
                              "mediumorchid2", "mediumorchid4","mediumorchid4"))+
  geom_point()+   
  ylab("ATP Production")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
#plot.Seahorse.HCT116.noAZ.ATP_Merge_KOclonesOnly.pdf

t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "WT")$`ATP Production`),
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "OE")$`ATP Production`)) # p= lower
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "WT")$`ATP Production`),
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO1")$`ATP Production`) ) # p=0.0008841
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "WT")$`ATP Production`),
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO2")$`ATP Production`) ) # p=1.432e-11
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "WT")$`ATP Production`),
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO3")$`ATP Production`) ) # p=3.889e-07
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO1")$`ATP Production`),
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO1.OE")$`ATP Production`) ) # p=lower
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO1")$`ATP Production`),
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO1.OE1")$`ATP Production`) ) # p=lower
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO2")$`ATP Production`),
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO2.OE")$`ATP Production`) ) # p=lower
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO2")$`ATP Production`),
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO2.OE1")$`ATP Production`) ) # p= OE higher
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO2")$`ATP Production`),
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO2.OE2")$`ATP Production`) ) # p= OE higher
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO3")$`ATP Production`),
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO3.OE")$`ATP Production`) ) # p=lower
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO3")$`ATP Production`),
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO3.OE2")$`ATP Production`) ) # p= lower
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO3")$`ATP Production`),
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO3.OE5")$`ATP Production`) ) # p=OE higher 



# Non-Mitochondrial Oxygen Consumption NS (UBE2H decreases,  AZ increases)
ggplot(SH.HCT.UBE2H_Merge, aes(x=`Group Name`, y=`Non-Mitochondrial Oxygen Consumption`, fill=`Group Name`, shape= AZ3146 ))+
  geom_boxplot()+
  scale_fill_manual(values= c("turquoise2","turquoise2", "turquoise4",  "turquoise4", 
                              "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4", "mediumorchid4", "mediumorchid4",
                              "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4", "mediumorchid4", "mediumorchid4","mediumorchid4", "mediumorchid4",
                              "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4", "mediumorchid4", "mediumorchid4", "mediumorchid4", "mediumorchid4"))+
  geom_point()+    
  ylab("Non-Mitochondrial\nOxygen Consumption")+ 
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
#plot.Seahorse.HCT116.NMOC_Merge.pdf


# Just look at AZ3146 (not AZ): Non-Mitochondrial Oxygen Consumption
ggplot(subset(SH.HCT.UBE2H_Merge, AZ3146=="0"), 
       aes(x=`Group Name`, y=`Non-Mitochondrial Oxygen Consumption`, fill=`Group Name`))+
  geom_boxplot()+
  scale_fill_manual(values= c("turquoise2", "turquoise4", 
                              "mediumorchid2", "mediumorchid4", "mediumorchid4", 
                              "mediumorchid2", "mediumorchid4", "mediumorchid4","mediumorchid4",
                              "mediumorchid2",  "mediumorchid4", "mediumorchid4", "mediumorchid4"))+
  geom_point()+   
  ylab("Non-Mitochondrial\nOxygen Consumption")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
# plot.Seahorse.HCT116.noAZ.NMOC_Merge.pdf

t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "WT")$`Non-Mitochondrial Oxygen Consumption`),
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "OE")$`Non-Mitochondrial Oxygen Consumption`)) # p= higher
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "WT")$`Non-Mitochondrial Oxygen Consumption`),
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO1")$`Non-Mitochondrial Oxygen Consumption`) ) # p=0.0109
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "WT")$`Non-Mitochondrial Oxygen Consumption`),
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO2")$`Non-Mitochondrial Oxygen Consumption`) ) # p=0.000141
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "WT")$`Non-Mitochondrial Oxygen Consumption`),
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO3")$`Non-Mitochondrial Oxygen Consumption`) ) # p=0.003167
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO1")$`Non-Mitochondrial Oxygen Consumption`),
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO1.OE")$`Non-Mitochondrial Oxygen Consumption`) ) # p=lower
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO1")$`Non-Mitochondrial Oxygen Consumption`),
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO1.OE1")$`Non-Mitochondrial Oxygen Consumption`) ) # p= OE higher
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO2")$`Non-Mitochondrial Oxygen Consumption`),
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO2.OE")$`Non-Mitochondrial Oxygen Consumption`) ) # p=lower
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO2")$`Non-Mitochondrial Oxygen Consumption`),
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO2.OE1")$`Non-Mitochondrial Oxygen Consumption`) ) # p= OE higher
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO2")$`Non-Mitochondrial Oxygen Consumption`),
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO2.OE2")$`Non-Mitochondrial Oxygen Consumption`) ) # p= OE higher
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO3")$`Non-Mitochondrial Oxygen Consumption`),
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO3.OE")$`Non-Mitochondrial Oxygen Consumption`) ) # p=lower
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO3")$`Non-Mitochondrial Oxygen Consumption`),
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO3.OE2")$`Non-Mitochondrial Oxygen Consumption`) ) # p= lower
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO3")$`Non-Mitochondrial Oxygen Consumption`),
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO3.OE5")$`Non-Mitochondrial Oxygen Consumption`) ) # p=0.04202 OE higher 









# Just look at AZ3146 (not AZ): Maxinal Respiration
ggplot(subset(SH.HCT.UBE2H_Merge, AZ3146=="0"), 
       aes(x=`Group Name`, y=`Maximal Respiration`, fill=`Group Name`))+
  geom_boxplot()+
  scale_fill_manual(values= c("turquoise2", "turquoise4", 
                              "mediumorchid2", "mediumorchid4", "mediumorchid4", 
                              "mediumorchid2", "mediumorchid4", "mediumorchid4","mediumorchid4",
                              "mediumorchid2",  "mediumorchid4", "mediumorchid4", "mediumorchid4"))+
  geom_point()+   
  ylab("Maxinal Respiration")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
# plot.Seahorse.HCT116.noAZ.MAX_Merge.pdf










    #### HCT116 Merge 1 & 2 & 3 KO.OEc (no bulk OE) ####

SH.HCT.UBE2H_Merge<- rbind(SH.HCT.UBE2H3, SH.HCT.UBE2H2) 
SH.HCT.UBE2H_Merge<- rbind(SH.HCT.UBE2H_Merge, SH.HCT.UBE2H[, 1:15])

SH.HCT.UBE2H_Merge$`Group Name` <- sub("H.", "", SH.HCT.UBE2H_Merge$`Group Name`) 
SH.HCT.UBE2H_Merge$`Group Name` <- sub(" ", ".", SH.HCT.UBE2H_Merge$`Group Name`) 
SH.HCT.UBE2H_Merge$`Group Name` <- sub(" \\+ AZ", ".AZ", SH.HCT.UBE2H_Merge$`Group Name`) 
SH.HCT.UBE2H_Merge$`Group Name` <- sub("\\+ AZ", ".AZ", SH.HCT.UBE2H_Merge$`Group Name`) 
SH.HCT.UBE2H_Merge$`Group Name` <- sub("\\.\\.", ".", SH.HCT.UBE2H_Merge$`Group Name`) 
SH.HCT.UBE2H_Merge$`Group Name` <- factor(SH.HCT.UBE2H_Merge$`Group Name`, 
                                          levels = c("WT", "WT.AZ", "OE", "OE.AZ", "KO1", "KO1.AZ", "KO1.OE", "KO1.OE.AZ", "KO1.OE1", "KO1.OE1.AZ",
                                                     "KO2", "KO2.AZ", "KO2.OE", "KO2.OE.AZ", "KO2.OE1", "KO2.OE1.AZ", "KO2.OE2", "KO2.OE2.AZ",
                                                     "KO3", "KO3.AZ", "KO3.OE", "KO3.OE.AZ", "KO3.OE2", "KO3.OE2.AZ", "KO3.OE5", "KO3.OE5.AZ"))

SH.HCT.UBE2H_Merge$Group2<- as.character(SH.HCT.UBE2H_Merge$`Group Name`)
for (i in 1:length(SH.HCT.UBE2H_Merge$Group2)) {
  if (SH.HCT.UBE2H_Merge$Group2[i] %in% c("KO1", "KO2", "KO3")){
    SH.HCT.UBE2H_Merge$Group2[i]<- "KO"
  } else if (SH.HCT.UBE2H_Merge$Group2[i] %in% c("KO1.OE", "KO2.OE", "KO3.OE", "KO1.OE1", "KO2.OE1", "KO3.OE2", "KO2.OE2", "KO3.OE5")){
    SH.HCT.UBE2H_Merge$Group2[i]<- "KO.OE"
  } else if (SH.HCT.UBE2H_Merge$Group2[i] %in% c("KO1.AZ", "KO2.AZ", "KO3.AZ")){
    SH.HCT.UBE2H_Merge$Group2[i]<- "KO.AZ"
  }else if (SH.HCT.UBE2H_Merge$Group2[i] %in% c("KO1.OE.AZ", "KO2.OE.AZ", "KO3.OE.AZ", "KO1.OE1.AZ", "KO2.OE1.AZ", "KO3.OE2.AZ", "KO2.OE2.AZ", "KO3.OE5.AZ")){
    SH.HCT.UBE2H_Merge$Group2[i]<- "KO.OE.AZ"
  }
}
SH.HCT.UBE2H_Merge$Group2 <- factor(SH.HCT.UBE2H_Merge$Group2, 
                                    levels = c("WT", "WT.AZ", "OE", "OE.AZ", "KO", "KO.AZ", "KO.OE", "KO.OE.AZ"))

SH.HCT.UBE2H_Merge<- subset(SH.HCT.UBE2H_Merge, `Group Name` %in% c("WT", "WT.AZ", "OE", "OE.AZ", 
                                                                     "KO1", "KO1.AZ","KO1.OE1", "KO1.OE1.AZ",
                                                                     "KO2", "KO2.AZ", "KO2.OE1", "KO2.OE1.AZ", "KO2.OE2", "KO2.OE2.AZ",
                                                                     "KO3", "KO3.AZ", "KO3.OE2", "KO3.OE2.AZ", "KO3.OE5", "KO3.OE5.AZ"))

SH.HCT.UBE2H_Merge$AZ3146<- factor(SH.HCT.UBE2H_Merge$AZ3146, levels = c(0,1))


# AZ3146 decreases Basal Respiration
ggplot(SH.HCT.UBE2H_Merge, aes(x=AZ3146, y=`Basal Respiration`, color= `Group Name`))+  
  geom_point()+  
  scale_color_manual(values= c("turquoise2","turquoise2", "turquoise4",  "turquoise4", 
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4",
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4","mediumorchid4", "mediumorchid4",
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4","mediumorchid4", "mediumorchid4"))+
  stat_summary(fun = "mean", size= 0.3, color= "red", geom = "crossbar")+
  ylab("Basal Respiration")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
# plot.Seahorse.AZ3146.BasalResp.pdf
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, AZ3146 == 0)$`Basal Respiration`),
       as.numeric(subset(SH.HCT.UBE2H_Merge, AZ3146 == 1)$`Basal Respiration`)  ,
       paired=TRUE) # p= 9.368e-06, difference = -0.2065278


# AZ3146 decreases Basal Respiration
ggplot(subset(SH.HCT.UBE2H_Merge, `Group Name` %in% c("WT", "WT.AZ")), 
       aes(x=AZ3146, y=`Basal Respiration`, color= `Group Name`))+  
  geom_point()+  
  scale_color_manual(values= c("turquoise2","turquoise2", "turquoise4",  "turquoise4", 
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4",
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4","mediumorchid4", "mediumorchid4",
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4","mediumorchid4", "mediumorchid4"))+
  stat_summary(fun = "mean", size= 0.3, color= "red", geom = "crossbar")+
  ylab("Basal Respiration")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))

ggplot(subset(SH.HCT.UBE2H_Merge, `Group Name` %in% c("KO1", "KO1.AZ", "KO2", "KO2.AZ", "KO3", "KO3.AZ")), 
       aes(x=AZ3146, y=`Basal Respiration`, color= `Group Name`))+  
  geom_point()+  
  scale_color_manual(values= c("turquoise2","turquoise2", "turquoise4",  "turquoise4", 
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4",
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4","mediumorchid4", "mediumorchid4",
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4","mediumorchid4", "mediumorchid4"))+
  stat_summary(fun = "mean", size= 0.3, color= "red", geom = "crossbar")+
  ylab("Basal Respiration")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))

ggplot(subset(SH.HCT.UBE2H_Merge, `Group Name` %in% c( "KO1.OE1", "KO1.OE1.AZ", 
                                                      "KO2.OE1", "KO2.OE1.AZ", "KO2.OE2", "KO2.OE2.AZ", 
                                                      "KO3.OE2", "KO3.OE2.AZ", "KO3.OE3", "KO3.OE5.AZ" )), 
       aes(x=AZ3146, y=`Basal Respiration`, color= `Group Name`))+  
  geom_point()+  
  scale_color_manual(values= c("turquoise2","turquoise2", "turquoise4",  "turquoise4", 
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4",
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4","mediumorchid4", "mediumorchid4",
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4","mediumorchid4", "mediumorchid4"))+
  stat_summary(fun = "mean", size= 0.3, color= "red", geom = "crossbar")+
  ylab("Basal Respiration")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
# plot.Seahorse.AZ3146.BasalResp.pdf




# AZ3146 decreases ATP production
ggplot(SH.HCT.UBE2H_Merge, aes(x=AZ3146, y=`ATP Production`, color= `Group Name`))+  
  geom_point()+  
  scale_color_manual(values= c("turquoise2","turquoise2", "turquoise4",  "turquoise4", 
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4",
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4","mediumorchid4", "mediumorchid4",
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4","mediumorchid4", "mediumorchid4"))+
  stat_summary(fun = "mean", size= 0.3, color= "black", geom = "crossbar")+
  ylab("ATP Production")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
# plot.Seahorse.AZ3146.ATP_run2.3.pdf
# Supplementary Figure 7B
# AZ3146 (paired all)
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, AZ3146 == "0")$`ATP Production`),
       as.numeric(subset(SH.HCT.UBE2H_Merge, AZ3146 == "1")$`ATP Production`)  ,
       paired=TRUE) # p= 2.203e-07, difference = -0.2208634



# AZ3146 Non-Mitochondrial Oxygen Consumption (NS)

# AZ3146  Maximal Respiration NS

# AZ3146  Protein Leak NS

#### 
ggplot(subset(SH.HCT.UBE2H_Merge, AZ3146=="0"),  
       aes(x=Group2, y=`Basal Respiration`, fill=Group2))+
  geom_boxplot()+
  scale_fill_manual(values= c("turquoise2", "turquoise4", "mediumorchid2", "mediumorchid4"))+
  geom_point()+   
  ylab("Basal Respiration")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
#plot.Seahorse.HCT116.BaselResp.pdf
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, Group2 == "WT")$`Basal Respiration`), 
       as.numeric(subset(SH.HCT.UBE2H_Merge, Group2 == "KO")$`Basal Respiration`) ) # p= 2.474e-09
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, Group2 == "WT")$`Basal Respiration`), 
       as.numeric(subset(SH.HCT.UBE2H_Merge, Group2 == "KO.OE")$`Basal Respiration`) ) # p= 3.04e-09
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, Group2 == "KO")$`Basal Respiration`), 
       as.numeric(subset(SH.HCT.UBE2H_Merge, Group2 == "KO.OE")$`Basal Respiration`) ) # p= NS



ggplot(subset(SH.HCT.UBE2H_Merge, AZ3146=="0" ),  
       aes(x=Group2, y=`ATP Production`, fill=Group2))+
  geom_boxplot()+
  scale_fill_manual(values= c("turquoise2", "turquoise4", "mediumorchid2", "mediumorchid4"))+
  geom_point()+   
  ylab("ATP production")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
#plot.Seahorse.HCT116.BaselResp.pdf
# Figure?
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, Group2 == "WT")$`ATP Production`), 
       as.numeric(subset(SH.HCT.UBE2H_Merge, Group2 == "KO")$`ATP Production`) ) # p= 3.49e-08
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, Group2 == "WT")$`ATP Production`), 
       as.numeric(subset(SH.HCT.UBE2H_Merge, Group2 == "KO.OE")$`ATP Production`) ) # p= 3.994e-08
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, Group2 == "KO")$`ATP Production`), 
       as.numeric(subset(SH.HCT.UBE2H_Merge, Group2 == "KO.OE")$`ATP Production`) ) # p= NS




ggplot(subset(SH.HCT.UBE2H_Merge, AZ3146=="0"),  aes(x=Group2, y=`Non-Mitochondrial Oxygen Consumption`, fill=Group2))+
  geom_boxplot()+
  scale_fill_manual(values= c("turquoise2", "turquoise4", "mediumorchid2", "mediumorchid4"))+
  geom_point()+   
  ylab("Non-Mitochondrial Oxygen Consumption")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
#plot.Seahorse.HCT116.NMOC.pdf
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, Group2 == "WT")$`Non-Mitochondrial Oxygen Consumption`), 
       as.numeric(subset(SH.HCT.UBE2H_Merge, Group2 == "KO")$`Non-Mitochondrial Oxygen Consumption`) ) # p= 0.0008607
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, Group2 == "WT")$`Non-Mitochondrial Oxygen Consumption`), 
       as.numeric(subset(SH.HCT.UBE2H_Merge, Group2 == "OE")$`Non-Mitochondrial Oxygen Consumption`) ) # p= NS
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, Group2 == "WT")$`Non-Mitochondrial Oxygen Consumption`), 
       as.numeric(subset(SH.HCT.UBE2H_Merge, Group2 == "KO.OE")$`Non-Mitochondrial Oxygen Consumption`) ) # p= 0.006833
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, Group2 == "KO")$`Non-Mitochondrial Oxygen Consumption`), 
       as.numeric(subset(SH.HCT.UBE2H_Merge, Group2 == "KO.OE")$`Non-Mitochondrial Oxygen Consumption`) ) # p= 0.14 higher


##
ggplot(subset(SH.HCT.UBE2H_Merge, AZ3146=="0"),  
       aes(x=Group2, y=`Basal Respiration`, fill=Group2))+
  geom_boxplot()+
  scale_fill_manual(values= c("turquoise2", "turquoise4", "mediumorchid2", "mediumorchid4"))+
  geom_point()+   
  ylab("Basal Respiration")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))


###
# Basal Respiration (UBE2H-KO decreases, AZ decreases)
ggplot(SH.HCT.UBE2H_Merge, aes(x=`Group Name`, y=`Basal Respiration`, fill=`Group Name`, shape= AZ3146 ))+
  geom_boxplot()+
  scale_fill_manual(values= c("turquoise2","turquoise2", "turquoise4",  "turquoise4", 
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4",
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4","mediumorchid4", "mediumorchid4",
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4","mediumorchid4", "mediumorchid4"))+
  geom_point()+   
  ylab("Basal Respiration")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
#plot.Seahorse.HCT116.BaselResp.pdf


# Just look at AZ3146 (not AZ): Basal Respiration: 
ggplot(subset(SH.HCT.UBE2H_Merge, AZ3146=="0"), 
       aes(x=`Group Name`, y=`Basal Respiration`, fill=`Group Name` ))+
  geom_boxplot()+
  scale_fill_manual(values= c("turquoise2", "turquoise4", 
                              "mediumorchid2", "mediumorchid4",
                              "mediumorchid2", "mediumorchid4","mediumorchid4",
                              "mediumorchid2", "mediumorchid4", "mediumorchid4"))+
  geom_point()+   
  ylab("Basal Respiration")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
# 5x4
#plot.Seahorse.HCT116.noAZ.BaselResp.pdf
# Supplementary Figure 7A
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "WT")$`Basal Respiration`),
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "OE")$`Basal Respiration`)) # p= lower
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "WT")$`Basal Respiration`),
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO1")$`Basal Respiration`) ) # p=4.691e-06
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "WT")$`Basal Respiration`),
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO2")$`Basal Respiration`) ) # p=4.005e-12
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "WT")$`Basal Respiration`),
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO3")$`Basal Respiration`) ) # p=1.625e-08
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO1")$`Basal Respiration`),
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO1.OE1")$`Basal Respiration`) ) # p= higher
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO2")$`Basal Respiration`),
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO2.OE1")$`Basal Respiration`) ) # p= 0.07 OE higher
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO2")$`Basal Respiration`),
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO2.OE2")$`Basal Respiration`) ) # p= OE higher
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO3")$`Basal Respiration`),
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO3.OE2")$`Basal Respiration`) ) # p= lower
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO3")$`Basal Respiration`), 
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO3.OE5")$`Basal Respiration`) ) # p= OE higher



# ATP Production (***) (UBE2H decreases, AZ decreases)
ggplot(SH.HCT.UBE2H_Merge, aes(x=`Group Name`, y=SH.HCT.UBE2H_Merge$`ATP Production`, fill=`Group Name`, shape= AZ3146 ))+
  geom_boxplot()+
  scale_fill_manual(values= c("turquoise2","turquoise2", "turquoise4",  "turquoise4", 
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4",
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4","mediumorchid4", "mediumorchid4",
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4","mediumorchid4", "mediumorchid4"))+
  geom_point()+   
  ylab("ATP Production")+  
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
#plot.Seahorse.HCT116.ATP_Merge.pdf

# Just look at AZ3146 (not AZ): ATP production: 
ggplot(subset(SH.HCT.UBE2H_Merge, AZ3146=="0"), 
       aes(x=`Group Name`, y=`ATP Production`, fill=`Group Name`))+
  geom_boxplot()+
  scale_fill_manual(values= c("turquoise2", "turquoise4", 
                              "mediumorchid2", "mediumorchid4",
                              "mediumorchid2", "mediumorchid4","mediumorchid4",
                              "mediumorchid2", "mediumorchid4", "mediumorchid4"))+
  geom_point()+   
  ylab("ATP Production")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
#plot.Seahorse.HCT116.noAZ.ATP_Merge_run2.3.pdf
table(SH.HCT.UBE2H_Merge$`Group Name`)

t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "WT")$`ATP Production`),
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "OE")$`ATP Production`)) # p= lower
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "WT")$`ATP Production`),
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO1")$`ATP Production`) ) # p=0.0008841
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "WT")$`ATP Production`),
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO2")$`ATP Production`) ) # p=1.432e-11
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "WT")$`ATP Production`),
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO3")$`ATP Production`) ) # p=3.889e-07
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO1")$`ATP Production`),
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO1.OE1")$`ATP Production`) ) # p=lower
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO2")$`ATP Production`),
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO2.OE1")$`ATP Production`) ) # p= OE higher
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO2")$`ATP Production`),
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO2.OE2")$`ATP Production`) ) # p= 0.06442OE higher
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO3")$`ATP Production`),
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO3.OE2")$`ATP Production`) ) # p= lower
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO3")$`ATP Production`),
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO3.OE5")$`ATP Production`) ) # p=OE higher 



# Non-Mitochondrial Oxygen Consumption NS (UBE2H decreases,  AZ increases)
ggplot(SH.HCT.UBE2H_Merge, aes(x=`Group Name`, y=`Non-Mitochondrial Oxygen Consumption`, fill=`Group Name`, shape= AZ3146 ))+
  geom_boxplot()+
  scale_fill_manual(values= c("turquoise2","turquoise2", "turquoise4",  "turquoise4", 
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4",
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4","mediumorchid4", "mediumorchid4",
                               "mediumorchid2","mediumorchid2", "mediumorchid4", "mediumorchid4","mediumorchid4", "mediumorchid4"))+
  geom_point()+    
  ylab("Non-Mitochondrial\nOxygen Consumption")+ 
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
#plot.Seahorse.HCT116.NMOC_Merge.pdf


# Just look at AZ3146 (not AZ): Non-Mitochondrial Oxygen Consumption
ggplot(subset(SH.HCT.UBE2H_Merge, AZ3146=="0"), 
       aes(x=`Group Name`, y=`Non-Mitochondrial Oxygen Consumption`, fill=`Group Name`))+
  geom_boxplot()+
  scale_fill_manual(values= c("turquoise2", "turquoise4", 
                              "mediumorchid2", "mediumorchid4",
                              "mediumorchid2", "mediumorchid4","mediumorchid4",
                              "mediumorchid2", "mediumorchid4", "mediumorchid4"))+
  geom_point()+   
  ylab("Non-Mitochondrial\nOxygen Consumption")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
# plot.Seahorse.HCT116.noAZ.NMOC_Merge.pdf

t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "WT")$`Non-Mitochondrial Oxygen Consumption`),
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "OE")$`Non-Mitochondrial Oxygen Consumption`)) # p= higher
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "WT")$`Non-Mitochondrial Oxygen Consumption`),
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO1")$`Non-Mitochondrial Oxygen Consumption`) ) # p=0.0109
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "WT")$`Non-Mitochondrial Oxygen Consumption`),
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO2")$`Non-Mitochondrial Oxygen Consumption`) ) # p=0.000141
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "WT")$`Non-Mitochondrial Oxygen Consumption`),
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO3")$`Non-Mitochondrial Oxygen Consumption`) ) # p=0.003167
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO1")$`Non-Mitochondrial Oxygen Consumption`),
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO1.OE")$`Non-Mitochondrial Oxygen Consumption`) ) # p=lower
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO1")$`Non-Mitochondrial Oxygen Consumption`),
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO1.OE1")$`Non-Mitochondrial Oxygen Consumption`) ) # p= OE higher
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO2")$`Non-Mitochondrial Oxygen Consumption`),
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO2.OE")$`Non-Mitochondrial Oxygen Consumption`) ) # p=lower
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO2")$`Non-Mitochondrial Oxygen Consumption`),
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO2.OE1")$`Non-Mitochondrial Oxygen Consumption`) ) # p= OE higher
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO2")$`Non-Mitochondrial Oxygen Consumption`),
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO2.OE2")$`Non-Mitochondrial Oxygen Consumption`) ) # p= OE higher
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO3")$`Non-Mitochondrial Oxygen Consumption`),
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO3.OE")$`Non-Mitochondrial Oxygen Consumption`) ) # p=lower
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO3")$`Non-Mitochondrial Oxygen Consumption`),
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO3.OE2")$`Non-Mitochondrial Oxygen Consumption`) ) # p= lower
t.test(as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO3")$`Non-Mitochondrial Oxygen Consumption`),
       as.numeric(subset(SH.HCT.UBE2H_Merge, `Group Name` == "KO3.OE5")$`Non-Mitochondrial Oxygen Consumption`) ) # p=0.04202 OE higher 



# Just look at non-AZ3146: Maximal Respiration
ggplot(subset(SH.HCT.UBE2H_Merge, AZ3146=="0"), 
       aes(x=`Group Name`, y=`Maximal Respiration`, fill=`Group Name`))+
  geom_boxplot()+
  scale_fill_manual(values= c("turquoise2", "turquoise4", 
                              "mediumorchid2", "mediumorchid4",
                              "mediumorchid2", "mediumorchid4","mediumorchid4",
                              "mediumorchid2", "mediumorchid4", "mediumorchid4"))+
  geom_point()+   
  ylab("Maxinal Respiration")+
  xlab("")+ 
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
# plot.Seahorse.HCT116.noAZ.MAX_Merge.pdf










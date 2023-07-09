# This is a copy of the R script used to generate the figures for, and perform the statistical analysis on, 
# the data for the perceived position of the target relative to the flanker ring that was reported by the 
# subject's in experiments 1-4 in Smithers et al.- "Large differences in target-flanker depth can increase 
# crowding- evidence from a multi-depth plane display". For an explanation of the statistical analysis conducted, 
# please refer to the statistical analysis section of the methods in the main manuscript.

# Script written by Dr Samuel P. Smithers, Northeastern University, 2022
# Last edited July 2023

# Corresponding authors SPS (s.smithers@northeastern.edu) and PJB (p.bex@northeastern.edu)

# Included in this script:
# - Code used to generate Figures 8 and 9 from the main manuscript.
# - Code used for the statistical analysis of perceived target position that is reported in the main manuscript.
# - Code used to generate Figure 8-figure supplements 1 and 2, and Figure 9-figure supplement 1 and 2.

# Session info:
# R version 4.2.2 (2022-10-31 ucrt)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19045)
# 
# Matrix products: default
# 
# locale:
#   [1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8    LC_MONETARY=English_United States.utf8
# [4] LC_NUMERIC=C                           LC_TIME=English_United States.utf8    
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] bannerCommenter_1.0.0 DHARMa_0.4.6          lme4_1.1-31           Matrix_1.5-1          cowplot_1.1.1         gridExtra_2.3        
# [7] svglite_2.1.0         ggplot2_3.4.0        
# 
# loaded via a namespace (and not attached):
#   [1] Rcpp_1.0.9        pillar_1.8.1      compiler_4.2.2    nloptr_2.0.3      tools_4.2.2       boot_1.3-28       lifecycle_1.0.3   tibble_3.2.1     
# [9] nlme_3.1-160      gtable_0.3.1      lattice_0.20-45   pkgconfig_2.0.3   rlang_1.1.0       cli_3.6.1         rstudioapi_0.14   withr_2.5.0      
# [17] dplyr_1.1.1       generics_0.1.3    vctrs_0.6.1       systemfonts_1.0.4 grid_4.2.2        tidyselect_1.2.0  glue_1.6.2        R6_2.5.1         
# [25] fansi_1.0.3       minqa_1.2.5       magrittr_2.0.3    scales_1.2.1      MASS_7.3-58.1     splines_4.2.2     colorspace_2.0-3  utf8_1.2.2       
# [33] munsell_0.5.0       

##Required dependencies/packages. 
#For plots
library(ggplot2) # For graphs
library(svglite) # Required to save .svg files (optional)
library(gridExtra) # For arranging/combining plots (optional)
library(cowplot) # for ggsave and arranging plots (optional)
#For stats
library(lme4) # Used for mixed effects models
library(DHARMa) # For residual diagnostics
#Other
library(bannerCommenter) # For the section banners (optional)

rm(list=ls(all=TRUE))

setwd("***Pathway to folder containing this R script and all '_accumulatedPerErr.csv' files***")

#############################################################################################################
#############################################################################################################
###                                                                                                       ###
###  PERCEIVED TARGET POSITION RELATIVE TO FLANKER RING GRAPHS FOR THE MAIN MANUSCRIPT (FIGURES 8 AND 9)  ###
###                                                                                                       ###
#############################################################################################################
#############################################################################################################
banner("Perceived Target Position Relative to Flanker Ring Graphs for the Main Manuscript (Figures 8 and 9)", emph = TRUE)

#Colour blind friendly palette
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#D55E00", "#0072B2", "#CC79A7","#F0E442")

#Create graph theme
Graph.theme<-theme_bw()+
  theme(axis.text=element_text(size=25,colour="black"),
        axis.title=element_text(size=30, colour="black"),
        axis.title.y=element_text(vjust=1),axis.title.x=element_text(vjust=-2),
        plot.margin = unit(c(1,1,1,1), "cm"),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20),
        legend.key.size=unit(1, "cm"),
        legend.background=element_rect(),
        legend.key = element_blank(),
        legend.position= "right",
        axis.line = element_line(colour = "black", linewidth=1),
        axis.ticks = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.background = element_rect(fill="white"),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        strip.background =element_rect(fill="white", colour= "black", linewidth=1),
        strip.text=element_text(face="bold", size=25,colour = "black"),
        plot.title = element_text(size=30, face="bold"))

AddExtraCols <- function(Inputdata){
  # Make column for flanker depth in m
  Inputdata["FlankerDepth_in_m"] <- Inputdata$Flanker_Screen
  Inputdata$FlankerDepth_in_m[Inputdata$FlankerDepth_in_m == "far"] <- "4"
  Inputdata$FlankerDepth_in_m[Inputdata$FlankerDepth_in_m == "mid"] <- "1.26"
  Inputdata$FlankerDepth_in_m[Inputdata$FlankerDepth_in_m == "near"] <- "0.4"
  # Make column for target depth in m
  Inputdata["TargetDepth_in_m"] <- Inputdata$Target_Screen
  Inputdata$TargetDepth_in_m[Inputdata$TargetDepth_in_m == "far"] <- "4"
  Inputdata$TargetDepth_in_m[Inputdata$TargetDepth_in_m == "mid"] <- "1.26"
  Inputdata$TargetDepth_in_m[Inputdata$TargetDepth_in_m == "near"] <- "0.4"
  # Make column for fixation depth in m
  Inputdata["FixationDepth_in_m"] <- Inputdata$Fixation_Screen
  Inputdata$FixationDepth_in_m[Inputdata$FixationDepth_in_m == "far"] <- "4"
  Inputdata$FixationDepth_in_m[Inputdata$FixationDepth_in_m == "mid"] <- "1.26"
  Inputdata$FixationDepth_in_m[Inputdata$FixationDepth_in_m == "near"] <- "0.4"
  return(Inputdata)
}

#Function to calculate the proportion of trials in which the subject reported seeing the target inside the flanker ring 
ProportionInsideT <- function(Inputdata, repeatno.){
  # Turn counts into proportion 
  Inputdata["PerceivedPosCount_insideT"] <- Inputdata$PerceivedPosCount_TCenter + Inputdata$PerceivedPosCount_TOffCenter
  Inputdata["PerceivedPosCount_outsideT"] <- Inputdata$PerceivedPosCount_TObstructed + Inputdata$PerceivedPosCount_TOutside + Inputdata$PerceivedPosCount_noRingOrUnsure
  Inputdata["Proportion_insideT"] <- Inputdata$PerceivedPosCount_insideT/repeatno.
  return(Inputdata)
}

##================================================================================================================
##  GRAPHS: Experiment 1: Target and fixation at 1.26 m with flanker ring presented at each depth  (Figure 8A)   =
##================================================================================================================
boxup("GRAPHS: Experiment 1: Target and fixation at 1.26 m with flanker ring presented at each depth  (Figure 8A)", bandChar = "=")

Exp1_data <-read.csv("Exp1_accumulatedPerErr.csv", header=T) # Load data
Exp1_data <- AddExtraCols(Exp1_data)
N <- nrow(unique(Exp1_data["ID"])) 
print(N) # n=22
repeatno. = 10 
Exp1_data <- ProportionInsideT(Exp1_data,repeatno.)

Exp1_PropInsideT <- ggplot(Exp1_data, aes(x=factor(Target_Flanker_Spacing_in_deg), y= Proportion_insideT, fill= FlankerDepth_in_m)) + 
  geom_boxplot(lwd=0.75, fatten =2, outlier.shape=NA) + xlab("Target-flanker spacing (deg)") + 
  geom_point(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 1), 
             shape =21,  size=1.5, alpha=0.5, colour = 'black') +
  facet_wrap(~FixationDepth_in_m,ncol = 3) +
  scale_y_continuous(breaks= seq(0,1,0.2), limits=c(-0.02,1.02), expand=c(0,0),  name="Proportion of trials") + 
  guides(fill=guide_legend(title="Flanker depth (m)")) + Graph.theme + scale_fill_manual(values=cbPalette) +
  scale_x_discrete(labels=c("0.625" = "0.625", "1.25" = "1.25", "2.5" = "2.5", "5" = "5", "Inf" = "no flankers")) 
Exp1_PropInsideT

##=====================================================================================================================
##  GRAPHS: Experiment 2: Target and fixation at 0.4 m or 4 m with flanker ring presented at each depth (Figure 8B)   =
##=====================================================================================================================
boxup("GRAPHS: Experiment 2: Target and fixation at 0.4 m or 4 m with flanker ring presented at each depth (Figure 8B)", bandChar = "=")

Exp2_data <-read.csv("Exp2_accumulatedPerErr.csv", header=T) # Load data
Exp2_data <- AddExtraCols(Exp2_data)
N <- nrow(unique(Exp2_data["ID"])) 
print(N) # n=19
repeatno. = 10 
Exp2_data <- ProportionInsideT(Exp2_data,repeatno.)

Exp2_PropInsideT <- ggplot(Exp2_data, aes(x=factor(Target_Flanker_Spacing_in_deg), y= Proportion_insideT, fill= FlankerDepth_in_m)) + 
  geom_boxplot(lwd=0.75, fatten =2, outlier.shape=NA) + xlab("Target-flanker spacing (deg)") + 
  geom_point(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 1), 
             shape =21,  size=1.5, alpha=0.5, colour = 'black') +
  facet_wrap(~FixationDepth_in_m,ncol = 3) +
  scale_y_continuous(breaks= seq(0,1,0.2), limits=c(-0.02,1.02), expand=c(0,0),  name="Proportion of trials") + 
  guides(fill=guide_legend(title="Flanker depth (m)")) + Graph.theme + scale_fill_manual(values=cbPalette) +
  scale_x_discrete(labels=c("0.625" = "0.625", "1.25" = "1.25", "2.5" = "2.5", "5" = "5", "Inf" = "no flankers")) 
Exp2_PropInsideT

##===============================================================================================================
##  GRAPHS: Experiment 3: Flanker ring and fixation at 1.26 m with target presented at each depth (Figure 9A)   =
##===============================================================================================================
boxup("GRAPHS: Experiment 3: Flanker ring and fixation at 1.26 m with target presented at each depth (Figure 9A)", bandChar = "=")

Exp3_data <-read.csv("Exp3_accumulatedPerErr.csv", header=T)
Exp3_data <- AddExtraCols(Exp3_data)
N <- nrow(unique(Exp3_data["ID"]))
print(N) # n=21
repeatno. = 10 
Exp3_data <- ProportionInsideT(Exp3_data,repeatno.)

Exp3_PropInsideT <- ggplot(Exp3_data, aes(x=factor(Target_Flanker_Spacing_in_deg), y= Proportion_insideT, fill= TargetDepth_in_m)) + 
  geom_boxplot(lwd=0.75, fatten =2, outlier.shape=NA) + xlab("Target-flanker spacing (deg)") + 
  geom_point(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 1), 
             shape =21,  size=1.5, alpha=0.5, colour = 'black') +
  facet_wrap(~FixationDepth_in_m,ncol = 3) +
  scale_y_continuous(breaks= seq(0,1,0.2), limits=c(-0.02,1.02), expand=c(0,0),  name="Proportion of trials") + 
  guides(fill=guide_legend(title="Target depth (m)")) + Graph.theme + scale_fill_manual(values=cbPalette) +
  scale_x_discrete(labels=c("0.625" = "0.625", "1.25" = "1.25", "2.5" = "2.5", "5" = "5", "Inf" = "no flankers")) 
Exp3_PropInsideT

##=====================================================================================================================
##  GRAPHS: Experiment 4: Flanker ring and fixation at 0.4 m or 4 m with target presented at each depth (Figure 9B)   =
##=====================================================================================================================
boxup("GRAPHS: Experiment 4: Flanker ring and fixation at 0.4 m or 4 m with target presented at each depth (Figure 9B)", bandChar = "=")

Exp4_data <-read.csv("Exp4_accumulatedPerErr.csv", header=T)
Exp4_data <- AddExtraCols(Exp4_data)
N <- nrow(unique(Exp4_data["ID"]))
print(N) # n=21
repeatno. = 8
Exp4_data <- ProportionInsideT(Exp4_data,repeatno.)

Exp4_PropInsideT <- ggplot(Exp4_data, aes(x=factor(Target_Flanker_Spacing_in_deg), y= Proportion_insideT, fill= TargetDepth_in_m)) + 
  geom_boxplot(lwd=0.75, fatten =2, outlier.shape=NA) + xlab("Target-flanker spacing (deg)") + 
  geom_point(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 1), 
             shape =21,  size=1.5, alpha=0.5, colour = 'black') +
  facet_wrap(~FixationDepth_in_m,ncol = 3) + 
  scale_y_continuous(breaks= seq(0,1,0.2), limits=c(-0.02,1.02), expand=c(0,0),  name="Proportion of trials") + 
  guides(fill=guide_legend(title="Target depth (m)")) + Graph.theme + scale_fill_manual(values=cbPalette) +
  scale_x_discrete(labels=c("0.625" = "0.625", "1.25" = "1.25", "2.5" = "2.5", "5" = "5", "Inf" = "no flankers")) 
Exp4_PropInsideT

##---------------------------------------------------------------
##              GRAPHS: Save graphs as .svg files               -
##---------------------------------------------------------------
boxup("GRAPHS: Save graphs as .svg files", bandChar = "-")
# The graphs are then saved as a .svg file to enable formatting, labeling, and other 
# edits required to prepare the figures for the paper.

CombinedPlots <- plot_grid(Exp1_PropInsideT, Exp2_PropInsideT,Exp3_PropInsideT,Exp4_PropInsideT, ncol = 2, rel_widths = c(5/8, 1))
# Combined plots will be saved to your working directory
ggsave(file ="Smithers et al- figures 8-9 (PropInsideRing).svg", device = 'svg', plot = CombinedPlots, width = 25, height = 15)

############################################################################################################
############################################################################################################
###                                                                                                      ###
###  STATISTICAL ANALYSIS OF PERCEIVED TARGET POSITION RELATIVE TO FLANKER RING FOR THE MAIN MANUSCRIPT  ###
###                                                                                                      ###
############################################################################################################
############################################################################################################
banner("Statistical Analysis of Perceived Target Position Relative to Flanker Ring for the Main Manuscript", emph = TRUE)
rm(list=ls(all=TRUE))

#Function to convert choice to a binary response variable  
BinaryChoiceCol <- function(Inputdata){
  # Make column for flanker depth in m
  Inputdata["BinaryChoice_TinsideF"] <- Inputdata$Perceived_TargetRing_Position
  Inputdata$BinaryChoice_TinsideF[Inputdata$BinaryChoice_TinsideF == "target in center of ring"] <- 1
  Inputdata$BinaryChoice_TinsideF[Inputdata$BinaryChoice_TinsideF == "target not in center of ring"] <- 1
  Inputdata$BinaryChoice_TinsideF[Inputdata$BinaryChoice_TinsideF == "ring obstructs target"] <- 0
  Inputdata$BinaryChoice_TinsideF[Inputdata$BinaryChoice_TinsideF == "target outside ring"] <- 0
  Inputdata$BinaryChoice_TinsideF[Inputdata$BinaryChoice_TinsideF == "unsure or no ring"] <- 0
  return(Inputdata)
}

##==================================================================================================
##  STATS: Experiment 1: Target and fixation at 1.26 m with flanker ring presented at each depth   =
##==================================================================================================
boxup("STATS: Experiment 1: Target and fixation at 1.26 m with flanker ring presented at each depth", bandChar = "=")

#Set pathway to folder containing the "..._rawData.csv" files for Experiment 1 (the folder must contain ONLY the raw data for Exp1) 
dir <- "***Pathway to folder containing all of the '_rawData.csv' files for Experiment 1 (and only Experiment 1)***"

all_files <- list.files(dir) # List all files in directory
Collated_Raw_Data <- do.call(rbind, lapply(paste(dir,all_files,sep ="/"), read.csv)) # Collate raw data from all subjects
df <- BinaryChoiceCol(Collated_Raw_Data)
df$Target_Screen <- factor(df$Target_Screen)
df$Flanker_Screen <- factor(df$Flanker_Screen)
df$ID <- factor(df$ID)
df$BinaryChoice_TinsideF <- factor(df$BinaryChoice_TinsideF)
df <- subset(df, Target_Flanker_Spacing_in_deg!=Inf) # Remove control (i.e. no flanker condition)

#Check data format is correct
str(df)

#Construct full mixed effects binary logistic regression model 
m0 <- glmer(BinaryChoice_TinsideF ~ Target_Flanker_Spacing_in_deg+Flanker_Screen +
              Target_Flanker_Spacing_in_deg:Flanker_Screen + (1|ID) , family="binomial", data=df)

#Set seed used to simulate residuals 
set.seed(2022) # Here we use the year this script was written  
# DHARMa uses a simulation-based approach for residual diagnostics. We therefore use a fixed seed which should ensure that DHARMa 
# produces reproducible results when we run "simulateResiduals()".

#Plot model residual and check residual diagnostics
simulationOutput <- simulateResiduals(fittedModel = m0, plot = T, seed = NULL)
# QQ plot looks good and there is no sign of a pattern in the residuals. Diagnostic tests are all non-significant.

#**************************************************************
#m0 is the min model. These are the stats reported in the paper.
m1 <-update(m0,~.-Target_Flanker_Spacing_in_deg:Flanker_Screen)
anova(m0, m1, test = "LRT") # Chisq= 36.119, Df=2, p = 1.435e-08 ***
#**************************************************************
rm(dir,all_files,Collated_Raw_Data,df,m0,m1)

##=========================================================================================================
##  STATS: Experiment 2: Target and fixation at 0.4 m or 4 m with flanker ring presented at each depth.   =
##=========================================================================================================
boxup("STATS: Experiment 2: Target and fixation at 0.4 m or 4 m with flanker ring presented at each depth.", bandChar = "=")

#Set pathway to folder containing the "..._rawData.csv" files for Experiment 2 (the folder must contain ONLY the raw data for Exp2) 
dir <- "***Pathway to folder containing all of the '_rawData.csv' files for Experiment 2 (and only Experiment 2)***"

all_files <- list.files(dir) # List all files in directory
Collated_Raw_Data <- do.call(rbind, lapply(paste(dir,all_files,sep ="/"), read.csv)) # Collate raw data from all subjects
df <- BinaryChoiceCol(Collated_Raw_Data)
df$Target_Screen <- factor(df$Target_Screen)
df$Flanker_Screen <- factor(df$Flanker_Screen)
df$ID <- factor(df$ID)
df$BinaryChoice_TinsideF <- factor(df$BinaryChoice_TinsideF)
df <- subset(df, Target_Flanker_Spacing_in_deg!=Inf) # Remove control (i.e. no flanker condition)

#Check data format is correct
str(df)

#Construct full mixed effects binary logistic regression model 
m0 <- glmer(BinaryChoice_TinsideF ~ Target_Flanker_Spacing_in_deg+Target_Screen+Flanker_Screen +
              Target_Flanker_Spacing_in_deg:Target_Screen +
              Target_Flanker_Spacing_in_deg:Flanker_Screen + 
              Target_Screen:Flanker_Screen + (1|ID) , family="binomial", data=df)

#Set seed used to simulate residuals 
set.seed(2022) # Here we use the year this script was written. 
# DHARMa uses a simulation-based approach for residual diagnostics. We therefore use a fixed seed which should ensure that DHARMa 
# produces reproducible results when we run "simulateResiduals()".

#Plot model residual and check residual diagnostics
simulationOutput <- simulateResiduals(fittedModel = m0, plot = T, seed = NULL)
# QQ plot looks good and there is no sign of a pattern in the residuals. Diagnostic tests are all non-significant. 

m1a <-update(m0,~.-Target_Screen:Flanker_Screen)
anova(m0, m1a, test = "LRT") # p< 2.2e-16 ***
m1b <-update(m0,~.-Target_Flanker_Spacing_in_deg:Flanker_Screen)
anova(m0, m1b, test = "LRT") #  Chi = 0.2412, df = 2, p =  0.8864 # Reported in paper
m1c <-update(m0,~.-Target_Flanker_Spacing_in_deg:Target_Screen)
# We get a convergence warning. While this is not necessarily something to worry about, we can try adding more iterations by 
# restarting from previous fit and using a different optimizer.
ss <- getME(m1c,c("theta","fixef"))
m1c <- update(m1c,start=ss,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
# That did the trick. 
anova(m0, m1c, test = "LRT") # p = 0.3589
# Drop Target_Flanker_Spacing_in_deg:Flanker_Screen

m2a <-update(m1b,~.-Target_Screen:Flanker_Screen)
anova(m1b, m2a, test = "LRT") # p < 2.2e-16 ***
m2b <-update(m1b,~.-Target_Flanker_Spacing_in_deg:Target_Screen)
anova(m1b, m2b, test = "LRT") # Chi = 1.64, df = 1, p =  0.2003 # Reported in paper 
# Drop Target_Flanker_Spacing_in_deg:Target_Screen

#**************************************************************
#m2b is the min model. These are the stats reported in the paper.
m3a <-update(m2b,~.-Target_Screen:Flanker_Screen)
anova(m2b, m3a, test = "LRT") # Chisq= 374.26, Df=2, p < 2.2e-16 ***
m3b <-update(m2b,~.-Target_Flanker_Spacing_in_deg)
anova(m2b, m3b, test = "LRT") # Chisq= 178.65, Df=1, p < 2.2e-16 ***
#**************************************************************
rm(dir,all_files,Collated_Raw_Data,df,m0,m1a,m1b,m1c,m2a,m2b,m3a,m3b,ss)

##==================================================================================================
##  STATS: Experiment 3: Flanker ring and fixation at 1.26 m with target presented at each depth   =
##==================================================================================================
boxup("STATS: Experiment 3: Flanker ring and fixation at 1.26 m with target presented at each depth", bandChar = "=")

#Set pathway to folder containing the "..._rawData.csv" files for Experiment 3 (the folder must contain ONLY the raw data for Exp3) 
dir <- "***Pathway to folder containing all of the '_rawData.csv' files for Experiment 3 (and only Experiment 3)***"

all_files <- list.files(dir) # List all files in directory
Collated_Raw_Data <- do.call(rbind, lapply(paste(dir,all_files,sep ="/"), read.csv)) # Collate raw data from all subjects
df <- BinaryChoiceCol(Collated_Raw_Data)
df$Target_Screen <- factor(df$Target_Screen)
df$Flanker_Screen <- factor(df$Flanker_Screen)
df$ID <- factor(df$ID)
df$BinaryChoice_TinsideF <- factor(df$BinaryChoice_TinsideF)
df <- subset(df, Target_Flanker_Spacing_in_deg!=Inf) # Remove control (i.e. no flanker condition)

#Check data format is correct
str(df)

#Construct full mixed effects binary logistic regression model 
m0 <- glmer(BinaryChoice_TinsideF ~ Target_Flanker_Spacing_in_deg+Target_Screen +
              Target_Flanker_Spacing_in_deg:Target_Screen + (1|ID) , family="binomial", data=df)

#Set seed used to simulate residuals 
set.seed(2022) # Here we use the year this script was written. 
# DHARMa uses a simulation-based approach for residual diagnostics. We therefore use a fixed seed which should ensure that DHARMa 
# produces reproducible results when we run "simulateResiduals()".

#Plot model residual and check residual diagnostics
simulationOutput <- simulateResiduals(fittedModel = m0, plot = T, seed = NULL)
# QQ plot looks good and there is no sign of a pattern in the residuals. Diagnostic tests are all non-significant. 

#**************************************************************
#m0 is the min model. These are the stats reported in the paper.
m1 <-update(m0,~.-Target_Flanker_Spacing_in_deg:Target_Screen)
anova(m0, m1, test = "LRT") # Chisq= 100.38, Df=2, p < 2.2e-16 ***
#**************************************************************
rm(dir,all_files,Collated_Raw_Data,df,m0,m1)

##========================================================================================================
##  STATS: Experiment 4: Flanker ring and fixation at 0.4 m or 4 m with target presented at each depth   =
##========================================================================================================
boxup("STATS: Experiment 4: Flanker ring and fixation at 0.4 m or 4 m with target presented at each depth", bandChar = "=")

#Set pathway to folder containing the "..._rawData.csv" files for Experiment 4 (the folder must contain ONLY the raw data for Exp4) 
dir <- "***Pathway to folder containing all of the '_rawData.csv' files for Experiment 4 (and only Experiment 4)***"

all_files <- list.files(dir) # List all files in directory
Collated_Raw_Data <- do.call(rbind, lapply(paste(dir,all_files,sep ="/"), read.csv)) # Collate raw data from all subjects
df <- BinaryChoiceCol(Collated_Raw_Data)
df$Target_Screen <- factor(df$Target_Screen)
df$Flanker_Screen <- factor(df$Flanker_Screen)
df$ID <- factor(df$ID)
df$BinaryChoice_TinsideF <- factor(df$BinaryChoice_TinsideF)
df <- subset(df, Target_Flanker_Spacing_in_deg!=Inf) # Remove control (i.e. no flanker condition)

#Check data format is correct
str(df)

#Construct full mixed effects binary logistic regression model 
m0 <- glmer(BinaryChoice_TinsideF ~ Target_Flanker_Spacing_in_deg+Target_Screen+Flanker_Screen +
              Target_Flanker_Spacing_in_deg:Target_Screen+
              Target_Flanker_Spacing_in_deg:Flanker_Screen +
              Target_Screen:Flanker_Screen + (1|ID) , family="binomial", data=df)
# We get a convergence warning. While this is not necessarily something to worry about we can try adding more iterations by 
# restarting from previous fit and using a different optimizer.
ss <- getME(m0,c("theta","fixef"))
m0 <- update(m0,start=ss,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
# That did the trick. 

#Set seed used to simulate residuals 
set.seed(2022) #Here we use the year this script was written.  
# DHARMa uses a simulation-based approach for residual diagnostics. We therefore use a fixed seed which should ensure that DHARMa 
# produces reproducible results when we run "simulateResiduals()".

#Plot model residual and check residual diagnostics
simulationOutput <- simulateResiduals(fittedModel = m0, plot = T, seed = NULL)
# QQ plot looks good and there is no sign of a pattern in the residuals. Diagnostic tests are all non-significant.

#**************************************************************
#m0 is the min model. These are the stats reported in the paper.
m1a <-update(m0,~.-Target_Screen:Flanker_Screen)
anova(m0, m1a, test = "LRT") # Chisq= 779.49, Df=2, p < 2.2e-16 ***
m1b <-update(m0,~.-Target_Flanker_Spacing_in_deg:Flanker_Screen)
anova(m0, m1b, test = "LRT") # Chisq= 56.229, Df=1, p = 6.45e-14 ***
m1c <-update(m0,~.-Target_Flanker_Spacing_in_deg:Target_Screen)
anova(m0, m1c, test = "LRT") # Chisq= 23.267, Df=2, p = 8.862e-06 ***
#**************************************************************
rm(dir,all_files,Collated_Raw_Data,df,m0,m1a,m1b,m1c,ss)


##################################################################################################################################################
##################################################################################################################################################
###                                                                                                                                            ###
###  PERCEIVED TARGET POSITION RELATIVE TO FLANKER RING GRAPHS FOR FIGURE 8-FIGURE SUPPLEMENT 1 AND 2, AND FIGURE 9-FIGURE SUPPLEMENT 1 AND 2  ###
###                                                                                                                                            ###
##################################################################################################################################################
##################################################################################################################################################
banner("Perceived Target Position Relative to Flanker Ring Graphs for Figure 8-figure supplement 1 and 2, and Figure 9-figure supplement 1 and 2", emph = TRUE)

# Colour blind friendly palette
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#D55E00", "#0072B2", "#CC79A7","#F0E442")

#Create graph theme
Graph.theme<-theme_bw()+
  theme(axis.text=element_text(size=25,colour="black"),
        axis.title=element_text(size=30, colour="black"),
        axis.title.y=element_text(vjust=1),axis.title.x=element_text(vjust=-2),
        plot.margin = unit(c(1,1,1,1), "cm"),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20),
        legend.key.size=unit(1, "cm"),
        legend.background=element_rect(),
        legend.key = element_blank(),
        legend.position= "right",
        axis.line = element_line(colour = "black", linewidth=1),
        axis.ticks = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.background = element_rect(fill="white"),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        strip.background =element_rect(fill="white", colour= "black", linewidth=1),
        strip.text=element_text(face="bold", size=25,colour = "black"),
        plot.title = element_text(size=30, face="bold"))

AddExtraCols <- function(Inputdata){
  # Make column for flanker depth in m
  Inputdata["FlankerDepth_in_m"] <- Inputdata$Flanker_Screen
  Inputdata$FlankerDepth_in_m[Inputdata$FlankerDepth_in_m == "far"] <- "4"
  Inputdata$FlankerDepth_in_m[Inputdata$FlankerDepth_in_m == "mid"] <- "1.26"
  Inputdata$FlankerDepth_in_m[Inputdata$FlankerDepth_in_m == "near"] <- "0.4"
  # Make column for target depth in m
  Inputdata["TargetDepth_in_m"] <- Inputdata$Target_Screen
  Inputdata$TargetDepth_in_m[Inputdata$TargetDepth_in_m == "far"] <- "4"
  Inputdata$TargetDepth_in_m[Inputdata$TargetDepth_in_m == "mid"] <- "1.26"
  Inputdata$TargetDepth_in_m[Inputdata$TargetDepth_in_m == "near"] <- "0.4"
  # Make column for fixation depth in m
  Inputdata["FixationDepth_in_m"] <- Inputdata$Fixation_Screen
  Inputdata$FixationDepth_in_m[Inputdata$FixationDepth_in_m == "far"] <- "4"
  Inputdata$FixationDepth_in_m[Inputdata$FixationDepth_in_m == "mid"] <- "1.26"
  Inputdata$FixationDepth_in_m[Inputdata$FixationDepth_in_m == "near"] <- "0.4"
  return(Inputdata)
}

#Function to make count table for plotting 
CountDataFunction <- function(Inputdata){
  Count <- c()
  Perception <- c()
  TargetFlankerSpacing <-c()
  FlankerDepth <- c()
  TargetDepth <- c()
  FixationDepth <- c()
  Count <-append(Count,Inputdata$PerceivedPosCount_TCenter)
  Perception <- append(Perception, rep("Target in center of ring",length(Inputdata$PerceivedPosCount_TCenter)))
  Count <-append(Count,Inputdata$PerceivedPosCount_TOffCenter)
  Perception <- append(Perception, rep("Target not in center of ring",length(Inputdata$PerceivedPosCount_TOffCenter)))
  Count <-append(Count,Inputdata$PerceivedPosCount_TObstructed)
  Perception <- append(Perception, rep("Ring obstructs target",length(Inputdata$PerceivedPosCount_TObstructed)))
  Count <-append(Count,Inputdata$PerceivedPosCount_TOutside)
  Perception <- append(Perception, rep("Target outside ring",length(Inputdata$PerceivedPosCount_TOutside)))
  Count <-append(Count,Inputdata$PerceivedPosCount_noRingOrUnsure)
  Perception <- append(Perception, rep("Unsure or no ring",length(Inputdata$PerceivedPosCount_noRingOrUnsure)))
  TargetFlankerSpacing <- append(TargetFlankerSpacing, rep(Inputdata$Target_Flanker_Spacing_in_deg,5))
  FlankerDepth <- append(FlankerDepth,rep(Inputdata$FlankerDepth_in_m,5))
  TargetDepth <- append(TargetDepth,rep(Inputdata$TargetDepth_in_m,5))
  FixationDepth <- append(FixationDepth,rep(Inputdata$FixationDepth_in_m,5))
  CountData <- data.frame(Count,Perception, TargetFlankerSpacing, FlankerDepth, TargetDepth, FixationDepth)
  return(CountData) 
}

##==================================================================================================================================================
##  SUPPLEMENTARY GRAPHS: Experiment 1: Target and fixation at 1.26 m with flanker ring presented at each depth  (Figure 8-figure supplement 1).   =
##==================================================================================================================================================
boxup("SUPPLEMENTARY GRAPHS: Experiment 1: Target and fixation at 1.26 m with flanker ring presented at each depth  (Figure 8-figure supplement 1).", bandChar = "=")

Exp1_data <-read.csv("Exp1_accumulatedPerErr.csv", header=T) # Load data
Exp1_data <- AddExtraCols(Exp1_data)

CountData <- CountDataFunction(Exp1_data)

repeatno. = 10  

PlotCountData <- subset(CountData, Perception == "Target in center of ring")
Exp1_plot_TCenter <- ggplot(PlotCountData, aes(x=as.factor(TargetFlankerSpacing), y= (Count/repeatno.), fill= FlankerDepth)) + 
  geom_boxplot(lwd=0.75, fatten =2, outlier.shape=NA) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 1), shape =21,  size=1.5, alpha=0.5, colour = 'black') +
  xlab("Target-flanker spacing (deg)") + Graph.theme + scale_fill_manual(values=cbPalette) +
  facet_wrap(~FixationDepth,ncol = 3) + ggtitle("Target in center of flanker ring") +
  scale_y_continuous(breaks= seq(0,1,0.2), limits=c(-0.02,1.02), expand=c(0,0),  name="Proportion of trials") + 
  guides(fill=guide_legend(title="Flanker depth (m)")) +
  scale_x_discrete(labels=c("0.625" = "0.625", "1.25" = "1.25", "2.5" = "2.5", "5" = "5", "Inf" = "no flankers")) 
Exp1_plot_TCenter

PlotCountData <- subset(CountData, Perception == "Target not in center of ring")
Exp1_plot_TOffCenter <- ggplot(PlotCountData, aes(x=as.factor(TargetFlankerSpacing), y= (Count/repeatno.), fill= FlankerDepth)) + 
  geom_boxplot(lwd=0.75, fatten =2, outlier.shape=NA) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 1), shape =21,  size=1.5, alpha=0.5, colour = 'black') +
  xlab("Target-flanker spacing (deg)") + Graph.theme + scale_fill_manual(values=cbPalette) +
  facet_wrap(~FixationDepth,ncol = 3) + ggtitle("Target inside, but not in center, of flanker ring") +
  scale_y_continuous(breaks= seq(0,1,0.2), limits=c(-0.02,1.02), expand=c(0,0),  name="Proportion of trials") + 
  guides(fill=guide_legend(title="Flanker depth (m)")) +
  scale_x_discrete(labels=c("0.625" = "0.625", "1.25" = "1.25", "2.5" = "2.5", "5" = "5", "Inf" = "no flankers")) 
Exp1_plot_TOffCenter

PlotCountData <- subset(CountData, Perception == "Ring obstructs target")
Exp1_plot_TObstructed <- ggplot(PlotCountData, aes(x=as.factor(TargetFlankerSpacing), y= (Count/repeatno.), fill= FlankerDepth)) + 
  geom_boxplot(lwd=0.75, fatten =2, outlier.shape=NA) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 1), shape =21,  size=1.5, alpha=0.5, colour = 'black') +
  xlab("Target-flanker spacing (deg)") + Graph.theme + scale_fill_manual(values=cbPalette) +
  facet_wrap(~FixationDepth,ncol = 3) + ggtitle("Flanker ring obstructs target") +
  scale_y_continuous(breaks= seq(0,1,0.2), limits=c(-0.02,1.02), expand=c(0,0),  name="Proportion of trials") + 
  guides(fill=guide_legend(title="Flanker depth (m)")) +
  scale_x_discrete(labels=c("0.625" = "0.625", "1.25" = "1.25", "2.5" = "2.5", "5" = "5", "Inf" = "no flankers")) 
Exp1_plot_TObstructed

PlotCountData <- subset(CountData, Perception == "Target outside ring")
Exp1_plot_TOutside <- ggplot(PlotCountData, aes(x=as.factor(TargetFlankerSpacing), y= (Count/repeatno.), fill= FlankerDepth)) + 
  geom_boxplot(lwd=0.75, fatten =2, outlier.shape=NA) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 1), shape =21,  size=1.5, alpha=0.5, colour = 'black') +
  xlab("Target-flanker spacing (deg)") + Graph.theme + scale_fill_manual(values=cbPalette) +
  facet_wrap(~FixationDepth,ncol = 3) + ggtitle("Target outside flanker ring") +
  scale_y_continuous(breaks= seq(0,1,0.2), limits=c(-0.02,1.02), expand=c(0,0),  name="Proportion of trials") + 
  guides(fill=guide_legend(title="Flanker depth (m)")) +
  scale_x_discrete(labels=c("0.625" = "0.625", "1.25" = "1.25", "2.5" = "2.5", "5" = "5", "Inf" = "no flankers")) 
Exp1_plot_TOutside

PlotCountData <- subset(CountData, Perception == "Unsure or no ring")
Exp1_plot_noRingOrUnsure <- ggplot(PlotCountData, aes(x=as.factor(TargetFlankerSpacing), y= (Count/repeatno.), fill= FlankerDepth)) + 
  geom_boxplot(lwd=0.75, fatten =2, outlier.shape=NA) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 1), shape =21,  size=1.5, alpha=0.5, colour = 'black') +
  xlab("Target-flanker spacing (deg)") + Graph.theme + scale_fill_manual(values=cbPalette) +
  facet_wrap(~FixationDepth,ncol = 3) + ggtitle("No flanker ring or unsure") +
  scale_y_continuous(breaks= seq(0,1,0.2), limits=c(-0.02,1.02), expand=c(0,0),  name="Proportion of trials") + 
  guides(fill=guide_legend(title="Flanker depth (m)")) +
  scale_x_discrete(labels=c("0.625" = "0.625", "1.25" = "1.25", "2.5" = "2.5", "5" = "5", "Inf" = "no flankers")) 
Exp1_plot_noRingOrUnsure

##---------------------------------------------------------------
##        SUPPLEMENTARY GRAPHS: Save figure S1 as .svg file       -
##---------------------------------------------------------------  
boxup("SUPPLEMENTARY GRAPHS: Save figure S1 as .svg files", bandChar = "-")
# The graphs are saved as a .svg file to enable formatting, labeling, and other 
# edits required to prepare the figures for the supplementary material.

#Save plots
JoinedPlots <- grid.arrange(Exp1_plot_TCenter, Exp1_plot_TOffCenter, Exp1_plot_TObstructed, 
                            Exp1_plot_TOutside, Exp1_plot_noRingOrUnsure, ncol = 2)
ggsave(file ="Smithers et al- Figure 8-figure supplement 1 (PerceivedTargetPos).svg", device = 'svg', plot = JoinedPlots, width = 23, height = 22)


##=======================================================================================================================================================
##  SUPPLEMENTARY GRAPHS: Experiment 2: Target and fixation at 0.4 m or 4 m with flanker ring presented at each depth (Figure 8-figure supplement 2).   =
##=======================================================================================================================================================
boxup("SUPPLEMENTARY GRAPHS: Experiment 2: Target and fixation at 0.4 m or 4 m with flanker ring presented at each depth (Figure 8-figure supplement 2).", bandChar = "=")

Exp2_data <-read.csv("Exp2_accumulatedPerErr.csv", header=T) # Load data
Exp2_data <- AddExtraCols(Exp2_data)

CountData <- CountDataFunction(Exp2_data)

repeatno. = 10

PlotCountData <- subset(CountData, Perception == "Target in center of ring")
Exp2_plot_TCenter <- ggplot(PlotCountData, aes(x=as.factor(TargetFlankerSpacing), y= (Count/repeatno.), fill= FlankerDepth)) + 
  geom_boxplot(lwd=0.75, fatten =2, outlier.shape=NA) + xlab("Target-flanker spacing (deg)") + 
  geom_point(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 1), shape =21,  size=1.5, alpha=0.5, colour = 'black') +
  facet_wrap(~FixationDepth,ncol = 3) + ggtitle("Target in center of flanker ring") +
  scale_y_continuous(breaks= seq(0,1,0.2), limits=c(-0.02,1.02), expand=c(0,0),  name="Proportion of trials") + 
  guides(fill=guide_legend(title="Flanker depth (m)")) + Graph.theme + scale_fill_manual(values=cbPalette) +
  scale_x_discrete(labels=c("0.625" = "0.625", "1.25" = "1.25", "2.5" = "2.5", "5" = "5", "Inf" = "no flankers")) 
Exp2_plot_TCenter

PlotCountData <- subset(CountData, Perception == "Target not in center of ring")
Exp2_plot_TOffCenter <- ggplot(PlotCountData, aes(x=as.factor(TargetFlankerSpacing), y= (Count/repeatno.), fill= FlankerDepth)) + 
  geom_boxplot(lwd=0.75, fatten =2, outlier.shape=NA) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 1), shape =21,  size=1.5, alpha=0.5, colour = 'black') +
  xlab("Target-flanker spacing (deg)") + Graph.theme + scale_fill_manual(values=cbPalette) +
  facet_wrap(~FixationDepth,ncol = 3) + ggtitle("Target inside, but not in center, of flanker ring") +
  scale_y_continuous(breaks= seq(0,1,0.2), limits=c(-0.02,1.02), expand=c(0,0),  name="Proportion of trials") + 
  guides(fill=guide_legend(title="Flanker depth (m)")) + 
  scale_x_discrete(labels=c("0.625" = "0.625", "1.25" = "1.25", "2.5" = "2.5", "5" = "5", "Inf" = "no flankers")) 
Exp2_plot_TOffCenter

PlotCountData <- subset(CountData, Perception == "Ring obstructs target")
Exp2_plot_TObstructed <- ggplot(PlotCountData, aes(x=as.factor(TargetFlankerSpacing), y= (Count/repeatno.), fill= FlankerDepth)) + 
  geom_boxplot(lwd=0.75, fatten =2, outlier.shape=NA) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 1), shape =21,  size=1.5, alpha=0.5, colour = 'black') +
  xlab("Target-flanker spacing (deg)") + Graph.theme + scale_fill_manual(values=cbPalette) +
  facet_wrap(~FixationDepth,ncol = 3) + ggtitle("Flanker ring obstructs target") +
  scale_y_continuous(breaks= seq(0,1,0.2), limits=c(-0.02,1.02), expand=c(0,0),  name="Proportion of trials") + 
  guides(fill=guide_legend(title="Flanker depth (m)")) + 
  scale_x_discrete(labels=c("0.625" = "0.625", "1.25" = "1.25", "2.5" = "2.5", "5" = "5", "Inf" = "no flankers")) 
Exp2_plot_TObstructed

PlotCountData <- subset(CountData, Perception == "Target outside ring")
Exp2_plot_TOutside <- ggplot(PlotCountData, aes(x=as.factor(TargetFlankerSpacing), y= (Count/repeatno.), fill= FlankerDepth)) + 
  geom_boxplot(lwd=0.75, fatten =2, outlier.shape=NA) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 1), shape =21,  size=1.5, alpha=0.5, colour = 'black') +
  xlab("Target-flanker spacing (deg)") + Graph.theme + scale_fill_manual(values=cbPalette) +
  facet_wrap(~FixationDepth,ncol = 3) + ggtitle("Target outside flanker ring") +
  scale_y_continuous(breaks= seq(0,1,0.2), limits=c(-0.02,1.02), expand=c(0,0),  name="Proportion of trials") + 
  guides(fill=guide_legend(title="Flanker depth (m)")) +
  scale_x_discrete(labels=c("0.625" = "0.625", "1.25" = "1.25", "2.5" = "2.5", "5" = "5", "Inf" = "no flankers")) 
Exp2_plot_TOutside

PlotCountData <- subset(CountData, Perception == "Unsure or no ring")
Exp2_plot_noRingOrUnsure <- ggplot(PlotCountData, aes(x=as.factor(TargetFlankerSpacing), y= (Count/repeatno.), fill= FlankerDepth)) + 
  geom_boxplot(lwd=0.75, fatten =2, outlier.shape=NA) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 1), shape =21,  size=1.5, alpha=0.5, colour = 'black') +
  xlab("Target-flanker spacing (deg)") + Graph.theme + scale_fill_manual(values=cbPalette) +
  facet_wrap(~FixationDepth,ncol = 3) + ggtitle("No flanker ring or unsure") +
  scale_y_continuous(breaks= seq(0,1,0.2), limits=c(-0.02,1.02), expand=c(0,0),  name="Proportion of trials") + 
  guides(fill=guide_legend(title="Flanker depth (m)")) +
  scale_x_discrete(labels=c("0.625" = "0.625", "1.25" = "1.25", "2.5" = "2.5", "5" = "5", "Inf" = "no flankers")) 
Exp2_plot_noRingOrUnsure

##---------------------------------------------------------------
##        SUPPLEMENTARY GRAPHS: Save figure S2 as .svg file       -
##---------------------------------------------------------------  
boxup("SUPPLEMENTARY GRAPHS: Save figure S2 as .svg files", bandChar = "-")
# The graphs are saved as a .svg file to enable formatting, labeling, and other 
# edits required to prepare the figures for the supplementary material.

#Save plots
JoinedPlots <- grid.arrange(Exp2_plot_TCenter, Exp2_plot_TOffCenter, Exp2_plot_TObstructed, 
                            Exp2_plot_TOutside, Exp2_plot_noRingOrUnsure, ncol = 2)
ggsave(file ="Smithers et al- Figure 8-figure supplement 2 (PerceivedTargetPos).svg", device = 'svg', plot = JoinedPlots, width = 23, height = 22)

##================================================================================================================================================
##  SUPPLEMENTARY GRAPHS: Experiment 3: Flanker ring and fixation at 1.26 m with target presented at each depth (Figure 9-figure supplement 1)   =
##================================================================================================================================================
boxup("SUPPLEMENTARY GRAPHS: Experiment 3: Flanker ring and fixation at 1.26 m with target presented at each depth (Figure 9-figure supplement 1)", bandChar = "=")

Exp3_data <-read.csv("Exp3_accumulatedPerErr.csv", header=T) # Load data
Exp3_data <- AddExtraCols(Exp3_data)

CountData <- CountDataFunction(Exp3_data)

repeatno. = 10  

PlotCountData <- subset(CountData, Perception == "Target in center of ring")
Exp3_plot_TCenter <- ggplot(PlotCountData, aes(x=as.factor(TargetFlankerSpacing), y= (Count/repeatno.), fill= TargetDepth)) + 
  geom_boxplot(lwd=0.75, fatten =2, outlier.shape=NA) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 1), shape =21,  size=1.5, alpha=0.5, colour = 'black') +
  xlab("Target-flanker spacing (deg)") + Graph.theme + scale_fill_manual(values=cbPalette) +
  facet_wrap(~FixationDepth,ncol = 3) +  ggtitle("Target in center of flanker ring") +
  scale_y_continuous(breaks= seq(0,1,0.2), limits=c(-0.02,1.02), expand=c(0,0),  name="Proportion of trials") + 
  guides(fill=guide_legend(title="Target depth (m)")) +
  scale_x_discrete(labels=c("0.625" = "0.625", "1.25" = "1.25", "2.5" = "2.5", "5" = "5", "Inf" = "no flankers")) 
Exp3_plot_TCenter

PlotCountData <- subset(CountData, Perception == "Target not in center of ring")
Exp3_plot_TOffCenter <- ggplot(PlotCountData, aes(x=as.factor(TargetFlankerSpacing), y= (Count/repeatno.), fill= TargetDepth)) + 
  geom_boxplot(lwd=0.75, fatten =2, outlier.shape=NA) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 1), shape =21,  size=1.5, alpha=0.5, colour = 'black') +
  xlab("Target-flanker spacing (deg)") + Graph.theme + scale_fill_manual(values=cbPalette) +
  facet_wrap(~FixationDepth,ncol = 3) + ggtitle("Target inside, but not in center, of flanker ring") +
  scale_y_continuous(breaks= seq(0,1,0.2), limits=c(-0.02,1.02), expand=c(0,0),  name="Proportion of trials") + 
  guides(fill=guide_legend(title="Flanker depth (m)")) + 
  scale_x_discrete(labels=c("0.625" = "0.625", "1.25" = "1.25", "2.5" = "2.5", "5" = "5", "Inf" = "no flankers")) 
Exp3_plot_TOffCenter

PlotCountData <- subset(CountData, Perception == "Ring obstructs target")
Exp3_plot_TObstructed <- ggplot(PlotCountData, aes(x=as.factor(TargetFlankerSpacing), y= (Count/repeatno.), fill= TargetDepth)) + 
  geom_boxplot(lwd=0.75, fatten =2, outlier.shape=NA) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 1), shape =21,  size=1.5, alpha=0.5, colour = 'black') +
  xlab("Target-flanker spacing (deg)") + Graph.theme + scale_fill_manual(values=cbPalette) +
  facet_wrap(~FixationDepth,ncol = 3) + ggtitle("Flanker ring obstructs target") +
  scale_y_continuous(breaks= seq(0,1,0.2), limits=c(-0.02,1.02), expand=c(0,0),  name="Proportion of trials") + 
  guides(fill=guide_legend(title="Flanker depth (m)")) + 
  scale_x_discrete(labels=c("0.625" = "0.625", "1.25" = "1.25", "2.5" = "2.5", "5" = "5", "Inf" = "no flankers")) 
Exp3_plot_TObstructed

PlotCountData <- subset(CountData, Perception == "Target outside ring")
Exp3_plot_TOutside <- ggplot(PlotCountData, aes(x=as.factor(TargetFlankerSpacing), y= (Count/repeatno.), fill= TargetDepth)) + 
  geom_boxplot(lwd=0.75, fatten =2, outlier.shape=NA) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 1), shape =21,  size=1.5, alpha=0.5, colour = 'black') +
  xlab("Target-flanker spacing (deg)") + Graph.theme + scale_fill_manual(values=cbPalette) +
  facet_wrap(~FixationDepth,ncol = 3) + ggtitle("Target outside flanker ring") +
  scale_y_continuous(breaks= seq(0,1,0.2), limits=c(-0.02,1.02), expand=c(0,0),  name="Proportion of trials") + 
  guides(fill=guide_legend(title="Flanker depth (m)")) + 
  scale_x_discrete(labels=c("0.625" = "0.625", "1.25" = "1.25", "2.5" = "2.5", "5" = "5", "Inf" = "no flankers")) 
Exp3_plot_TOutside

PlotCountData <- subset(CountData, Perception == "Unsure or no ring")
Exp3_plot_noRingOrUnsure <- ggplot(PlotCountData, aes(x=as.factor(TargetFlankerSpacing), y= (Count/repeatno.), fill= TargetDepth)) + 
  geom_boxplot(lwd=0.75, fatten =2, outlier.shape=NA) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 1), shape =21,  size=1.5, alpha=0.5, colour = 'black') +
  xlab("Target-flanker spacing (deg)") + Graph.theme + scale_fill_manual(values=cbPalette) +
  facet_wrap(~FixationDepth,ncol = 3) + ggtitle("No flanker ring or unsure") +
  scale_y_continuous(breaks= seq(0,1,0.2), limits=c(-0.02,1.02), expand=c(0,0),  name="Proportion of trials") + 
  guides(fill=guide_legend(title="Flanker depth (m)")) +
  scale_x_discrete(labels=c("0.625" = "0.625", "1.25" = "1.25", "2.5" = "2.5", "5" = "5", "Inf" = "no flankers")) 
Exp3_plot_noRingOrUnsure

##---------------------------------------------------------------
##        SUPPLEMENTARY GRAPHS: Save figure S3 as .svg file       -
##---------------------------------------------------------------  
boxup("SUPPLEMENTARY GRAPHS: Save figure S3 as .svg files", bandChar = "-")
# The graphs are saved as a .svg file to enable formatting, labeling, and other 
# edits required to prepare the figures for the supplementary material.

#Save plots
JoinedPlots <- grid.arrange(Exp3_plot_TCenter, Exp3_plot_TOffCenter, Exp3_plot_TObstructed, 
                            Exp3_plot_TOutside, Exp3_plot_noRingOrUnsure, ncol = 2)
ggsave(file ="Smithers et al- Figure 9-figure supplement 1 (PerceivedTargetPos).svg", device = 'svg', plot = JoinedPlots, width = 23, height = 22)


##======================================================================================================================================================
##  SUPPLEMENTARY GRAPHS: Experiment 4: Flanker ring and fixation at 0.4 m or 4 m with target presented at each depth (Figure 9-figure supplement 2)   =
##======================================================================================================================================================
boxup("SUPPLEMENTARY GRAPHS: Experiment 4: Flanker ring and fixation at 0.4 m or 4 m with target presented at each depth (Figure 9-figure supplement 2)", bandChar = "=")

Exp4_data <-read.csv("Exp4_accumulatedPerErr.csv", header=T) # Load data
Exp4_data <- AddExtraCols(Exp4_data)

CountData <- CountDataFunction(Exp4_data)

repeatno. = 8

PlotCountData <- subset(CountData, Perception == "Target in center of ring")
Exp4_plot_TCenter <- ggplot(PlotCountData, aes(x=as.factor(TargetFlankerSpacing), y= (Count/repeatno.), fill= TargetDepth)) + 
  geom_boxplot(lwd=0.75, fatten =2, outlier.shape=NA) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 1), shape =21,  size=1.5, alpha=0.5, colour = 'black') +
  xlab("Target-flanker spacing (deg)") + Graph.theme + scale_fill_manual(values=cbPalette) +
  facet_wrap(~FixationDepth,ncol = 3) + 
  ggtitle("Target in center of flanker ring") +
  scale_y_continuous(breaks= seq(0,1,0.2), limits=c(-0.02,1.02), expand=c(0,0),  name="Proportion of trials") + 
  guides(fill=guide_legend(title="Target depth (m)")) +
  scale_x_discrete(labels=c("0.625" = "0.625", "1.25" = "1.25", "2.5" = "2.5", "5" = "5", "Inf" = "no flankers")) 
Exp4_plot_TCenter

PlotCountData <- subset(CountData, Perception == "Target not in center of ring")
Exp4_plot_TOffCenter <- ggplot(PlotCountData, aes(x=as.factor(TargetFlankerSpacing), y= (Count/repeatno.), fill= TargetDepth)) + 
  geom_boxplot(lwd=0.75, fatten =2, outlier.shape=NA) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 1), shape =21,  size=1.5, alpha=0.5, colour = 'black') +
  xlab("Target-flanker spacing (deg)") + Graph.theme + scale_fill_manual(values=cbPalette) +
  facet_wrap(~FixationDepth,ncol = 3) + 
  ggtitle("Target inside, but not in center, of flanker ring") +
  scale_y_continuous(breaks= seq(0,1,0.2), limits=c(-0.02,1.02), expand=c(0,0),  name="Proportion of trials") + 
  guides(fill=guide_legend(title="Flanker depth (m)")) +
  scale_x_discrete(labels=c("0.625" = "0.625", "1.25" = "1.25", "2.5" = "2.5", "5" = "5", "Inf" = "no flankers")) 
Exp4_plot_TOffCenter

PlotCountData <- subset(CountData, Perception == "Ring obstructs target")
Exp4_plot_TObstructed <- ggplot(PlotCountData, aes(x=as.factor(TargetFlankerSpacing), y= (Count/repeatno.), fill= TargetDepth)) + 
  geom_boxplot(lwd=0.75, fatten =2, outlier.shape=NA) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 1), shape =21,  size=1.5, alpha=0.5, colour = 'black') +
  xlab("Target-flanker spacing (deg)") + Graph.theme + scale_fill_manual(values=cbPalette) +
  facet_wrap(~FixationDepth,ncol = 3) + 
  ggtitle("Flanker ring obstructs target") +
  scale_y_continuous(breaks= seq(0,1,0.2), limits=c(-0.02,1.02), expand=c(0,0),  name="Proportion of trials") + 
  guides(fill=guide_legend(title="Flanker depth (m)")) +
  scale_x_discrete(labels=c("0.625" = "0.625", "1.25" = "1.25", "2.5" = "2.5", "5" = "5", "Inf" = "no flankers")) 
Exp4_plot_TObstructed

PlotCountData <- subset(CountData, Perception == "Target outside ring")
Exp4_plot_TOutside <- ggplot(PlotCountData, aes(x=as.factor(TargetFlankerSpacing), y= (Count/repeatno.), fill= TargetDepth)) + 
  geom_boxplot(lwd=0.75, fatten =2, outlier.shape=NA) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 1), shape =21,  size=1.5, alpha=0.5, colour = 'black') +
  xlab("Target-flanker spacing (deg)") + Graph.theme + scale_fill_manual(values=cbPalette) +
  facet_wrap(~FixationDepth,ncol = 3) + 
  ggtitle("Target outside flanker ring") +
  scale_y_continuous(breaks= seq(0,1,0.2), limits=c(-0.02,1.02), expand=c(0,0),  name="Proportion of trials") + 
  guides(fill=guide_legend(title="Flanker depth (m)")) +
  scale_x_discrete(labels=c("0.625" = "0.625", "1.25" = "1.25", "2.5" = "2.5", "5" = "5", "Inf" = "no flankers")) 
Exp4_plot_TOutside

PlotCountData <- subset(CountData, Perception == "Unsure or no ring")
Exp4_plot_noRingOrUnsure <- ggplot(PlotCountData, aes(x=as.factor(TargetFlankerSpacing), y= (Count/repeatno.), fill= TargetDepth)) + 
  geom_boxplot(lwd=0.75, fatten =2, outlier.shape=NA) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 1), shape =21,  size=1.5, alpha=0.5, colour = 'black') +
  xlab("Target-flanker spacing (deg)") + Graph.theme + scale_fill_manual(values=cbPalette) +
  facet_wrap(~FixationDepth,ncol = 3) + 
  ggtitle("No flanker ring or unsure") +
  scale_y_continuous(breaks= seq(0,1,0.2), limits=c(-0.02,1.02), expand=c(0,0),  name="Proportion of trials") + 
  guides(fill=guide_legend(title="Flanker depth (m)")) +
  scale_x_discrete(labels=c("0.625" = "0.625", "1.25" = "1.25", "2.5" = "2.5", "5" = "5", "Inf" = "no flankers")) 
Exp4_plot_noRingOrUnsure

##---------------------------------------------------------------
##        SUPPLEMENTARY GRAPHS: Save figure S4 as .svg file       -
##---------------------------------------------------------------  
boxup("SUPPLEMENTARY GRAPHS: Save figure S4 as .svg files", bandChar = "-")
# The graphs are saved as a .svg file to enable formatting, labeling, and other 
# edits required to prepare the figures for the supplementary material.

#Save plots
JoinedPlots <- grid.arrange(Exp4_plot_TCenter, Exp4_plot_TOffCenter, Exp4_plot_TObstructed, 
                            Exp4_plot_TOutside, Exp4_plot_noRingOrUnsure, ncol = 2)
ggsave(file ="Smithers et al- Figure 9-figure supplement 2 (PerceivedTargetPos).svg", device = 'svg', plot = JoinedPlots, width = 23, height = 22)

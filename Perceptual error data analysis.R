# This is a copy of the R script used to generate the figures for, and perform the statistical analysis on, 
# the perceptual error data in Smithers et al. (2023)- "Large depth differences between target and flankers 
# can increase crowding: Evidence from a multi-depth plane display", eLife. For an explanation of the 
# statistical analysis conducted, please refer to the statistical analysis section of the methods in the 
# main manuscript.

# Script written by Dr Samuel P. Smithers, Northeastern University, 2022-2023
# Last edited July 2023

# Corresponding authors SPS (s.smithers@northeastern.edu) and PJB (p.bex@northeastern.edu)

# Included in this script:
# - Code used to generate Figures 5, 6 and 7 from the main manuscript.
# - Code used for the statistical analysis of perceptual error that is reported in the main manuscript.
# - Code used to generate Figure 5-figure supplement 1 and Figure 6-figure supplement 1.
# - Code used for the supporting statistical analysis of perceptual error that is reported in the legend 
# for Figure 5-figure supplement 1 and Figure 6-figure supplement 1.

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
#   [1] bannerCommenter_1.0.0 lme4_1.1-31           Matrix_1.5-1          cowplot_1.1.1         gridExtra_2.3         svglite_2.1.0        
# [7] ggplot2_3.4.0        
# 
# loaded via a namespace (and not attached):
#   [1] Rcpp_1.0.9        pillar_1.8.1      compiler_4.2.2    nloptr_2.0.3      tools_4.2.2       boot_1.3-28       lifecycle_1.0.3   tibble_3.2.1     
# [9] nlme_3.1-160      gtable_0.3.1      lattice_0.20-45   pkgconfig_2.0.3   rlang_1.1.0       cli_3.6.1         rstudioapi_0.14   withr_2.5.0      
# [17] dplyr_1.1.1       generics_0.1.3    vctrs_0.6.1       systemfonts_1.0.4 grid_4.2.2        tidyselect_1.2.0  glue_1.6.2        R6_2.5.1         
# [25] fansi_1.0.3       minqa_1.2.5       magrittr_2.0.3    scales_1.2.1      MASS_7.3-58.1     splines_4.2.2     colorspace_2.0-3  utf8_1.2.2       
# [33] munsell_0.5.0    

## Required dependencies/packages. 
#For plots
library(ggplot2) # For graphs
library(svglite) # Required to save .svg files (optional)
library(gridExtra) # For arranging/combining plots (optional)
library(cowplot) # for ggsave and arranging plots (optional)
#For stats
library(lme4) # Used for mixed effects models (also requires Matrix)
#Other
library(bannerCommenter) # For the section banners (optional)

rm(list=ls(all=TRUE))

#Set up Working directory 
setwd("***Pathway to folder containing this R script and all '_accumulatedPerErr.csv' and '_accumulatedPerErr_TargetInsideRing.csv' files***")

##############################################################################
##############################################################################
###                                                                        ###
###  PERCEPTUAL ERROR GRAPHS FOR THE MAIN MANUSCRIPT (FIGURES 5, 6 AND 7)  ###
###                                                                        ###
##############################################################################
##############################################################################
banner("Perceptual Error Graphs for the Main Manuscript (Figures 5, 6 and 7)", emph = TRUE)

#Colour blind friendly palette
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#D55E00", "#0072B2", "#CC79A7","#F0E442")

#Create graph theme
Graph.theme<-theme_bw()+
  theme(axis.text=element_text(size=25,colour="black"),
        axis.title=element_text(size=30, colour="black"),
        axis.title.y=element_text(vjust=1),axis.title.x=element_text(vjust=-3),
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
        strip.text=element_text(face="bold", size=25,colour = "black"))

#Make function to add extra column in which depth is given in meters.
AddExtraCols <- function(Inputdata){
  # Make column for flanker depth in m
  Inputdata["Flanker_depth_in_m"] <- Inputdata$Flanker_Screen
  Inputdata$Flanker_depth_in_m[Inputdata$Flanker_depth_in_m == "far"] <- "4"
  Inputdata$Flanker_depth_in_m[Inputdata$Flanker_depth_in_m == "mid"] <- "1.26"
  Inputdata$Flanker_depth_in_m[Inputdata$Flanker_depth_in_m == "near"] <- "0.4"
  # Make column for target depth in m
  Inputdata["Target_depth_in_m"] <- Inputdata$Target_Screen
  Inputdata$Target_depth_in_m[Inputdata$Target_depth_in_m == "far"] <- "4"
  Inputdata$Target_depth_in_m[Inputdata$Target_depth_in_m == "mid"] <- "1.26"
  Inputdata$Target_depth_in_m[Inputdata$Target_depth_in_m == "near"] <- "0.4"
  # Make column for fixation depth in m
  Inputdata["Fixation_depth_in_m"] <- Inputdata$Fixation_Screen
  Inputdata$Fixation_depth_in_m[Inputdata$Fixation_depth_in_m == "far"] <- "4"
  Inputdata$Fixation_depth_in_m[Inputdata$Fixation_depth_in_m == "mid"] <- "1.26"
  Inputdata$Fixation_depth_in_m[Inputdata$Fixation_depth_in_m == "near"] <- "0.4"
  return(Inputdata)
}

##=================================================================================================================
##  GRAPHS: Experiment 1: Target and fixation at 1.26 m with flanker ring presented at each depth  (Figure 5A)    =
##=================================================================================================================
boxup("GRAPHS: Experiment 1: Target and fixation at 1.26 m with flanker ring presented at each depth  (Figure 5A)", bandChar = "=")

Exp1_data <-read.csv("Exp1_accumulatedPerErr.csv", header=T) # Load data
Exp1_data <- AddExtraCols(Exp1_data)
N <- nrow(unique(Exp1_data["ID"])) 
print(N) # N=22

Exp_1_boxplot <- ggplot(Exp1_data, aes(x=factor(Target_Flanker_Spacing_in_deg), y= CircSD_deg, fill= Flanker_depth_in_m)) + 
  geom_boxplot(lwd=0.75, fatten =2, outlier.shape=NA) + xlab("Target-flanker spacing (deg)") + Graph.theme + scale_fill_manual(values=cbPalette) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15), shape = 21, size=2, alpha=0.6, colour = 'black') +
  facet_wrap(~Fixation_depth_in_m,ncol = 3) +
  scale_y_continuous(breaks= seq(0,180,20), limits=c(0,130), expand=c(0,0),  name="Perceptual error (deg)") + 
  guides(fill=guide_legend(title="Flanker depth (m)")) +
  scale_x_discrete(labels=c("0.625" = "0.625", "1.25" = "1.25", "2.5" = "2.5", "5" = "5", "Inf" = "no flankers"))
Exp_1_boxplot 

##======================================================================================================================
##  GRAPHS: Experiment 2: Target and fixation at 0.4 m or 4 m with flanker ring presented at each depth (Figure 5B)    =
##======================================================================================================================
boxup("GRAPHS: Experiment 2: Target and fixation at 0.4 m or 4 m with flanker ring presented at each depth (Figure 5B)", bandChar = "=")

Exp2_data <-read.csv("Exp2_accumulatedPerErr.csv", header=T)
Exp2_data <- AddExtraCols(Exp2_data)
N <- nrow(unique(Exp2_data["ID"]))
print(N) # N=19

Exp_2_boxplot <- ggplot(Exp2_data, aes(x=factor(Target_Flanker_Spacing_in_deg), y= CircSD_deg, fill= Flanker_depth_in_m)) + 
  geom_boxplot(lwd=0.75, fatten =2, outlier.shape=NA) + xlab("Target-flanker spacing (deg)") + Graph.theme  + scale_fill_manual(values=cbPalette) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15), shape = 21, size=2, alpha=0.6, colour = 'black') +
  facet_wrap(~Fixation_depth_in_m,ncol = 2) +
  scale_y_continuous(breaks= seq(0,180,20), limits=c(0,120), expand=c(0,0), name="Perceptual error (deg)") + 
  guides(fill=guide_legend(title="Flanker depth (m)")) +
  scale_x_discrete(labels=c("0.625" = "0.625", "1.25" = "1.25", "2.5" = "2.5", "5" = "5", "Inf" = "no flankers"))
Exp_2_boxplot


##===============================================================================================================
##  GRAPHS: Experiment 3: Flanker ring and fixation at 1.26 m with target presented at each depth (Figure 6A)   =
##===============================================================================================================
boxup("GRAPHS: Experiment 3: Flanker ring and fixation at 1.26 m with target presented at each depth (Figure 6A)", bandChar = "=")

Exp3_data <-read.csv("Exp3_accumulatedPerErr.csv", header=T)
Exp3_data <- AddExtraCols(Exp3_data)
N <- nrow(unique(Exp3_data["ID"]))
print(N) # N=21

Exp_3_boxplot <- ggplot(Exp3_data, aes(x=factor(Target_Flanker_Spacing_in_deg), y= CircSD_deg, fill= Target_depth_in_m)) + 
  geom_boxplot(lwd=0.75, fatten =2, outlier.shape=NA) + xlab("Target-flanker spacing (deg)") + Graph.theme + scale_fill_manual(values=cbPalette) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15), shape = 21, size=2, alpha=0.6, colour = 'black') +
  facet_wrap(~Fixation_depth_in_m,ncol = 2, scales = "free_x") + force_panelsizes(cols = c(1, 1/4)) +
  scale_y_continuous(breaks= seq(0,180,20), limits=c(0,150), expand=c(0,0),  name="Perceptual error (deg)") + 
  guides(fill=guide_legend(title="Target depth (m)"))+
  scale_x_discrete(labels=c("0.625" = "0.625", "1.25" = "1.25", "2.5" = "2.5", "5" = "5", "Inf" = "no flankers"))
Exp_3_boxplot

##=====================================================================================================================
##  GRAPHS: Experiment 4: Flanker ring and fixation at 0.4 m or 4 m with target presented at each depth (Figure 6B)   =
##=====================================================================================================================
boxup("GRAPHS: Experiment 4: Flanker ring and fixation at 0.4 m or 4 m with target presented at each depth (Figure 6B)", bandChar = "=")

Exp4_data <-read.csv("Exp4_accumulatedPerErr.csv", header=T)
Exp4_data <- AddExtraCols(Exp4_data)
N <- nrow(unique(Exp4_data["ID"]))
print(N) # N=21

Exp_4_boxplot <- ggplot(Exp4_data, aes(x=factor(Target_Flanker_Spacing_in_deg), y= CircSD_deg, fill= Target_depth_in_m)) + 
  geom_boxplot(lwd=0.75, fatten =2, outlier.shape=NA) + xlab("Target-flanker spacing (deg)") + Graph.theme  + scale_fill_manual(values=cbPalette) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15), shape = 21, size=2, alpha=0.6, colour = 'black') +
  facet_wrap(~Fixation_depth_in_m,ncol = 3) + force_panelsizes(cols = c(1, 1, 1/3)) +
  scale_y_continuous(breaks= seq(0,180,20),limits=c(0,140), expand=c(0,0), name="Perceptual error (deg)") + 
  guides(fill=guide_legend(title="Target depth (m)"))+
  scale_x_discrete(labels=c("0.625" = "0.625", "1.25" = "1.25", "2.5" = "2.5", "5" = "5", "Inf" = "no flankers"))
Exp_4_boxplot 

##=======================================================================================
##  GRAPHS: Experiment 5: Target and flanker ring always at fixation depth (Figure 7)   =
##=======================================================================================
boxup("GRAPHS: Experiment 5: Target and flanker ring always at fixation depth (Figure 7)", bandChar = "=")

Exp5_data <-read.csv("Exp5_accumulatedPerErr.csv", header=T)
Exp5_data <- AddExtraCols(Exp5_data)
N <- nrow(unique(Exp5_data["ID"]))
print(N) # N=15

Exp_5_boxplot <- ggplot(Exp5_data, aes(x=as.factor(Target_Flanker_Spacing_in_deg), y= CircSD_deg, fill= Target_depth_in_m)) + 
  geom_boxplot(lwd=0.75, fatten =2, outlier.shape=NA) + xlab("Target-flanker spacing (deg)") + Graph.theme + scale_fill_manual(values=cbPalette) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15), shape = 21, size=2, alpha=0.6, colour = 'black') +
  scale_y_continuous(breaks= seq(0,180,20), limits=c(0,120), expand=c(0,0), name="Perceptual error (deg)") + 
  guides(fill=guide_legend(title="Depth (m)"))+ 
  scale_x_discrete(labels=c("0.625" = "0.625", "1.25" = "1.25", "2.5" = "2.5", "5" = "5", "Inf" = "no flankers"))
Exp_5_boxplot

##---------------------------------------------------------------
##              GRAPHS: Save graphs as .svg files               -
##---------------------------------------------------------------
boxup("GRAPHS: Save graphs as .svg files", bandChar = "-")
# The graphs are then saved as a .svg file to enable formatting, labeling, and other 
# edits required to prepare the figures for the paper.
CombinedPlots <- plot_grid(Exp_1_boxplot, Exp_2_boxplot, Exp_3_boxplot, Exp_4_boxplot, Exp_5_boxplot, ncol = 2, rel_widths = c(5/8, 1))
#Combined plots will be saved to your working directory 
ggsave(file ="Smithers et al- figures 5-7 (PerErr).svg", device = 'svg', plot = CombinedPlots, width = 26, height = 22)

############################################################################
############################################################################
###                                                                      ###
###   STATISTICAL ANALYSIS OF PERCEPTUAL ERROR FOR THE MAIN MANUSCRIPT   ###
###                                                                      ###
############################################################################
############################################################################
banner("Statistical Analysis of Perceptual Error for the Main Manuscript", emph = TRUE)
rm(list=ls(all=TRUE))
# For the stats we do not include the control in the mixed effects models. An explanation for this decision 
# is provided in the statistical analysis section of the methods in the main manuscript.

##===================================================================================================
##  STATS: Experiment 1: Target and fixation at 1.26 m with flanker ring presented at each depth.   =
##===================================================================================================
boxup("STATS: Experiment 1: Target and fixation at 1.26 m with flanker ring presented at each depth.", bandChar = "=")
Exp1_data <-read.csv("Exp1_accumulatedPerErr.csv", header=T) # Load data
Exp1_data<- subset(Exp1_data, Target_Flanker_Spacing_in_deg!=Inf) # Remove control (i.e. no flanker condition)
#The explanation for why the control is removed is provided in the statistical analysis section of the methods.

#Fit full mixed model containing fixed effects and all possible pairwise interactions of the fixed effects 
Exp1_M0<-lmer(log(CircSD_deg)~ Target_Flanker_Spacing_in_deg+Flanker_Screen + 
                Target_Flanker_Spacing_in_deg:Flanker_Screen + (1|ID), data= Exp1_data, REML=TRUE) 
#The response variable is natural log transformed to correct for positive skew

#Check normality and homoscedasticity
sresid <- resid(Exp1_M0)  # Extract the standardized residuals
hist(sresid) # Plot histogram and QQ plot to look at distribution of residuals
qqnorm(sresid)
qqline(sresid)
fitted.glmm <- fitted(Exp1_M0, level=1) # Extract the fitted (predicted) values
plot(sresid ~ fitted.glmm) # Check for homoscedasticity 
plot(Exp1_M0)
#These look fine and acceptable. 

# Test significance of the interaction using a likelihood ratio rest (LRT) to compare the full model 
# with the same model but with the interaction removed. 
m1 <-update(Exp1_M0,~.-Target_Flanker_Spacing_in_deg:Flanker_Screen)
anova(Exp1_M0, m1, test = "LRT") # Chisq=5.5642, Df=2, p= 0.06191 # Reported in paper
# Drop 2-way interaction from the full model as it is not significant 

# We then use the minimum model (i.e. model containing all fixed effects of interest plus any 
# significant two-way interactions) to report the significance of each main effect or important 
# interaction(s) (where applicable) by using a LRT to compare the minimum model with same model but with the 
# effect/interaction of interest removed. 
#***************************************************************
#Stats reported in paper. 
m2a <-update(m1,~.-Target_Flanker_Spacing_in_deg)
anova(m1, m2a, test = "LRT") # Chisq=31.713, Df=1, p = 1.787e-08 ***
m2b <-update(m1,~.-Flanker_Screen)
anova(m1, m2b, test = "LRT") # Chisq=20.673, Df=2, p = 3.243e-05 ***
#***************************************************************
rm(list=ls(all=TRUE))

##=========================================================================================================
##  STATS: Experiment 2: Target and fixation at 0.4 m or 4 m with flanker ring presented at each depth.   =
##=========================================================================================================
boxup("STATS: Experiment 2: Target and fixation at 0.4 m or 4 m with flanker ring presented at each depth.", bandChar = "=")

Exp2_data <-read.csv("Exp2_accumulatedPerErr.csv", header=T)
Exp2_data<- subset(Exp2_data, Target_Flanker_Spacing_in_deg!=Inf) # Remove control (i.e. no flanker condition)

#Fit full mixed model containing fixed effects and all possible pairwise interactions of the fixed effects
Exp2_M0<-lmer(log(CircSD_deg)~ Target_Flanker_Spacing_in_deg+Target_Screen+Flanker_Screen +
                 Target_Flanker_Spacing_in_deg:Target_Screen + 
                 Target_Flanker_Spacing_in_deg:Flanker_Screen +
                 Target_Screen:Flanker_Screen + (1|ID), data= Exp2_data, REML=TRUE) 
#The response variable is natural log transformed to correct for positive skew

#Check normality and homoscedasticity
sresid <- resid(Exp2_M0) # Extract the standardized residuals
hist(sresid) # Plot histogram and QQ plot to look at distribution of residuals
qqnorm(sresid)
qqline(sresid)
fitted.glmm <- fitted(Exp2_M0, level=1) # Extract the fitted (predicted) values
plot(sresid ~ fitted.glmm) # Check for homoscedasticity 
plot(Exp2_M0)
#These look fine and acceptable. 

m1a <-update(Exp2_M0,~.-Target_Screen:Flanker_Screen)
anova(Exp2_M0, m1a, test = "LRT") # Chisq=17.357, Df=2, p = 0.0001702 ***
m1b <-update(Exp2_M0,~.-Target_Flanker_Spacing_in_deg:Flanker_Screen)
anova(Exp2_M0, m1b, test = "LRT") # Chisq= 4.5239, Df=2, p = 0.1041 # Reported in paper
m1C <-update(Exp2_M0,~.-Target_Flanker_Spacing_in_deg:Target_Screen)
anova(Exp2_M0, m1C, test = "LRT") # Chisq=9.2907, Df=1, p = 0.002303 **
#Drop the Target_Flanker_Spacing_in_deg:Flanker_Screen interaction

#***************************************************************
#m1b is the min model. These are the stats report in the paper. 
m2a <-update(m1b,~.-Target_Screen:Flanker_Screen)
anova(m1b, m2a, test = "LRT") # Chisq=17.181, Df=2, p = 0.0001858 ***
m2b <-update(m1b,~.-Target_Flanker_Spacing_in_deg:Target_Screen)
anova(m1b, m2b, test = "LRT") # Chisq=9.196, Df=1, p = 0.002425 **
#*****************************************************************
rm(list=ls(all=TRUE))

##==================================================================================================
##  STATS: Experiment 3: Flanker ring and fixation at 1.26 m with target presented at each depth   =
##==================================================================================================
boxup("STATS: Experiment 3: Flanker ring and fixation at 1.26 m with target presented at each depth", bandChar = "=")

Exp3_data <-read.csv("Exp3_accumulatedPerErr.csv", header=T)
Exp3_data<- subset(Exp3_data, Target_Flanker_Spacing_in_deg!=Inf) # Remove control (i.e. no flanker condition)

#Fit full mixed model containing fixed effects and all possible pairwise interactions of the fixed effects 
Exp3_M0<-lmer(log(CircSD_deg)~ Target_Flanker_Spacing_in_deg + Target_Screen +
                Target_Flanker_Spacing_in_deg:Target_Screen + (1|ID), data= Exp3_data, REML=TRUE) 
#The response variable is natural log transformed to correct for positive skew

#Check normality and homoscedasticity
sresid <- resid(Exp3_M0) # Extract the standardized residuals
hist(sresid) # Plot histogram and QQ plot to look at distribution of residuals
qqnorm(sresid)
qqline(sresid)
fitted.glmm <- fitted(Exp3_M0, level=1) # Extract the fitted (predicted) values
plot(sresid ~ fitted.glmm) # Check for homoscedasticity 
plot(Exp3_M0)
#These look fine and acceptable. 

m1 <-update(Exp3_M0,~.-Target_Flanker_Spacing_in_deg:Target_Screen)
anova(Exp3_M0, m1, test = "LRT") #  Chisq=3.5483 , Df=2, p = 0.1696. # Reported in paper
#Drop the Target_Flanker_Spacing_in_deg:Target_Screen interaction

#*#***************************************************************
#m1 is the min model. These are the stats reported in the paper. 
m2a <-update(m1,~.-Target_Flanker_Spacing_in_deg)
anova(m1, m2a, test = "LRT") # Chisq=45.988, Df=1, p = 1.19e-11 ***
m2a <-update(m1,~.-Target_Screen)
anova(m1, m2a, test = "LRT") # Chisq=36.931, Df=2, p = 9.563e-09 ***
#***************************************************************
rm(list=ls(all=TRUE))

##========================================================================================================
##  STATS: Experiment 4: Flanker ring and fixation at 0.4 m or 4 m with target presented at each depth   =
##========================================================================================================
boxup("STATS: Experiment 4: Flanker ring and fixation at 0.4 m or 4 m with target presented at each depth", bandChar = "=")

Exp4_data <-read.csv("Exp4_accumulatedPerErr.csv", header=T)
Exp4_data <- subset(Exp4_data, Target_Flanker_Spacing_in_deg!=Inf)

#Fit full mixed model containing fixed effects and all possible pairwise interactions of the fixed effects
Exp4_M0<-lmer(log(CircSD_deg)~ Target_Flanker_Spacing_in_deg + Target_Screen + Flanker_Screen + 
                Target_Flanker_Spacing_in_deg:Target_Screen + 
                Target_Flanker_Spacing_in_deg:Flanker_Screen +
                Target_Screen:Flanker_Screen + (1|ID), data= Exp4_data, REML=TRUE) 
# The response variable is natural log transformed to correct for positive skew.

#Check normality and homoscedasticity
sresid <- resid(Exp4_M0) # Extract the standardized residuals
hist(sresid) # Plot histogram and QQ plot to look at distribution of residuals
qqnorm(sresid)
qqline(sresid)

fitted.glmm <- fitted(Exp4_M0, level=1) # Extract the fitted (predicted) values
plot(sresid ~ fitted.glmm) # Check for homoscedasticity 
plot(Exp4_M0)
#These look fine and acceptable. 

m1a <-update(Exp4_M0,~.-Target_Screen:Flanker_Screen)
anova(Exp4_M0, m1a, test = "LRT") # Chisq=62.063  Df=2, p = 3.336e-14 ***
m1b <-update(Exp4_M0,~.-Target_Flanker_Spacing_in_deg:Flanker_Screen)
anova(Exp4_M0, m1b, test = "LRT") # Chisq= 0.1387, Df=1, p = 0.7096 # Reported in paper
m1C <-update(Exp4_M0,~.-Target_Flanker_Spacing_in_deg:Target_Screen)
anova(Exp4_M0, m1C, test = "LRT") # Chisq=1.4339, Df=2, p = 0.4882
#Drop the Target_Flanker_Spacing_in_deg:Flanker_Screen interaction

m2a <-update(m1b,~.-Target_Screen:Flanker_Screen)
anova(m1b, m2a, test = "LRT") # Chisq=62.046 , Df=2, p = 3.364e-14 ***
m2b <-update(m1b,~.-Target_Flanker_Spacing_in_deg:Target_Screen)
anova(m1b, m2b, test = "LRT") # Chisq=1.4335, Df=2, p = 0.4883 # Reported in paper
#Drop Target_Flanker_Spacing_in_deg:Target_Screen

#***************************************************************
#m2b is the min model. These are the stats reported in the paper.
m3a <-update(m2b,~.-Target_Screen:Flanker_Screen)
anova(m2b, m3a, test = "LRT") # Chisq=61.874, Df=2, p = 3.667e-14 ***
m3b <-update(m2b,~.-Target_Flanker_Spacing_in_deg)
anova(m2b, m3b, test = "LRT") # Chisq=29.917, Df=1, p = 4.51e-08 ****
#*****************************************************************
rm(list=ls(all=TRUE))

##===========================================================================
##  STATS: Experiment 5: Target and flanker ring always at fixation depth   =
##===========================================================================
boxup("STATS: Experiment 5: Target and flanker ring always at fixation depth", bandChar = "=")

Exp5_data <-read.csv("Exp5_accumulatedPerErr.csv", header=T)
Exp5_data <- subset(Exp5_data, Target_Flanker_Spacing_in_deg!=Inf) # Remove control

#Fit full mixed model containing fixed effects and all possible pairwise interactions of the fixed effects
Exp5_M0<-lmer(log(CircSD_deg)~ Target_Flanker_Spacing_in_deg + Target_Screen + 
                Target_Flanker_Spacing_in_deg:Target_Screen + (1|ID), data= Exp5_data, REML=TRUE) 
# The response variable is natural log transformed to correct for positive skew

#Check normality and homoscedasticity
sresid <- resid(Exp5_M0) # Extract the standardized residuals
hist(sresid) # Plot histogram and QQ plot to look at distribution of residuals
qqnorm(sresid)
qqline(sresid)
fitted.glmm <- fitted(Exp5_M0, level=1) # Extract the fitted (predicted) values
plot(sresid ~ fitted.glmm) # Check for homoscedasticity 
plot(Exp5_M0)
#These look fine and acceptable. 

m1 <-update(Exp5_M0,~.-Target_Flanker_Spacing_in_deg:Target_Screen)
anova(Exp5_M0, m1, test = "LRT") # Chisq=1.7161, Df=2, p = 0.424 # Reported in paper
#Drop Target_Flanker_Spacing_in_deg:Target_Screen interaction

#***************************************************************
#m1 is the min model. These are the stats reported in the paper
m2a <-update(m1,~.-Target_Flanker_Spacing_in_deg)
anova(m2a, m1, test = "LRT") # Chisq=64.802, Df=1, p = 8.28e-16 ***
m2b <-update(m1,~.-Target_Screen)
anova(m2b, m1, test = "LRT") # Chisq=8.7476, Df=2, p = 0.0126 *
#***************************************************************
rm(list=ls(all=TRUE))

###################################################################################################################
###################################################################################################################
###                                                                                                             ###
###  PERCEPTUAL ERROR GRAPHS FOR SUPPLEMENTARY (FIGURE 5-FIGURE SUPPLEMENT 1 AND FIGURE 6-FIGURE SUPPLEMENT 1)  ###
###                                                                                                             ###
###################################################################################################################
###################################################################################################################
banner("Perceptual Error Graphs for supplementary (Figure 5-figure supplement 1 and Figure 6-figure supplement 1)", emph = TRUE)
#Perceptual error calculated using only trials in which the observer reported seeing the target inside the flanker ring

# Colour blind friendly palette
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#D55E00", "#0072B2", "#CC79A7","#F0E442")

#Create graph theme
Graph.theme<-theme_bw()+
  theme(axis.text=element_text(size=25,colour="black"),
        axis.title=element_text(size=30, colour="black"),
        axis.title.y=element_text(vjust=1),axis.title.x=element_text(vjust=-3),
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
        strip.text=element_text(face="bold", size=25,colour = "black"))

#Make function to add extra column in which depth is given in meters
AddExtraCols <- function(Inputdata){
  # Make column for flanker depth in m
  Inputdata["Flanker_depth_in_m"] <- Inputdata$Flanker_Screen
  Inputdata$Flanker_depth_in_m[Inputdata$Flanker_depth_in_m == "far"] <- "4"
  Inputdata$Flanker_depth_in_m[Inputdata$Flanker_depth_in_m == "mid"] <- "1.26"
  Inputdata$Flanker_depth_in_m[Inputdata$Flanker_depth_in_m == "near"] <- "0.4"
  # Make column for target depth in m
  Inputdata["Target_depth_in_m"] <- Inputdata$Target_Screen
  Inputdata$Target_depth_in_m[Inputdata$Target_depth_in_m == "far"] <- "4"
  Inputdata$Target_depth_in_m[Inputdata$Target_depth_in_m == "mid"] <- "1.26"
  Inputdata$Target_depth_in_m[Inputdata$Target_depth_in_m == "near"] <- "0.4"
  # Make column for fixation depth in m
  Inputdata["Fixation_depth_in_m"] <- Inputdata$Fixation_Screen
  Inputdata$Fixation_depth_in_m[Inputdata$Fixation_depth_in_m == "far"] <- "4"
  Inputdata$Fixation_depth_in_m[Inputdata$Fixation_depth_in_m == "mid"] <- "1.26"
  Inputdata$Fixation_depth_in_m[Inputdata$Fixation_depth_in_m == "near"] <- "0.4"
  return(Inputdata)
}


##===================================================================================================================================================
##  SUPPLEMENTARY GRAPHS: Experiment 1: Target and fixation at 1.26 m with flanker ring presented at each depth  (Figure 5-figure supplement 1A)    =
##===================================================================================================================================================
boxup("SUPPLEMENTARY GRAPHS: Experiment 1: Target and fixation at 1.26 m with flanker ring presented at each depth  (Figure 5-figure supplement 1A)", bandChar = "=")

Exp1_data <-read.csv("Exp1_accumulatedPerErr_TargetInsideRing.csv", header=T) # Load data
Exp1_data <- AddExtraCols(Exp1_data)
Exp1_data <- na.omit(Exp1_data)

Exp_1_boxplot <- ggplot(Exp1_data, aes(x=factor(Target_Flanker_Spacing_in_deg), y= CircSD_deg, fill= Flanker_depth_in_m)) + 
  geom_boxplot(lwd=0.75, fatten =2, outlier.shape=NA) + xlab("Target-flanker spacing (deg)") + Graph.theme + scale_fill_manual(values=cbPalette) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15), shape = 21, size=2, alpha=0.6, colour = 'black') +
  facet_wrap(~Fixation_depth_in_m,ncol = 3) +
  scale_y_continuous(breaks= seq(0,180,20), limits=c(0,130), expand=c(0,0),  name="Perceptual error (deg)") + 
  guides(fill=guide_legend(title="Flanker depth (m)")) +
  scale_x_discrete(labels=c("0.625" = "0.625", "1.25" = "1.25", "2.5" = "2.5", "5" = "5", "Inf" = "no flankers"))
Exp_1_boxplot 


##========================================================================================================================================================
##  SUPPLEMENTARY GRAPHS: Experiment 2: Target and fixation at 0.4 m or 4 m with flanker ring presented at each depth (Figure 5-figure supplement 1B)    =
##========================================================================================================================================================
boxup("SUPPLEMENTARY GRAPHS: Experiment 2: Target and fixation at 0.4 m or 4 m with flanker ring presented at each depth (Figure 5-figure supplement 1B)", bandChar = "=")

Exp2_data <-read.csv("Exp2_accumulatedPerErr_TargetInsideRing.csv", header=T)
Exp2_data <- AddExtraCols(Exp2_data)
Exp2_data <- na.omit(Exp2_data)

Exp_2_boxplot <- ggplot(Exp2_data, aes(x=factor(Target_Flanker_Spacing_in_deg), y= CircSD_deg, fill= Flanker_depth_in_m)) + 
  geom_boxplot(lwd=0.75, fatten =2, outlier.shape=NA) + xlab("Target-flanker spacing (deg)") + Graph.theme  + scale_fill_manual(values=cbPalette) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15), shape = 21, size=2, alpha=0.6, colour = 'black') +
  facet_wrap(~Fixation_depth_in_m,ncol = 2) +
  scale_y_continuous(breaks= seq(0,180,20), limits=c(0,130), expand=c(0,0), name="Perceptual error (deg)") + 
  guides(fill=guide_legend(title="Flanker depth (m)")) +
  scale_x_discrete(labels=c("0.625" = "0.625", "1.25" = "1.25", "2.5" = "2.5", "5" = "5", "Inf" = "no flankers"))
Exp_2_boxplot


##=================================================================================================================================================
##  SUPPLEMENTARY GRAPHS: Experiment 3: Flanker ring and fixation at 1.26 m with target presented at each depth (Figure 6-figure supplement 1A)   =
##=================================================================================================================================================
boxup("SUPPLEMENTARY GRAPHS: Experiment 3: Flanker ring and fixation at 1.26 m with target presented at each depth (Figure 6-figure supplement 1A)", bandChar = "=")

Exp3_data <-read.csv("Exp3_accumulatedPerErr_TargetInsideRing.csv", header=T)
Exp3_data <- AddExtraCols(Exp3_data)
Exp3_data <- na.omit(Exp3_data)

Exp_3_boxplot <- ggplot(Exp3_data, aes(x=factor(Target_Flanker_Spacing_in_deg), y= CircSD_deg, fill= Target_depth_in_m)) + 
  geom_boxplot(lwd=0.75, fatten =2, outlier.shape=NA) + xlab("Target-flanker spacing (deg)") + Graph.theme + scale_fill_manual(values=cbPalette) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15), shape = 21, size=2, alpha=0.6, colour = 'black') +
  facet_wrap(~Fixation_depth_in_m,ncol = 2) + force_panelsizes(cols = c(1, 1/4)) +
  scale_y_continuous(breaks= seq(0,180,20), limits=c(0,150), expand=c(0,0),  name="Perceptual error (deg)") + 
  guides(fill=guide_legend(title="Target depth (m)"))+
  scale_x_discrete(labels=c("0.625" = "0.625", "1.25" = "1.25", "2.5" = "2.5", "5" = "5", "Inf" = "no flankers"))
Exp_3_boxplot

##=======================================================================================================================================================
##  SUPPLEMENTARY GRAPHS: Experiment 4: Flanker ring and fixation at 0.4 m or 4 m with target presented at each depth (Figure 6-figure supplement 1B)   =
##=======================================================================================================================================================
boxup("SUPPLEMENTARY GRAPHS: Experiment 4: Flanker ring and fixation at 0.4 m or 4 m with target presented at each depth (Figure 6-figure supplement 1B)", bandChar = "=")

Exp4_data <-read.csv("Exp4_accumulatedPerErr_TargetInsideRing.csv", header=T)
Exp4_data <- AddExtraCols(Exp4_data)
Exp4_data <- na.omit(Exp4_data)

Exp_4_boxplot <- ggplot(Exp4_data, aes(x=factor(Target_Flanker_Spacing_in_deg), y= CircSD_deg, fill= Target_depth_in_m)) + 
  geom_boxplot(lwd=0.75, fatten =2, outlier.shape=NA) + xlab("Target-flanker spacing (deg)") + Graph.theme  + scale_fill_manual(values=cbPalette) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15), shape = 21, size=2, alpha=0.6, colour = 'black') +
  facet_wrap(~Fixation_depth_in_m,ncol = 3) + force_panelsizes(cols = c(1, 1, 1/3)) +
  scale_y_continuous(breaks= seq(0,180,20),limits=c(0,140), expand=c(0,0), name="Perceptual error (deg)") + 
  guides(fill=guide_legend(title="Target depth (m)"))+
  scale_x_discrete(labels=c("0.625" = "0.625", "1.25" = "1.25", "2.5" = "2.5", "5" = "5", "Inf" = "no flankers"))
Exp_4_boxplot 

##---------------------------------------------------------------
##        SUPPLEMENTARY GRAPHS: Save graphs as .svg files       -
##---------------------------------------------------------------  
boxup("SUPPLEMENTARY GRAPHS: Save graphs as .svg files", bandChar = "-")
# The graphs are then saved as a .svg file to enable formatting, labeling, and other 
# edits required to prepare the figures for the supplementary material.
CombinedPlots <- plot_grid(Exp_1_boxplot, Exp_2_boxplot, Exp_3_boxplot, Exp_4_boxplot, ncol = 2, rel_widths = c(5/8, 1))
#Combined plots will be saved to your working directory 
ggsave(file ="Smithers et al- Figure 5 and 6-figure supplement 1 (PerErr-TargetInsideRingOnly).svg", device = 'svg', plot = CombinedPlots, width = 26, height = 17)


############################################################################################
############################################################################################
###                                                                                      ###
###  SUPPORTING STATISTICAL ANALYSIS OF PERCEPTUAL ERROR FOR THE SUPPLEMENTARY MATERIAL  ###
###                                                                                      ###
############################################################################################
############################################################################################
banner("Supporting Statistical Analysis of Perceptual Error for the Supplementary Material", emph = TRUE)
#Perceptual error calculated using only trials in which the observer reported seeing the target inside the flanker ring

rm(list=ls(all=TRUE))

##================================================================================================================
##  SUPPLEMENTARY STATS: Experiment 1: Target and fixation at 1.26 m with flanker ring presented at each depth   =
##================================================================================================================
boxup("SUPPLEMENTARY STATS: Experiment 1: Target and fixation at 1.26 m with flanker ring presented at each depth", bandChar = "=")

Exp1_data <-read.csv("Exp1_accumulatedPerErr_TargetInsideRing.csv", header=T) # Load data
Exp1_data <- subset(Exp1_data, Target_Flanker_Spacing_in_deg!=Inf) # Remove control

#Fit full mixed model containing fixed effects and all possible pairwise interactions of the fixed effects 
Exp1_M0<-lmer(log(CircSD_deg)~ Target_Flanker_Spacing_in_deg+Flanker_Screen + 
                Target_Flanker_Spacing_in_deg:Flanker_Screen + (1|ID), data= Exp1_data, REML=TRUE) 
# The response variable is natural log transformed to correct for positive skew

#Check normality and homoscedasticity
sresid <- resid(Exp1_M0) # Extract the standardized residuals
hist(sresid) # Plot histogram and QQ plot to look at residuals
qqnorm(sresid)
qqline(sresid)
fitted.glmm <- fitted(Exp1_M0, level=1) # Extract the fitted (predicted) values
plot(sresid ~ fitted.glmm) # Check for homoscedasticity 
plot(Exp1_M0)
#These look acceptable. 

m1 <-update(Exp1_M0,~.-Target_Flanker_Spacing_in_deg:Flanker_Screen)
anova(Exp1_M0, m1, test = "LRT") # Chisq=1.9162, Df=2, p = 0.3836 # Reported in supplementary material. 
# 2-way interaction is not significant so it is dropped

#***************************************************************
#m1 is the min model. These are the stats reported in the supplementary material
m2a <-update(m1,~.-Target_Flanker_Spacing_in_deg)
anova(m1, m2a, test = "LRT") # Chisq=20.026, Df=1, p = 7.641e-06 ***
m2b <-update(m1,~.-Flanker_Screen)
anova(m1, m2b, test = "LRT") # Chisq=15.686, Df=2, p = 0.0003925 ***
#***************************************************************
rm(list=ls(all=TRUE))

##=======================================================================================================================
##  SUPPLEMENTARY STATS: Experiment 2: Target and fixation at 0.4 m or 4 m with flanker ring presented at each depth.   =
##=======================================================================================================================
boxup("SUPPLEMENTARY STATS: Experiment 2: Target and fixation at 0.4 m or 4 m with flanker ring presented at each depth.", bandChar = "=")

Exp2_data <-read.csv("Exp2_accumulatedPerErr_TargetInsideRing.csv", header=T)
Exp2_data<- subset(Exp2_data, Target_Flanker_Spacing_in_deg!=Inf) # Remove control (i.e. no flanker condition)

#Fit full mixed model containing fixed effects and all possible pairwise interactions of the fixed effects 
Exp2_M0<-lmer(log(CircSD_deg)~ Target_Flanker_Spacing_in_deg+Target_Screen+Flanker_Screen +
                Target_Flanker_Spacing_in_deg:Target_Screen + 
                Target_Flanker_Spacing_in_deg:Flanker_Screen +
                Target_Screen:Flanker_Screen + (1|ID), data= Exp2_data, REML=TRUE) 
# The response variable is natural log transformed to correct for positive skew

#Check normality and homoscedasticity
sresid <- resid(Exp2_M0) # Extract the standardized residuals
hist(sresid) # Plot histogram and QQ plot to look at residuals
qqnorm(sresid)
qqline(sresid)
fitted.glmm <- fitted(Exp2_M0, level=1) # Extract the fitted (predicted) values
plot(sresid ~ fitted.glmm) # Check for homoscedasticity 
plot(Exp2_M0)
#These look fine and acceptable. 

m1a <-update(Exp2_M0,~.-Target_Screen:Flanker_Screen)
anova(Exp2_M0, m1a, test = "LRT") # Chisq=3.4421, Df=2, p = 0.1789 # Reported in supplementary material
m1b <-update(Exp2_M0,~.-Target_Flanker_Spacing_in_deg:Flanker_Screen)
anova(Exp2_M0, m1b, test = "LRT") # Chisq= 5.1059, Df=2, p = 0.07785
m1C <-update(Exp2_M0,~.-Target_Flanker_Spacing_in_deg:Target_Screen)
anova(Exp2_M0, m1C, test = "LRT") # Chisq=7.0476, Df=1, p = 0.007937 **
#Drop the Target_Screen:Flanker_Screen interaction

m2a <-update(m1a,~.-Target_Flanker_Spacing_in_deg:Flanker_Screen)
anova(m1a, m2a, test = "LRT") # Chisq= 4.8548, Df=2, p = 0.08827 # Reported in supplementary material
m2b <-update(m1a,~.-Target_Flanker_Spacing_in_deg:Target_Screen)
anova(m1a, m2b, test = "LRT") # Chisq=7.1249, Df=1, p = 0.007602 **
#Drop the Target_Flanker_Spacing_in_deg:Flanker_Screen interaction

#***************************************************************
#m2b is the min model. These are the stats reported in the supplementary material
m3a <-update(m2a,~.-Flanker_Screen)
anova(m2a, m3a, test = "LRT") # Chisq=5.7007, Df=2, p = 0.05782
m3b <-update(m2a,~.-Target_Flanker_Spacing_in_deg:Target_Screen)
anova(m2a, m3b, test = "LRT") # Chisq=7.32, Df=1, p = 0.006819 **
#*****************************************************************
rm(list=ls(all=TRUE))

##================================================================================================================
##  SUPPLEMENTARY STATS: Experiment 3: Flanker ring and fixation at 1.26 m with target presented at each depth   =
##================================================================================================================
boxup("SUPPLEMENTARY STATS: Experiment 3: Flanker ring and fixation at 1.26 m with target presented at each depth", bandChar = "=")

Exp3_data <-read.csv("Exp3_accumulatedPerErr_TargetInsideRing.csv", header=T)
Exp3_data<- subset(Exp3_data, Target_Flanker_Spacing_in_deg!=Inf) # Remove control (i.e. no flanker condition)

#Fit full mixed model containing fixed effects and all possible pairwise interactions of the fixed effects 
Exp3_M0<-lmer(log(CircSD_deg)~ Target_Flanker_Spacing_in_deg + Target_Screen +
                Target_Flanker_Spacing_in_deg:Target_Screen + (1|ID), data= Exp3_data, REML=TRUE) 
# The response variable is natural log transformed to correct for positive skew.

#Check normality and homoscedasticity
sresid <- resid(Exp3_M0) # Extract the standardized residuals
hist(sresid) # Plot histogram and QQ plot to look at residuals
qqnorm(sresid)
qqline(sresid)
fitted.glmm <- fitted(Exp3_M0, level=1) # Extract the fitted (predicted) values
plot(sresid ~ fitted.glmm) # Check for homoscedasticity 
plot(Exp3_M0)
#These look fine and acceptable. 

m1 <-update(Exp3_M0,~.-Target_Flanker_Spacing_in_deg:Target_Screen)
anova(Exp3_M0, m1, test = "LRT") #  Chisq=3.5874, Df=2, p = 0.1663 # Reported in supplementary material

#****************************************************************
#m1 is the min model. These are the stats reported in the supplementary material
m2a <-update(m1,~.-Target_Flanker_Spacing_in_deg)
anova(m1, m2a, test = "LRT") # Chisq=50.629, Df=1, p = 1.116e-12 ***
m2b <-update(m1,~.-Target_Screen)
anova(m1, m2b, test = "LRT") # Chisq=36.978, Df=2, p = 9.34e-09 ***
#***************************************************************
rm(list=ls(all=TRUE))

##======================================================================================================================
##  SUPPLEMENTARY STATS: Experiment 4: Flanker ring and fixation at 0.4 m or 4 m with target presented at each depth   =
##======================================================================================================================
boxup("SUPPLEMENTARY STATS: Experiment 4: Flanker ring and fixation at 0.4 m or 4 m with target presented at each depth", bandChar = "=")

Exp4_data <-read.csv("Exp4_accumulatedPerErr_TargetInsideRing.csv", header=T)
Exp4_data <- subset(Exp4_data, Target_Flanker_Spacing_in_deg!=Inf)

#Fit full mixed model containing fixed effects and all possible pairwise interactions of the fixed effects
Exp4_M0<-lmer(log(CircSD_deg)~ Target_Flanker_Spacing_in_deg + Target_Screen + Flanker_Screen + 
                Target_Flanker_Spacing_in_deg:Target_Screen + 
                Target_Flanker_Spacing_in_deg:Flanker_Screen +
                Target_Screen:Flanker_Screen + (1|ID), data= Exp4_data, REML=TRUE) 
# The response variable is natural log transformed to correct for positive skew

#Check normality and homoscedasticity
sresid <- resid(Exp4_M0) # Extract the standardized residuals
hist(sresid) # Plot histogram and QQ plot to look at residuals
qqnorm(sresid)
qqline(sresid)
fitted.glmm <- fitted(Exp4_M0, level=1) # Extract the fitted (predicted) values
plot(sresid ~ fitted.glmm) # Check for homoscedasticity 
plot(Exp4_M0)
#These look acceptable enough 

m1a <-update(Exp4_M0,~.-Target_Screen:Flanker_Screen)
anova(Exp4_M0, m1a, test = "LRT") # Chisq=8.4162  Df=2, p = 0.01487 *
m1b <-update(Exp4_M0,~.-Target_Flanker_Spacing_in_deg:Flanker_Screen)
anova(Exp4_M0, m1b, test = "LRT") # Chisq= 0.79, Df=1, p = 0.3741
m1c <-update(Exp4_M0,~.-Target_Flanker_Spacing_in_deg:Target_Screen)
anova(Exp4_M0, m1c, test = "LRT") # Chisq=0.981, Df=2, p = 0.6123 # Reported in supplementary material
#Drop the Target_Flanker_Spacing_in_deg:Target_Screen interaction

m2a <-update(m1c,~.-Target_Screen:Flanker_Screen)
anova(m1c, m2a, test = "LRT") # Chisq=8.5752, Df=2, p = 0.01374 *
m2b <-update(m1c,~.-Target_Flanker_Spacing_in_deg:Flanker_Screen)
anova(m1c, m2b, test = "LRT") # Chisq=0.6364, Df=1, p = 0.425 # Reported in supplementary material
#Drop Target_Flanker_Spacing_in_deg:Flanker_Screen

#***************************************************************
#m2b is the min model. These are the stats reported in the supplementary material
m3a <-update(m2b,~.-Target_Screen:Flanker_Screen)
anova(m2b, m3a, test = "LRT") # Chisq=8.697, Df=2, p = 0.01293 *
m3b <-update(m2b,~.-Target_Flanker_Spacing_in_deg)
anova(m2b, m3b, test = "LRT") # Chisq=7.5783, Df=1, p = 0.005908 **
#*****************************************************************
rm(list=ls(all=TRUE))

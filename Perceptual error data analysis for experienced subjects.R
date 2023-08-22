# This is a copy of the R script used to generate figures for the the perceptual error data from the experienced 
# subjects who did all five experiments. These figures form part of Smithers et al. (2023)- "Large depth differences 
# between target and flankers can increase crowding: Evidence from a multi-depth plane display", eLife.

# Script written by Dr Samuel P. Smithers, Northeastern University, 2023
# Last edited July 2023

# Corresponding authors SPS (s.smithers@northeastern.edu) and PJB (p.bex@northeastern.edu)

# Included in this script:
# - Code used to generate Figures 10, 11 and 12 from the main manuscript.

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
#   [1] bannerCommenter_1.0.0 svglite_2.1.0         ggplot2_3.4.0        
# 
# loaded via a namespace (and not attached):
#   [1] fansi_1.0.3       withr_2.5.0       dplyr_1.0.10      utf8_1.2.2        grid_4.2.2        R6_2.5.1         
# [7] lifecycle_1.0.3   gtable_0.3.1      magrittr_2.0.3    scales_1.2.1      pillar_1.8.1      rlang_1.0.6      
# [13] cli_3.4.1         generics_0.1.3    vctrs_0.5.0       tools_4.2.2       glue_1.6.2        munsell_0.5.0    
# [19] compiler_4.2.2    systemfonts_1.0.4 pkgconfig_2.0.3   colorspace_2.0-3  tidyselect_1.2.0  tibble_3.1.8   

## Required dependencies/packages. 
#For plots
library(ggplot2) # For graphs
#library(Rmisc)
library(svglite) # Required to save .svg files (optional)
#Other
library(bannerCommenter) # For the section banners (optional)

rm(list=ls(all=TRUE))

#Set up Working directory 
setwd("***Pathway to folder containing this R script and all '_accumulatedPerErr_experiencedSubjects.csv' files***")

#######################################################################################################
#######################################################################################################
###                                                                                                 ###
###  PERCEPTUAL ERROR GRAPHS FOR REPEAT OF STUDY WITH EXPERIENCED SUBJECTS (FIGURES 10, 11 AND 12)  ###
###                                                                                                 ###
#######################################################################################################
#######################################################################################################
banner("Perceptual Error Graphs For Repeat of Study With Experienced Subjects (Figures 10, 11 and 12)", emph = TRUE)

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
        #panel.border = element_blank(),
        panel.background = element_rect(fill="white"),
        panel.grid.major.x=element_blank(),
        panel.grid.minor=element_blank(),
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


##==========================================================================================================================================================
##  GRAPHS: Experiments 1 and 2 with experienced subjects: Target at fixation with flanker ring presented in front, at, and behind fixation  (Figure 10)   =
##==========================================================================================================================================================
boxup("GRAPHS: Experiments 1 and 2 with experienced subjects: Target at fixation with flanker ring presented in front, at, and behind fixation  (Figure 10)", bandChar = "=")

Exp1_data <-read.csv("Exp1_accumulatedPerErr_experiencedSubjects.csv", header=T) # Load data
Exp2_data <-read.csv("Exp2_accumulatedPerErr_experiencedSubjects.csv", header=T)
TargetAtFixationData <- rbind(Exp1_data,Exp2_data)
TargetAtFixationData <- AddExtraCols(TargetAtFixationData)
N <- nrow(unique(TargetAtFixationData["ID"])) 
print(N) # N=4

# Plots for each individual subject. 
Exp_1and2_IndPerErr <- ggplot(TargetAtFixationData, aes(x=factor(Target_Flanker_Spacing_in_deg), y= CircSD_deg, color =  Flanker_depth_in_m)) + 
  geom_point(position=position_dodge(width = 0.5), size=5) +
  geom_errorbar(aes(ymin= pmax(CircSD_deg-BootStrappedCI,0), ymax=CircSD_deg+BootStrappedCI), width=0,
                position=position_dodge(width = 0.5)) +
  xlab("Target-flanker spacing (deg)") + Graph.theme + scale_color_manual(values=cbPalette) +
  facet_grid(ID~Fixation_depth_in_m) +
  scale_y_continuous(breaks= seq(0,180,40), limits=c(0,160), expand=c(0.01,0),  name="Perceptual error (deg)") + 
  guides(color=guide_legend(title="Flanker depth (m)")) +
  scale_x_discrete(labels=c("0.625" = "0.625", "1.25" = "1.25", "2.5" = "2.5", "5" = "5", "Inf" = "no\nflankers")) +
  theme(panel.spacing.y = unit(2, "lines"), legend.text=element_text(size=25),legend.title=element_text(size=25))
Exp_1and2_IndPerErr

##=========================================================================================================================================================
##  GRAPHS: Experiment 3 and 4 with experienced subjects: Flanker ring at fixation with target presented in front, at, and behind fixation  (Figure 11)   =
##=========================================================================================================================================================
boxup("GRAPHS: Experiment 3 and 4 with experienced subjects: Flanker ring at fixation with target presented in front, at, and behind fixation  (Figure 11)", bandChar = "=")

Exp3_data <-read.csv("Exp3_accumulatedPerErr_experiencedSubjects.csv", header=T)
Exp4_data <-read.csv("Exp4_accumulatedPerErr_experiencedSubjects.csv", header=T)
FlankerAtFixationData <- rbind(Exp3_data,Exp4_data)
FlankerAtFixationData <- AddExtraCols(FlankerAtFixationData)
N <- nrow(unique(FlankerAtFixationData["ID"])) 
print(N) # N=4

# Plots for each individual subject. 
Exp_3and4_IndPerErr <- ggplot(FlankerAtFixationData, aes(x=factor(Target_Flanker_Spacing_in_deg), y= CircSD_deg, color =  Target_depth_in_m)) + 
  geom_point(position=position_dodge(width = 0.5), size=5) +
  geom_errorbar(aes(ymin= pmax(CircSD_deg-BootStrappedCI,0), ymax=CircSD_deg+BootStrappedCI), width=0,
                position=position_dodge(width = 0.5)) +
  xlab("Target-flanker spacing (deg)") + Graph.theme + scale_color_manual(values=cbPalette) +
  facet_grid(ID~Fixation_depth_in_m) +
  scale_y_continuous(breaks= seq(0,180,40), limits=c(0,170), expand=c(0.01,0),  name="Perceptual error (deg)") + 
  guides(color=guide_legend(title="Target depth (m)")) +
  scale_x_discrete(labels=c("0.625" = "0.625", "1.25" = "1.25", "2.5" = "2.5", "5" = "5", "Inf" = "no\nflankers")) +
  theme(panel.spacing.y = unit(2, "lines"), legend.text=element_text(size=25),legend.title=element_text(size=25))
Exp_3and4_IndPerErr

##==================================================================================================================
##  GRAPHS: Experiment 5 with experienced subjects: Target and flanker ring always at fixation depth (Figure 12)   =
##==================================================================================================================
boxup("GRAPHS: Experiment 5 with experienced subjects: Target and flanker ring always at fixation depth (Figure 12)", bandChar = "=")

Exp5_data <-read.csv("Exp5_accumulatedPerErr_experiencedSubjects.csv", header=T)
Exp5_data <- AddExtraCols(Exp5_data)
N <- nrow(unique(Exp5_data["ID"]))
print(N) # N=4

# Plots for each individual subject. 
Exp_5_IndPerErr <- ggplot(Exp5_data, aes(x=factor(Target_Flanker_Spacing_in_deg), y= CircSD_deg, color =  Target_depth_in_m)) + 
  geom_point(position=position_dodge(width = 0.5), size=5) +
  geom_errorbar(aes(ymin= pmax(CircSD_deg-BootStrappedCI,0), ymax=CircSD_deg+BootStrappedCI), width=0,
                position=position_dodge(width = 0.5)) +
  xlab("Target-flanker spacing (deg)") + Graph.theme + scale_color_manual(values=cbPalette) +
  facet_grid(ID~1) +
  scale_y_continuous(breaks= seq(0,180,40), limits=c(0,170), expand=c(0.01,0),  name="Perceptual error (deg)") + 
  guides(color=guide_legend(title="Depth (m)")) +
  scale_x_discrete(labels=c("0.625" = "0.625", "1.25" = "1.25", "2.5" = "2.5", "5" = "5", "Inf" = "no flankers")) +
  theme(panel.spacing.y = unit(2, "lines"), legend.text=element_text(size=25),legend.title=element_text(size=25))
Exp_5_IndPerErr 

##---------------------------------------------------------------
##              GRAPHS: Save graphs as .svg files               -
##---------------------------------------------------------------
boxup("GRAPHS: Save graphs as .svg files", bandChar = "-")
# The graphs are then saved as a .svg file to enable formatting, labeling, and other 
# edits required to prepare the figures for the paper.

# Save plots for individual subjects. 
ggsave(file ="Smithers et al- Figure 10.svg", device = 'svg', plot = Exp_1and2_IndPerErr, width = 18, height = 15)
ggsave(file ="Smithers et al- Figure 11.svg", device = 'svg', plot = Exp_3and4_IndPerErr, width = 18, height = 15)
ggsave(file ="Smithers et al- Figure 12.svg", device = 'svg', plot = Exp_5_IndPerErr, width = 10, height = 14)

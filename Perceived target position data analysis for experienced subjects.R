# This is a copy of the R script used to generate supplementary figures for the perceived position of the 
# target relative to the flanker ring that was reported by the experienced subjects who did all five experiments.
# These figures form part of Smithers et al.- "Large differences in target-flanker depth can increase 
# crowding- evidence from a multi-depth plane display". 

# Script written by Dr Samuel P. Smithers, Northeastern University, 2023
# Last edited July 2023

# Corresponding authors SPS (s.smithers@northeastern.edu) and PJB (p.bex@northeastern.edu)

# Included in this script:
# - Code used to generate Figure 10-figure supplements 1 and 2, and Figure 11-figure supplement 1 and 2.

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
#   [1] bannerCommenter_1.0.0 Hmisc_4.8-0           Formula_1.2-4         survival_3.4-0        lattice_0.20-45       gridExtra_2.3        
# [7] svglite_2.1.0         ggplot2_3.4.0        
# 
# loaded via a namespace (and not attached):
#   [1] tidyselect_1.2.0    xfun_0.34           splines_4.2.2       colorspace_2.0-3    vctrs_0.6.1         generics_0.1.3      htmltools_0.5.3    
# [8] base64enc_0.1-3     utf8_1.2.2          rlang_1.1.0         pillar_1.8.1        foreign_0.8-83      glue_1.6.2          withr_2.5.0        
# [15] RColorBrewer_1.1-3  jpeg_0.1-9          lifecycle_1.0.3     stringr_1.4.1       munsell_0.5.0       gtable_0.3.1        htmlwidgets_1.5.4  
# [22] latticeExtra_0.6-30 knitr_1.40          fastmap_1.1.0       fansi_1.0.3         htmlTable_2.4.1     Rcpp_1.0.9          scales_1.2.1       
# [29] backports_1.4.1     checkmate_2.1.0     deldir_1.0-6        systemfonts_1.0.4   interp_1.1-3        png_0.1-7           digest_0.6.30      
# [36] stringi_1.7.8       dplyr_1.1.1         grid_4.2.2          cli_3.6.1           tools_4.2.2         magrittr_2.0.3      tibble_3.2.1       
# [43] cluster_2.1.4       pkgconfig_2.0.3     Matrix_1.5-1        data.table_1.14.8   rstudioapi_0.14     R6_2.5.1            rpart_4.1.19       
# [50] nnet_7.3-18         compiler_4.2.2   

##Required dependencies/packages. 
#For plots
library(ggplot2) # For graphs
library(svglite) # Required to save .svg files (optional)
library(gridExtra) # For arranging/combining plots (optional)
library(Hmisc) #For calculating confidence intervals
#Other
library(bannerCommenter) # For the section banners (optional)

rm(list=ls(all=TRUE))

setwd("***Pathway to folder containing this R script and all '_accumulatedPerErr_experiencedSubjects.csv' files***")

#######################################################################################################################################
#######################################################################################################################################
###                                                                                                                                 ###
###  PERCEIVED TARGET POSITION RELATIVE TO FLANKER RING GRAPHS FOR FIGURE 10-FIGURE SUPPLEMENT 1 AND FIGURE 11-FIGURE SUPPLEMENT 1  ###
###                                                                                                                                 ###
#######################################################################################################################################
#######################################################################################################################################
banner("Perceived Target Position Relative to Flanker Ring Graphs for Figure 10-figure supplement 1 and Figure 11-figure supplement 1", emph = TRUE)

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
        #panel.border = element_blank(),
        panel.background = element_rect(fill="white"),
        panel.grid.major.x=element_blank(),
        panel.grid.minor=element_blank(),
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
  # Calulate Wilson score intervals. 
  Output <- binconf(Inputdata[,"PerceivedPosCount_insideT"],repeatno., alpha=0.05,
                    method=c("wilson"),
                    include.x=TRUE, include.n=TRUE, return.df=FALSE)
  Inputdata["Proportion_insideT"] <- Output[,"PointEst"]
  Inputdata["Lower_WilsonCI"] <- Output[,"Lower"]
  Inputdata["Upper_WilsonCI"] <- Output[,"Upper"]
  return(Inputdata)
}



##============================================================================================================================================
##  SUPPLEMENTARY GRAPHS: Experiment 1 and 2: Target at fixation with flanker ring presented at each depth (Figure 10-figure supplement 1)   =
##============================================================================================================================================
boxup("SUPPLEMENTARY GRAPHS: Experiment 1 and 2: Target at fixation with flanker ring presented at each depth (Figure 10-figure supplement 1)", bandChar = "=")

Exp1_data <-read.csv("Exp1_accumulatedPerErr_experiencedSubjects.csv", header=T) # Load data
Exp2_data <-read.csv("Exp2_accumulatedPerErr_experiencedSubjects.csv", header=T)
TargetAtFixationData <- rbind(Exp1_data,Exp2_data)
TargetAtFixationData <- AddExtraCols(TargetAtFixationData)
N <- nrow(unique(TargetAtFixationData["ID"])) 
print(N) # n=4
repeatno. = 10 
TargetAtFixationData <- ProportionInsideT(TargetAtFixationData,repeatno.)

# Plots for each individual subject. 
Exp1and2_IndPropInside <- ggplot(TargetAtFixationData, aes(x=factor(Target_Flanker_Spacing_in_deg), y= Proportion_insideT, color =  FlankerDepth_in_m)) + 
  geom_point(position=position_dodge(width = 0.5), size=5) +
  geom_errorbar(aes(ymin= Lower_WilsonCI, ymax= Upper_WilsonCI), width=0,
                position=position_dodge(width = 0.5)) +
  xlab("Target-flanker spacing (deg)") + Graph.theme + scale_color_manual(values=cbPalette) +
  facet_grid(ID~FixationDepth_in_m) +
  scale_y_continuous(breaks= seq(0,1,0.2), limits=c(-0.02,1.02), expand=c(0.1,0),  name="Proportion of trials") + 
  guides(color=guide_legend(title="Flanker depth (m)")) +
  scale_x_discrete(drop=FALSE, labels=c("0.625" = "0.625", "1.25" = "1.25", "2.5" = "2.5", "5" = "5", "Inf" = "no\nflankers")) +
  theme(panel.spacing.y = unit(2, "lines"), legend.text=element_text(size=25),legend.title=element_text(size=25))
Exp1and2_IndPropInside



##=============================================================================================================================================
##  SUPPLEMENTARY GRAPHS: Experiment 3 and 4: Flanker ring at fixation with target presented at each depth (Figure 11-figure supplement 1)    =
##=============================================================================================================================================
boxup("SUPPLEMENTARY GRAPHS: Experiment 3 and 4: Flanker ring at fixation with target presented at each depth (Figure 11-figure supplement 1) ", bandChar = "=")

Exp3_data <-read.csv("Exp3_accumulatedPerErr_experiencedSubjects.csv", header=T) # Load data
repeatno. = 10 
Exp3_data <- ProportionInsideT(Exp3_data,repeatno.)
Exp4_data <-read.csv("Exp4_accumulatedPerErr_experiencedSubjects.csv", header=T)
repeatno. = 8
Exp4_data <- ProportionInsideT(Exp4_data,repeatno.)
FlankerAtFixationData <- rbind(Exp3_data,Exp4_data)
FlankerAtFixationData <- AddExtraCols(FlankerAtFixationData)
N <- nrow(unique(FlankerAtFixationData["ID"])) 
print(N) # n=4

# Plots for each individual subject. 
Exp3and4_IndPropInside <- ggplot(FlankerAtFixationData, aes(x=factor(Target_Flanker_Spacing_in_deg), y= Proportion_insideT, color =  TargetDepth_in_m)) + 
  geom_point(position=position_dodge(width = 0.5), size=5) +
  geom_errorbar(aes(ymin= Lower_WilsonCI, ymax= Upper_WilsonCI), width=0,
                position=position_dodge(width = 0.5)) +
  xlab("Target-flanker spacing (deg)") + Graph.theme + scale_color_manual(values=cbPalette) +
  facet_grid(ID~FixationDepth_in_m) +
  scale_y_continuous(breaks= seq(0,1,0.2), limits=c(-0.02,1.02), expand=c(0.1,0),  name="Proportion of trials") + 
  guides(color=guide_legend(title="Target depth (m)")) +
  scale_x_discrete(drop=FALSE, labels=c("0.625" = "0.625", "1.25" = "1.25", "2.5" = "2.5", "5" = "5", "Inf" = "no\nflankers")) +
  theme(panel.spacing.y = unit(2, "lines"), legend.text=element_text(size=25),legend.title=element_text(size=25))
Exp3and4_IndPropInside
  

##---------------------------------------------------------------
##              GRAPHS: Save graphs as .svg files               -
##---------------------------------------------------------------
boxup("GRAPHS: Save graphs as .svg files", bandChar = "-")
# The graphs are then saved as a .svg file to enable formatting, labeling, and other 
# edits required to prepare the figures for the paper.

# Save plots for individual subjects. 
ggsave(file ="Smithers et al- Figure 10-figure supplement 1.svg", device = 'svg', plot = Exp1and2_IndPropInside, width = 18, height = 15)
ggsave(file ="Smithers et al- Figure 11-figure supplement 1.svg", device = 'svg', plot = Exp3and4_IndPropInside, width = 18, height = 15)


#######################################################################################################################################
#######################################################################################################################################
###                                                                                                                                 ###
###  PERCEIVED TARGET POSITION RELATIVE TO FLANKER RING GRAPHS FOR FIGURE 10-FIGURE SUPPLEMENT 2 AND FIGURE 11-FIGURE SUPPLEMENT 2  ###
###                                                                                                                                 ###
#######################################################################################################################################
#######################################################################################################################################
banner("Perceived Target Position Relative to Flanker Ring Graphs for Figure 10-figure supplement 2 and Figure 11-figure supplement 2", emph = TRUE)

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
        #panel.border = element_blank(),
        panel.background = element_rect(fill="white"),
        panel.grid.major.x=element_blank(),
        panel.grid.minor=element_blank(),
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
CountDataFunction <- function(Inputdata,repeatno.){
  Count <- c()
  Perception <- c()
  TargetFlankerSpacing <-c()
  FlankerDepth <- c()
  TargetDepth <- c()
  FixationDepth <- c()
  ID <-c()
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
  ID <- append(ID,rep(Inputdata$ID,5))
  CountData <- data.frame(Count,Perception, TargetFlankerSpacing, FlankerDepth, TargetDepth, FixationDepth, ID)
  CountData["Proportion"] <- CountData$Count/repeatno.
  return(CountData) 
}

##===============================================================================================================================================
##  SUPPLEMENTARY GRAPHS: Experiments 1 and 2: Target at fixation with flanker ring presented at each depth  (Figure 10-figure supplement 2).   =
##===============================================================================================================================================
boxup("SUPPLEMENTARY GRAPHS: Experiments 1 and 2: Target at fixation with flanker ring presented at each depth  (Figure 10-figure supplement 2).", bandChar = "=")

Exp1_data <-read.csv("Exp1_accumulatedPerErr_experiencedSubjects.csv", header=T) # Load data
Exp2_data <-read.csv("Exp2_accumulatedPerErr_experiencedSubjects.csv", header=T)
TargetAtFixationData <- rbind(Exp1_data,Exp2_data)
TargetAtFixationData <- AddExtraCols(TargetAtFixationData)
repeatno. = 10 
CountData <- CountDataFunction(TargetAtFixationData,repeatno.)


PlotCountData <- subset(CountData, Perception == "Target in center of ring")
Exp1and2_plot_TCenter <- ggplot(PlotCountData, aes(x=factor(TargetFlankerSpacing), y= Proportion, color= FlankerDepth)) + 
  geom_point(position=position_dodge(width = 0.75), size=5) +
  xlab("Target-flanker spacing (deg)") + Graph.theme + scale_color_manual(values=cbPalette) +
  ggtitle("Target in center of flanker ring") +
  facet_grid(ID~FixationDepth) +
  scale_y_continuous(breaks= seq(0,1,0.2), limits=c(-0.02,1.02), expand=c(0.1,0),  name="Proportion of trials") + 
  guides(color=guide_legend(title="Flanker depth (m)")) +
  scale_x_discrete(labels=c("0.625" = "0.625", "1.25" = "1.25", "2.5" = "2.5", "5" = "5", "Inf" = "no\nflankers")) +
  theme(legend.text=element_text(size=25),legend.title=element_text(size=25))
Exp1and2_plot_TCenter

PlotCountData <- subset(CountData, Perception == "Target not in center of ring")
Exp1and2_plot_TOffCenter <- ggplot(PlotCountData, aes(x=factor(TargetFlankerSpacing), y= Proportion, color= FlankerDepth)) + 
  geom_point(position=position_dodge(width = 0.75), size=5) +
  xlab("Target-flanker spacing (deg)") + Graph.theme + scale_color_manual(values=cbPalette) +
  ggtitle("Target inside, but not in center, of flanker ring") +
  facet_grid(ID~FixationDepth) +
  scale_y_continuous(breaks= seq(0,1,0.2), limits=c(-0.02,1.02), expand=c(0.1,0),  name="Proportion of trials") + 
  guides(color=guide_legend(title="Flanker depth (m)")) +
  scale_x_discrete(labels=c("0.625" = "0.625", "1.25" = "1.25", "2.5" = "2.5", "5" = "5", "Inf" = "no\nflankers")) +
  theme(legend.text=element_text(size=25),legend.title=element_text(size=25))
Exp1and2_plot_TOffCenter

PlotCountData <- subset(CountData, Perception == "Ring obstructs target")
Exp1and2_plot_TObstructed <- ggplot(PlotCountData, aes(x=factor(TargetFlankerSpacing), y= Proportion, color= FlankerDepth)) + 
  geom_point(position=position_dodge(width = 0.75), size=5) +
  xlab("Target-flanker spacing (deg)") + Graph.theme + scale_color_manual(values=cbPalette) +
  ggtitle("Flanker ring obstructs target") +
  facet_grid(ID~FixationDepth) +
  scale_y_continuous(breaks= seq(0,1,0.2), limits=c(-0.02,1.02), expand=c(0.1,0),  name="Proportion of trials") + 
  guides(color=guide_legend(title="Flanker depth (m)")) +
  scale_x_discrete(labels=c("0.625" = "0.625", "1.25" = "1.25", "2.5" = "2.5", "5" = "5", "Inf" = "no\nflankers")) +
  theme(legend.text=element_text(size=25),legend.title=element_text(size=25))
Exp1and2_plot_TObstructed

PlotCountData <- subset(CountData, Perception == "Target outside ring")
Exp1and2_plot_TOutside <- ggplot(PlotCountData, aes(x=factor(TargetFlankerSpacing), y= Proportion, color= FlankerDepth)) + 
  geom_point(position=position_dodge(width = 0.75), size=5) +
  xlab("Target-flanker spacing (deg)") + Graph.theme + scale_color_manual(values=cbPalette) +
  ggtitle("Target outside flanker ring") +
  facet_grid(ID~FixationDepth) +
  scale_y_continuous(breaks= seq(0,1,0.2), limits=c(-0.02,1.02), expand=c(0.1,0),  name="Proportion of trials") + 
  guides(color=guide_legend(title="Flanker depth (m)")) +
  scale_x_discrete(labels=c("0.625" = "0.625", "1.25" = "1.25", "2.5" = "2.5", "5" = "5", "Inf" = "no\nflankers")) +
  theme(legend.text=element_text(size=25),legend.title=element_text(size=25))
Exp1and2_plot_TOutside
  
PlotCountData <- subset(CountData, Perception == "Unsure or no ring")
Exp1and2_plot_noRingOrUnsure <- ggplot(PlotCountData, aes(x=factor(TargetFlankerSpacing), y= Proportion, color= FlankerDepth)) + 
  geom_point(position=position_dodge(width = 0.75), size=5) +
  xlab("Target-flanker spacing (deg)") + Graph.theme + scale_color_manual(values=cbPalette) +
  ggtitle("No flanker ring or unsure") +
  facet_grid(ID~FixationDepth) +
  scale_y_continuous(breaks= seq(0,1,0.2), limits=c(-0.02,1.02), expand=c(0.1,0),  name="Proportion of trials") + 
  guides(color=guide_legend(title="Flanker depth (m)")) +
  scale_x_discrete(labels=c("0.625" = "0.625", "1.25" = "1.25", "2.5" = "2.5", "5" = "5", "Inf" = "no\nflankers")) +
  theme(legend.text=element_text(size=25),legend.title=element_text(size=25)) 
Exp1and2_plot_noRingOrUnsure

##----------------------------------------------------------------------------
##  SUPPLEMENTARY GRAPHS: Save Figure 10-figure supplement 2 as .svg files   -
##---------------------------------------------------------------------------- 
boxup("SUPPLEMENTARY GRAPHS: Save Figure 10-figure supplement 2 as .svg files", bandChar = "-")
# The graphs are saved as a .svg file to enable formatting, labeling, and other 
# edits required to prepare the figures for the supplementary material.

#Save plots
JoinedPlots <- grid.arrange(Exp1and2_plot_TCenter, Exp1and2_plot_TOffCenter, Exp1and2_plot_TObstructed, 
                            Exp1and2_plot_TOutside, Exp1and2_plot_noRingOrUnsure, ncol = 2)
ggsave(file ="Smithers et al- Figure 10-figure supplement 2.svg", device = 'svg', plot = JoinedPlots, width = 25, height = 32)

##===============================================================================================================================================
##  SUPPLEMENTARY GRAPHS: Experiments 3 and 4: Flanker ring at fixation with target presented at each depth  (Figure 11-figure supplement 2).   =
##===============================================================================================================================================
boxup("SUPPLEMENTARY GRAPHS: Experiments 3 and 4: Flanker ring at fixation with target presented at each depth  (Figure 11-figure supplement 2).", bandChar = "=")

Exp3_data <-read.csv("Exp3_accumulatedPerErr_experiencedSubjects.csv", header=T) # Load data
Exp3_data <- AddExtraCols(Exp3_data)
repeatno. = 10 
CountData_exp3 <- CountDataFunction(Exp3_data,repeatno.)
Exp4_data <-read.csv("Exp4_accumulatedPerErr_experiencedSubjects.csv", header=T)
Exp4_data <- AddExtraCols(Exp4_data)
repeatno. = 8
CountData_exp4 <- CountDataFunction(Exp4_data,repeatno.)
CountData<- rbind(CountData_exp3,CountData_exp4)

PlotCountData <- subset(CountData, Perception == "Target in center of ring")
Exp3and4_plot_TCenter <- ggplot(PlotCountData, aes(x=factor(TargetFlankerSpacing), y= Proportion, color= TargetDepth)) + 
  geom_point(position=position_dodge(width = 0.75), size=5) +
  xlab("Target-flanker spacing (deg)") + Graph.theme + scale_color_manual(values=cbPalette) + 
  ggtitle("Target in center of flanker ring") +
  facet_grid(ID~FixationDepth) +
  scale_y_continuous(breaks= seq(0,1,0.2), limits=c(-0.02,1.02), expand=c(0.1,0),  name="Proportion of trials") + 
  guides(color=guide_legend(title="Target depth (m)")) +
  scale_x_discrete(labels=c("0.625" = "0.625", "1.25" = "1.25", "2.5" = "2.5", "5" = "5", "Inf" = "no\nflankers")) +
  theme(legend.text=element_text(size=25),legend.title=element_text(size=25))
Exp3and4_plot_TCenter

PlotCountData <- subset(CountData, Perception == "Target not in center of ring")
Exp3and4_plot_TOffCenter <- ggplot(PlotCountData, aes(x=factor(TargetFlankerSpacing), y= Proportion, color= TargetDepth)) + 
  geom_point(position=position_dodge(width = 0.75), size=5) +
  xlab("Target-flanker spacing (deg)") + Graph.theme + scale_color_manual(values=cbPalette) + 
  ggtitle("Target inside, but not in center, of flanker ring") +
  facet_grid(ID~FixationDepth) +
  scale_y_continuous(breaks= seq(0,1,0.2), limits=c(-0.02,1.02), expand=c(0.1,0),  name="Proportion of trials") + 
  guides(color=guide_legend(title="Target depth (m)")) +
  scale_x_discrete(labels=c("0.625" = "0.625", "1.25" = "1.25", "2.5" = "2.5", "5" = "5", "Inf" = "no\nflankers")) +
  theme(legend.text=element_text(size=25),legend.title=element_text(size=25))
Exp3and4_plot_TOffCenter

PlotCountData <- subset(CountData, Perception == "Ring obstructs target")
Exp3and4_plot_TObstructed <- ggplot(PlotCountData, aes(x=factor(TargetFlankerSpacing), y= Proportion, color= TargetDepth)) + 
  geom_point(position=position_dodge(width = 0.75), size=5) +
  xlab("Target-flanker spacing (deg)") + Graph.theme + scale_color_manual(values=cbPalette) + 
  ggtitle("Flanker ring obstructs target") +
  facet_grid(ID~FixationDepth) +
  scale_y_continuous(breaks= seq(0,1,0.2), limits=c(-0.02,1.02), expand=c(0.1,0),  name="Proportion of trials") + 
  guides(color=guide_legend(title="Target depth (m)")) +
  scale_x_discrete(labels=c("0.625" = "0.625", "1.25" = "1.25", "2.5" = "2.5", "5" = "5", "Inf" = "no\nflankers")) +
  theme(legend.text=element_text(size=25),legend.title=element_text(size=25))
Exp3and4_plot_TObstructed

PlotCountData <- subset(CountData, Perception == "Target outside ring")
Exp3and4_plot_TOutside <- ggplot(PlotCountData, aes(x=factor(TargetFlankerSpacing), y= Proportion, color= TargetDepth)) + 
  geom_point(position=position_dodge(width = 0.75), size=5) +
  xlab("Target-flanker spacing (deg)") + Graph.theme + scale_color_manual(values=cbPalette) + 
  ggtitle("Target outside flanker ring") +
  facet_grid(ID~FixationDepth) +
  scale_y_continuous(breaks= seq(0,1,0.2), limits=c(-0.02,1.02), expand=c(0.1,0),  name="Proportion of trials") + 
  guides(color=guide_legend(title="Target depth (m)")) +
  scale_x_discrete(labels=c("0.625" = "0.625", "1.25" = "1.25", "2.5" = "2.5", "5" = "5", "Inf" = "no\nflankers")) +
  theme(legend.text=element_text(size=25),legend.title=element_text(size=25))
Exp3and4_plot_TOutside

PlotCountData <- subset(CountData, Perception == "Unsure or no ring")
Exp3and4_plot_noRingOrUnsure <- ggplot(PlotCountData, aes(x=factor(TargetFlankerSpacing), y= Proportion, color= TargetDepth)) + 
  geom_point(position=position_dodge(width = 0.75), size=5) +
  xlab("Target-flanker spacing (deg)") + Graph.theme + scale_color_manual(values=cbPalette) + 
  ggtitle("No flanker ring or unsure") +
  facet_grid(ID~FixationDepth) +
  scale_y_continuous(breaks= seq(0,1,0.2), limits=c(-0.02,1.02), expand=c(0.1,0),  name="Proportion of trials") + 
  guides(color=guide_legend(title="Target depth (m)")) +
  scale_x_discrete(labels=c("0.625" = "0.625", "1.25" = "1.25", "2.5" = "2.5", "5" = "5", "Inf" = "no\nflankers")) +
  theme(legend.text=element_text(size=25),legend.title=element_text(size=25))
Exp3and4_plot_noRingOrUnsure

##----------------------------------------------------------------------------
##  SUPPLEMENTARY GRAPHS: Save Figure 11-figure supplement 2 as .svg files   -
##----------------------------------------------------------------------------
boxup("SUPPLEMENTARY GRAPHS: Save Figure 11-figure supplement 2 as .svg files", bandChar = "-")
# The graphs are saved as a .svg file to enable formatting, labeling, and other 
# edits required to prepare the figures for the supplementary material.

#Save plots
JoinedPlots <- grid.arrange(Exp3and4_plot_TCenter, Exp3and4_plot_TOffCenter, Exp3and4_plot_TObstructed, 
                            Exp3and4_plot_TOutside, Exp3and4_plot_noRingOrUnsure, ncol = 2)
ggsave(file ="Smithers et al- Figure 11-figure supplement 2.svg", device = 'svg', plot = JoinedPlots, width = 25, height = 32)

# Supplementary materials: Large depth differences between target and flankers can increase crowding: Evidence from a multi-depth plane display
### Supplementary information for: Samuel P. Smithers, Yulong Shao, James Altham, and Peter J. Bex. (eLife, 2023) *Large depth differences between target and flankers can increase crowding: Evidence from a multi-depth plane display.* DOI: https://doi.org/10.7554/eLife.85143

## Repository contents
- Raw date from Experiments 1-5 (files ending "*_rawData.csv*").
- ```Cal_CircSD.m```: Matlab code used to accumulate the raw data from each experiment and calculate perceptual error. Outputs processed data as "*_accumulatedPerErr.csv*" or "*_accumulatedPerErr_TargetInsideRing.csv*".
- Accumulated perceptual error data for Experiments 1-5 generated by ```Cal_CircSD.m``` (files ending "*_accumulatedPerErr.csv*").
- Accumulated perceptual error data for Experiments 1-4 based only on trials in which the subject reported seeing the target inside the flanker ring generated by ```Cal_CircSD.m``` (files ending "*_accumulatedPerErr_TargetInsideRing.csv*").
- Accumulated perceptual error data for repeat of Experiments 1-5 with small number of experienced subjects generated by ```Cal_CircSD.m``` (files ending "*_accumulatedPerErr_experiencedSubjects.csv*").
- ```Perceptual error data analysis.R```: R code used to generate the figures for, and perform the statistical analysis on, the perceptual error data from the main study reported in the main manuscript and supplementary figures.
- ```Perceived target position data analysis.R```: R code used to generate the figures for, and perform the statistical analysis on, the data for the perceived position of the target relative to the flanker ring that was reported by the naive observers in Experiments 1-4 of the main study.
- ```Perceptual error data analysis for experienced subjects.R```: R code used to generate the figures for the perceptual error data from the small number of experienced observers.
- ```Perceived target position data analysis for experienced subjects.R```: R code used to generate the figures for the data for the perceived position of the target relative to the flanker ring that was reported by the small number of experienced observers in the repeat of Experiments 1-4. 

## Calculating perceptual error (i.e. circular standard deviation) & accumulating raw data (MATLAB)
### Requirements & dependencies
- MATLAB R2022a (other versions should work but this one is tested)
- [MATLAB Circular Statistics Toolbox](https://www.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics)

#### Running instructions
1. Manually download ```Cal_CircSD.m``` and raw data files from this repository. The raw data from each experiment should be kept in separate folders for the script to work. Alternatively, you can clone this repository ([instructions here for MATLAB](https://www.mathworks.com/help/simulink/ug/clone-git-repository.html)). 
2. Install any uninstalled dependencies, above.
3. Open ```Cal_CircSD.m``` and follow instructions within the script. 

## Statistical analysis & plotting (R)
### Requirements & dependencies
- R version >= 4
- The following packages must be installed: 
  - ggplot2
  - ggh4x
  - lme4
  - Matrix
  - DHARMa
  - svglite (optional to save plots)
  - gridExtra (optional to arrange/combine plots for saving)
  - cowplot (optional to arrange/combine plots for saving)
  - bannerCommenter (optional)

#### Running instructions
1. Manually download ```Perceptual error data analysis.R``` and ```Perceived target position data analysis.R``` from this repository. You will also need to download the "*_accumulatedPerErr.csv*" and "*_accumulatedPerErr_TargetInsideRing.csv*" data files (required by both R scripts) along with the raw data files (required for ```Perceived target position data analysis.R``` only). The raw data from each experiment should be kept in separate folders for the script to work. Alternatively, you can clone this repository ([instructions for Rstudio here](https://datacarpentry.org/rr-version-control/03-git-in-rstudio/index.html)). 
2. Install any uninstalled dependencies, above.
3. Open ```Perceptual error data analysis.R``` and/or ```Perceived target position data analysis.R``` and follow instructions within the script.
4. The process is the same for running the R scripts for the data from the experienced subjects (```Perceptual error data analysis for experienced subjects.R``` and ```Perceived target position data analysis for experienced subjects.R```), except you will need to download the accumulated and raw data for the experienced subjects instead.

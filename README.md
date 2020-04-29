# Detect Marine Heatwaves in the Arctic region 

The codes were written for BSc Marine Science 4th year student dissertation project on mapping marine heatwaves in the Arctic region. The codes written in MATLAB were aimed to study the characteristics of MHWs in the region north of 60°N latitude from 1983 to 2012 (30 years). The MHW event is defined as an anomalously warm-water event lasting for at least 5 consecutive days with daily SST higher than the 90th percentile over a 30-year historical baseline period, and two successive events with a temporal gap of less than 3 days are to be treated as a single joint event (Hobday et al. 2016, pp. 227-238). 

Previous appraoches to detecting and analysing global MHWs had been done by Eric C. J. Oliver (2015) in Python, Schlegel and Smit (2018) in R, and Zhao and Marin (2019) in MATLAB. 

The National Oceanic and Atmospheric Administration (NOAA) Optimum Interpolation Sea Surface Temperature (OISST) V2.0 high-resolution gridded data set was used for this study and was provided by NOAA/OAR/ESRL PSD (Boulder, CO, USA) from their website (https://psl.noaa.gov/). Additional MATLAB functions, which are 'findND' and 'trend' functions, were applied in this code and can be downloaded from https://uk.mathworks.com/matlabcentral/fileexchange/64383-findnd and https://uk.mathworks.com/matlabcentral/fileexchange/46363-trend, respectively (Greene 2020; Rik 2020).


# Contents

The codes are seperated into 3 MATLAB scripts to generate the results of MHWs detection and analysis. These 3 scripts are:

|File                 |Description|
|---------------------|-----------|
|MHW_detect.m         |Detection of MHWs, which includes loading OISST data, calculations of the climatological mean and threshold, recording every MHW in each location and calculate the duration, mean intensity, maximum intensity, cumulative intensity and intensity variability (standard deviation) of each MHW event|
|MHW_annual.m         |Documentation folder|
|MHW_analysis.m       |Software license information|


# Acknowledgement

The scripts was 


# Reference

Greene, C. (2020) trend. MATLAB Central File Exchange. Available at: https://www.mathworks.com/matlabcentral/fileexchange/46363-trend

Hobday, A.J. et al. (2016) A hierarchical approach to defining marine
heatwaves. Progress in Oceanography, 141, pp. 227-238 https://doi.org/10.1016/j.pocean.2015.12.014

Oliver E. C. J. (2015) Marine Heatwaves detection code. Available at: https://github.com/ecjoliver/marineHeatWaves

Rik (2020). findND. MATLAB Central File Exchange. Available at: https://www.mathworks.com/matlabcentral/fileexchange/64383-findnd

Schlegel, R. W. and Smit, A. J, (2018) heatwaveR: A central algorithm for the detection of heatwaves and cold-spells. The Journal of Open Source Software, 3, p.821  https://doi.org/10.21105/joss.00821

Zhao, Z. and Marin, M. (2019) A MATLAB toolbox to detect and analyze marine heatwaves. Journal of Open Source Software, 4(33), 1124, https://doi.org/10.21105/joss.01124

# Contact

Shi Quan Ooi

The Scottish Association for Marine Science,

University of the Highlands and Islands

<17002380@uhi.ac.uk> 


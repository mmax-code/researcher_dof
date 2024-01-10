# Addressing researcher degrees of freedom through minP adjustment

This repository contains material to reproduce the results of the article "Addressing researcher degrees of freedom through minP adjustment" by Maximilian Mandl, Andrea Becker-Pennrich, Christian Hinske, Sabine Hoffmann and Anne-Laure Boulesteix. 

Session Info:
R version 4.2.2 Patched (2022-11-10 r83330)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Debian GNU/Linux 12 (bookworm)

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] doParallel_1.0.17 iterators_1.0.14  foreach_1.5.2     viridis_0.6.2     viridisLite_0.4.1
 [6] reshape2_1.4.4    dplyr_1.0.10      broom_1.0.3       fastDummies_1.7.3 mvtnorm_1.1-3    
[11] ggplot2_3.4.1     rlist_0.4.6.2    

The repository contains of four R-scripts and one results folder: functions.R, analysis_minP.R, analysis_fwer.R. To reproduce the results, please set the working directory to the repository and the number of cores to use
in both files. analysis_fwer.R produces Figure 2 and analysis_nimP.R Figure 3 of the manuscript.

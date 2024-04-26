This folder contains the metrics following flow cytometry for each biological replicate:
311020 = replicate 1
091120 = replicate 2 
161220 = replicate 3

Within each '...c.xlsx' file, the mean, standard deviation and number of cells per subpopulation has been calculated using the scripts in the 'YTK-heterogeneity-analysis' repository.
- mu1-3 = mean of subpopulation 1-3
- sd1-3 = standard deviation of subpopulation 1-3
- lam1-3 = number of cells in subpopulation 1-3
- n = total number of cells sampled 

In some '...c.xlsx' files, a row has been deleted. This is the case if the calculated number of subpopulations did not agree with the other 2 replicates. This is to avoid skewing the mean of individual subpopulations. 


The 'long_map_labels' file contains metrics such as the sfGFP intensities, max OD600, doubling time and time to midlog as an average from all three replicates. pTDH3 in YEPG was deleted due to lack of detectable growth. 


The plots generated here, were used to make Figure 1 and Supplementary Figures S1 and S2

This folder contains the metrics following flow cytometry for each biological replicate. 

Within each '...c.xlsx' file, the mean, standard deviation and number of cells per subpopulation has been calculated using the scripts in the 'YTK-heterogeneity-analysis' folder. 
- mu1-3 = mean of subpopulation 1-3
- sd1-3 = standard deviation of subpopulation 1-3
- lam1-3 = number of cells in subpopulation 1-3
- n = total number of cells sampled 

In some '...c.xlsx' files, a row has been deleted. This is the case if the calculated number of subpopulations did not agree with the other 2 replicates. This is to avoid skewing the mean of individual subpopulations. 
Complete tables containing the calculated number of subpopulations can be found in YTK_promoter_flow>...run2>...summary_stats


The 'replicate_metrics_for_long_map' folder contains the metrics for each biological replicate that was used to create the 'long_map_labels' file - there is a README in that folder for more information. 


The 'savedplots' folder contains the figures generated from the 'means_2023_2.R' script:

expplots3 - figure showing exponential phase expression of sfGFP by YTK promoters
staplots - figure showing stationary phase expression of sfGFP by YTK promoters
x1plots - figure showing exponential phase expression of sfGFP by YTK promoters in YPD10 media. The figure used in Fig.1 of the paper.

All files named '...annotated' are labelled versions of their corresponding file. 
The plots are representative of the average expression across 3 biological replicates.



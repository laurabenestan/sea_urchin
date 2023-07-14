# Sea_urchin

Here, we investigate population structure in sea urchin population genomics. 
First a sampling map showing all sampling locations in relation to oceanic currents was created with the R script.


Figure 1. Sampling design

We use the function `genet.dist` available in [hierfstat](https://rdrr.io/cran/hierfstat/man/genet.dist.html) R package.
This function estimate the pairwise FST values according to Weir & Cockerham (1984).
FST measures the amount of genetic differentiation among populations (and simultaneously the extent to which individuals within populations are similar to one another).
A total of 55 pairwise comparisons (i.e. (n*(n-1))/2; 11*10)/2) ranges from 0 to 0.055.

![Figure 4a - Index of genetic differentiation](FST-matrix.png){weight=60%}

A total of 48 P-values were non significant.
The P-values were tighly associated to the boostrap tree obtained on FST.

![Figure 4b - UPGMA tree with bootstraps values on genetic differentiation index](tree-bootstap.png){weight=60%}


|Variation	 %var  |F-stat|F-value|c.i.2.5% |c.i.97.5%|	P-value	   | 
|------------------|------|-------|---------|---------|------------|
|Within individuals|0.669 |	Fit	  |0.331	  |0.327	  |0.334	--	 |
|Within location	 |0.310 |	Fis	  |0.317	  |0.313	  |0.320	0.001|
|Among location	   |0.003 |	Fsc	  |0.004	  |0.003	  |0.005	0.001|	
|Among regions	   |0.017 |	Fct	  |0.017	  |0.016	  |0.018	0.001|	



Figure 6. IBD nucléaire versus mtDNA / davantage perméabilité nucléaire mt
Du point 1 au point 6 pour le chapitre d'Indra

Figure 7. Outlier detection - description 


Figure 8. RDA tous SNPs - description - tableau variables significativité 
- RDA partielle


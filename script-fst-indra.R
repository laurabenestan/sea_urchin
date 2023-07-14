' --------------------------------------------------------------------------   @Header
#'
#' @title Fst calculations 
#'
#' @description
#' This R script
#'
#' @author  Laura Benestan
#'
#' @date 2023/07/13
#'
#' --------------------------------------------------------------------------   @libraries

### Download libraries
library(ggplot2)
library(reshape2)
library(dplyr)
library(vcfR)
library(tidyr)
library(hierfstat)
library(adegenet)
library(poppr)
library(stats)
library(RColorBrewer)
library(phangorn)

###### STEP 1 : DOWNLOAD DATASET  ####

### Download vcf
vcf <- vcfR::read.vcfR("../00-Data/02-filtering/01-21597snps-87ind/21585snps_87ind.recode.vcf")

### Download pop info
pop <- read.table("../00-Data/02-filtering/01-21597snps-87ind/population_map_87individuals.txt", header=TRUE, sep="\t", stringsAsFactors = TRUE)

### Transform vcf to genind
genind <- vcfR::vcfR2genind(vcf)

### Define population map 
genind@pop <- as.factor(pop$CODE) #### Calculate fst on this vcf data

###### STEP 2 : FILTER DATASET  ####

### Remove loci with too many missing
loci <- poppr::missingno(genind, type = "loci", cutoff = 0.3, quiet = FALSE, freq = FALSE)
loci

### Remove individuals with too many missing
loci.individuals <- poppr::missingno(loci, type = "geno", cutoff = 0.3, quiet = FALSE, freq = FALSE)
loci.individuals

### Apply a MAF criteria
loci.individuals.maf <- poppr::informloci(loci.individuals, cutoff = 2/nInd(loci.individuals), MAF = 0.05, quiet = FALSE)
loci.individuals.maf

### Define population map 
loci.individuals.maf@pop <- pop$CODE

###### STEP 3 : CALCULATE FST  ####

### Estimate fst
urchin_fst = genet.dist(loci.individuals.maf, method = "WC84")
save.image("fst_urchins.Rdata")

### Upload fst data
load("03-fst/fst_urchins.Rdata")

# Convert dist object to data.frame
fst.matrix = as.matrix(urchin_fst)
ind = which( upper.tri(fst.matrix), arr.ind = TRUE)
fst.df = data.frame(Site1 = dimnames(fst.matrix)[[2]][ind[,2]],
                    Site2 = dimnames(fst.matrix)[[1]][ind[,1]], Fst = fst.matrix[ ind ] %>% round(digits = 3))
heatmap(fst.matrix)

### Convert minus values to zero
fst.df$Fst[fst.df$Fst < 0] = 0

### Print data.frame summary
fst.df %>% str

### Fst italic label
fst.label = expression(italic("F")[ST])

### Extract middle Fst value for gradient argument
mid = max(fst.df$Fst) / 2

### Check the range of Fst values
summary(fst.df)

### Plot heatmap
fst.df %>%
ggplot(aes(x = Site1, y = Site2, fill = Fst))+
  geom_tile(colour = "black")+
  geom_text(aes(label = Fst), color="black", size = 3)+
  scale_fill_gradient2(low = "white", mid = "orange", high = "red", midpoint = mid, name = fst.label,limits=c(0,0.05))+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0), position = "right")+
  theme(axis.text = element_text(colour = "black", size = 10, face = "bold"),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 10)
  )
ggsave("FST_urchins.pdf", dpi=600, width=6, height=6)

########### VISUALIZE REGARDNG REGIONS

### Order by Longitude
geo <- read.table("../00-Data/02-filtering/01-21597snps-87ind/population_map_87individuals.txt", header=TRUE)
pop.geo.fst <- merge(fst.df, geo, by.x="Site1", by.y="CODE")

### Rearrange complete fst matrix
fst.complete <- melt(fst.matrix)
colnames(fst.complete) <- c("Site1","Site2","Fst")
fst.complete$Fst <- round(fst.complete$Fst,4)

# Plot heatmap
fst.complete %>%
  mutate(Site1 = fct_relevel(Site1,"SAB","BAT","PAR","MAL","BAL","LOM","LAN","LUW","SIK","BAN")) %>%
  mutate(Site2 = fct_relevel(Site2,"AMB", "BAN","SIK","LUW","LAN","LOM","BAL","MAL","PAR","BAT","SAB")) %>%
  ggplot(aes(x = Site1, y = Site2, fill = Fst))+
  geom_tile(colour = "black")+
  geom_text(aes(label = Fst), color="black", size = 3)+
  scale_fill_gradient2(low = "white", mid = "orange", high = "red", midpoint = mid, name = fst.label,limits=c(0,0.05))+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0), position = "right")+
  theme(axis.text = element_text(colour = "black", size = 10, face = "bold"),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 10)
  )
ggsave("FST_matrix.pdf", dpi=600, width=8, height=8)

### Perform bootstrap
tree <- nj(fst.matrix)
tree <- ape::as.phylo(stats::hclust(stats::dist(fst.matrix), method = "average"))
plot(tree, type="phylogram", show.tip=TRUE)

pdf("PlotBS.pdf")
bootstrap.value <- ape::boot.phylo(phy = tree, x = fst.matrix, FUN = function(xx) ape::as.phylo(stats::hclust(stats::dist(xx), method = "average")) , block = 1, B = 10000, trees = FALSE, rooted = TRUE) 
bootstrap.value <- round((bootstrap.value/10000)*100, 0)
tree$node.label <- bootstrap.value
plotBS(tree)
dev.off()

#### ANALYSIS P-VALUES ####

p.values <- read.table("02-fst/pvalues-indra.txt", header=TRUE, dec=".", stringsAsFactors = TRUE)
p.values.matrix <- as.matrix(p.values[,2:12])
row.names(p.values) <- p.values$p.values
heatmap(as.matrix(p.values[, -1]))

### Melt P-values
pvalues.df <- melt(p.values)

### Subset the non significant P-values
pvalues.non.significant <- subset(pvalues.df, value >0.05)

### Subset the non significant P-values
dim(pvalues.non.significant)

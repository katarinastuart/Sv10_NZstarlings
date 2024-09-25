#!/usr/bin/env Rscript
# Define input arguments --------------------------------------------#####
args <- commandArgs(trailingOnly = TRUE)
infileprefix = args[1]
outfile = args[2]
# infileprefix <- "data/vcftools_sum_stats/STACKS/populations.ALL.filt_bial.snps"
# outfile <- "results/ALL.filt_20bial.all.stats.pdf"
lqual <- paste(infileprefix, ".lqual", sep = "")
ldepthmean <- paste(infileprefix, ".ldepth.mean", sep = "")
lmiss <- paste(infileprefix, ".lmiss", sep = "")
frq <- paste(infileprefix, ".frq", sep = "")
idepth <- paste(infileprefix, ".idepth", sep = "")
imiss <- paste(infileprefix, ".imiss", sep = "")
het <- paste(infileprefix, ".het", sep = "")
# Following https://speciationgenomics.github.io/filtering_vcfs/ ----#####
# NOTE THAT THE VCF FILE IN THE LINK IS CREATED FROM VARIANT CALLING
# USING GATK AND NOT STACKS. STACKS DOES NOT PRODUCE VALUES FOR THE 
# QUAL COLUMN

# load tidyverse package
library(tidyverse)
# pdf("results/ALL.filt_20bial.all.stats.pdf")
pdf(outfile)
# Read in site quality scores ---------------------------------------#####
var_qual <- read_delim(lqual, delim = "\t",
                       col_names = c("chr", "pos", "qual"), skip = 1)
# Note that STACKS does not provide site quality scores and hence an 
# empty plot
if (all(var_qual$qual == -1)){
  print("No QUAL values")
} else {
  a <- ggplot(var_qual, aes(qual)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
  print(a + theme_light())
  
  print(a + theme_light() + xlim(c(0,50)) + ggtitle("QUAL per site"))
  hist(var_qual$qual, breaks = c(seq(0,200, by = 1), Inf), xlim = c(0, 200))
  
}

# Read in mean depths per site --------------------------------------#####
var_depth <- read_delim(ldepthmean, 
                        delim = "\t",
                        col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
a <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light() + ggtitle("depth per site")
# Summary of the mean depths per site
summary(var_depth$mean_depth)

q <- quantile(var_depth$mean_depth, probs = c(0, 0.01, 0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975, 0.99))
print(q)

a + theme_light() + xlim(0, 50)

# a + theme_light() + xlim(0, 25)
quantile(var_depth$mean_depth, probs = c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1))
hist(var_depth$mean_depth, 
     breaks = c(seq(0,25, by = 1), Inf), 
     xlim = c(0, 25),
     main = "Mean depth per site")

hist(var_depth$mean_depth, 
     breaks = c(seq(0,150, by = 1), Inf), 
     xlim = c(0, 150),
     main = "Mean depth per site")

# Read in site missingness ------------------------------------------#####
var_miss <- read_delim(lmiss, 
                       delim = "\t",
                       col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
a <- ggplot(var_miss, aes(fmiss)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light() + ggtitle("missingness per site")
a <- ggplot(var_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3, binwidth = 0.01, boundary = 0 )
a + theme_light() + ggtitle("missingness per site")
a + theme_light() + ggtitle("missingness per site") + xlim(c(0,0.5))

# a + theme_light() + xlim(c(0,1.0))
summary(var_miss$fmiss)

# Read in minor allele frequencies ----------------------------------#####
var_freq <- read_delim(frq, delim = "\t",
                       col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)

# find minor allele frequency for each site
var_freq$maf <- var_freq %>% select(a1, a2) %>% apply(1, function(z) min(z))

a <- ggplot(var_freq, aes(maf)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light() + ggtitle("MAF")

summary(var_freq$maf)

# Individual-based statistics ---------------------------------------#####
## Read in mean depth per individual --------------------------------#####
ind_depth <- read_delim(idepth, delim = "\t",
                        col_names = c("ind", "nsites", "depth"), skip = 1)

a <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light() + ggtitle("Mean depth per individual")
## Read in proporation of missing data per individual ---------------#####
ind_miss  <- read_delim(imiss, delim = "\t",
                        col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
a <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light() + ggtitle("Missingness per individual")
## Read in heterozygosity and inbreeding coefficient per individual -#####
ind_het <- read_delim(het, delim = "\t",
                      col_names = c("ind","ho", "he", "nsites", "f"), skip = 1)
a <- ggplot(ind_het, aes(f)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light() + ggtitle("Inbreeding coefficient per individual")
dev.off()

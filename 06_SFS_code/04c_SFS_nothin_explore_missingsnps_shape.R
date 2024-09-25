# This script makes plots SFS for publication for the ALL dataset
#
# This script is for the BCFtools dataset
#
# Libraries ---------------------------------------------------------
library(dartR)
library(tidyverse)
library(vcfR)
source("scripts/functions.R")

# Define input files ------------------------------------------------
vcffileALL <- "data/filtered_VCF/starling_noduprel_qual_miss_filt_bial.vcf.gz"
# vcffileALL is the path to the vcf file that has not been filtered for singletons and doubeltons. This path will have to be changed to match the file path 
metadtfile <- "data/Metadata_NZ_AU_UK_BE_ReplicatesSibRemoved2.csv"
# Create folders to store some data/outputs -------------------------
# dir.create("results/", recursive = T, showWarnings = F)
# Read in VCF file --------------------------------------------------
# glALL <- gl.read.vcf(vcffileALL)
vcf_in <- read.vcfR(vcffileALL, verbose = FALSE )
glALL <- vcfR2genlight(vcf_in)
vcfFIX <- getFIX(vcf_in)
vcfFIXmerged <- cbind(vcfFIX, glALL$other$loc.metrics)
vcfFIXmerged <- data.frame(vcfFIXmerged)
glALL$other$loc.metrics <- vcfFIXmerged
glALL@other$loc.metrics.flags$monomorphs <- F # Has to assign this to stop some errors downstream
# Read in metadata --------------------------------------------------
metadt <- read.csv(metadtfile)
# Subset metadt for the genlight object
# indiv_names <- gsub('.sorted.$', '', indNames(glALL))
indiv_names <- indNames(glALL)
metadtsubALL <- metadt %>% 
  dplyr::rename(ID = id) %>% 
  filter(ID %in% indiv_names) %>% 
  arrange(factor(ID, levels = indiv_names)) 
# Attach individual metadata to genlight object ---------------------
glALL@other$ind.metrics <- metadtsubALL
# Assign individual ID to individual metadata
glALL@other$ind.metrics$ID <- indNames(glALL)
# Attach population information to genlight object
pop(glALL) <- glALL$other$ind.metrics$pop2
# Order of SFS plots for main text ----------------------------------
# poporder <- c("IND: Other", "IND: Maharashtra subpopulation A",
#               "AUS: Melbourne", "Fiji", "NZ: Napier", "NZ: Other", 
#               "AUS: Sydney", "AUS: Gold Coast",
#               "Hawaii", "South Africa")
poporder <- c("ANT (BE)", "MKW (UK)", "NWC (UK)", "ORG (AU)", "MLV (AU)", 
              "CAN (NZ)", "BLN (NZ)", "UHT (NZ)", "PLM (NZ)", "AUK (NZ)")
# poporder <- unique(glALL$other$ind.metrics$pop2)
# quick summary
glALL$other$ind.metrics %>%
  dplyr::count(pop2)

glnomissing <- gl.filter.callrate(glALL, method = "loc", threshold = 1)
lockeep <- setdiff(locNames(glALL), locNames(glnomissing))
glonlymissingsnps <- gl.keep.loc(glALL, lockeep)
sfsdt <- plot_dt_SFS_compare_same_n_conf(glonlymissingsnps,
                                        popname_ls = poporder,
                                        iterations = 100,
                                        indmin = 10)

# sfsdt$median_quantile_plot_nLoc
save(sfsdt, file = "data/SFS_nothin_100_iterations_only_missingsnps_explore.Rdata")
# load("data/SFS_nothin_100_iterations_only_missingsnps_explore.Rdata")

# Define colours ----------------------------------------------------
lab_order <- poporder
lab_order_renamed <- poporder

dtcolshapes <- data.frame(pop = lab_order,
                          pop_lab = lab_order_renamed,
                          col = c("#a6cee3", "#b2df8a", "#cab2d6", 
                                  "#ff7f00", "#33a02c", "#1f78b4",
                                  "#e31a1c", "#fdbf6f", "#fb9a99",
                                  "#6a3d9a"),
                          shape = c(16, 3, 12, 18, 8, 14, 11, 17,15,7))
## Figure 6 ---------------------------------------------------------
pdf("results/SFS_ALL_nothin_only_missingsnps.pdf", width = 8.3, height = 3.5)
sfsdt$median_quantile_plot_nLoc + 
  facet_grid(.~popn) +
  scale_y_continuous(expand = expansion(mult = c(0, .02))) +
  theme(axis.text.x = element_text(angle = 90, colour = "black", size = 8, vjust = 0.5),
        axis.title.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(size = 10, face = "bold"),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.justification = "center",
        panel.grid = element_blank(),
        panel.spacing = unit(0.2, "lines"),
        strip.text.x = element_blank()) +
  scale_fill_manual(values = dtcolshapes$col,
                    breaks = dtcolshapes$pop,
                    labels = dtcolshapes$pop_lab) +
  scale_x_continuous(breaks = c(0, 25, 50),
                     labels = c(0, 25, 50)) +
  xlab("Minor allele frequency (%)") +
  ylab("Median number of SNPs") +
  ggtitle("SFS plot, only SNPs with missingness per pop (explore biases), subset n = 10, 100 iterations")
dev.off()


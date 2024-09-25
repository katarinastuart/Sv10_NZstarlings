# load libraries
library(radiator)
library(StAMPP)
library(tibble)
library(stringr)
library(ggrepel)

# set working directory
setwd("C:/Users/htan626/OneDrive - The University of Auckland/Documents/Starlings")

#### calculate genetic distance ####
## import genetic data
# original vcf file is in /nesi/nobackup/uoa02613/kstuart_projects/Sv10_NZstarlings/data/processing_rawdata/BCFtools/final_filtering/starling_noduprel_qual_miss_filt_NZ.recode.vcf
# convert vcf > plink > bed
# downloaded bed files
svulgaris <- genomic_converter(data = "./data/starling_noduprel_qual_miss_filt_NZ_removeAKL.recode.bed",
                           strata = "./data/starling_noduprel_qual_miss_filt_NZ_removeAKL.recode-pop.txt",
                           output = c("hierfstat", "genind", "genepop", "genlight"),
                           filter.common.markers = FALSE)

pop_file <- read.delim("./data/starling_noduprel_qual_miss_filt_NZ_removeAKL.recode-pop.txt", header = TRUE)

# query genlight
svulgaris$genlight 

# update @pop slot
adegenet::pop(svulgaris$genlight) <- pop_file$STRATA 

# calculate genetic distance
dist.sample.pair <- stamppNeisD(svulgaris$genlight, pop = FALSE)
dist.pop.pair <- stamppNeisD(svulgaris$genlight, pop = TRUE)

# formatting matrix for output - sample
#dist.sample.pair.df <- as.data.frame(dist.sample.pair)
#list <- c("Sample", row.names(dist.sample.pair.df))
#dist.sample.pair.df <- tibble::rownames_to_column(dist.sample.pair.df, "Population")
#colnames(dist.sample.pair.df) <- list
#write.table(dist.sample.pair.df, "./outputs/starling_geneDist-sample.txt", sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)

# formatting matrix for output - population
#dist.pop.pair.df <- as.data.frame(dist.pop.pair)
#list <- c("Population", row.names(dist.pop.pair.df))
#dist.pop.pair.df <- tibble::rownames_to_column(dist.pop.pair.df, "Population")
#colnames(dist.pop.pair.df) <- list
#write.table(dist.pop.pair.df, "./outputs/starling_geneDist-population.txt", sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)


#### calculate geographical distance ####
library(raster)
library(geodist)
library(geodata)
library(terra)

# import coordinates
sample.coord <-read.delim("./data/starling_noduprel_qual_miss_filt_NZ_removeAKL.recode-pop-wCoord.txt", header = TRUE)
colnames(sample.coord) <- c("individual", "strata", "longitude", "latitude")
head(sample.coord)

# calculate geographical distance
geodist.svulgaris <- geodist(sample.coord[,3:4], measure="geodesic")

# formatting matrix for output
#geodist.svulgaris.df <- as.data.frame(geodist.svulgaris)
#geodist.svulgaris.df2 <- cbind(sample.coord$individual, geodist.svulgaris.df)
#colnames(geodist.svulgaris.df2) <- c("sample", sample.coord$individual)
#write.table(geodist.svulgaris.df2, "./outputs/starling_geoDist-sample.txt", sep = '\t', col.names = TRUE, row.names =FALSE, quote = FALSE)

#### calculate environmental distance ####

## get climatic data
# Set boundaries for New Zealand
climdata <- getData('worldclim',download=TRUE,var='bio',res=.5, lon=c(160, 180), lat=c(-50, -30)) # resolution=30s
#climdata <- getData('worldclim',download=TRUE,var='bio',res=5) # resolution=5 minutes
#climdataNZ <- worldclim_country("New Zealand", path = tempdir(), var='bio',res=5) # did not get it to work

# extract climatic data for points
projection <- CRS("+proj=longlat +datum=WGS84 +no_defs")
points_D3 <- SpatialPoints(sample.coord[,3:4], proj4string=projection) # coords need to be in log, lat format
values_D3 <- extract(climdata,points_D3)
clim.points <- cbind.data.frame(sample.coord,values_D3)
head(clim.points)
#write.table(clim.points, "./outputs/starling_climate_metadata-res30s.txt", sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)

# environmental distance
library(rdist)
env.all <- clim.points[,5:23]
envdist.svulgaris <- pdist(env.all)

# formatting matrix for output
#envdist.svulgaris.df <- as.data.frame(envdist.svulgaris)
#envdist.svulgaris.df2 <- cbind(sample.coord$individual, envdist.svulgaris.df)
#colnames(envdist.svulgaris.df2) <- c("sample", sample.coord$individual)
#write.table(envdist.svulgaris.df2, "./outputs/starling_envDist_res30s-sample.txt", sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)

#### save distance matrices to file ####
save(dist.sample.pair, geodist.svulgaris, envdist.svulgaris, file = "./outputs/starling_mmrr_allDistMatrices_removeAKL-res30s.RData")

#### run MMRR ####

# MMRR performs Multiple Matrix Regression with Randomization analysis
# Y is a dependent distance matrix
# X is a list of independent distance matrices (with optional names)

MMRR<-function(Y,X,nperm=999){
  #compute regression coefficients and test statistics
  nrowsY<-nrow(Y)
  y<-unfold(Y)
  if(is.null(names(X)))names(X)<-paste("X",1:length(X),sep="")
  Xmats<-sapply(X,unfold)
  fit<-lm(y~Xmats)
  coeffs<-fit$coefficients
  summ<-summary(fit)
  r.squared<-summ$r.squared
  tstat<-summ$coefficients[,"t value"]
  Fstat<-summ$fstatistic[1]
  tprob<-rep(1,length(tstat))
  Fprob<-1
  
  #perform permutations
  for(i in 1:nperm){
    rand<-sample(1:nrowsY)
    Yperm<-Y[rand,rand]
    yperm<-unfold(Yperm)
    fit<-lm(yperm~Xmats)
    summ<-summary(fit)
    Fprob<-Fprob+as.numeric(summ$fstatistic[1]>=Fstat)
    tprob<-tprob+as.numeric(abs(summ$coefficients[,"t value"])>=abs(tstat))
  }
  
  #return values
  tp<-tprob/(nperm+1)
  Fp<-Fprob/(nperm+1)
  names(r.squared)<-"r.squared"
  names(coeffs)<-c("Intercept",names(X))
  names(tstat)<-paste(c("Intercept",names(X)),"(t)",sep="")
  names(tp)<-paste(c("Intercept",names(X)),"(p)",sep="")
  names(Fstat)<-"F-statistic"
  names(Fp)<-"F p-value"
  return(list(r.squared=r.squared,
              coefficients=coeffs,
              tstatistic=tstat,
              tpvalue=tp,
              Fstatistic=Fstat,
              Fpvalue=Fp))
}

# unfold converts the lower diagonal elements of a matrix into a vector
# unfold is called by MMRR

unfold<-function(X){
  x<-vector()
  for(i in 2:nrow(X)) x<-c(x,X[i,1:i-1])
  x<-scale(x, center=TRUE, scale=TRUE)  # Comment this line out if you wish to perform the analysis without standardizing the distance matrices! 
  return(x)
}

# Tutorial for data files gendist.txt, geodist.txt, and ecodist.txt

# Read the matrices from files.
# The read.matrix function requires {tseries} package to be installed and loaded.
# If the files have a row as a header (e.g. column names), then specify 'header=TRUE', default is 'header=FALSE'.
library(tseries)

# set working directory
setwd("C:/Users/htan626/OneDrive - The University of Auckland/Documents/Starlings/outputs")

# read in distance matrices
load("./starling_mmrr_allDistMatrices_removeAKL-res30s.RData")

geneMat <- dist.sample.pair
geoMat <- geodist.svulgaris
envMat <- envdist.svulgaris

# Make a list of the explanatory (X) matrices.
# Names are optional.  Order doesn't matter.
# Can include more than two matrices, if desired.
Xmats <- list(geographical=geoMat,environmental=envMat)
XmatsGeo <- list(geographical=geoMat)
XmatsEnv <- list(environmental=envMat)

# Run MMRR function using genMat as the response variable and Xmats as the explanatory variables.
# nperm does not need to be specified, default is nperm=999)
results <- MMRR(geneMat,Xmats,nperm=9999)
save(results, file = "starling_mmrr_results_nperm9999_removeAKL-res30s.RData")
load("./starling_mmrr_results_nperm9999_removeAKL-res30s.RData")
# view results
results

#### visualisation ####
library(ggplot2)
library(ggpubr)

setwd("C:/Users/htan626/OneDrive - The University of Auckland/Documents/Starlings/outputs")

## load matrices and results (if not yet loaded above)
# read in distance matrices
#load("./starling_mmrr_allDistMatrices.RData")
#geneMat <- dist.sample.pair
#geoMat <- geodist.svulgaris
#envMat <- envdist.svulgaris
# read in results
#load("./starling_mmrr_results_nperm9999.RData")

# see Fig. 2 of Wang 2013 (paper describing mmrr)

## merge matrices into dataframe, with corresponding pairwise comparisons on the same row
## does not include diagonal of the matrices (pairwise comparison of same sample)
# genetic
geneMat1 <- geneMat
geneMat1[upper.tri(geneMat1, diag = TRUE)] <- NA
geneMat2 <- data.frame(
  row = c(t(row(geneMat1))),
  col = c(t(col(geneMat1))),
  geneticDist = c(t(geneMat1))
)

# geographical
geoMat1 <- geoMat
geoMat1[upper.tri(geoMat1, diag = TRUE)] <- NA
geoMat2 <- data.frame(
  row = c(t(row(geoMat1))),
  col = c(t(col(geoMat1))),
  geographicDist = c(t(geoMat1))
)

# environmental
envMat1 <- envMat
envMat1[upper.tri(envMat1, diag = TRUE)] <- NA
envMat2 <- data.frame(
  row = c(t(row(envMat1))),
  col = c(t(col(envMat1))),
  environmentalDist = c(t(envMat1))
)

# merge
df <- cbind(geneMat2, geoMat2, envMat2)
df2 <- na.omit(df)
df3 <- df2[,-c(4,5,7,8)]
df3$combinedEffect <- (0.123*df3$geographicDist)+(0.203*df3$environmentalDist) # update coefficients from results

## get identity of pairwise comparisons for potential labeling and grouping purposes
# get row names
rname <- as.data.frame(rownames(geneMat))
df3$row2 <- NA
df3$col2 <- NA
# get identity of first sample in pair
for (i in 1:nrow(df3)){
  df3[i,7] <- rname[df3[i,1],1]
}
# get identity of second sample in pair
for (i in 1:nrow(df3)){
  df3[i,8] <- rname[df3[i,2],1]
}
# remove "-R.sorted.bam"
df3$row2 <- gsub("-R.sorted.bam","", as.character(df3$row2))
df3$col2 <- gsub("-R.sorted.bam","", as.character(df3$col2))
df3$row2 <- gsub(".sorted.bam","", as.character(df3$row2))
df3$col2 <- gsub(".sorted.bam","", as.character(df3$col2))
# get unique pairwise name
df3$pair <- str_c(df3$row2,'_',df3$col2)
# add column labelling samples across certain threshold for a certain factor
df4 <- df3 %>%
  mutate(label = ifelse(geographicDist < 400000 & environmentalDist > 600,
                           df3$pair, ""))

## plots
# genetic distance vs geographic distance
plotGenGeo <- ggplot(df4, aes(x = geographicDist, y = geneticDist)) +
  geom_point() +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  #geom_text_repel(aes(label = label), box.padding = 2, size = 2) +
  labs(x = "geographic distance", y = "genetic distance")

# genetic distance vs environmental distance
plotGenEnvt <- ggplot(df4, aes(x = environmentalDist, y = geneticDist)) +
  geom_point() +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  #geom_text_repel(aes(label = label), box.padding = 2, size = 2) +
  labs(x = "environmental distance", y = "genetic distance")

# genetic distance vs combined effects
plotGenComb <- ggplot(df4, aes(x = combinedEffect, y = geneticDist)) +
  geom_point() +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  #geom_text_repel(aes(label = label), box.padding = 2, size = 2) +
  labs(x = "0.123(geo. dist.) + 0.203(envt. dist.)", y = "genetic distance") # update coefficients from results

# geographic distance vs environmental distance
plotEnvtGeo <- ggplot(df4, aes(x = geographicDist, y = environmentalDist)) +
  geom_point() +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  #geom_text_repel(aes(label = label), box.padding = 2, max.overlaps = 500, size = 2) +
  labs(x = "geographic distance", y = "environmental distance")

# combine plots
png("mmrr-removeAKL-res30s.png", width = 8, height = 8, units = "in", res = 300)
plotAll <- ggarrange(plotGenGeo, plotGenEnvt, plotGenComb, plotEnvtGeo,
          nrow = 2, ncol = 2)
annotate_figure(plotAll, top = text_grob("Pairwise distance amongst samples excluding AKL (n=58)"))
dev.off()

#### visualisation - TEST - investigating select data points ####
df5 <- df4
df5[c("row-pop", "row-indiv")] <- str_split_fixed(df5$row2, "-" , 2)
df5[c("col-pop", "col-indiv")] <- str_split_fixed(df5$col2, "-" , 2)
df5$rowPopNew <- ifelse(df5$`row-pop` == "BAD" | df5$`row-pop` == "BRA", "BLN", df5$`row-pop`)
df5$colPopNew <- ifelse(df5$`col-pop` == "BAD" | df5$`col-pop` == "BRA", "BLN", df5$`col-pop`)
df5$pair2 <- paste0(df5$rowPopNew,"-",df5$colPopNew)

# colour grouping 1 (highlighting any pairwise involving Auckland)
df5$colour1 <- ifelse(df5$rowPopNew == df5$colPopNew, "within pop", "different pop")
df5$colour1 <- ifelse(grepl("AUK", df5$pair2), "Auckland", df5$colour1)

# colour group 2 (highlighting any pairwise involving Upper Hutt)
df5$colour2 <- ifelse(df5$rowPopNew == df5$colPopNew, "within pop", "different pop")
df5$colour2 <- ifelse(grepl("UHT", df5$pair2), "Upper Hutt", df5$colour2)

## plots
# genetic distance vs geographic distance
plotGenGeo <- ggplot(df5, aes(x = geographicDist, y = geneticDist, colour = colour1)) +
  geom_point(alpha=0.3) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_text_repel(aes(label = label), box.padding = 2, size = 2) +
  labs(x = "geographic distance", y = "genetic distance")

# genetic distance vs environmental distance
plotGenEnvt <- ggplot(df5, aes(x = environmentalDist, y = geneticDist, colour = colour1)) +
  geom_point(alpha=0.3) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_text_repel(aes(label = label), box.padding = 2, size = 2) +
  labs(x = "environmental distance", y = "genetic distance")

# genetic distance vs combined effects
plotGenComb <- ggplot(df5, aes(x = combinedEffect, y = geneticDist, colour = colour1)) +
  geom_point(alpha=0.3) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_text_repel(aes(label = label), box.padding = 2, size = 2) +
  labs(x = "0.11354(geo. dist.) + 0.20274(envt. dist.)", y = "genetic distance") # update coefficients from results

# geographic distance vs environmental distance
plotEnvtGeo <- ggplot(df5, aes(x = geographicDist, y = environmentalDist, colour=colour2)) +
  geom_point(alpha=0.3) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  #geom_text_repel(aes(label = label), box.padding = 2, max.overlaps = 500, size = 2) +
  labs(x = "geographic distance", y = "environmental distance")

# combine plots
png("mmrr-allSamples-coloured.png", width = 10, height = 8, units = "in", res = 300)
plotAll <- ggarrange(plotGenGeo, plotGenEnvt, plotGenComb, plotEnvtGeo,
                     nrow = 2, ncol = 2)
annotate_figure(plotAll, top = text_grob("Pairwise distance amongst samples excluding AKL (n=58)"))
dev.off()

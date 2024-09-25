# This script contains utility functions that are reused
library(ggplot2)
### Define function for quick plotting of individual vs indiv. call rates
gl_indvsCR <- function(gl){
  require(ggplot2)
  # Add population info to individual metrics table
  gl$other$ind.metrics$pop_name <- pop(gl)
  # 1. Convert to matrix; individuals as rows and loci as columns
  mat_all <- as.matrix(gl)
  # 2. Calculate the number of "NA" per row, referring to the number of
  # missing loci per individual
  dt_ind_all <- data.frame(ind = rownames(mat_all), NA_loci = apply(is.na(mat_all), 1, sum))
  # 3. Calculate call rate from the number of NA per row and the total 
  # number of loci
  dt_ind_all$call_rate <- 1 - dt_ind_all$NA_loci/nLoc(gl)
  # 4. Merge the summarised table of individuals ID and call rates to the
  # table of individual metrics to link individuals to populations
  dt_ind_all <- merge(dt_ind_all, gl$other$ind.metrics,
                      by.x = 'ind', by.y = 'DART_ID')
  
  # 5. Make plot of individuals vs call rates, coloured by the populations
  p_ind_callrate_all <- ggplot(data = dt_ind_all, aes(x = ind, y = call_rate, colour = COUNTRY.OF.ORIGIN, label = pop_name)) + 
    geom_point() +
    theme_bw() +
    theme(legend.position = c(0.95,0.05), legend.justification = c(1,0)) +
    labs(colour = "Populations") +
    ggtitle("Indiv. vs call rate. check table/interactive plot for indiv name.")
  return(list(dt_indvsCR = dt_ind_all, plot_indvsCR = p_ind_callrate_all))
}


### Define function for quick filter to explore filters

gl_filt_CR_rep_MAC <- function(gl, CR = 1, rep = 1, MAC_threshold = 1, SNP_in_single_samp_rm = F){
  # NOTE that this is Minor Allele Count filter over Minor Allele Frequency
  # filter. This is because missing data might mess up with MAF, but not MAC.
  #
  # The way this function calculates MAC also allows us to filter out loci with
  # one lone sample with a "single homozygous alternate or reference" (I.e. All 
  # samples are homozygous reference except for one sample) more easily.
  # How this is done is commented out at the end of this function. 
  glsub <- gl.filter.callrate(gl, threshold = CR, mono.rm = T)
  glsub <- gl.filter.reproducibility(glsub, threshold = rep)
  glsub <- gl.filter.monomorphs(glsub)
  # 1. Convert to matrix; individuals as rows and loci as columns
  matb <- as.matrix(glsub)
  # 2. Calculate the number of 
  dtb <- data.frame(loci = colnames(matb), loci_MAC_1 = apply(matb == 1, 2, sum, na.rm = T), loci_homo_alt = apply(matb == 2, 2, sum, na.rm = T), loci_homo_ref = apply(matb == 0, 2, sum, na.rm = T))
  dtb$Alt_Alnct <- 1*dtb$loci_MAC_1 + 2*dtb$loci_homo_alt
  dtb$Ref_Alnct <- 1*dtb$loci_MAC_1 + 2*dtb$loci_homo_ref
  dtb$MAC <- dtb$Alt_Alnct
  dtb$MAC[dtb$Alt_Alnct > dtb$Ref_Alnct] <- dtb$Ref_Alnct[dtb$Alt_Alnct > dtb$Ref_Alnct]
  dtbfilt <- subset(dtb, MAC <= MAC_threshold)
  glsub <- gl.drop.loc(glsub, loc.list = c(dtbfilt$loci))
  if (SNP_in_single_samp_rm == T){
    dtbfilt1 <- subset(dtb, MAC <= 1)
    dtbfilt2 <- subset(dtb, loci_homo_alt == 1 & MAC == 2)
    dtbfilt3 <- subset(dtb, loci_homo_ref == 1 & MAC == 2)
    glsub <- gl.drop.loc(glsub, loc.list = c(dtbfilt1$loci))
    glsub <- gl.drop.loc(glsub, loc.list = c(dtbfilt2$loci))
    glsub <- gl.drop.loc(glsub, loc.list = c(dtbfilt3$loci))
  } 
  # singletons
  # dtb1 <- subset(dtb, loci_MAC_1 == 1 & loci_homo_alt == 0)
  # dtb2 <- subset(dtb, loci_MAC_1 == 1 & loci_homo_ref == 0)
  # single homozygous calls
  # dtb3 <- subset(dtb, loci_MAC_1 == 0 & loci_homo_alt == 1)
  # dtb4 <- subset(dtb, loci_MAC_1 == 0 & loci_homo_ref == 1)
  # Subset to remove singletons and single homozygous calls
  # glsub <- gl.drop.loc(glsub, loc.list = c(dtb1$loci, dtb2$loci, dtb3$loci, dtb4$loci))
  return(glsub)
}



gl.oneSNPpermarker <- function(gl){
  require(dplyr)
  # This function samples one SNP per locus (marker)
  ## 1. Count number of SNPs per marker for subsampling
  dt_SNPperloci <- gl$other$loc.metrics %>% 
    group_by(CloneID) %>%
    dplyr::summarise(nSNP = length(AlleleID))
  ## dt_SNPperloci <- ddply(gl$other$loc.metrics, 
  ##                       .(CloneID), summarise, 
  ##                       nSNP = length(AlleleID))
  ## 1.b. Select markers with multiple SNPs
  dt_SNPperloci1a <- subset(dt_SNPperloci, nSNP > 1)
  ## 1.c. Select markers with only one SNP
  dt_SNPperloci1b <- subset(dt_SNPperloci, nSNP == 1)
  ## 2. Sample one SNP per marker for markers with multiple SNPs
  dt_SNPperloci1c <- gl$other$loc.metrics %>%
    filter(CloneID %in% dt_SNPperloci1a$CloneID) %>%
    group_by(CloneID) %>% 
    dplyr::summarise(AlleleID = sample(AlleleID, 1))
  ## dt_SNPperloci1c <- ddply(subset(gl$other$loc.metrics, CloneID %in% dt_SNPperloci1a$CloneID), 
  ##                         .(CloneID), summarise, 
  ##                         AlleleID = sample(AlleleID, 1))
  ## 3. Create a table with all the subsampled loci. 1 SNP per loci
  dt_1SNPperloci <- rbind(subset(gl$other$loc.metrics, CloneID %in% dt_SNPperloci1b$CloneID), subset(gl$other$loc.metrics, AlleleID %in% dt_SNPperloci1c$AlleleID))
  ## 4. Make a list of loci to keep based on locNames function
  loc2keep <- locNames(gl)[gl$other$loc.metrics$AlleleID %in% dt_1SNPperloci$AlleleID]
  ## 5. Retain the targeted SNPs
  glsub <- gl.keep.loc(gl, loc.list = loc2keep)
  return(glsub)
}



sort_sample <- function(x, thresholdbp = 50000){
  # This function is used in gl.prune.SNP.dist to return SNP
  # position every thresholdbp (default to 50kb)
  threshold <- thresholdbp 
  x <- sort(x)
  for(i in 1:length(x)){
    if(i == 1){
      ls <- x[i]
    } else {
      if ((x[i] - ls[length(ls)]) > threshold){
        ls <- c(ls, x[i])
      }
    }
  }
  return(ls)
}

gl.prune.SNP.dist <- function(gl, bp_threshold = 100000){
  require(dplyr)
  # This function retains 1 SNP every bp_threshold (default = 10kb)
  # Returns a
  dt <- gl$other$loc.metrics
  dt$ChromPos_plus_SnpPos <- dt$ChromPos_Indian_Myna_v2.1 + dt$SnpPosition
  dt_LDpass1 <- dt %>%
    group_by(Chrom_Indian_Myna_v2.1) %>%
    dplyr::summarise(ChromPos_plus_SnpPos = sort_sample(ChromPos_plus_SnpPos, thresholdbp = bp_threshold))
  # dt_LDpass1 <- ddply(dt,
  #                    .(Chrom_Indian_Myna_v2.1), plyr::here(summarise),
  #                    ChromPos_plus_SnpPos = sort_sample(ChromPos_plus_SnpPos, thresholdbp = bp_threshold))
  dt_LDpass <- merge(dt_LDpass1, dt)
  
  ## 4. Make a list of loci to keep based on locNames function
  LDfiltloc2keep <- locNames(gl)[gl$other$loc.metrics$AlleleID %in% dt_LDpass$AlleleID]
  ## 5. Retain the targeted SNPs
  glsub <- gl.keep.loc(gl, loc.list = LDfiltloc2keep)
  return(glsub)
}

gl.prune.SNP.dist.FROM.VCF <- function(gl, bp_threshold = 100000){
  require(dplyr)
  # This function retains 1 SNP every bp_threshold (default = 10kb)
  # Returns a
  dt <- gl$other$loc.metrics
  dt$POS <- as.numeric(dt$POS)
  dt_LDpass1 <- dt %>%
    group_by(CHROM) %>%
    dplyr::summarise(POS = sort_sample(POS, thresholdbp = bp_threshold))
  # dt_LDpass1 <- ddply(dt,
  #                    .(Chrom_Indian_Myna_v2.1), plyr::here(summarise),
  #                    ChromPos_plus_SnpPos = sort_sample(ChromPos_plus_SnpPos, thresholdbp = bp_threshold))
  dt_LDpass <- merge(dt_LDpass1, dt)
  
  ## 4. Make a list of loci to keep based on locNames function
  LDfiltloc2keep <- locNames(gl)[gl$other$loc.metrics$ID %in% dt_LDpass$ID]
  ## 5. Retain the targeted SNPs
  glsub <- gl.keep.loc(gl, loc.list = LDfiltloc2keep)
  return(glsub)
}

# PCA ---------------------------------------------------------------#####
PCA_subset_summary <- function(gl, ind_no, output.pdf, MAF_threshold = 0, reclasspop = NULL){
  # This function performs PCA on subsetted genlight dataset and plots it
  # to a PDF, and allows for MAF filters
  require(dartR)
  require(plyr)
  gl$other$ind.metrics$pop_subsampled <- pop(gl)
  dt_sub <- ddply(gl$other$ind.metrics, .(pop_subsampled), plyr::here(summarise), subsampled_ID = sample(ID, ind_no))
  gl_sub <- gl.keep.ind(gl, ind.list = dt_sub$subsampled_ID, recalc = T)
  if (MAF_threshold != 0){
    gl_sub <- gl.filter.maf(gl_sub, threshold = MAF_threshold)
  }
  if (!is.null(reclasspop)){
    pop(gl_sub) <- gl_sub$other$ind.metrics[[reclasspop]]
  }
  pdf(output.pdf)
  pc_sub <- gl.pcoa(gl_sub, nfactors=10, parallel = T, n.cores = 8)
  gl.pcoa.plot(pc_sub, gl_sub, labels="pop", xaxis=1, yaxis=2)
  gl.pcoa.plot(pc_sub, gl_sub, labels="pop", xaxis=2, yaxis=3)
  gl.pcoa.plot(pc_sub, gl_sub, labels="pop", xaxis=3, yaxis=4)
  dev.off()
}

PCA_screeplot <- function(pc_obj, nPC = 15){
  require(ggplot2)
  screeplotdt <- data.frame(pc = seq(1, length(pc_obj$eig)), 
                            variance_explained = (pc_obj$eig/sum(pc_obj$eig))*100)
  p2 <- ggplot(data = subset(screeplotdt, pc <= nPC), 
               aes(x = pc, y = variance_explained)) +
    geom_bar(stat = "identity", colour = "black") +
    theme_classic() + 
    theme(legend.position = "none",
          # axis.title = element_blank(),
          # axis.text = element_blank(),
          # axis.ticks = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=1)) +
    ylab("Variance explained (%)") +
    xlab("PC") 
  # scale_y_log10() +
  # annotate(geom = "text", x = 8, y = 6, label = "Eigenvalues")
  return(p2)
}

PCA_var_explained <- function(pc_obj){
  screeplotdt <- data.frame(pc = seq(1, length(pc_obj$eig)), 
                            variance_explained = (pc_obj$eig/sum(pc_obj$eig))*100)
  return(screeplotdt)
}

plotPCASI <- function(dtpc, varls, dtcolshapes, xaxis = 1, yaxis = 2, popcol = "pop"){
  # dtpc: dataframe with pc and population definition
  # varls: vector of variance explained for each axis (output from PCA_var_explained)
  # xaxis: pc axis for x axis
  # yaxis: 
  xcolname <- paste("PC", xaxis, sep = "")
  ycolname <- paste("PC", yaxis, sep = "")
  dtpc$xaxis <- dtpc[[xcolname]]
  dtpc$yaxis <- dtpc[[ycolname]]
  dtpc$pop <- dtpc[[popcol]]
  xlabel <- paste("PCA ", xaxis, " (", varls[xaxis], "%)", sep = "")
  ylabel <- paste("PCA ", yaxis, " (", varls[yaxis], "%)", sep = "")
  p <- ggplot() + geom_point(data = dtpc, 
                             aes(x = xaxis, y = yaxis, 
                                 colour = pop, 
                                 shape = pop)) +  
    scale_colour_manual(breaks = levels(dtpc$pop),
                        labels = dtcolshapes$pop_lab[match(levels(dtpc$pop), dtcolshapes$pop)], 
                        values = dtcolshapes$col[match(levels(dtpc$pop), dtcolshapes$pop)]) + 
    scale_shape_manual(breaks = levels(dtpc$pop),
                       labels = dtcolshapes$pop_lab[match(levels(dtpc$pop), dtcolshapes$pop)], 
                       values = dtcolshapes$shape[match(levels(dtpc$pop), dtcolshapes$pop)]) +
    xlab(xlabel) + ylab (ylabel) 
  return(p)
}
PCA_theme <- theme_classic() + theme(axis.title = element_text(size = 12, face = "bold"),
                                     axis.text = element_text(size = 10, face = "bold"),
                                     legend.title = element_blank(),
                                     legend.key.size = unit(0.75, 'lines'),
                                     legend.text = element_text(size = 8, face = "bold"),
                                     legend.background = element_rect(color = "black", size = 0.01, linetype = "solid"),
                                     legend.position = "bottom",
                                     aspect.ratio=1)


# PCAdapt -----------------------------------------------------------#####
pcadapt_plots <- function(gl, path2filesuffix, outputfile.pdf, MAF_threshold = 0.05, K_no = 40){
  require(LEA)
  require(dartR)
  require(plyr)
  genofile <- paste(path2filesuffix, ".geno", sep = "")
  lfmmfile <- paste(path2filesuffix, ".lfmm", sep = "")
  pdf(outputfile.pdf)
  mat <- t(as.matrix(gl))
  write.table(mat, file = genofile,
              quote = F, sep = "", na = "9", row.names = F, 
              col.names = F)
  geno2lfmm(genofile, lfmmfile)
  pcadapt_all <- read.pcadapt(lfmmfile, type = "lfmm")
  x1 <- pcadapt(input = pcadapt_all, K = K_no, min.maf = MAF_threshold)
  plot(x1, option = "screeplot")
  poplist.names <- gl@pop
  plot(x1, option = "scores", pop = poplist.names)
  plot(x1 , option = "manhattan")
  plot(x1, option = "qqplot") 
  hist(x1$pvalues, xlab = "p-values", main = "min maf = 0.05", breaks = 50, col = "orange") 
  plot(x1, option = "stat.distribution")
  dev.off()
  return(x1)
}

pcadapt_plots_lfmm <- function(gl, lfmmfile, outputfile.pdf, MAF_threshold = 0.05, K_no = 40){
  require(LEA)
  require(dartR)
  require(plyr)
  pdf(outputfile.pdf)
  pcadapt_all <- read.pcadapt(lfmmfile, type = "lfmm")
  x1 <- pcadapt(input = pcadapt_all, K = K_no, min.maf = MAF_threshold)
  plot(x1, option = "screeplot")
  poplist.names <- gl@pop
  plot(x1, option = "scores", pop = poplist.names)
  plot(x1 , option = "manhattan")
  plot(x1, option = "qqplot") 
  hist(x1$pvalues, xlab = "p-values", breaks = 50, col = "orange", main = paste("min maf = ", MAF_threshold, sep = "")) 
  plot(x1, option = "stat.distribution")
  dev.off()
  return(x1)
}


pcadapt_plot_subset <- function(gl, ind_no, path2filesuffix, outputfile.pdf, iteration_no, MAF_threshold = 0.05, K_no = 40, plot_type){
  require(LEA)
  require(dartR)
  require(plyr)
  require(pcadapt)
  pdf(outputfile.pdf)
  for (i in seq(1, iteration_no)){
    genofile <- paste(path2filesuffix, "i", i, ".geno", sep = "")
    lfmmfile <- paste(path2filesuffix, "i", i, ".lfmm", sep = "")
    dt_sub <- ddply(gl$other$ind.metrics, .(pop), plyr::here(summarise), subsampled_ID = sample(ID, ind_no))
    gl_sub <- gl.keep.ind(gl, ind.list = dt_sub$subsampled_ID, recalc = T)
    gl_sub <- gl.filter.monomorphs(gl_sub)
    mat_sub <- t(as.matrix(gl_sub))
    write.table(mat_sub, file = genofile,
                quote = F, sep = "", na = "9", row.names = F, 
                col.names = F)
    geno2lfmm(genofile, lfmmfile)
    pcadapt_filt5 <- read.pcadapt(lfmmfile, type = "lfmm")
    x1 <- pcadapt(input = pcadapt_filt5, K = K_no, min.maf = MAF_threshold)
    if (plot_type == "PCA_and_screeplot"){
      p1 <- plot(x1, option = "screeplot")
      p1 + scale_y_log10()
      poplist.names <- gl_sub@pop
      plot(x1, option = "scores", pop = poplist.names)
    }
    if (plot_type == "screeplot"){
      p1 <- plot(x1, option = "screeplot")
      p1 + scale_y_log10()
    } 
    if (plot_type == "PCAscores"){
      poplist.names <- gl_sub@pop
      plot(x1, option = "scores", pop = poplist.names)
    } 
    if (plot_type == "manhattan"){
      plot(x1 , option = "manhattan")
    } 
    if (plot_type == "qqplot"){
      plot(x1, option = "qqplot") 
    } 
    if (plot_type == "pvalue_hist"){
      hist(x1$pvalues, xlab = "p-values", breaks = 50, col = "orange") 
    } 
    if (plot_type == "stat.distribution"){
      plot(x1, option = "stat.distribution")
    } 
  }
  dev.off()
}

pcadapt_plot_subset_lfmm <- function(gl, path2filesuffix, outputfile.pdf, iteration_no, MAF_threshold = 0.05, K_no = 40, plot_type){
  require(LEA)
  require(dartR)
  require(plyr)
  require(ggplot2)
  require(pcadapt)
  pdf(outputfile.pdf)
  for (i in seq(1, iteration_no)){
    lfmmfile <- paste(path2filesuffix, "i", i, ".lfmm", sep = "")
    print(lfmmfile)
    pcadapt_filt5 <- read.pcadapt(lfmmfile, type = "lfmm")
    x1 <- pcadapt(input = pcadapt_filt5, K = K_no, min.maf = MAF_threshold)
    if (plot_type == "PCA_and_screeplot"){
      p1 <- plot(x1, option = "screeplot")
      p1 + scale_y_log10()
      poplist.names <- gl@pop
      plot(x1, option = "scores", pop = poplist.names)
    }
    if (plot_type == "screeplot"){
      p1 <- plot(x1, option = "screeplot")
      p1 + scale_y_log10()
    } 
    if (plot_type == "manhattan"){
      plot(x1 , option = "manhattan")
    } 
    if (plot_type == "qqplot"){
      plot(x1, option = "qqplot") 
    } 
    if (plot_type == "pvalue_hist"){
      fdradj <- p.adjust(x1$pvalues, method = "BH")
      outliers_a <- which(fdradj < 0.1)
      outliers_b <- which(fdradj < 0.05)
      dt0 <- data.frame(p_values = x1$pvalues)
      dt0$outliers <- "Not outliers"
      dt0$outliers[outliers_a] <- "FDR < 0.1"
      dt0$outliers[outliers_b] <- "FDR < 0.05"
      p1 <- ggplot(dt0, aes(x = p_values, col = outliers, fill = outliers)) + geom_histogram()
      print(p1)
      # hist(x1$pvalues, xlab = "p-values", breaks = 50, col = "orange") 
    } 
    if (plot_type == "pvalue_uncorr_hist"){
      x1$pvalues_gif1 <- as.numeric(pchisq(x1$chi2.stat*x1$gif/1, df = Kval, lower.tail = FALSE))
      fdradj <- p.adjust(x1$pvalues_gif1, method = "BH")
      outliers_a <- which(fdradj < 0.1)
      outliers_b <- which(fdradj < 0.05)
      dt0 <- data.frame(p_values_uncorr = x1$pvalues_gif1)
      dt0$outliers <- "Not outliers"
      dt0$outliers[outliers_a] <- "FDR < 0.1"
      dt0$outliers[outliers_b] <- "FDR < 0.05"
      p1 <- ggplot(dt0, aes(x = p_values_uncorr, col = outliers, fill = outliers)) + geom_histogram()
      print(p1)
      # hist(x1$pvalues, xlab = "p-values", breaks = 50, col = "orange") 
    } 
    if (plot_type == "stat.distribution"){
      plot(x1, option = "stat.distribution")
    } 
  }
  dev.off()
}

gl_PCAdapt_iterations <- function(gl, n_per_pop, output, K_no, minmaf_threshold = 0.05, #FDR, 
                                  iterations, padj_method = "fdr", new_gif = NULL){
  require(plyr)
  require(dartR)
  require(pcadapt)
  require(LEA)
  mat1 <- matrix(nrow = nLoc(gl), ncol = iterations)
  mat2 <- matrix(nrow = nLoc(gl), ncol = iterations)
  mat3 <- matrix(nrow = nLoc(gl), ncol = iterations)
  mat4 <- matrix(nrow = nLoc(gl), ncol = iterations)
  matx <- matrix(nrow = nLoc(gl), ncol = iterations)
  row.names(mat1) <- locNames(gl)
  row.names(mat2) <- locNames(gl)
  row.names(mat3) <- locNames(gl)
  row.names(mat4) <- locNames(gl)
  row.names(matx) <- locNames(gl)
  if (is.null(new_gif)){
    new_gif = NULL
    mat5 <- NULL
    mat6 <- NULL
  } else {
    mat5 <- matrix(nrow = nLoc(gl), ncol = iterations)
    mat6 <- matrix(nrow = nLoc(gl), ncol = iterations)
    row.names(mat5) <- locNames(gl)
    row.names(mat6) <- locNames(gl)
  }
  v <- rep(NA, iterations)
  for (i in seq(1, iterations)){
    print(i)
    dt_subsample_pop <- ddply(gl$other$ind.metrics, .(pop), plyr::here(summarise), subsampled_ID = sample(ID, n_per_pop))
    glsub <- gl.keep.ind(gl, ind.list = dt_subsample_pop$subsampled_ID, recalc = T)
    glsub <- gl.filter.monomorphs(glsub)
    matsub <- t(as.matrix(glsub))
    genofile <- paste(output, "i", sprintf("%05d", i), ".geno", sep = "")
    lfmmfile <- paste(output, "i", sprintf("%05d", i), ".lfmm", sep = "")
    locnamefile <- paste(output, "i", sprintf("%05d", i), ".loci.list.csv", sep = "")
    write.table(matsub, file = genofile,
                quote = F, sep = "", na = "9", row.names = F, col.names = F)
    geno2lfmm(genofile, lfmmfile)
    pcadapt_gl <- read.pcadapt(lfmmfile, type = "lfmm")
    x1 <- pcadapt(input = pcadapt_gl, K = K_no, min.maf = minmaf_threshold)
    chi2_stat_gif1 <- x1$chi2.stat*x1$gif
    x1$pvalues_gif1 <- as.numeric(pchisq(x1$chi2.stat*x1$gif/1, df = Kval, lower.tail = FALSE))
    padj_uncorr_gif <- p.adjust(x1$pvalues_gif1,method=padj_method)
    padj <- p.adjust(x1$pvalues,method=padj_method)
    # alpha <- FDR
    # outliers <- which(padj < alpha)
    glsubloc <- locNames(glsub)
    write.table(glsubloc, file = locnamefile, row.names = FALSE, col.names = FALSE, sep=',')
    mat1[glsubloc,i] <- x1$pvalues
    mat2[glsubloc,i] <- x1$pvalues_gif1
    mat3[glsubloc,i] <- padj_uncorr_gif
    mat4[glsubloc,i] <- padj
    matx[glsubloc,i] <- chi2_stat_gif1
    v[i] <- x1$gif
    if (!is.null(new_gif)){
      x1$pvalues_newgif <- as.numeric(pchisq(x1$chi2.stat*x1$gif/new_gif, df = Kval, lower.tail = FALSE))
      padj_newgif <- p.adjust(x1$pvalues_newgif,method=padj_method)
      mat5[glsubloc,i] <- x1$pvalues_newgif
      mat6[glsubloc,i] <- padj_newgif
    }
  }
  colnames(mat1) <- paste("pvalues_i_", sprintf("%05d", seq(1, iterations)), sep = "")
  colnames(mat2) <- paste("pvalues_uncorrected_i_", padj_method, sprintf("%05d", seq(1, iterations)), sep = "")
  colnames(mat3) <- paste("padj_uncorrected_i_", sprintf("%05d", seq(1, iterations)), sep = "")
  colnames(mat4) <- paste("padj_i_", padj_method, sprintf("%05d", seq(1, iterations)), sep = "")
  return(list(pval = mat1, pval_uncorr_gif = mat2, GIF = v, padj_uncorr_gif = mat3, padj = mat4, pval_newgif = mat5, padj_newgif = mat6, newgif = new_gif, chi2stat_gif1 = matx))
}

gl_PCAdapt_iterations_lfmm <- function(gl, n_per_pop, lfmmfile_suffix, K_no, minmaf_threshold = 0.05, #FDR, 
                                       iterations, padj_method = "fdr", new_gif = NULL){
  require(plyr)
  require(dartR)
  require(pcadapt)
  require(LEA)
  mat1 <- matrix(nrow = nLoc(gl), ncol = iterations)
  mat2 <- matrix(nrow = nLoc(gl), ncol = iterations)
  mat3 <- matrix(nrow = nLoc(gl), ncol = iterations)
  mat4 <- matrix(nrow = nLoc(gl), ncol = iterations)
  matx <- matrix(nrow = nLoc(gl), ncol = iterations)
  row.names(mat1) <- locNames(gl)
  row.names(mat2) <- locNames(gl)
  row.names(mat3) <- locNames(gl)
  row.names(mat4) <- locNames(gl)
  row.names(matx) <- locNames(gl)
  v <- rep(NA, iterations)
  if (is.null(new_gif)){
    new_gif = NULL
    mat5 <- NULL
    mat6 <- NULL
  } else {
    mat5 <- matrix(nrow = nLoc(gl), ncol = iterations)
    mat6 <- matrix(nrow = nLoc(gl), ncol = iterations)
    row.names(mat5) <- locNames(gl)
    row.names(mat6) <- locNames(gl)
  }
  for (i in seq(1, iterations)){
    print(i)
    lfmmfile <- paste(lfmmfile_suffix, "i", sprintf("%05d", i), ".lfmm", sep = "")
    locnamefile <- paste(lfmmfile_suffix, "i", sprintf("%05d", i), ".loci.list.csv", sep = "")
    pcadapt_gl <- read.pcadapt(lfmmfile, type = "lfmm")
    pcadapth_sub <- pcadapt(input = pcadapt_gl, K = K_no, min.maf = minmaf_threshold)
    chi2_stat_gif1 <- pcadapth_sub$chi2.stat*pcadapth_sub$gif
    pcadapth_sub$pvalues_gif1 <- as.numeric(pchisq(pcadapth_sub$chi2.stat*pcadapth_sub$gif/1, df = Kval, lower.tail = FALSE))
    padj_uncorr_gif <- p.adjust(pcadapth_sub$pvalues_gif1,method=padj_method)
    padj <- p.adjust(pcadapth_sub$pvalues,method=padj_method)
    # alpha <- FDR
    # outliers <- which(padj < alpha)
    # glsubloc <- locNames(glsub)
    glsubloc <- read.csv(locnamefile, header = F)[[1]]
    mat1[glsubloc,i] <- pcadapth_sub$pvalues
    mat2[glsubloc,i] <- pcadapth_sub$pvalues_gif1
    mat3[glsubloc,i] <- padj_uncorr_gif
    mat4[glsubloc,i] <- padj
    matx[glsubloc,i] <- chi2_stat_gif1
    v[i] <- pcadapth_sub$gif
    if (!is.null(new_gif)){
      # pcadapth_sub$pvalues_newgif <- as.numeric(pchisq(pcadapth_sub$chi2.stat, df = Kval, lower.tail = FALSE))
      pcadapth_sub$pvalues_newgif <- as.numeric(pchisq(pcadapth_sub$chi2.stat*pcadapth_sub$gif/new_gif, df = Kval, lower.tail = FALSE))
      padj_newgif <- p.adjust(pcadapth_sub$pvalues_newgif,method=padj_method)
      mat5[glsubloc,i] <- pcadapth_sub$pvalues_newgif
      mat6[glsubloc,i] <- padj_newgif
    }
  }
  colnames(mat1) <- paste("pvalues_i_", sprintf("%05d", seq(1, iterations)), sep = "")
  colnames(mat2) <- paste("padj_i_", padj_method, sprintf("%05d", seq(1, iterations)), sep = "")
  return(list(pval = mat1, pval_uncorr_gif = mat2, GIF = v, padj_uncorr_gif = mat3, padj = mat4, pval_newgif = mat5, padj_newgif = mat6, newgif = new_gif, chi2stat_gif1 = matx))
}

# SFS functions -----------------------------------------------------#####
plot_dt_SFS_compare_same_n_conf <- function(gl, popname_ls, iterations = 10, indmin = NULL, nomissingness = F){
  # Use: Make SFS plot of allele frequency vs counts
  # Require:  gl - genlight object
  #           popname_ls - list/vector of names of population of interest
  #           iterations - number of iterations (default = 10)
  #
  # Return:   list of dataframe of raw data, dataframe of summarised data,
  #           number of individuals, and a ggplot object with the plot for SFS MAF vs counts 
  # 
  require(plyr)
  require(dartR)
  require(ggplot2)
  # subset the genlight object for only the target populations
  glsub <- gl.keep.pop(gl, pop.list = popname_ls)
  # Determine the number of individuals in the smallest population
  if (is.null(indmin)){
    indmin <- Inf
    for (pop in popname_ls){
      if(nInd(gl.keep.pop(glsub, pop)) < indmin){
        indmin <- nInd(gl.keep.pop(glsub, pop))
      }
    }
  } 
  # Determine binsize from the sample size
  binsize <- 50/indmin
  breaks <- seq(0+binsize/2,50+binsize/2, by = binsize)
  breakssimp <- round(breaks, digits = 2)
  brklab <- paste(breakssimp[1:length(breakssimp)-1], breakssimp[2:length(breakssimp)], sep = " - ")
  brkmid <- (breaks[2:length(breaks)] + breaks[1:length(breaks)-1])/2
  # Making sure that the population is found in the ind.metrics table
  glsub$other$ind.metrics$POP_SUB <- glsub@pop
  # Iterate through the iterations
  for (i in seq(1, iterations)){
    print(paste("Iteration: ", i, sep = ""))
    # dt_subsample_pop <- ddply(glsub$other$ind.metrics, .(POP_SUB), summarise, subsampled_ID = sample(ID, size = eval(indmin), replace = F))
    # dt_subsample_pop <- do.call("ddply",list(glsub$other$ind.metrics, .(POP_SUB), summarize, subsampled_ID = call("sample", as.symbol(id), indmin)))
    dt_subsample_pop <- ddply(glsub$other$ind.metrics, .(POP_SUB), here(summarize),
                              subsampled_ID = sample(ID, eval(indmin)))
    glsubsampled <- gl.keep.ind(glsub, ind.list = dt_subsample_pop$subsampled_ID, mono.rm = T, recalc = T)
    # glsubsampled <- gl.filter.callrate(glsubsampled, threshold = 1, recalc = T, plot.out = F)
    subsetdt <- gl.percent.freq(glsubsampled)
    if (nomissingness == T){
      subsetdt <- subsetdt[subsetdt$nmissing == 0,]
    }
    subsetdt$MAF <- subsetdt$frequency
    # Fold the SFS
    subsetdt$MAF[subsetdt$MAF > 50 & !is.na(subsetdt$MAF)] <- 100 - subsetdt$MAF[subsetdt$MAF > 50 & !is.na(subsetdt$MAF)]
    # Remove monomorphic loci
    subsetdt <- subset(subsetdt, MAF != 0)
    # Initiate counter
    ctr <- 1
    for (pop in popname_ls){
      subsetdtx <- subset(subsetdt, popn == pop)
      nLoc <- dim(subsetdtx)[1]
      print(nLoc)
      dtx <- data.frame(table(cut(subsetdtx$MAF, breaks=breaks, right=TRUE, labels = brklab)))
      # dtx$popn <- paste(pop, nLoc, sep = "; nLoc = ")
      dtx$popn <- pop
      dtx$prop <- dtx$Freq/sum(dtx$Freq)
      dtx$brkmid <- brkmid
      dtx$iterations_nLoc <- nLoc
      if (ctr == 1){
        sumdt <- dtx
      } else {
        sumdt <- rbind(sumdt, dtx)
      }
      ctr <- ctr + 1
    }
    sumdt$iterations <- i
    if (i == 1){
      # print("Assign sumdt to dtcomb")
      dtcomb <- sumdt
    } else {
      # print("rbind sumdt to dtcomb")
      dtcomb <- rbind(dtcomb, sumdt)
    }
  }
  dtcombsum <- ddply(dtcomb, .(popn, Var1, brkmid), summarise, min_nLoc = min(iterations_nLoc), max_nLoc = max(iterations_nLoc), mean_nLoc = mean(iterations_nLoc),
                     prop_min = min(prop), prop_max = max(prop), prop_mean = mean(prop), prop_sd = sd(prop), 
                     prop_quantile2.5 = quantile(prop, probs = c(0.025)), 
                     prop_quantile25 = quantile(prop, probs = c(0.25)), 
                     prop_median = median(prop), 
                     prop_quantile75 = quantile(prop, probs = c(0.75)), 
                     prop_quantile97.5 = quantile(prop, probs = c(0.975)),
                     nLoc_min = min(Freq), nLoc_max = max(Freq), nLoc_mean = mean(Freq), nLoc_sd = sd(Freq), 
                     nLoc_quantile2.5 = quantile(Freq, probs = c(0.025)), 
                     nLoc_quantile25 = quantile(Freq, probs = c(0.25)), 
                     nLoc_median = median(Freq), 
                     nLoc_quantile75 = quantile(Freq, probs = c(0.75)), 
                     nLoc_quantile97.5 = quantile(Freq, probs = c(0.975)),
                     number_of_iterations = max(iterations))
  dtcombsum$poplab <- paste(dtcombsum$popn, "; nLoc = ", round(dtcombsum$mean_nLoc), " (", dtcombsum$min_nLoc, " - ", dtcombsum$max_nLoc, ")", sep = "")
  dtcombsum$popn <- factor(dtcombsum$popn, 
                           levels = popname_ls)
  p1 <- ggplot(data = dtcombsum, mapping = aes(x = brkmid, y = prop_median, ymin = prop_min, 
                                               ymax = prop_max, lower = prop_quantile25,
                                               middle = prop_median, upper = prop_quantile75, 
                                               fill = popn))+
    geom_histogram(stat = "identity", position = "dodge") +
    geom_boxplot(aes(group = interaction(brkmid, popn)),stat = "identity", alpha = 0) +
    xlab("MAF (%)") + ylab("Median proportion of SNPs (25%, 75% and min/max quantiles)") +
    ggtitle(paste("SFS MAF vs frequency; ", paste(popname_ls, collapse = " vs "), "; n = ", indmin, "; i = ", iterations, sep = "")) +
    theme_bw() +
    labs(fill = "Population") +
    theme(legend.background = element_rect(color = "black"),
          legend.position = c(0.95,0.95), 
          legend.justification = c(1,1))
  p2 <- ggplot(data = dtcombsum, mapping = aes(x = brkmid, y = prop_mean, fill = popn))+
    geom_histogram(stat = "identity", position = "dodge") + 
    geom_errorbar(aes(ymin = prop_mean - prop_sd, ymax = prop_mean + prop_sd), width = 2, 
                  position = position_dodge(width = 0.9*binsize)) +
    # xlim(c(0,50)) +
    xlab("MAF (%)") + ylab("Mean proportion of SNPs (+/- SD)") +
    ggtitle(paste("SFS MAF vs frequency; ", paste(popname_ls, collapse = " vs "), "; n = ", indmin, "; i = ", iterations, sep = "")) +
    theme_bw() +
    labs(fill = "Population") +
    theme(legend.background = element_rect(color = "black"),
          legend.position = c(0.95,0.95), 
          legend.justification = c(1,1))
  p3 <- ggplot(data = dtcombsum, mapping = aes(x = brkmid, y = nLoc_median, ymin = nLoc_min, 
                                               ymax = nLoc_max, lower = nLoc_quantile25,
                                               middle = nLoc_median, upper = nLoc_quantile75, 
                                               fill = popn))+
    geom_histogram(stat = "identity", position = "dodge") +
    geom_boxplot(aes(group = interaction(brkmid, popn)),stat = "identity", alpha = 0) +
    xlab("MAF (%)") + ylab("Median number of SNPs (25%, 75% and min/max quantiles)") +
    ggtitle(paste("SFS MAF vs frequency; ", paste(popname_ls, collapse = " vs "), "; n = ", indmin, "; i = ", iterations, sep = "")) +
    theme_bw() +
    labs(fill = "Population") +
    theme(legend.background = element_rect(color = "black"),
          legend.position = c(0.95,0.95), 
          legend.justification = c(1,1))
  return(list(data = dtcomb, datasum = dtcombsum, n = indmin, median_quantile_plot = p1, mean_SD_plot = p2, median_quantile_plot_nLoc = p3))
}

plot_dt_SFS_compare_same_n_conf_only_missingsnps <- function(gl, popname_ls, iterations = 10, indmin = NULL){
  # Use: Make SFS plot of allele frequency vs counts
  # Require:  gl - genlight object
  #           popname_ls - list/vector of names of population of interest
  #           iterations - number of iterations (default = 10)
  #
  # Return:   list of dataframe of raw data, dataframe of summarised data,
  #           number of individuals, and a ggplot object with the plot for SFS MAF vs counts 
  # 
  require(plyr)
  require(dartR)
  require(ggplot2)
  # subset the genlight object for only the target populations
  glsub <- gl.keep.pop(gl, pop.list = popname_ls)
  # Determine the number of individuals in the smallest population
  if (is.null(indmin)){
    indmin <- Inf
    for (pop in popname_ls){
      if(nInd(gl.keep.pop(glsub, pop)) < indmin){
        indmin <- nInd(gl.keep.pop(glsub, pop))
      }
    }
  } 
  # Determine binsize from the sample size
  binsize <- 50/indmin
  breaks <- seq(0+binsize/2,50+binsize/2, by = binsize)
  breakssimp <- round(breaks, digits = 2)
  brklab <- paste(breakssimp[1:length(breakssimp)-1], breakssimp[2:length(breakssimp)], sep = " - ")
  brkmid <- (breaks[2:length(breaks)] + breaks[1:length(breaks)-1])/2
  # Making sure that the population is found in the ind.metrics table
  glsub$other$ind.metrics$POP_SUB <- glsub@pop
  # Iterate through the iterations
  for (i in seq(1, iterations)){
    print(paste("Iteration: ", i, sep = ""))
    # dt_subsample_pop <- ddply(glsub$other$ind.metrics, .(POP_SUB), summarise, subsampled_ID = sample(ID, size = eval(indmin), replace = F))
    # dt_subsample_pop <- do.call("ddply",list(glsub$other$ind.metrics, .(POP_SUB), summarize, subsampled_ID = call("sample", as.symbol(id), indmin)))
    dt_subsample_pop <- ddply(glsub$other$ind.metrics, .(POP_SUB), here(summarize),
                              subsampled_ID = sample(ID, eval(indmin)))
    glsubsampled <- gl.keep.ind(glsub, ind.list = dt_subsample_pop$subsampled_ID, mono.rm = T, recalc = T)
    # glsubsampled <- gl.filter.callrate(glsubsampled, threshold = 1, recalc = T, plot.out = F)
    subsetdt <- gl.percent.freq(glsubsampled)
    subsetdt <- subsetdt[subsetdt$nmissing > 0,] # Only SNPs with missingness
    subsetdt$MAF <- subsetdt$frequency
    # Fold the SFS
    subsetdt$MAF[subsetdt$MAF > 50 & !is.na(subsetdt$MAF)] <- 100 - subsetdt$MAF[subsetdt$MAF > 50 & !is.na(subsetdt$MAF)]
    # Remove monomorphic loci
    subsetdt <- subset(subsetdt, MAF != 0)
    # Initiate counter
    ctr <- 1
    for (pop in popname_ls){
      subsetdtx <- subset(subsetdt, popn == pop)
      nLoc <- dim(subsetdtx)[1]
      print(nLoc)
      dtx <- data.frame(table(cut(subsetdtx$MAF, breaks=breaks, right=TRUE, labels = brklab)))
      # dtx$popn <- paste(pop, nLoc, sep = "; nLoc = ")
      dtx$popn <- pop
      dtx$prop <- dtx$Freq/sum(dtx$Freq)
      dtx$brkmid <- brkmid
      dtx$iterations_nLoc <- nLoc
      if (ctr == 1){
        sumdt <- dtx
      } else {
        sumdt <- rbind(sumdt, dtx)
      }
      ctr <- ctr + 1
    }
    sumdt$iterations <- i
    if (i == 1){
      # print("Assign sumdt to dtcomb")
      dtcomb <- sumdt
    } else {
      # print("rbind sumdt to dtcomb")
      dtcomb <- rbind(dtcomb, sumdt)
    }
  }
  dtcombsum <- ddply(dtcomb, .(popn, Var1, brkmid), summarise, min_nLoc = min(iterations_nLoc), max_nLoc = max(iterations_nLoc), mean_nLoc = mean(iterations_nLoc),
                     prop_min = min(prop), prop_max = max(prop), prop_mean = mean(prop), prop_sd = sd(prop), 
                     prop_quantile2.5 = quantile(prop, probs = c(0.025)), 
                     prop_quantile25 = quantile(prop, probs = c(0.25)), 
                     prop_median = median(prop), 
                     prop_quantile75 = quantile(prop, probs = c(0.75)), 
                     prop_quantile97.5 = quantile(prop, probs = c(0.975)),
                     nLoc_min = min(Freq), nLoc_max = max(Freq), nLoc_mean = mean(Freq), nLoc_sd = sd(Freq), 
                     nLoc_quantile2.5 = quantile(Freq, probs = c(0.025)), 
                     nLoc_quantile25 = quantile(Freq, probs = c(0.25)), 
                     nLoc_median = median(Freq), 
                     nLoc_quantile75 = quantile(Freq, probs = c(0.75)), 
                     nLoc_quantile97.5 = quantile(Freq, probs = c(0.975)),
                     number_of_iterations = max(iterations))
  dtcombsum$poplab <- paste(dtcombsum$popn, "; nLoc = ", round(dtcombsum$mean_nLoc), " (", dtcombsum$min_nLoc, " - ", dtcombsum$max_nLoc, ")", sep = "")
  dtcombsum$popn <- factor(dtcombsum$popn, 
                           levels = popname_ls)
  p1 <- ggplot(data = dtcombsum, mapping = aes(x = brkmid, y = prop_median, ymin = prop_min, 
                                               ymax = prop_max, lower = prop_quantile25,
                                               middle = prop_median, upper = prop_quantile75, 
                                               fill = popn))+
    geom_histogram(stat = "identity", position = "dodge") +
    geom_boxplot(aes(group = interaction(brkmid, popn)),stat = "identity", alpha = 0) +
    xlab("MAF (%)") + ylab("Median proportion of SNPs (25%, 75% and min/max quantiles)") +
    ggtitle(paste("SFS MAF vs frequency; ", paste(popname_ls, collapse = " vs "), "; n = ", indmin, "; i = ", iterations, sep = "")) +
    theme_bw() +
    labs(fill = "Population") +
    theme(legend.background = element_rect(color = "black"),
          legend.position = c(0.95,0.95), 
          legend.justification = c(1,1))
  p2 <- ggplot(data = dtcombsum, mapping = aes(x = brkmid, y = prop_mean, fill = popn))+
    geom_histogram(stat = "identity", position = "dodge") + 
    geom_errorbar(aes(ymin = prop_mean - prop_sd, ymax = prop_mean + prop_sd), width = 2, 
                  position = position_dodge(width = 0.9*binsize)) +
    # xlim(c(0,50)) +
    xlab("MAF (%)") + ylab("Mean proportion of SNPs (+/- SD)") +
    ggtitle(paste("SFS MAF vs frequency; ", paste(popname_ls, collapse = " vs "), "; n = ", indmin, "; i = ", iterations, sep = "")) +
    theme_bw() +
    labs(fill = "Population") +
    theme(legend.background = element_rect(color = "black"),
          legend.position = c(0.95,0.95), 
          legend.justification = c(1,1))
  p3 <- ggplot(data = dtcombsum, mapping = aes(x = brkmid, y = nLoc_median, ymin = nLoc_min, 
                                               ymax = nLoc_max, lower = nLoc_quantile25,
                                               middle = nLoc_median, upper = nLoc_quantile75, 
                                               fill = popn))+
    geom_histogram(stat = "identity", position = "dodge") +
    geom_boxplot(aes(group = interaction(brkmid, popn)),stat = "identity", alpha = 0) +
    xlab("MAF (%)") + ylab("Median number of SNPs (25%, 75% and min/max quantiles)") +
    ggtitle(paste("SFS MAF vs frequency; ", paste(popname_ls, collapse = " vs "), "; n = ", indmin, "; i = ", iterations, sep = "")) +
    theme_bw() +
    labs(fill = "Population") +
    theme(legend.background = element_rect(color = "black"),
          legend.position = c(0.95,0.95), 
          legend.justification = c(1,1))
  return(list(data = dtcomb, datasum = dtcombsum, n = indmin, median_quantile_plot = p1, mean_SD_plot = p2, median_quantile_plot_nLoc = p3))
}

plot_dt_SFS_compare_same_n_conf_nomissingsnps <- function(gl, popname_ls, iterations = 10, indmin = NULL){
  # Use: Make SFS plot of allele frequency vs counts
  # Require:  gl - genlight object
  #           popname_ls - list/vector of names of population of interest
  #           iterations - number of iterations (default = 10)
  #
  # Return:   list of dataframe of raw data, dataframe of summarised data,
  #           number of individuals, and a ggplot object with the plot for SFS MAF vs counts 
  # 
  require(plyr)
  require(dartR)
  require(ggplot2)
  # subset the genlight object for only the target populations
  glsub <- gl.keep.pop(gl, pop.list = popname_ls)
  # Determine the number of individuals in the smallest population
  if (is.null(indmin)){
    indmin <- Inf
    for (pop in popname_ls){
      if(nInd(gl.keep.pop(glsub, pop)) < indmin){
        indmin <- nInd(gl.keep.pop(glsub, pop))
      }
    }
  } 
  # Determine binsize from the sample size
  binsize <- 50/indmin
  breaks <- seq(0+binsize/2,50+binsize/2, by = binsize)
  breakssimp <- round(breaks, digits = 2)
  brklab <- paste(breakssimp[1:length(breakssimp)-1], breakssimp[2:length(breakssimp)], sep = " - ")
  brkmid <- (breaks[2:length(breaks)] + breaks[1:length(breaks)-1])/2
  # Making sure that the population is found in the ind.metrics table
  glsub$other$ind.metrics$POP_SUB <- glsub@pop
  # Iterate through the iterations
  for (i in seq(1, iterations)){
    print(paste("Iteration: ", i, sep = ""))
    # dt_subsample_pop <- ddply(glsub$other$ind.metrics, .(POP_SUB), summarise, subsampled_ID = sample(ID, size = eval(indmin), replace = F))
    # dt_subsample_pop <- do.call("ddply",list(glsub$other$ind.metrics, .(POP_SUB), summarize, subsampled_ID = call("sample", as.symbol(id), indmin)))
    dt_subsample_pop <- ddply(glsub$other$ind.metrics, .(POP_SUB), here(summarize),
                              subsampled_ID = sample(ID, eval(indmin)))
    glsubsampled <- gl.keep.ind(glsub, ind.list = dt_subsample_pop$subsampled_ID, mono.rm = T, recalc = T)
    # glsubsampled <- gl.filter.callrate(glsubsampled, threshold = 1, recalc = T, plot.out = F)
    subsetdt <- gl.percent.freq(glsubsampled)
    subsetdt <- subsetdt[subsetdt$nmissing == 0,] # Only SNPs with no missingness
    subsetdt$MAF <- subsetdt$frequency
    # Fold the SFS
    subsetdt$MAF[subsetdt$MAF > 50 & !is.na(subsetdt$MAF)] <- 100 - subsetdt$MAF[subsetdt$MAF > 50 & !is.na(subsetdt$MAF)]
    # Remove monomorphic loci
    subsetdt <- subset(subsetdt, MAF != 0)
    # Initiate counter
    ctr <- 1
    for (pop in popname_ls){
      subsetdtx <- subset(subsetdt, popn == pop)
      nLoc <- dim(subsetdtx)[1]
      print(nLoc)
      dtx <- data.frame(table(cut(subsetdtx$MAF, breaks=breaks, right=TRUE, labels = brklab)))
      # dtx$popn <- paste(pop, nLoc, sep = "; nLoc = ")
      dtx$popn <- pop
      dtx$prop <- dtx$Freq/sum(dtx$Freq)
      dtx$brkmid <- brkmid
      dtx$iterations_nLoc <- nLoc
      if (ctr == 1){
        sumdt <- dtx
      } else {
        sumdt <- rbind(sumdt, dtx)
      }
      ctr <- ctr + 1
    }
    sumdt$iterations <- i
    if (i == 1){
      # print("Assign sumdt to dtcomb")
      dtcomb <- sumdt
    } else {
      # print("rbind sumdt to dtcomb")
      dtcomb <- rbind(dtcomb, sumdt)
    }
  }
  dtcombsum <- ddply(dtcomb, .(popn, Var1, brkmid), summarise, min_nLoc = min(iterations_nLoc), max_nLoc = max(iterations_nLoc), mean_nLoc = mean(iterations_nLoc),
                     prop_min = min(prop), prop_max = max(prop), prop_mean = mean(prop), prop_sd = sd(prop), 
                     prop_quantile2.5 = quantile(prop, probs = c(0.025)), 
                     prop_quantile25 = quantile(prop, probs = c(0.25)), 
                     prop_median = median(prop), 
                     prop_quantile75 = quantile(prop, probs = c(0.75)), 
                     prop_quantile97.5 = quantile(prop, probs = c(0.975)),
                     nLoc_min = min(Freq), nLoc_max = max(Freq), nLoc_mean = mean(Freq), nLoc_sd = sd(Freq), 
                     nLoc_quantile2.5 = quantile(Freq, probs = c(0.025)), 
                     nLoc_quantile25 = quantile(Freq, probs = c(0.25)), 
                     nLoc_median = median(Freq), 
                     nLoc_quantile75 = quantile(Freq, probs = c(0.75)), 
                     nLoc_quantile97.5 = quantile(Freq, probs = c(0.975)),
                     number_of_iterations = max(iterations))
  dtcombsum$poplab <- paste(dtcombsum$popn, "; nLoc = ", round(dtcombsum$mean_nLoc), " (", dtcombsum$min_nLoc, " - ", dtcombsum$max_nLoc, ")", sep = "")
  dtcombsum$popn <- factor(dtcombsum$popn, 
                           levels = popname_ls)
  p1 <- ggplot(data = dtcombsum, mapping = aes(x = brkmid, y = prop_median, ymin = prop_min, 
                                               ymax = prop_max, lower = prop_quantile25,
                                               middle = prop_median, upper = prop_quantile75, 
                                               fill = popn))+
    geom_histogram(stat = "identity", position = "dodge") +
    geom_boxplot(aes(group = interaction(brkmid, popn)),stat = "identity", alpha = 0) +
    xlab("MAF (%)") + ylab("Median proportion of SNPs (25%, 75% and min/max quantiles)") +
    ggtitle(paste("SFS MAF vs frequency; ", paste(popname_ls, collapse = " vs "), "; n = ", indmin, "; i = ", iterations, sep = "")) +
    theme_bw() +
    labs(fill = "Population") +
    theme(legend.background = element_rect(color = "black"),
          legend.position = c(0.95,0.95), 
          legend.justification = c(1,1))
  p2 <- ggplot(data = dtcombsum, mapping = aes(x = brkmid, y = prop_mean, fill = popn))+
    geom_histogram(stat = "identity", position = "dodge") + 
    geom_errorbar(aes(ymin = prop_mean - prop_sd, ymax = prop_mean + prop_sd), width = 2, 
                  position = position_dodge(width = 0.9*binsize)) +
    # xlim(c(0,50)) +
    xlab("MAF (%)") + ylab("Mean proportion of SNPs (+/- SD)") +
    ggtitle(paste("SFS MAF vs frequency; ", paste(popname_ls, collapse = " vs "), "; n = ", indmin, "; i = ", iterations, sep = "")) +
    theme_bw() +
    labs(fill = "Population") +
    theme(legend.background = element_rect(color = "black"),
          legend.position = c(0.95,0.95), 
          legend.justification = c(1,1))
  p3 <- ggplot(data = dtcombsum, mapping = aes(x = brkmid, y = nLoc_median, ymin = nLoc_min, 
                                               ymax = nLoc_max, lower = nLoc_quantile25,
                                               middle = nLoc_median, upper = nLoc_quantile75, 
                                               fill = popn))+
    geom_histogram(stat = "identity", position = "dodge") +
    geom_boxplot(aes(group = interaction(brkmid, popn)),stat = "identity", alpha = 0) +
    xlab("MAF (%)") + ylab("Median number of SNPs (25%, 75% and min/max quantiles)") +
    ggtitle(paste("SFS MAF vs frequency; ", paste(popname_ls, collapse = " vs "), "; n = ", indmin, "; i = ", iterations, sep = "")) +
    theme_bw() +
    labs(fill = "Population") +
    theme(legend.background = element_rect(color = "black"),
          legend.position = c(0.95,0.95), 
          legend.justification = c(1,1))
  return(list(data = dtcomb, datasum = dtcombsum, n = indmin, median_quantile_plot = p1, mean_SD_plot = p2, median_quantile_plot_nLoc = p3))
}

plot_SFS_from_gl <- function(gl){
  # Use: Make SFS plot of allele frequency vs counts
  # Require:  gl - genlight object
  #           popname_ls - list/vector of names of population of interest
  #           iterations - number of iterations (default = 10)
  #
  # Return:   list of dataframe of raw data, dataframe of summarised data,
  #           number of individuals, and a ggplot object with the plot for SFS MAF vs counts 
  # 
  require(plyr)
  require(dartR)
  require(ggplot2)
  # subset the genlight object for only the target populations
  af_dt <- gl.percent.freq(gl)
  
  # Determine binsize from the sample size
  binsize <- 50/indmin
  breaks <- seq(0+binsize/2,50+binsize/2, by = binsize)
  breakssimp <- round(breaks, digits = 2)
  brklab <- paste(breakssimp[1:length(breakssimp)-1], breakssimp[2:length(breakssimp)], sep = " - ")
  brkmid <- (breaks[2:length(breaks)] + breaks[1:length(breaks)-1])/2
  # Making sure that the population is found in the ind.metrics table
  glsub$other$ind.metrics$POP_SUB <- glsub@pop
  # Iterate through the iterations
  for (i in seq(1, iterations)){
    print(paste("Iteration: ", i, sep = ""))
    # dt_subsample_pop <- ddply(glsub$other$ind.metrics, .(POP_SUB), summarise, subsampled_ID = sample(ID, size = eval(indmin), replace = F))
    # dt_subsample_pop <- do.call("ddply",list(glsub$other$ind.metrics, .(POP_SUB), summarize, subsampled_ID = call("sample", as.symbol(id), indmin)))
    dt_subsample_pop <- ddply(glsub$other$ind.metrics, .(POP_SUB), here(summarize),
                              subsampled_ID = sample(ID, eval(indmin)))
    glsubsampled <- gl.keep.ind(glsub, ind.list = dt_subsample_pop$subsampled_ID, mono.rm = T, recalc = T)
    glsubsampled <- gl.filter.callrate(glsubsampled, threshold = 1, recalc = T, plot.out = F)
    subsetdt <- gl.percent.freq(glsubsampled)
    subsetdt$MAF <- subsetdt$frequency
    # Fold the SFS
    subsetdt$MAF[subsetdt$MAF > 50 & !is.na(subsetdt$MAF)] <- 100 - subsetdt$MAF[subsetdt$MAF > 50 & !is.na(subsetdt$MAF)]
    # Remove monomorphic loci
    subsetdt <- subset(subsetdt, MAF != 0)
    # Initiate counter
    ctr <- 1
    for (pop in popname_ls){
      subsetdtx <- subset(subsetdt, popn == pop)
      nLoc <- dim(subsetdtx)[1]
      print(nLoc)
      dtx <- data.frame(table(cut(subsetdtx$MAF, breaks=breaks, right=TRUE, labels = brklab)))
      # dtx$popn <- paste(pop, nLoc, sep = "; nLoc = ")
      dtx$popn <- pop
      dtx$prop <- dtx$Freq/sum(dtx$Freq)
      dtx$brkmid <- brkmid
      dtx$iterations_nLoc <- nLoc
      if (ctr == 1){
        sumdt <- dtx
      } else {
        sumdt <- rbind(sumdt, dtx)
      }
      ctr <- ctr + 1
    }
    sumdt$iterations <- i
    if (i == 1){
      # print("Assign sumdt to dtcomb")
      dtcomb <- sumdt
    } else {
      # print("rbind sumdt to dtcomb")
      dtcomb <- rbind(dtcomb, sumdt)
    }
  }
  dtcombsum <- ddply(dtcomb, .(popn, Var1, brkmid), summarise, min_nLoc = min(iterations_nLoc), max_nLoc = max(iterations_nLoc), mean_nLoc = mean(iterations_nLoc),
                     prop_min = min(prop), prop_max = max(prop), prop_mean = mean(prop), prop_sd = sd(prop), 
                     prop_quantile2.5 = quantile(prop, probs = c(0.025)), 
                     prop_quantile25 = quantile(prop, probs = c(0.25)), 
                     prop_median = median(prop), 
                     prop_quantile75 = quantile(prop, probs = c(0.75)), 
                     prop_quantile97.5 = quantile(prop, probs = c(0.975)),
                     nLoc_min = min(Freq), nLoc_max = max(Freq), nLoc_mean = mean(Freq), nLoc_sd = sd(Freq), 
                     nLoc_quantile2.5 = quantile(Freq, probs = c(0.025)), 
                     nLoc_quantile25 = quantile(Freq, probs = c(0.25)), 
                     nLoc_median = median(Freq), 
                     nLoc_quantile75 = quantile(Freq, probs = c(0.75)), 
                     nLoc_quantile97.5 = quantile(Freq, probs = c(0.975)),
                     number_of_iterations = max(iterations))
  dtcombsum$poplab <- paste(dtcombsum$popn, "; nLoc = ", round(dtcombsum$mean_nLoc), " (", dtcombsum$min_nLoc, " - ", dtcombsum$max_nLoc, ")", sep = "")
  dtcombsum$popn <- factor(dtcombsum$popn, 
                           levels = popname_ls)
  p1 <- ggplot(data = dtcombsum, mapping = aes(x = brkmid, y = prop_median, ymin = prop_min, 
                                               ymax = prop_max, lower = prop_quantile25,
                                               middle = prop_median, upper = prop_quantile75, 
                                               fill = popn))+
    geom_histogram(stat = "identity", position = "dodge") +
    geom_boxplot(aes(group = interaction(brkmid, popn)),stat = "identity", alpha = 0) +
    xlab("MAF (%)") + ylab("Median proportion of SNPs (25%, 75% and min/max quantiles)") +
    ggtitle(paste("SFS MAF vs frequency; ", paste(popname_ls, collapse = " vs "), "; n = ", indmin, "; i = ", iterations, sep = "")) +
    theme_bw() +
    labs(fill = "Population") +
    theme(legend.background = element_rect(color = "black"),
          legend.position = c(0.95,0.95), 
          legend.justification = c(1,1))
  p2 <- ggplot(data = dtcombsum, mapping = aes(x = brkmid, y = prop_mean, fill = popn))+
    geom_histogram(stat = "identity", position = "dodge") + 
    geom_errorbar(aes(ymin = prop_mean - prop_sd, ymax = prop_mean + prop_sd), width = 2, 
                  position = position_dodge(width = 0.9*binsize)) +
    # xlim(c(0,50)) +
    xlab("MAF (%)") + ylab("Mean proportion of SNPs (+/- SD)") +
    ggtitle(paste("SFS MAF vs frequency; ", paste(popname_ls, collapse = " vs "), "; n = ", indmin, "; i = ", iterations, sep = "")) +
    theme_bw() +
    labs(fill = "Population") +
    theme(legend.background = element_rect(color = "black"),
          legend.position = c(0.95,0.95), 
          legend.justification = c(1,1))
  p3 <- ggplot(data = dtcombsum, mapping = aes(x = brkmid, y = nLoc_median, ymin = nLoc_min, 
                                               ymax = nLoc_max, lower = nLoc_quantile25,
                                               middle = nLoc_median, upper = nLoc_quantile75, 
                                               fill = popn))+
    geom_histogram(stat = "identity", position = "dodge") +
    geom_boxplot(aes(group = interaction(brkmid, popn)),stat = "identity", alpha = 0) +
    xlab("MAF (%)") + ylab("Median number of SNPs (25%, 75% and min/max quantiles)") +
    ggtitle(paste("SFS MAF vs frequency; ", paste(popname_ls, collapse = " vs "), "; n = ", indmin, "; i = ", iterations, sep = "")) +
    theme_bw() +
    labs(fill = "Population") +
    theme(legend.background = element_rect(color = "black"),
          legend.position = c(0.95,0.95), 
          legend.justification = c(1,1))
  return(list(data = dtcomb, datasum = dtcombsum, n = indmin, median_quantile_plot = p1, mean_SD_plot = p2, median_quantile_plot_nLoc = p3))
}


# Diversity indices -------------------------------------------------#####

genalex2genepop <- function(genalex_csv, genepop_file, genepop_title = ""){
  # This function converts genalex csv file (exported through gl2genalex) to
  # genepop file.
  #
  # Arguments:  genalex_csv = path to csv file exported through gl2genalex function in dartR
  #             genepop_file = path to output genepop file
  #             genepop_title = information to go into the first line of the genepop file
  dtint1 <- read.csv(genalex_csv, skip = 2, check.names = F)
  for (i in 3:dim(dtint1)[2]){
    dtint1[,i] <- sprintf("%02d",dtint1[,i]) 
  }
  dtint2  <- dtint1
  for (i in seq(3, dim(dtint2)[2], 2)){
    dtint2[,i] <- paste(dtint2[,i], dtint2[,i+1], sep = "")
  }
  dtint2 <- dtint2[,c(1,2,seq(3, dim(dtint2)[2], 2))]
  dtint3 <- dtint2[,-1]
  dtint3[,1] <- paste(dtint3[,1], ",", sep = '')
  locnames <- colnames(dtint3)[seq(2, dim(dtint3)[2], 1)]
  filename <- genepop_file
  write.table(paste('Title:', genepop_title), file = filename, col.names = F, row.names = F, quote = F)
  write.table(locnames, file = filename, append = T, row.names = F, col.names = F, quote = F)
  for (i in unique(dtint3$Pop)){
    write.table("Pop", file = filename, append = T, row.names = F, col.names = F, quote = F)
    write.table(subset(dtint3, Pop == i), file = filename, append = T, row.names = F, col.names = F, quote = F)
  }
}
# 
calc_prop_polymorphic_loci_perpop <- function(gl){
  require(dartR)
  library(tidyverse)
  poplist <- seppop(gl)
  count <- 0
  for (i in poplist) {
    count <- count + 1
    matb <- as.matrix(i)
    #####
    dtb <- data.frame(loci = colnames(matb), loci_hetero = apply(matb == 1, 2, sum, na.rm = T), loci_homo_alt = apply(matb == 2, 2, sum, na.rm = T), loci_homo_ref = apply(matb == 0, 2, sum, na.rm = T))
    dtb$Alt_Alnct <- 1*dtb$loci_hetero + 2*dtb$loci_homo_alt
    dtb$Ref_Alnct <- 1*dtb$loci_hetero + 2*dtb$loci_homo_ref
    dtb$polymorph_loci <- 0
    dtb$polymorph_loci[dtb$Alt_Alnct > 0 & dtb$Ref_Alnct > 0] <- 1
    npoly <- sum(dtb$polymorph_loci, na.rm = T)
    ntot <- dim(dtb)[1]
    if (count == 1) {
      Population <- names(poplist)[count]
      result <- data.frame(pop = Population,
                           n_loci_poly = npoly,
                           n_loci_tot = ntot,
                           prop_poly_loci = npoly/ntot)
    } else {
      Population <- names(poplist)[count]
      r <- data.frame(pop = Population,
                      n_loci_poly = npoly,
                      n_loci_tot = ntot,
                      prop_poly_loci = npoly/ntot)
      result <- rbind(result, r)
    }
  }
  return(result)
}

calc_prop_polymorphic_loci_perpop_rarefact <- function(gl, n_per_pop, n_iter){
  require(dartR)
  require(plyr)
  # library(tidyverse)
  count <- 0
  gl$other$ind.metrics$popsub <- pop(gl)
  for (i in seq(1, n_iter)){
    print(i)
    count <- count + 1
    dt_subsample_pop <- ddply(gl$other$ind.metrics, .(popsub), plyr::here(summarise), subsampled_ID = sample(ID, n_per_pop))
    glsub <- gl.keep.ind(gl, ind.list = dt_subsample_pop$subsampled_ID)
    if (i == 1){
      dtsub <- calc_prop_polymorphic_loci_perpop(glsub)
      dtsub$iteration_number <- i
    } else {
      dtsubx <- calc_prop_polymorphic_loci_perpop(glsub)
      dtsubx$iteration_number <- i
      dtsub <- rbind(dtsub, dtsubx)
    }
  }
  dtsum <- ddply(dtsub, .(pop), plyr::here(summarise), 
                 mean_prop_poly_loci = mean(prop_poly_loci), 
                 sd_prop_poly_loci = sd(prop_poly_loci), 
                 median_prop_poly_loci = median(prop_poly_loci), 
                 prop_poly_loci_0_025 = quantile(prop_poly_loci, probs = 0.025),
                 prop_poly_loci_0_975 = quantile(prop_poly_loci, probs = 0.975),
                 mean_poly_loci_n = mean(n_loci_poly),
                 sd_poly_loci_n = sd(n_loci_poly), 
                 median_poly_loci_n = median(n_loci_poly), 
                 prop_loci_n_0_025 = quantile(n_loci_poly, probs = 0.025),
                 prop_loci_n_0_975 = quantile(n_loci_poly, probs = 0.975))
  return(list(n_per_pop = n_per_pop, dt_full = dtsub, dt_summary = dtsum))
}


# LEA ---------------------------------------------------------------#####

LEA_snmf_multiplot <- function(glsub, snmf_sub, qmatsuffix, pdfoutput, k2plot, indperline_plot, nrow_plot, MAF_threshold = 0){
  require(gridExtra)
  require(grid)
  require(LEA)
  require(dartR)
  require(plyr)
  require(pophelper)
  for (i in k2plot){
    ce = cross.entropy(snmf_sub, K = i)
    best = which.min(ce)
    qmatrix.all.ki = Q(snmf_sub, K = i, run = best)
    write.table(qmatrix.all.ki, file = paste(qmatsuffix, sprintf("%02d", i), ".out", sep = ""), row.names = F, col.names = F)
  }
  pathqmat <- sub("/[^/]+$", "", qmatsuffix)
  qmatpattern <- paste(sub(".*/", "", qmatsuffix), "*", sep = "")
  print(c(pathqmat, qmatpattern))
  matk_ALL <- readQBasic(files = list.files(path = pathqmat, pattern = qmatpattern, full.names = T))
  pop_label_ALL <- data.frame(pop_names = as.character(pop(glsub)), 
                              individual_id = indNames(glsub))
  
  pdf(pdfoutput, width = 10, height = 5)
  plot(snmf_sub, col = "blue4", cex = 1.4, pch = 19, main = "cross entropy plot")
  for (i in seq(1, length(matk_ALL))){
    p2ki <- plotQMultiline(matk_ALL[i], spl = indperline_plot, lpp = nrow_plot, returnplot=T,exportplot=F,basesize=11, grplab = pop_label_ALL[, c("pop_names"), drop = F], ordergrp = T)
    grid.arrange(p2ki$plot[[1]][[1]], top=textGrob(paste("K =", k2plot[i])))
  }
  dev.off()
}

# Merge sNMF repetitions 
merge_sNMF_qmatrix_multiK <- function(snmf_obj, gl, k2plot){
  # This function reads in the snmf obj and return lists of qmatrix which can be
  # used for plotQ or plotQmultiline
  # Loop through each K value
  require(LEA)
  require(pophelper)
  samplenames <- indNames(gl)
  q.df.ls <- list()
  for (i in k2plot){
    ce = cross.entropy(snmf_obj, K = i)
    rep = 1:length(ce)
    q.df <- list()
    # Loop through each repetition
    for (j in rep){
      q.matrix <- Q(snmf_obj, K = i, run = j)
      rownames(q.matrix) <- samplenames
      colnames(q.matrix) <- paste0("Cluster",1:ncol(q.matrix))
      # Add q.matrix to list
      q.df[[paste("replicate_", sep = "", j)]] <- as.data.frame(q.matrix)
    }
    # merge each K value
    q.df <- alignK(as.qlist(q.df))
    q.merged <- mergeQ(q.df)
    # store in list 
    q.merged <- as.data.frame(q.merged)
    colnames(q.merged) <- paste0("Cluster",1:ncol(q.merged))
    q.df.ls[[paste("K_", sep = "", i)]] <- as.data.frame(q.merged)
  }
  q.df.ls <- as.qlist(q.df.ls)
  return(q.df.ls)
}

# HWE filters -------------------------------------------------------#####
HWExact_multipop <- function(gl, pval_type = "selome"){
  require(dartR)
  require(HardyWeinberg)
  poplist <- seppop(gl)
  count <- 0
  npops2plot <- nPop(gl)
  for (i in poplist) {
    count <- count + 1
    ii <- gl.filter.monomorphs(i,verbose =0)
    matINDoth <- as.matrix(ii)
    matHWExact <- matrix(nrow = dim(matINDoth)[2], ncol = 3)
    rownames(matHWExact) <- colnames(matINDoth)
    colnames(matHWExact) <- c("AA", "AB", "BB")
    matHWExact[,1] <- apply(matINDoth == 0, 2, sum, na.rm = T)
    matHWExact[,2] <- apply(matINDoth == 1, 2, sum, na.rm = T)
    matHWExact[,3] <- apply(matINDoth == 2, 2, sum, na.rm = T)
    HW.test_pval <- HWExactMat(matHWExact, pvaluetype = pval_type)
    if (count == 1) {
      Population <- rep(names(poplist)[count],length(HW.test_pval))
      result <- data.frame(Pop = Population,
                           loci_name = rownames(matHWExact),
                           pval = HW.test_pval)
    } else {
      Population <- rep(names(poplist)[count],length(HW.test_pval))
      r <- data.frame(Pop = Population,
                      loci_name = rownames(matHWExact),
                      pval = HW.test_pval)
      result <- rbind(result, r)
    }
  }
  return(result)
}


# Define functions
compare_rare_alle_f2 <- function(gl, popname_ls, ref_popname, rare_alle_f = 10, major_alle_f = 50){
  # Function determines rare loci with low allele frequencies in reference
  # population and check what happens to them in each of the population in 
  # popname_ls 
  require(plyr)
  require(dartR)
  require(reshape2)
  poptot <- c(popname_ls, ref_popname)
  glsub <- gl.keep.pop(gl, pop.list = poptot)
  glsub$other$ind.metrics$POP_SUB <- glsub@pop
  pop.ls <- popname_ls
  # Calculate allele frequency
  dt <- gl.percent.freq(gl)
  dt1 <- subset(dt, popn == ref_popname)
  nref <- max(dt1$n)
  dt1a <- subset(dt1, frequency <= rare_alle_f)
  nrefloc <- dim(dt1a)[1]
  pop_ref_loci <- dt1a$locus
  dt2 <- subset(dt, (locus %in% dt1a$locus) & (popn %in% popname_ls))
  # pop.ls <- unique(dt2$popn)
  colnames <- c("reference_pop", "reference_pop_n", "reference_pop_MAF_loc_n", 
                "popnames", "n_ind", 
                "n_rare_allele_lost", "n_rare_allele_to_major_allele")
  pop.mat <- matrix(nrow = length(pop.ls), ncol = length(colnames), dimnames = list(pop.ls, colnames))
  pop.dt <- as.data.frame(pop.mat, index = NULL)
  pop.dt$reference_pop <- ref_popname
  pop.dt$reference_pop_n <- nref
  pop.dt$reference_pop_MAF_loc_n <- nrefloc
  ls1 <- list()
  ls2 <- list()
  ls3 <- list()
  # Initiate counter
  ctr <- 1
  for (pop in pop.ls){
    dt3 <- subset(dt2, popn == pop)
    ncomp <- max(dt3$n)
    dt3a <- subset(dt3, frequency == 0) # lost
    dt3b <- subset(dt3, frequency >= major_alle_f) # become major alle f
    dt3c <- subset(dt3, (frequency > 0) & (frequency < major_alle_f)) # lost test
    pop.dt$popnames[ctr] <- pop
    pop.dt$n_ind[ctr] <- ncomp
    pop.dt$n_rare_allele_lost[ctr] <- dim(dt3a)[1]
    pop.dt$n_rare_allele_to_major_allele[ctr] <- dim(dt3b)[1]
    sumdt <- pop.dt
    ls1[[pop]] <- dt3a$locus
    ls2[[pop]] <- dt3b$locus
    ls3[[pop]] <- dt3c$locus
    ctr <- ctr + 1
  }
  # Make matrix of 0,1,2; 0 = lost, 1 = the same, 2 = became Major allele
  mat1 <- matrix(nrow = length(pop_ref_loci), ncol = length(poptot))
  colnames(mat1) <- poptot
  rownames(mat1) <- pop_ref_loci
  mat1[as.character(pop_ref_loci), ref_popname] <- "MAF"
  for (pop in popname_ls){
    mat1[as.character(ls1[[pop]]),pop] <- "lost"
    mat1[as.character(ls2[[pop]]),pop] <- "Major allele"
    mat1[as.character(ls3[[pop]]),pop] <- "MAF"
  }
  # Convert matrix to long format with some form of ordering
  dt1 <- data.frame(mat1, check.names = F)
  dt1$loci <- as.character(rownames(mat1))
  ctr <- 0
  for (pop in popname_ls){
    ctr <- ctr + 1
    dt1 <- dt1[order(dt1[[pop]]),]
    dt1[[paste("loci_order", ctr, sep = "")]] <- seq(1, dim(dt1)[1])
  }
  dt1_col <- colnames(dt1)
  dt1_id.vars <- dt1_col[!(dt1_col %in% poptot)] 
  dt2 <- melt(dt1, id.vars = dt1_id.vars,
              variable_name = "pops")
  # output
  return(list(dt = sumdt, 
              mat_sum = mat1,
              dt_sum_wide = dt1,
              dt_sum_long = dt2,
              n_ref = nref, 
              n_loc_ref = nrefloc, 
              pop_ref_MAF_loci = pop_ref_loci,
              loci_lost = ls1, 
              loci_to_MAF = ls2, 
              MAF_loci_retained_as_is = ls3,
              rare_alle_f = rare_alle_f, 
              major_alle_f = major_alle_f))
}

# Reich et al. (2009) FST estimates


reich.fst <- function(gl, bootstrap=FALSE, plot=FALSE, verbose=TRUE) { 
  # The function below has been downloaded from:
  # https://github.com/jessicarick/reich-fst
  #
  # The script was used in Junker et al. 2020, Molecular Ecology (https://doi.org/10.1111/mec.15559).
  #
  ## reich fst estimator
  ## vectorized version
  ## input=genlight object
  ## FST will be calculated between pops in genlight object
  ## specify number of bootstraps using "bootstrap=100"
  if (!require("matrixStats",character.only=T, quietly=T)) {
    install.packages("matrixStats")
    library(matrixStats, character.only=T)
  }
  if (!require("dplyr",character.only=T, quietly=T)) {
    install.packages("dplyr")
    library(dplyr, character.only=T)
  }
  
  nloc <- gl@n.loc
  npop <- length(levels(gl@pop))
  
  fsts <- matrix(nrow=npop,
                 ncol=npop,
                 dimnames=list(levels(gl@pop),levels(gl@pop)))
  
  if (bootstrap != FALSE){
    n.bs <- bootstrap
    bs <- data.frame(matrix(nrow=nrow(combinat::combn2(levels(gl@pop))),
                            ncol=n.bs+5))
  }
  
  k <- 0
  
  for (p1 in levels(gl@pop)){
    for (p2 in levels(gl@pop)){
      if (which(levels(gl@pop) == p1) < which(levels(gl@pop) == p2)) {
        k <- 1+k
        
        pop1 <- gl.keep.pop(gl, p1, mono.rm=FALSE, v=0)
        pop2 <- gl.keep.pop(gl, p2, mono.rm=FALSE, v=0)
        
        a1 <- colSums2(as.matrix(pop1),na.rm=T)
        a2 <- colSums2(as.matrix(pop2),na.rm=T)
        n1 <- apply(as.matrix(pop1),2,function(x) 2*sum(!is.na(x)))
        n2 <- apply(as.matrix(pop2),2,function(x) 2*sum(!is.na(x)))
        
        h1 <- (a1*(n1-a1))/(n1*(n1-1))
        h2 <- (a2*(n2-a2))/(n2*(n2-1))
        
        N <- (a1/n1 - a2/n2)^2 - h1/n1 - h2/n2
        D <- N + h1 + h2
        
        F <- sum(N, na.rm=T)/sum(D, na.rm=T)
        fsts[p2,p1] <- F
        if (verbose == TRUE) {
          print(paste("Pop1: ",p1,", Pop2: ",p2,", Reich FST: ",F,sep=""))
        }
        
        if (bootstrap != FALSE) {
          if (verbose == TRUE) {
            print("beginning bootstrapping")
          }
          
          bs[k,1:3] <- c(p2,p1,as.numeric(F))
          
          for (i in 1:n.bs){
            loci <- sample((1:nloc), nloc, replace=TRUE)
            
            pop1.bs <- matrix(as.matrix(pop1)[,loci],
                              ncol=length(loci))
            pop2.bs <- matrix(as.matrix(pop2)[,loci],
                              ncol=length(loci))
            
            a1 <- colSums2(as.matrix(pop1.bs),na.rm=T)
            a2 <- colSums2(as.matrix(pop2.bs),na.rm=T)
            n1 <- apply(as.matrix(pop1.bs),2,function(x) 2*sum(!is.na(x)))
            n2 <- apply(as.matrix(pop2.bs),2,function(x) 2*sum(!is.na(x)))
            
            h1 <- (a1*(n1-a1))/(n1*(n1-1))
            h2 <- (a2*(n2-a2))/(n2*(n2-1))
            
            N <- (a1/n1 - a2/n2)^2 - h1/n1 - h2/n2
            D <- N + h1 + h2
            
            F.bs <- sum(N, na.rm=T)/sum(D, na.rm=T)
            bs[k,i+5] <- F.bs
          }
          if (verbose == TRUE){
            print(paste("bootstrapping 95% CI: ",
                        quantile(bs[k,6:(n.bs+5)],0.025,na.rm=T),"-",
                        quantile(bs[k,6:(n.bs+5)],0.975,na.rm=T)))
          }
          
          bs[k,4:5] <- c(quantile(bs[k,6:(n.bs+5)],0.025,na.rm=T),
                         quantile(bs[k,6:(n.bs+5)],0.975,na.rm=T))
        }
        
      }
    }
  }
  
  fsts[fsts < 0] <- 0
  
  if (bootstrap != FALSE){
    colnames(bs)[1:5] <- c("pop1","pop2","fst_estimate","min_CI","max_CI")
    fst.list <- list(fsts,bs)
    names(fst.list) <- c("fsts","bootstraps")
    
    if (plot == TRUE){
      print("drawing plot with bootstraps")
      
      if (!require("ggplot2",character.only=T, quietly=T)) {
        install.packages("ggplot2")
        library(ggplot2, character.only=T)
      }
      
      plot.data <- bs[,1:5]
      plot.data$fst_estimate <- as.numeric(plot.data$fst_estimate)
      plot.data$min_CI <- as.numeric(plot.data$min_CI)
      plot.data$max_CI <- as.numeric(plot.data$max_CI)
      plot.data$pop_pair <- paste(plot.data$pop1,plot.data$pop2,sep="_")
      plot.data$signif <- case_when(plot.data$min_CI > 0 ~ TRUE,
                                    TRUE ~ FALSE)
      
      
      bs.plot <- ggplot(plot.data, aes(x=pop_pair,y=fst_estimate,col=signif)) + 
        geom_point(size=2) + 
        coord_flip() + 
        geom_errorbar(aes(ymin=min_CI,ymax=max_CI),width=0.1,size=1) + 
        geom_hline(yintercept=0, lty=2, lwd=1, col="gray50") + 
        theme_minimal() + 
        theme(legend.position="none")
      
      print(bs.plot)
    }
  } else {
    fst.list <- list(fsts)
    names(fst.list) <- "fsts"
    
    if (plot == TRUE){
      print("drawing plot without bootstraps")
      
      if (!require("ggplot2",character.only=T, quietly=T)) {
        install.packages("ggplot2")
        library(ggplot2, character.only=T)
      }
      
      plot.data <- data.frame(combinat::combn2(row.names(fsts)),
                              fst_estimate=fsts[lower.tri(fsts)])
      plot.data$pop_pair <- paste(plot.data$X1,plot.data$X2,sep="_")
      
      fst.plot <- ggplot(plot.data, aes(x=pop_pair,y=fst_estimate)) + 
        geom_point(size=2) + 
        coord_flip() + 
        geom_hline(yintercept=0, lty=2, lwd=1, col="gray50") + 
        theme_minimal() + 
        theme(legend.position="none")
      
      print(fst.plot)
    }
  }
  
  return(fst.list)
  beepr::beep()
}

# Plot FST matrix ---------------------------------------------------#####
plot_lowertri_FST_mat_from_xlsx <- function(FSTxlsx, poporder,
                                            dt1spreadsheetname = "FST_100bs",
                                            dt2spreadsheetname = "pval_100bs",
                                            fontsize = 2.5){
  # This function plots FST matrix with just the values in the lower
  # triangle. p-value > 0.05 is italicised
  require(tidyverse)
  require(xlsx)
  dt1 <- read.xlsx(FSTxlsx,
                   sheetName = dt1spreadsheetname)
  dt2 <- read.xlsx(FSTxlsx,
                   sheetName = dt2spreadsheetname)
  # Clean up dt1 ------------------------------------------------------#####
  colnames(dt1)[1] <- "population"
  colnames(dt1)[-1] <- dt1$population
  rownames(dt1) <- dt1$population
  dt1$population <- NULL
  mat1 <- round(as.matrix(dt1),3)
  # Create symmetric matrix for plotting fst heat map colour/fill -----#####
  # Store values for the lower triangle (FST values)
  val <- mat1[lower.tri(mat1)]
  # transpose matrix
  mat1a <- t(mat1)
  # Add previously stored values from the lower triangle to the current
  # lower triangle which is all NA
  mat1a[lower.tri(mat1a)] <- val
  # Create a list of order based on latitude and longitude from higher up
  # dtlatlon <- gl@other$ind.metrics %>%
  #   group_by(loc_precise_time) %>%
  #   summarise(lat = mean(latitude),
  #             lon = mean(longitude))
  print(poporder)
  mat1a <- mat1a[poporder, poporder]
  # Quick heat map plot 
  # pheatmap(mat1a, display_numbers = T, cluster_rows = F, cluster_cols = F)
  
  ## Convert mat1a (assymetric matrix with FST values in lower triangle)#####
  # Make matrix with just lower triangle 
  mat1b <- mat1a
  mat1b[upper.tri(mat1b)] <- NA
  # pheatmap(mat1b, display_numbers = T, cluster_rows = F, cluster_cols = F)
  dt1b <- data.frame(mat1b, check.names = F)
  # Add population column
  dt1b$population <- rownames(dt1b)
  # Convert data from wide to long
  dat <- dt1b %>%
    pivot_longer(cols = colnames(dt1b)[colnames(dt1b) != "population"],
                 names_to = "population2")
  # Order the populations accordingly by converting the population
  # columns to factors and assigning levels
  dat$population <- factor(dat$population, levels = rev(rownames(mat1b)))
  dat$population2 <- factor(dat$population2, levels = (rownames(mat1b)))
  dat$value <- sprintf("%.3f", dat$value)
  dat$value[dat$value == "NA"] <- NA
  
  # Make pvalue matrix ------------------------------------------------#####
  colnames(dt2)[1] <- "population"
  colnames(dt2)[-1] <- dt2$population
  rownames(dt2) <- dt2$population
  dt2$population <- NULL
  # Exclude populations with less than 5 samples
  dt2 <- dt2[poporder,poporder]
  mat2 <- as.matrix(dt2)
  # convert to symmetric and then assign NA to upper triangle
  val <- mat2[lower.tri(mat2)]
  # transpose matrix
  mat2 <- t(mat2)
  # Add previously stored values from the lower triangle to the current
  # lower triangle which is all NA
  mat2[lower.tri(mat2)] <- val
  mat2 <- mat2[poporder, poporder]
  mat2[upper.tri(mat2)] <- NA
  dt2a <- data.frame(mat2, check.names = F)
  # Add population column
  dt2a$population <- rownames(dt2a)
  # Convert data from wide to long
  datpval <- dt2a %>%
    pivot_longer(cols = colnames(dt2a)[colnames(dt2a) != "population"],
                 names_to = "population2")
  colnames(datpval)[3] <- "pvalue"
  
  # Order the populations accordingly by converting the population
  # columns to factors and assigning levels
  datpval$population <- factor(datpval$population, levels = rev(rownames(mat2)))
  datpval$population2 <- factor(datpval$population2, levels = (rownames(mat2)))
  
  # Merge datpval to dat
  dat <- merge(dat, datpval, by = c("population", "population2"))
  dat$label <- dat$value
  idx1 <- (dat$pvalue <= 0.05 & !is.na(dat$pvalue))
  idx2 <- (dat$pvalue <= 0.01 & !is.na(dat$pvalue))
  dat$label[idx1] <- paste(dat$value[idx1], "*", sep = "")
  dat$label[idx2] <- paste(dat$value[idx2], "**", sep = "")
  dat$italicise <- "plain"
  dat$italicise[dat$pvalue > 0.05 & !is.na(dat$pvalue)] <- "italic"
  # poporder <- poporder[poporder %in% rownames(mat1a)]
  # mat1a <- mat1a[poporder, poporder]
  dat$value <- as.numeric(dat$value)
  # Make plot
  p <- ggplot() +
    geom_tile(data = dat, aes(x = population2, y = population, fill = value)) +                                         # Add values to heatmap
    geom_text(data = dat, aes(x = population2, y = population, label = value, fontface = italicise),
              size = fontsize)+                                         # Add values to heatmap
    # geom_text(data = dat2a, aes(x = population2, y = population, label = value),
    #           size = fontsize) + 
    scale_fill_distiller(palette = "Spectral",
                         direction = 1,
                         # limits = c(-0.001, 0.32),
                         na.value = "white") + 
    theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=0),
          axis.title = element_blank()) +
    scale_y_discrete(position = "right", expand = c(0,0)) +
    scale_x_discrete(expand = c(0,0)) +
    coord_fixed() +
    guides(fill = guide_colourbar(barwidth = 20, barheight = 1)) +
    theme(legend.position = "bottom")
  return(p)
}


plot_FST_mat_from_xlsx <- function(FSTxlsx, poporder,
                                   dt1spreadsheetname = "FST_100bs",
                                   dt2spreadsheetname = "FST_100bs_val",
                                   dt3spreadsheetname = "pval_100bs",
                                   fontsize = 2.5){
  # This function creates ggplot object of FST matrix with 
  # pairwise FST values in the lower triangle, 
  # 95% confidence interval in the upper triangle
  # * = pvalue <= 0.05
  # ** = pvalue <= 0.01
  # Note that scale_fill_distiller may have to be updated to match the 
  # value ranges or kept constant across figures.
  #
  # The values here are -0.01 to 0.25.
  require(tidyverse)
  require(xlsx)
  dt1 <- read.xlsx(FSTxlsx,
                   sheetName = dt1spreadsheetname)
  dt2 <- read.xlsx(FSTxlsx,
                   sheetName = dt2spreadsheetname)
  dt3 <- read.xlsx(FSTxlsx,
                   sheetName = dt3spreadsheetname)
  # Clean up dt1 ------------------------------------------------------#####
  colnames(dt1)[1] <- "population"
  colnames(dt1)[-1] <- dt1$population
  rownames(dt1) <- dt1$population
  dt1$population <- NULL
  colnames(dt1) <- rownames(dt1)
  # dt1 <- dt1[poporder,poporder]
  mat1 <- round(as.matrix(dt1),3)
  
  # Create symmetric matrix for plotting fst heat map colour/fill -----#####
  # Store values for the lower triangle (FST values)
  val <- mat1[lower.tri(mat1)]
  # transpose matrix
  mat1a <- t(mat1)
  # Add previously stored values from the lower triangle to the current
  # lower triangle which is all NA
  mat1a[lower.tri(mat1a)] <- val
  # Create a list of order based on latitude and longitude from higher up
  # dtlatlon <- gl@other$ind.metrics %>%
  #   group_by(loc_precise_time) %>%
  #   summarise(lat = mean(latitude),
  #             lon = mean(longitude))
  poporder <- poporder[poporder %in% rownames(mat1a)]
  mat1a <- mat1a[poporder, poporder]
  # Quick heat map plot 
  # pheatmap(mat1a, display_numbers = T, cluster_rows = F, cluster_cols = F)
  
  # Clean up dt2 ------------------------------------------------------#####
  dt2$CI <- paste(sprintf("%.3f", round(dt2$Lower.bound.CI.limit,3)),
                  sprintf("%.3f", round(dt2$Upper.bound.CI.limit,3)), sep = "\n")
  ## Extract 95% confidence interval based on bootstrap ---------------#####
  dt2a <- dt2 %>% 
    dplyr::select(Population1, Population2, CI) %>%
    tidyr::pivot_wider(names_from = c("Population2"),
                       values_from = c("CI"))
  ### Clean up 95% confidence interval matrix -------------------------#####
  popls <- dt2a$Population1
  dt2a$Population1 <- NA
  colnames(dt2a)
  mat2a <- as.matrix(dt2a)
  mat2a <- rbind(mat2a, matrix(data=NA, ncol=ncol(mat2a), nrow=1))
  colnames(mat2a)[1] <- popls[1]
  # Convert to dataframe as you cannot add row names to tibble, and somehow
  # I cannot do it for this matrix
  mat2a <- as.data.frame(mat2a)
  # Add rownames
  rownames(mat2a) <- colnames(mat2a)
  # Convert to matrix
  mat2a <- as.matrix(mat2a)
  # Make matrix symmetric so that it can be reordered
  val <- mat2a[upper.tri(mat2a)]
  # transpose matrix
  mat2a <- t(mat2a)
  # Add previously stored values from the lower triangle to the current
  # lower triangle which is all NA
  mat2a[upper.tri(mat2a)] <- val
  # Reorder matrix
  mat2a <- mat2a[poporder, poporder]
  # Remove lower triangle
  mat2a[lower.tri(mat2a)] <- NA
  
  # Convert matrix to dataframe for plotting in ggplot ----------------#####
  # Convert dt1a (symmetric matrix with FST values) -------------------#####
  # This is for the "fill" in the matrix
  dt1a <- data.frame(mat1a, check.names = F)
  dt1a$population <- rownames(dt1a)
  # Convert data from wide to long
  dat1a <- dt1a %>%
    tidyr::pivot_longer(cols = colnames(dt1a)[colnames(dt1a) != "population"],
                        names_to = "population2")
  # Order the populatins accoridingly
  dat1a$population <- factor(dat1a$population, levels = rev(rownames(mat1a)))
  dat1a$population2 <- factor(dat1a$population2, levels = (rownames(mat1a)))
  
  ## Convert dt1b (assymetric matrix with FST values in lower triangle)#####
  # Make matrix with just lower triangle 
  mat1b <- mat1a
  mat1b[upper.tri(mat1b)] <- NA
  dt1b <- data.frame(mat1b, check.names = F)
  # Add population column
  dt1b$population <- rownames(dt1b)
  # Convert data from wide to long
  dat <- dt1b %>%
    tidyr::pivot_longer(cols = colnames(dt1b)[colnames(dt1b) != "population"],
                        names_to = "population2")
  # Order the populations accordingly by converting the population
  # columns to factors and assigning levels
  dat$population <- factor(dat$population, levels = rev(rownames(mat1b)))
  dat$population2 <- factor(dat$population2, levels = (rownames(mat1b)))
  # Make quick plot to see 
  # ggplot(dat, aes(x = population2, y = population)) +
  #   geom_tile(aes(fill = value)) +                                         # Add values to heatmap
  #   geom_text(aes(label = round(value, 3)))
  # Convert dt2a (assymetric matrix with FST CI in upper triangle) ----#####
  dt2a <- data.frame(mat2a, check.names = F)
  dt2a$population <- rownames(dt2a)
  dat2a <- dt2a %>%
    tidyr::pivot_longer(cols = colnames(dt2a)[colnames(dt2a) != "population"],
                        names_to = "population2")
  dat2a$population <- factor(dat2a$population, levels = rev(rownames(mat2a)))
  dat2a$population2 <- factor(dat2a$population2, levels = (rownames(mat2a)))
  
  dat$value <- sprintf("%.3f", dat$value)
  dat$value[dat$value == "NA"] <- NA
  
  # Make pvalue matrix ------------------------------------------------#####
  colnames(dt3)[1] <- "population"
  colnames(dt3)[-1] <- dt3$population
  rownames(dt3) <- dt3$population
  dt3$population <- NULL
  mat3 <- as.matrix(dt3)
  # convert to symmetric and then assign NA to upper triangle
  val <- mat3[lower.tri(mat3)]
  # transpose matrix
  mat3 <- t(mat3)
  # Add previously stored values from the lower triangle to the current
  # lower triangle which is all NA
  mat3[lower.tri(mat3)] <- val
  mat3 <- mat3[poporder, poporder]
  mat3[upper.tri(mat3)] <- NA
  dt3a <- data.frame(mat3, check.names = F)
  # Add population column
  dt3a$population <- rownames(dt3a)
  # Convert data from wide to long
  datpval <- dt3a %>%
    tidyr::pivot_longer(cols = colnames(dt3a)[colnames(dt3a) != "population"],
                        names_to = "population2")
  colnames(datpval)[3] <- "pvalue"
  
  # Order the populations accordingly by converting the population
  # columns to factors and assigning levels
  datpval$population <- factor(datpval$population, levels = rev(rownames(mat3)))
  datpval$population2 <- factor(datpval$population2, levels = (rownames(mat3)))
  
  # Merge datpval to dat
  dat <- merge(dat, datpval, by = c("population", "population2"))
  dat$label <- dat$value
  idx1 <- (dat$pvalue <= 0.05 & !is.na(dat$pvalue))
  idx2 <- (dat$pvalue <= 0.01 & !is.na(dat$pvalue))
  dat$label[idx1] <- paste(dat$value[idx1], "*", sep = "")
  dat$label[idx2] <- paste(dat$value[idx2], "**", sep = "")
  #
  # Create a list of order based on latitude and longitude from higher up
  # dtlatlon <- gl@other$ind.metrics %>%
  #   group_by(loc_precise_time) %>%
  #   summarise(lat = mean(latitude),
  #             lon = mean(longitude))
  # poporder <- dtlatlon$loc_precise_time[order(dtlatlon$lon)]
  poporder <- poporder[poporder %in% rownames(mat1a)]
  mat1a <- mat1a[poporder, poporder]
  # Make plot
  p <- ggplot() +
    geom_tile(data = dat1a, aes(x = population2, y = population, fill = value)) +                                         # Add values to heatmap
    geom_text(data = dat, aes(x = population2, y = population, label = label),
              size = fontsize)+                                         # Add values to heatmap
    geom_text(data = dat2a, aes(x = population2, y = population, label = value),
              size = fontsize) + 
    scale_fill_distiller(palette = "Spectral",
                         direction = 1,
                         limits = c(-0.001, 0.25)) + 
    theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=0),
          axis.title = element_blank()) +
    scale_y_discrete(position = "right", expand = c(0,0)) +
    scale_x_discrete(expand = c(0,0)) +
    coord_fixed() +
    guides(fill = guide_colourbar(barwidth = 20, barheight = 1)) +
    theme(legend.position = "bottom")
  return(p)
}

# Subsampling -------------------------------------------------------#####
subsample_or_max_ind <- function(x, n){
  # This function returns the entire vector if the subsampled number is
  # larger than the vector. If the number to subsample is less than the
  # length of the vector, subsample the vector
  if (n > length(x)){
    return(x)
  } else {
    return(sample(x, n))
  }
}

# Making colours ----------------------------------------------------
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
#!/bin/bash -e
#SBATCH --job-name=00_bial_SNPs
#SBATCH --output=00_errors/00_%j.out
#SBATCH --error=00_errors/00_%j.err
#SBATCH --mail-user=kats326@aucklanduni.ac.nz
#SBATCH --mail-type=END
#SBATCH --time=1:00:00
#SBATCH --mem=4G
#SBATCH --ntasks=1
#SBATCH --account=uoa02613

module purge
ml BCFtools/1.16-GCC-11.3.0
ml R/4.2.1-gimkl-2022a
ml VCFtools/0.1.15-GCC-9.2.0-Perl-5.30.1

inputa=variant_calls_annotate.vcf.gz
foldera=/nesi/nobackup/uoa02613/kstuart_projects/Sv10_NZstarlings/analysis/A_SFS/data/
outputa=starlings_bial.vcf.gz

cd ${foldera}

# A bit of an overkill
bcftools view -m2 -M2 -v snps ${inputa} | bcftools view -c1 | bcftools filter -e 'AC==0 || AC==AN' -Oz -o ${foldera}${outputa}

# view -m2 -M2 makes sure that we only keep biallelic SNPs
# view -c1 removes monomorphic SNPs in case the first filter did went through
# -e 'AC==0 || AC=AN" excludes cases when there are no alternate alleles (monomorphic SNPs) or if the number of alternate alleles 
# called are equal to the number of alleles genotypes; essentially another monomorphic SNP filter, in case the -c1 filter did not catch these.

bcftools index ${foldera}${outputa}

VCF=${outputa}
outfile="$(basename ${outputa} .vcf.gz)"
folderb=vcftools_sum_stats/
folderc=${outfile}/
mkdir -p ${folderb}${folderc}
outpath=${folderb}${folderc}
OUT=${outpath}${outfile}

# Allele frequency
vcftools --gzvcf $VCF --freq2 --out $OUT --max-alleles 2
# Mean depth per individual
vcftools --gzvcf $VCF --depth --out $OUT
# Mean depth per site
vcftools --gzvcf $VCF --site-mean-depth --out $OUT
# Mean site quality
vcftools --gzvcf $VCF --site-quality --out $OUT
# Proportion of missing data per individual
vcftools --gzvcf $VCF --missing-indv --out $OUT
# Proportion of missing data per site
vcftools --gzvcf $VCF --missing-site --out $OUT
# Heterozygousity and inbreeding coefficient per individual
vcftools --gzvcf $VCF --het --out $OUT
# Relatedness2
vcftools --gzvcf $VCF --relatedness2 --out $OUT

# Plot the results

cd $outpath
outputpdf=${outfile}.all.stats.pdf

Rscript --vanilla /nesi/nobackup/uoa02613/kstuart_projects/Sv10_NZstarlings/analysis/A_SFS/scripts/01_data_exploration_site_indiv_stats.R $outfile $outputpdf

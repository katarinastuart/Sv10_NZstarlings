#!/bin/bash -e
#SBATCH --job-name=02_bial_SNPs
#SBATCH --output=02_errors/02_%j.out
#SBATCH --error=02_errors/02_%j.err
#SBATCH --mail-user=kats326@aucklanduni.ac.nz
#SBATCH --mail-type=END
#SBATCH --time=6:00:00
#SBATCH --mem=4G
#SBATCH --ntasks=1
#SBATCH --account=uoa02613

module purge
ml R/4.2.1-gimkl-2022a

foldera=/nesi/nobackup/uoa02613/kstuart_projects/Sv10_NZstarlings/analysis/A_SFS/data/filtered_VCF/
outputa=starling_noduprel_qual_miss_filt_bial.vcf.gz

vcfa=${foldera}starling_noduprel_qual_miss_filt_bial.vcf.gz
PopLDfolder=PopLDdecay/
outputLD="$(basename ${vcfa} .vcf.gz)".LDdecay
outpath=${foldera}${PopLDfolder}${outputLD}

cd ${foldera}
mkdir -p ${PopLDfolder}

cd /nesi/nobackup/uoa02613/kstuart_projects/Sv10_NZstarlings/analysis/A_SFS/software/PopLDdecay/bin/

./PopLDdecay -InVCF ${vcfa} -OutStat ${outpath}
perl Plot_OnePop.pl -inFile ${outpath}.stat.gz -output ${outpath}.figure

# in case it does not produce a figure
Rscript ${outpath}.figure.r

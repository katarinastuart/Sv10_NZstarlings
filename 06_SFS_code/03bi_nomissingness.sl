#!/bin/bash -e
#SBATCH --job-name=03bi
#SBATCH --output=03_errors/03bi_%j.out
#SBATCH --error=03_errors/03bi_%j.err
#SBATCH --mail-user=kats326@aucklanduni.ac.nz
#SBATCH --mail-type=END
#SBATCH --time=1:00:00
#SBATCH --mem=4G
#SBATCH --ntasks=1
#SBATCH --account=uoa02613

module purge
ml BCFtools/1.16-GCC-11.3.0
ml VCFtools/0.1.15-GCC-9.2.0-Perl-5.30.1

inputa=starling_noduprel_qual_miss_filt_bial.vcf.gz
foldera=/nesi/nobackup/uoa02613/kstuart_projects/Sv10_NZstarlings/analysis/A_SFS/data/filtered_VCF/
outputa=starling_noduprel_qual_miss_filt_bial_nomissingness.vcf.gz
thinbp=50000

cd ${foldera}

bcftools view -i 'F_MISSING<=0' ${inputa} | bcftools filter -e 'AC==0 || AC=AN' -Oz -o ${outputa}


cd /nesi/nobackup/uoa02613/kstuart_projects/Sv10_NZstarlings/analysis/A_SFS/
mkdir -p software
cd software/
#
# Files will be within a folder
git clone https://github.com/BGI-shenzhen/PopLDdecay.git

cd PopLDdecay/bin/

chmod +x *.pl
chmod +x PopLDdecay

cd /nesi/nobackup/uoa02613/kstuart_projects/Sv10_NZstarlings/analysis/A_SFS/scripts/
cp /nesi/nobackup/uoa02613/A_data/scripts/vcf_filter_and_summary_stats/01_data_exploration_site_indiv_stats.R ./

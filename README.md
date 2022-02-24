# Liver_rat_Multiome
R processing for 10x Genomics Multiome

2 replicates of 3 different conditions of rat liver: Control (CT), Cirrhotic (CH) and Regression (R)

FastQ files:

Replicate 1: ATAC data: /analysisdata/rawseq/fastq/SHARED/000118/AGING_scMultiome/HN00164268/HN00164268_10X_RawData_Outs
              RNA data: /analysisdata/rawseq/fastq/SHARED/000118/AGING_scMultiome/HN00164269/HN00164269_10X_RawData_Outs

Replicate 2 : /analysisdata/rawseq/fastq/SHARED/000121/Rat_liver_multiome_2

STEPS:
1- CREATE NEW RAT GENOME REFERENCE FOR MULTIOME
2- CREATE A CSV with the required data


3- RUN cellranger-arc:

Example: 
cd /work/marina.v/Multiome/
export PATH=/work/marina.v/Multiome/opt/cellranger-arc-2.0.1:$PATH

cellranger-arc count --id=Rat_liver_multiome_CT_1\
			--reference=/work/marina.v/Multiome/mRatBN7-2 \
			--libraries=/work/marina.v/Multiome/CSV/224-225_multiome.csv \
			--localcores=16	





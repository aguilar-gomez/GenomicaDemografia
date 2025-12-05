################### FASTQ ##########################
zcat Bede001_R1_001.fastq.gz |head -n 100 > individual1_R1_subset.fastq

#Download from NCBI
https://www.ncbi.nlm.nih.gov/datasets/taxonomy/9771/
https://www.ncbi.nlm.nih.gov/sra/SRX20092011[accn]
fasterq-dump -o blueWhale6165	 --outdir ./rawreads --mem 12G --split-files SRR24296622	

#fastqc
fastqc individual1_R1_subset.fastq 

#Challenge, trim the adapters and run fastqc again
#HINT: use trimmomatic

################### FASTA ##########################
#Downaload genome from NCBI
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/873/245/GCF_009873245.2_mBalMus1.pri.v3/
wget link to file

#Extract sequence from genome: 
samtools faidx GCF_009873245.2_mBalMus1.pri.v3_genomic.fna NC_045785.1:1-10000 > blueWhale_chr1.fasta

################### GENERATE BAM ##########################

################### CONVERT TO SAM ########################

################### GENERATE VCF ##########################



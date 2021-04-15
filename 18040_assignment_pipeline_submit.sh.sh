#!/bin/bash

# Advanced Bioinformatics Course Assignment 2021. 
# Candidate number: 18040
# Date: 15/05/21

# This pipeline processes NGS read data into a VCF. 
# See comments at each step for more details on the processing and analyses performed. 

##################################################################################################
### 1. Install required tools for the project ####################################################
##################################################################################################

cd ~/

# Install Anaconda 
wget https://repo.anaconda.com/archive/Anaconda3-2020.02-Linux-x86_64.sh
chmod +x ./Anaconda3-2020.02-Linux-x86_64.sh
bash ./Anaconda3-2020.02-Linux-x86_64.sh
source ~/.bashrc

# Install required packages with Anaconda 
conda install samtools bwa freebayes picard bedtools trimmomatic fastqc vcflib


##################################################################################################
### 2. Make initial directory structure ##########################################################
##################################################################################################

mkdir ngs_assignment
cd ngs_assignment
mkdir data meta results logs

# Create a simple README description of the project 
touch README.md 
echo "This project runs a standard Bioinformatics NGS pipeline to perform read \
alignment, variant discovery and annotation on some raw sequencing data." \
> README.md 


##################################################################################################
### 3. Download initial input files ##############################################################
##################################################################################################

cd data
mkdir untrimmed_fastq trimmed_fastq

# download FASTQ raw read data
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R1.fastq.qz
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R2.fastq.qz

mv *fastq.qz untrimmed_fastq

# download BED file
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/annotation.bed

# Convert .qz FASTQ files to more usable .gz files
cd untrimmed_fastq
mv NGS0001.R1.fastq.qz NGS0001.R1.fastq.gz
mv NGS0001.R2.fastq.qz NGS0001.R2.fastq.gz

# Download reference genome FASTA file  
mkdir ../reference
cd ../reference
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
# generate index of reference genome fasta file for alignment step 
bwa index hg19.fa.gz # run this step now as it takes ~45 mins, making troubleshooting of the steps later in the pipeline easier

# Download and setup Annovar database for variant annotation 
cd ~/annovar
tar -zxvf annovar.latest.tar.gz # N.B. this file must be downloaded from the Annovar website as it is not opensource
# Application form link = http://download.openbioinformatics.org/annovar_download_form.php

# setup databases 
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar knownGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar ensGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20180603 humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar exac03 humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp31a_interpro humandb/
# ./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar snp138 humandb/ 
# due to the size of the dbSNP 138 database, I had to download and run this step in my local environment

cd ~/ 


##################################################################################################
### 4. Run FastQC quality assessment on raw sequencing data ######################################
##################################################################################################

workdir=~/ngs_assignment # set workdir variable, enabled me to use my script from the module workshop more easily 

# run FastQC to generate quality metrics on untrimmed reads
fastqc -t 4 $workdir/data/untrimmed_fastq/*.fastq.gz  # note the parameter to run on all 4 CPUs

# Move FastQC results to a new directory 
mkdir $workdir/results/fastqc_untrimmed_reads
mv $workdir/data/untrimmed_fastq/*fastqc* $workdir/results/fastqc_untrimmed_reads/


##################################################################################################
### 5. Perform read trimming with Trimmomatic on raw sequencing data #############################
##################################################################################################

# use Trimmomatic to trim away adapters and filter out poor quality score reads
# Drops reads below 50 bp in length, and trims bases from the end of reads, if below a threshold quality (25).
# Also removes Illumina adapter sequences. 
trimmomatic PE  \
  -threads 4 \
  -phred33 \
  $workdir/data/untrimmed_fastq/NGS0001.R1.fastq.gz $workdir/data/untrimmed_fastq/NGS0001.R2.fastq.gz \
  -baseout $workdir/data/trimmed_fastq/NGS0001_trimmed_R \
  ILLUMINACLIP:/home/ubuntu/anaconda3/pkgs/trimmomatic-0.39-1/share/trimmomatic-0.39-1/adapters/NexteraPE-PE.fa:2:30:10 \
  TRAILING:25 MINLEN:50 

# remove unneeded data to save disk space
rm $workdir/data/untrimmed_fastq/NGS0001.R1.fastq.gz $workdir/data/untrimmed_fastq/NGS0001.R2.fastq.gz


##################################################################################################
### 6. FastQC quality assessment on trimmed sequencing data ######################################
##################################################################################################

# run FastQC to generate quality metrics on trimmed reads
fastqc -t 4 $workdir/data/trimmed_fastq/NGS0001_trimmed_R_1P \
  $workdir/data/trimmed_fastq/NGS0001_trimmed_R_2P

# Move FastQC results to a new directory   
mkdir $workdir/results/fastqc_trimmed_reads
mv $workdir/data/trimmed_fastq/*fastqc* $workdir/results/fastqc_trimmed_reads/


##################################################################################################
### 7.1 Alignment of reads to reference genome ###################################################
##################################################################################################

mkdir $workdir/data/aligned_data
alignment_dir=$workdir/data/aligned_data

# alignment of trimmed reads to  hg19 reference genome, using BWA-MEM algorithm. 
# read group info obtained from the raw read FASTQ files. 
bwa mem -t 4 -v 1 -R '@RG\tID:11V6WR1:111:D1375ACXX:1.NGS0001\tSM:NGS0001\tPL:ILLUMINA\' \
  -I 250,50  $workdir/data/reference/hg19.fa.gz $workdir/data/trimmed_fastq/NGS0001_trimmed_R_1P \
  $workdir/data/trimmed_fastq/NGS0001_trimmed_R_2P > $workdir/data/aligned_data/NGS0001.sam

rm $workdir/data/trimmed_fastq/* # remove uneeded data to save disk space


##################################################################################################
### 7.2 Convert, process and index SAM file ######################################################
##################################################################################################

samtools view -h -b $alignment_dir/NGS0001.sam > $alignment_dir/NGS0001.bam # covert SAM to BAM file
rm $alignment_dir/NGS0001.sam 	# remove unneeded data to save disk space
samtools sort $alignment_dir/NGS0001.bam > $alignment_dir/NGS0001_sorted.bam # sort the BAM file
samtools index $alignment_dir/NGS0001_sorted.bam 	# generate index 


##################################################################################################
### 7.3 Mark duplicate reads #####################################################################
##################################################################################################

# this step locates and marks duplicated reads in the BAM file 
picard MarkDuplicates I=$alignment_dir/NGS0001_sorted.bam O=$alignment_dir/NGS0001_sorted_marked.bam \
  M=$alignment_dir/marked_dup_metrics.txt
samtools index $alignment_dir/NGS0001_sorted_marked.bam # generate index 
rm $alignment_dir/NGS0001_sorted.bam # remove uneeded data to save disk space


##################################################################################################
### 7.4 Filter BAM based on mapping quality and bitwise flags using samtools #####################
##################################################################################################

# Filter reads below a minimum MAPQ score (-q 20). The samtools flag -F 1796 filters reads that are
# unmapped, not in the primary alignment, that fail platform/vendor quality checks, or that are 
# PCR/optical duplicates. For more info see: https://broadinstitute.github.io/picard/explain-flags.html 
samtools view -F 1796  -q 20 -o $alignment_dir/NGS0001_sorted_marked_filtered.bam \
  $alignment_dir/NGS0001_sorted_marked.bam

samtools index $alignment_dir/NGS0001_sorted_marked_filtered.bam # generate index 


##################################################################################################
### 7.5 Generate alignment statistics ############################################################
##################################################################################################

mkdir $alignment_dir/alignment_stats # make a new directory for these stats 
stats_dir=$alignment_dir/alignment_stats

# Generate flagstats 
samtools flagstats $alignment_dir/NGS0001_sorted_marked_filtered.bam \
  > $stats_dir/flagstats_output.txt

# view the BAM file in the command line 
# samtools view -h $alignment_dir/NGS0001_sorted_marked_filtered.bam | less

# Generate alignment statistics per chromosome with idxstats
samtools idxstats  $alignment_dir/NGS0001_sorted_marked_filtered.bam \
  > $stats_dir/idxstats_output.txt

# Determine the distribution of insert sizes between read pairs with Picard tools 
picard CollectInsertSizeMetrics -H $stats_dir/CollectInsertSizeMetrics_histogram.pdf \
	-I $alignment_dir/NGS0001_sorted_marked_filtered.bam -O $stats_dir/CollectInsertSizeMetrics_output.txt

# Calculate Depth of coverage. 
# First get the all BAM regions that overlap with the input, then use this output to calculate the coverage.
# 	(This saves memory usage)
bedtools intersect -bed -a NGS0001_sorted_marked_filtered.bam -b $workdir/data/annotation.bed \
  | bedtools coverage -d -a $workdir/data/annotation.bed -b - \
  > $stats_dir/coverageBed_output.txt

rm $alignment_dir/NGS0001_sorted_marked.bam # remove unneeded data to save disk space


##################################################################################################
# 8.1 Variant calling with Freebayes #############################################################
##################################################################################################

zcat $workdir/data/reference/hg19.fa.gz > $workdir/data/reference/hg19.fa # uncompress reference genome, as this is required by freebayes
samtools faidx $workdir/data/reference/hg19.fa # generate index 

# use FreeBayes to report putative variants in the sequencing data compared to the reference allele
freebayes --bam $alignment_dir/NGS0001_sorted_marked_filtered.bam \
  --fasta-reference $workdir/data/reference/hg19.fa --vcf $workdir/results/NGS0001.vcf

bgzip -f $workdir/results/NGS0001.vcf # compress VCF 
tabix -p vcf $workdir/results/NGS0001.vcf.gz # index VCF 

rm $workdir/data/reference/hg19.fa # remove unneeded data to save disk space


##################################################################################################
### 8.2 Filtering the VCF ########################################################################
##################################################################################################

# remove "bad" variant calls from the VCF. 
# "QUAL > 1" removes horrible sites, "QUAL / AO > 10" additional contribution of each obs should be 
# 	10 log units (~ Q10 per read), "SAF > 0 & SAR > 0" reads on both strands, "RPR > 1 & RPL > 1" at 
# 	least two reads “balanced” to each side of the site. 
vcffilter -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" \
  $workdir/results/NGS0001.vcf.gz > $workdir/results/NGS0001_filtered.vcf

# use bedtools to filter VCF for regions present in the annotation bed file. 
bedtools intersect -header -wa -a $workdir/results/NGS0001_filtered.vcf -b $workdir/data/annotation.bed \
  > $workdir/results/NGS0001_filtered_in_bedfile.vcf

bgzip -f $workdir/results/NGS0001_filtered_in_bedfile.vcf # compress filtered VCF
tabix -p vcf $workdir/results/NGS0001_filtered_in_bedfile.vcf.gz # index filtered VCF 


##################################################################################################
### 9. Variant Annotation and Prioritisation #####################################################
##################################################################################################

# Convert VCF to annovar input 
~/tools/annovar/convert2annovar.pl -format vcf4 $workdir/results/NGS0001_filtered_in_bedfile.vcf.gz  \
  > $workdir/results/NGS0001_filtered_in_bedfile.avinput

# Run Annovar to annotate variants with database frequencies/functional consequences. 
# 	Output is a .csv that can be opened in MS Excel. 
~/tools/annovar/table_annovar.pl $workdir/results/NGS0001_filtered_in_bedfile.avinput ~/tools/annovar/humandb/ -buildver hg19 \
  -out $workdir/results/NGS0001_filtered_in_bedfile -remove \
  -protocol refGene,ensGene,clinvar_20180603,exac03,dbnsfp31a_interpro -operation g,g,f,f,f -otherinfo -nastring . -csvout

# # Line for including dbSNP 138 database (due to size of database I had to do this in my local env)
# ~/tools/annovar/table_annovar.pl $workdir/results/NGS0001_filtered_in_bedfile.avinput ~/tools/annovar/humandb/ -buildver hg19 \
#   -out $workdir/results/NGS0001_filtered_in_bedfile -remove \
#   -protocol refGene,ensGene,clinvar_20180603,exac03,dbnsfp31a_interpro,snp138 -operation g,g,f,f,f,f -otherinfo -nastring . -csvout
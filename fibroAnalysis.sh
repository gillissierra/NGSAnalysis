#!/bin/bash
#PBS -l nodes=1:ppn=8
#PBS -l walltime=2:00:00
#PBS -N fibro
#PBS -o fibroAnalysis.out.log
#PBS -e fibroAnalysis.err.log

cd $PBS_O_WORKDIR

### BUILD DIRECTORY TREE ###########################################################################
# mkdir -p references/fasta references/bwa
mkdir -p fibro/fastq fibro/fastqc fibro/bwa

####################################################################################################
### SET UP REFERENCE FILES #########################################################################
### copy files
#cd references/fasta
#cp /scratch/d/danfldhs/mchan/public/inclass/hg19.fa .

### BUILD BWA INDEX ################################################################################
# link file in bw and bulid index. When aligning data you will point to the chr??.fa file located in the same directory as the index files
cd $SCRATCH/MajorProject/references/bwa
ln -s $SCRATCH/MajorProject/references/fasta/hg19.fa

module load bwakit
bwa index hg19.fa
# this will create multiple index files. you align you need to point to the chr??.fa file in the same directory as the index files

####################################################################################################
### PROCESS DATA ###################################################################################
cd $SCRATCH/MajorProject/fibro/fastq
cp $SCRATCH/MajorProject/sample1_fibro*gz .

### RUN FASTQC #####################################################################################
cd $SCRATCH/MajorProject/fibro/fastqc
# links any gz file in this directory
ln -s $SCRATCH/MajorProject/fibro/fastq/*gz .

# run fastqc
module load use.own 
module load fastqc
# runs fastqc on all gz files in this directory
fastqc *.gz

### RUN BWA ALIGNMENT ##############################################################################
cd $SCRATCH/MajorProject/fibro/bwa

# @RG is read group information about the sample. alter details per sample
# paths to references may need to be changed
bwa mem -M -R "@RG\tID:S1Fibro\tSM:S1Fibro\tPL:Illumina\tPU:Illumina1\tLB:library1" $SCRATCH/MajorProject/references/bwa/hg19.fa $SCRATCH/MajorProject/fibro/fastq/sample1_fibro_1.fastq.gz $SCRATCH/MajorProject/fibro/fastq/sample1_fibro_2.fastq.gz > $SCRATCH/MajorProject/fibro/bwa/S1Fibro.sam

### CONVERT SAM TO BAM FILE ########################################################################
module load java
module load picard-tools
java -jar $SCINET_PICARD_JAR SamFormatConverter I=$SCRATCH/MajorProject/fibro/bwa/S1Fibro.sam O=$SCRATCH/MajorProject/fibro/bwa/S1Fibro.bam

### COORDINATE SORT BAM ############################################################################
java -jar $SCINET_PICARD_JAR SortSam I=$SCRATCH/MajorProject/fibro/bwa/S1Fibro.bam O=$SCRATCH/MajorProject/fibro/bwa/S1Fibro_sorted.bam SO=coordinate

### INDEX BAM FILE #################################################################################
# This is for viewing in IGV but not mandatory for downstream analysis
java -jar $SCINET_PICARD_JAR BuildBamIndex I=$SCRATCH/MajorProject/fibro/bwa/S1Fibro_sorted.bam

### COLLAPSE DATA ##################################################################################
java -jar $SCINET_PICARD_JAR MarkDuplicates I=$SCRATCH/MajorProject/fibro/bwa/S1Fibro_sorted.bam O=$SCRATCH/MajorProject/fibro/bwa/S1Fibro_sorted_dedup.bam M=$SCRATCH/MajorProject/fibro/bwa/S1Fibro_sorted_dedup.metrics

### INDEX BAM FILE #################################################################################
### THIS STEP is required for next steps
java -jar $SCINET_PICARD_JAR BuildBamIndex I=$SCRATCH/MajorProject/fibro/bwa/S1Fibro_sorted_dedup.bam

### BASIC STATS FROM A BAM FILE ###
module load samtools
samtools flagstat $SCRATCH/MajorProject/fibro/bwa/S1Fibro_sorted_dedup.bam > $SCRATCH/MajorProject/fibro/bwa/S1Fibro_sorted_dedup.stats

####################################################################################################
### VARIANT CALLING ################################################################################

### MAKE GATK REFERENCES ###########################################################################
cd $SCRATCH/MajorProject/references  #this path may need to be changed
# mkdir gatk
cd $SCRATCH/MajorProject/references/gatk
ln -s $SCRATCH/MajorProject/references/fasta/hg19.fa

java -jar $SCINET_PICARD_JAR CreateSequenceDictionary R=$SCRATCH/MajorProject/references/gatk/hg19.fa O=$SCRATCH/MajorProject/references/gatk/hg19.dict

#module load samtools;
samtools faidx $SCRATCH/MajorProject/references/gatk/hg19.fa

# cp /scratch/d/danfldhs/mchan/public/reference_files/* .

### GATK ###########################################################################################
cd $SCRATCH/MajorProject/fibro/bwa
mkdir gatk
cd $SCRATCH/MajorProject/fibro/bwa/gatk
module load java
module load GATK

### BASE RECAL ###   The paths may need to be changed depending on your given path
# adjust qualities for systematic errors  
java -jar $SCINET_GATK_JAR -T BaseRecalibrator -R $SCRATCH/MajorProject/references/gatk/hg19.fa -I $SCRATCH/MajorProject/fibro/bwa/S1Fibro_sorted_dedup.bam -L hg19 -knownSites $SCRATCH/MajorProject/references/gatk/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -knownSites $SCRATCH/MajorProject/references/gatk/dbsnp_138.hg19.excluding_sites_after_129.vcf -knownSites $SCRATCH/MajorProject/references/gatk/dbsnp_138.hg19.vcf -o $SCRATCH/MajorProject/fibro/bwa/gatk/recal_data.table

java -jar $SCINET_GATK_JAR -T BaseRecalibrator -R $SCRATCH/MajorProject/references/gatk/hg19.fa -I $SCRATCH/MajorProject/fibro/bwa/S1Fibro_sorted_dedup.bam -L hg19 -BQSR $SCRATCH/MajorProject/fibro/bwa/gatk/recal_data.table -knownSites $SCRATCH/MajorProject/references/gatk/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -knownSites $SCRATCH/MajorProject/references/gatk/dbsnp_138.hg19.vcf -knownSites $SCRATCH/MajorProject/references/gatk/dbsnp_138.hg19.excluding_sites_after_129.vcf -o $SCRATCH/MajorProject/fibro/bwa/gatk/recal_data.table

### PRINT READS ###
java -jar $SCINET_GATK_JAR -T PrintReads -R $SCRATCH/MajorProject/references/gatk/hg19.fa -I $SCRATCH/MajorProject/fibro/bwa/S1Fibro_sorted_dedup.bam -L hg19 -BQSR $SCRATCH/MajorProject/fibro/bwa/gatk/recal_data.table -o $SCRATCH/MajorProject/fibro/bwa/gatkS1Fibro_sorted_dedup_recal.bam

### HAPLOTYPE CALLER ### The paths may need to be changed depending on your given path
# calls variants
java -jar $SCINET_GATK_JAR -T HaplotypeCaller -R $SCRATCH/MajorProject/references/gatk/hg19.fa -I $SCRATCH/MajorProject/fibro/bwa/gatk/S1Fibro_sorted_dedup_recal.bam -L hg19 --genotyping_mode DISCOVERY -o $SCRATCH/MajorProject/fibro/bwa/gatk/S1Fibro_sorted_dedup_recal_raw_variants.vcf

### filters snvs
java -jar $SCINET_GATK_JAR -T SelectVariants -R $SCRATCH/MajorProject/references/gatk/hg19.fa -V $SCRATCH/MajorProject/fibro/bwa/gatk/S1Fibro_sorted_dedup_recal_raw_variants.vcf -selectType SNP -o $SCRATCH/MajorProject/fibro/bwa/gatk/S1Fibro_sorted_dedup_recal_raw_variants_SNPs.vcf 

java -jar $SCINET_GATK_JAR -T VariantFiltration -R $SCRATCH/MajorProject/references/gatk/hg19.fa -V $SCRATCH/MajorProject/fibro/bwa/gatk/S1Fibro_sorted_dedup_recal_raw_variants_SNPs.vcf --filterExpression "QUAL < 100" --filterName "lowqual" -o $SCRATCH/MajorProject/fibro/bwa/gatk/S1Fibro_sorted_dedup_recal_raw_variants_SNPs_filtered.vcf

# filters indels
java -jar $SCINET_GATK_JAR -T SelectVariants -R $SCRATCH/MajorProject/references/gatk/hg19.fa -V $SCRATCH/MajorProject/fibro/bwa/gatk/S1Fibro_sorted_dedup_recal_raw_variants.vcf -selectType INDEL -o $SCRATCH/MajorProject/fibro/bwa/gatk/S1Fibro_sorted_dedup_recal_raw_variants_INDELs.vcf

java -jar $SCINET_GATK_JAR -T VariantFiltration -R $SCRATCH/MajorProject/references/gatk/hg19.fa -V $SCRATCH/MajorProject/fibro/bwa/gatk/S1Fibro_sorted_dedup_recal_raw_variants_INDELs.vcf --filterExpression "QUAL < 100" --filterName "lowqual" -o $SCRATCH/MajorProject/fibro/bwa/gatk/S1Fibro_sorted_dedup_recal_raw_variants_INDELs_filtered.vcf

### ANNOVAR #########################################################################################
# make sure you have annovar module in private modules
module load annovar

table_annovar.pl $SCRATCH/MajorProject/fibro/bwa/gatk/S1Fibro_sorted_dedup_recal_raw_variants_SNPs_filtered.vcf /scratch/d/danfldhs/mchan/src/annovar/annovar/humandb/  -buildver hg19  -out $SCRATCH/MajorProject/fibro/bwa/gatk/S1Fibro_sorted_dedup_recal_raw_variants_SNPs_filtered_annotated.vcf -protocol refGene -operation g -vcfinput

echo $PBS_JOBID > $SCRATCH/MajorProject/fibroAnalysis.out.log

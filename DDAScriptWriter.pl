#!/usr/bin/perl

use strict;
use warnings;

## When running this perl prrm
## go into the general folder you haev your fastq files in and that you want your 
## tree built in as well

#expects input and output file to be passed as arguments when the program is run
my $fastq1  = shift(@ARGV); # first read file
my $fastq2 = shift(@ARGV); # second read file
my $name = shift(@ARGV); #File name base eg S1N
my $dirPath = shift(@ARGV); # directory path; eg sample1/normal which would also have /sample1/refchrm/tumor
my $outFile = shift(@ARGV); # name of your bash script
my $refchrm = shift(@ARGV); # name of your ref chrm beforet the .fa


#opens the file you created for your .sh
open(OUT, ">> $outFile") || die("Error: Could not open OUTPUT file... '$outFile' $!\n"); #append

#the following variable holds the bash script

print OUT "### BUILD DIRECTORY TREE ###########################################################################
mkdir -p references/fasta references/bwa
mkdir -p $dirPath/fastq $dirPath/fastqc $dirPath/bwa

####################################################################################################
### SET UP REFERENCE FILES #########################################################################
### copy files
cd references/fasta
cp /scratch/d/danfldhs/mchan/public/inclass/$refchrm.fa .

### BUILD BWA INDEX ################################################################################
# link file in bw and bulid index. When aligning data you will point to the chr??.fa file located in the same directory as the index files
cd ../bwa
ln -s ../fasta/$refchrm.fa

module load bwakit
bwa index $refchrm.fa
# this will create multiple index files. you align you need to point to the chr??.fa file in the same directory as the index files

####################################################################################################
### PROCESS DATA ###################################################################################
cd ../../$dirPath/fastq
cp /scratch/d/danfldhs/mchan/public/inclass2/*gz .

### RUN FASTQC #####################################################################################
cd ../fastqc
# links any gz file in this directory
ln -s ../fastq/*gz .

# run fastqc
module load use.own 
module load fastqc
# runs fastqc on all gz files in this directory
fastqc *.gz

### RUN BWA ALIGNMENT ##############################################################################
cd ../bwa

# \@RG is read group information about the sample. alter details per sample
# paths to references may need to be changed
bwa mem -M -R \"\@RG\\tID:$name\\tSM:$name\\tPL:Illumina\\tPU:Illumina1\\tLB:library1\" ../../../references/bwa/$refchrm.fa ../fastq/$fastq1 ../fastq/$fastq2 > $name.sam

### CONVERT SAM TO BAM FILE ########################################################################
module load java
module load picard-tools
java -jar \$SCINET_PICARD_JAR SamFormatConverter I=$name.sam O=$name.bam

### COORDINATE SORT BAM ############################################################################
java -jar \$SCINET_PICARD_JAR SortSam I=$name.bam O=$name\_sorted.bam SO=coordinate

### INDEX BAM FILE #################################################################################
# This is for viewing in IGV but not mandatory for downstream analysis
java -jar \$SCINET_PICARD_JAR BuildBamIndex I=$name\_sorted.bam

### COLLAPSE DATA ##################################################################################
java -jar \$SCINET_PICARD_JAR MarkDuplicates I=$name\_sorted.bam O=$name\_sorted_dedup.bam M=$name\_sorted_dedup.metrics

### INDEX BAM FILE #################################################################################
### THIS STEP is required for next steps
java -jar \$SCINET_PICARD_JAR BuildBamIndex I=$name\_sorted_dedup.bam

### BASIC STATS FROM A BAM FILE ###
module load samtools
samtools flagstat $name\_sorted.bam > $name\_sorted.stats

####################################################################################################
### VARIANT CALLING ################################################################################

### MAKE GATK REFERENCES ###########################################################################
cd ../../../references/  #this path may need to be changed
mkdir gatk
cd gatk
ln -s ../fasta/$refchrm.fa

java -jar \$SCINET_PICARD_JAR CreateSequenceDictionary R=$refchrm.fa O=$refchrm.dict

module load samtools;
samtools faidx $refchrm.fa

cp /scratch/d/danfldhs/mchan/public/reference_files/* .

### GATK ###########################################################################################
cd ../../$dirPath/bwa/
mkdir gatk
cd gatk
module load java
module load GATK

### BASE RECAL ###   The paths may need to be changed depending on your given path
# adjust qualities for systematic errors
java -jar \$SCINET_GATK_JAR -T BaseRecalibrator -R ../../../../references/gatk/$refchrm.fa -I ../$name\_sorted_dedup.bam -L $refchrm -knownSites ../../../../references/gatk/Mills_and_1000G_gold_standard.indels.$refchrm.vcf -knownSites ../../../../references/gatk/dbsnp_138.$refchrm.vcf -o recal_data.table

java -jar \$SCINET_GATK_JAR -T BaseRecalibrator -R ../../../../references/gatk/$refchrm.fa -I ../$name\_sorted_dedup.bam -L $refchrm -BQSR recal_data.table -knownSites ../../../../references/gatk/Mills_and_1000G_gold_standard.indels.$refchrm.vcf -knownSites ../../../../references/gatk/dbsnp_138.$refchrm.vcf -o recal_data.table

### PRINT READS ###
java -jar \$SCINET_GATK_JAR -T PrintReads -R ../../../../references/gatk/$refchrm.fa -I ../$name\_sorted_dedup.bam -L $refchrm -BQSR recal_data.table -o $name\_sorted_dedup_recal.bam

### HAPLOTYPE CALLER ### The paths may need to be changed depending on your given path
# calls variants
java -jar \$SCINET_GATK_JAR -T HaplotypeCaller -R ../../../../references/gatk/$refchrm.fa -I $name\_sorted_dedup_recal.bam -L $refchrm --genotyping_mode DISCOVERY -o $name\_sorted_dedup_recal_raw_variants.vcf

### filters snvs
java -jar \$SCINET_GATK_JAR -T SelectVariants -R ../../../../references/gatk/$refchrm.fa -V $name\_sorted_dedup_recal_raw_variants.vcf -selectType SNP -o $name\_sorted_dedup_recal_raw_variants_SNPs.vcf 

java -jar \$SCINET_GATK_JAR -T VariantFiltration -R ../../../../references/gatk/$refchrm.fa -V $name\_sorted_dedup_recal_raw_variants_SNPs.vcf --filterExpression \"QUAL < 100\" --filterName \"lowqual\" -o $name\_sorted_dedup_recal_raw_variants_SNPs_filtered.vcf

# filters indels
java -jar \$SCINET_GATK_JAR -T SelectVariants -R ../../../../references/gatk/$refchrm.fa -V $name\_sorted_dedup_recal_raw_variants.vcf -selectType INDEL -o $name\_sorted_dedup_recal_raw_variants_INDELs.vcf

java -jar \$SCINET_GATK_JAR -T VariantFiltration -R ../../../../references/gatk/$refchrm.fa -V $name\_sorted_dedup_recal_raw_variants_INDELs.vcf --filterExpression \"QUAL < 100\" --filterName \"lowqual\" -o $name\_sorted_dedup_recal_raw_variants_INDELs_filtered.vcf

### ANNOVAR #########################################################################################
# make sure you have annovar module in private modules
module load annovar

table_annovar.pl $name\_sorted_dedup_recal_raw_variants_SNPs_filtered.vcf /scratch/d/danfldhs/mchan/src/annovar/annovar/humandb/  -buildver hg19  -out $name\_sorted_dedup_recal_raw_variants_SNPs_filtered_annotated.vcf -protocol refGene -operation g -vcfinput";

# Run MUTECT2 
#close script
close(OUT) || die "Could not close output file '$outFile' properly... $!";
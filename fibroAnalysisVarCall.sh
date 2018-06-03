#!/bin/bash
#PBS -l nodes=1:ppn=8
#PBS -l walltime=8:00:00
#PBS -N fibroVarCall
#PBS -o fibroAnalysisVarCall.out.log
#PBS -e fibroAnalysisVarCall.err.log

cd $PBS_O_WORKDIR

module load use.own 
####################################################################################################
### VARIANT CALLING ################################################################################

### MAKE GATK REFERENCES ###########################################################################
#cd $SCRATCH/MajorProject/references  #this path may need to be changed
# mkdir gatk
#cd $SCRATCH/MajorProject/references/gatk
#ln -s $SCRATCH/MajorProject/references/fasta/hg19.fa

#java -jar $SCINET_PICARD_JAR CreateSequenceDictionary R=$SCRATCH/MajorProject/references/gatk/hg19.fa O=$SCRATCH/MajorProject/references/gatk/hg19.dict

#module load samtools;
#samtools faidx $SCRATCH/MajorProject/references/gatk/hg19.fa

# cp /scratch/d/danfldhs/mchan/public/reference_files/* .

### GATK ###########################################################################################
cd $SCRATCH/MajorProject/fibro/bwa
mkdir gatk
cd $SCRATCH/MajorProject/fibro/bwa/gatk
module load java
module load GATK

### BASE RECAL ###   The paths may need to be changed depending on your given path
# adjust qualities for systematic errors  
java -jar $SCINET_GATK_JAR -T BaseRecalibrator -R $SCRATCH/MajorProject/references/gatk/hg19.fa -I $SCRATCH/MajorProject/fibro/bwa/S1Fibro_sorted_dedup.bam -knownSites $SCRATCH/MajorProject/references/gatk/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -knownSites $SCRATCH/MajorProject/references/gatk/dbsnp_138.hg19.excluding_sites_after_129.vcf -knownSites $SCRATCH/MajorProject/references/gatk/dbsnp_138.hg19.vcf -o $SCRATCH/MajorProject/fibro/bwa/gatk/recal_data.table

java -jar $SCINET_GATK_JAR -T BaseRecalibrator -R $SCRATCH/MajorProject/references/gatk/hg19.fa -I $SCRATCH/MajorProject/fibro/bwa/S1Fibro_sorted_dedup.bam -BQSR $SCRATCH/MajorProject/fibro/bwa/gatk/recal_data.table -knownSites $SCRATCH/MajorProject/references/gatk/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -knownSites $SCRATCH/MajorProject/references/gatk/dbsnp_138.hg19.vcf -knownSites $SCRATCH/MajorProject/references/gatk/dbsnp_138.hg19.excluding_sites_after_129.vcf -o $SCRATCH/MajorProject/fibro/bwa/gatk/recal_data.table

### PRINT READS ###
java -jar $SCINET_GATK_JAR -T PrintReads -R $SCRATCH/MajorProject/references/gatk/hg19.fa -I $SCRATCH/MajorProject/fibro/bwa/S1Fibro_sorted_dedup.bam -BQSR $SCRATCH/MajorProject/fibro/bwa/gatk/recal_data.table -et NO_ET -K $SCRATCH/MajorProject/gsamembers_broadinstitute.org.key -o $SCRATCH/MajorProject/fibro/bwa/gatk/S1Fibro_sorted_dedup_recal.bam

echo $PBS_JOBID > $SCRATCH/MajorProject/fibroAnalysisVarCall.out.log
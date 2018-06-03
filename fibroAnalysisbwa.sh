#!/bin/bash
#PBS -l nodes=1:ppn=8
#PBS -l walltime=8:00:00
#PBS -N fibrobwa
#PBS -o fibroAnalysisbwa.out.log
#PBS -e fibroAnalysisbwa.err.log

cd $PBS_O_WORKDIR

module load bwakit
module load use.own 
module load java
module load picard-tools

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
samtools flagstat $SCRATCH/MajorProject/fibro/bwa/S1Fibro_sorted.bam > $SCRATCH/MajorProject/fibro/bwa/S1Fibro_sorted.stats

echo $PBS_JOBID > $SCRATCH/MajorProject/fibroAnalysisbwa.out.log
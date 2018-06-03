#!/bin/bash
#PBS -l nodes=1:ppn=8
#PBS -l walltime=8:00:00
#PBS -N fibroSortSam
#PBS -o fibroAnalysisSortSam.out.log
#PBS -e fibroAnalysisSortSam.err.log

cd $PBS_O_WORKDIR

module load bwakit
module load use.own 
module load java
module load picard-tools

### COORDINATE SORT BAM ############################################################################
java -jar $SCINET_PICARD_JAR SortSam I=$SCRATCH/MajorProject/fibro/bwa/S1Fibro.bam O=$SCRATCH/MajorProject/fibro/bwa/S1Fibro_sorted.bam SO=coordinate


echo $PBS_JOBID > $SCRATCH/MajorProject/fibroAnalysisSortSam.out.log
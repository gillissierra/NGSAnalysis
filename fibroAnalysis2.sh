#!/bin/bash
#PBS -l nodes=1:ppn=8
#PBS -l walltime=8:00:00
#PBS -N fibro2
#PBS -o fibroAnalysis2.out.log
#PBS -e fibroAnalysis2.err.log

cd $PBS_O_WORKDIR

### BUILD DIRECTORY TREE ###########################################################################
#mkdir -p references/fasta references/bwa
#mkdir -p fibro/fastq fibro/fastqc fibro/bwa

####################################################################################################
### SET UP REFERENCE FILES #########################################################################
### copy files
#cd references/fasta
#cp /scratch/d/danfldhs/mchan/public/inclass/hg19.fa .

### BUILD BWA INDEX ################################################################################
# link file in bw and bulid index. When aligning data you will point to the chr??.fa file located in the same directory as the index files
#cd $SCRATCH/MajorProject/references/bwa
#ln -s $SCRATCH/MajorProject/references/fasta/hg19.fa

module load bwakit
#bwa index $SCRATCH/MajorProject/references/bwa/hg19.fa
# this will create multiple index files. you align you need to point to the chr??.fa file in the same directory as the index files

####################################################################################################
### PROCESS DATA ###################################################################################
#cd $SCRATCH/MajorProject/fibro/fastq
#cp $SCRATCH/MajorProject/sample1_fibro*gz .

### RUN FASTQC #####################################################################################
#cd $SCRATCH/MajorProject/fibro/fastqc
# links any gz file in this directory
#ln -s $SCRATCH/MajorProject/fibro/fastq/*gz .

# run fastqc
module load use.own 
module load fastqc
# runs fastqc on all gz files in this directory
#fastqc *.gz

### RUN BWA ALIGNMENT ##############################################################################
cd $SCRATCH/MajorProject/fibro/bwa

# @RG is read group information about the sample. alter details per sample
# paths to references may need to be changed
#bwa mem -M -R @RG\tID:S1Fibro\tSM:S1Fibro\tPL:Illumina\tPU:Illumina1\tLB:library1" $SCRATCH/MajorProject/references/bwa/hg19.fa $SCRATCH/MajorProject/fibro/fastq/sample1_fibro_1.fastq.gz $SCRATCH/MajorProject/fibro/fastq/sample1_fibro_2.fastq.gz > $SCRATCH/MajorProject/fibro/bwa/S1Fibro.sam

### CONVERT SAM TO BAM FILE ########################################################################
module load java
module load picard-tools
java -Djava.io.tmpdir=. -jar $SCINET_PICARD_JAR SamFormatConverter I=$SCRATCH/MajorProject/fibro/bwa/S1Fibro.sam O=$SCRATCH/MajorProject/fibro/bwa/S1Fibro.bam


echo $PBS_JOBID > $SCRATCH/MajorProject/fibroAnalysis2.out.log

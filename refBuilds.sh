#!/bin/bash
#PBS -l nodes=1:ppn=8
#PBS -l walltime=6:00:00
#PBS -N refBuild
#PBS -o refBuild.out.log
#PBS -e refBuild.err.log

cd $PBS_O_WORKDIR

### BUILD BWA INDEX 
################################################################################
# link file in bw and bulid index. When aligning data you will point to the chr??.fa file located in the same directory as the index files
#cd $SCRATCH/MajorProject/references/bwa
#ln -s $SCRATCH/MajorProject/references/fasta/hg19.fa

module load bwakit
#bwa index hg19.fa

module load use.own 
module load fastqc
module load java
module load picard-tools

####################################################################################################
### VARIANT CALLING ################################################################################

### MAKE GATK REFERENCES ###########################################################################
cd $SCRATCH/MajorProject/references  #this path may need to be changed
# mkdir gatk
cd $SCRATCH/MajorProject/references/gatk
ln -s $SCRATCH/MajorProject/references/fasta/hg19.fa

java -jar $SCINET_PICARD_JAR CreateSequenceDictionary R=$SCRATCH/MajorProject/references/gatk/hg19.fa O=$SCRATCH/MajorProject/references/gatk/hg19.dict

module load samtools
samtools faidx $SCRATCH/MajorProject/references/gatk/hg19.fa

echo $PBS_JOBID > refBuild.out.log

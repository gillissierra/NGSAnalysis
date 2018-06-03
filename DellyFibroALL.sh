#!/bin/bash
#PBS -l nodes=1:ppn=8
#PBS -l walltime=20:00:00
#PBS -N DellyFibroALL
#PBS -o DellyFibroALL.out.log
#PBS -e DellyFibroALL.err.log

cd $PBS_O_WORKDIR

module load use.own 
module load java
module load picard-tools
module load GATK
module load delly

mkdir $SCRATCH/MajorProject/delly;
cd $SCRATCH/MajorProject/delly;

echo "ALL	tumor" > samplesFibroALL;
echo "S1Fibro	control" >> samplesFibroALL;

delly call -t DEL -o $SCRATCH/MajorProject/delly/FibroALLDEL.bcf -g $SCRATCH/MajorProject/references/fasta/hg19.fa $SCRATCH/MajorProject/ALL_T2/bwa/ALL_sorted_dedup.bam $SCRATCH/MajorProject/fibro/bwa/S1Fibro_sorted_dedup.bam
delly filter -t DEL -f somatic -o $SCRATCH/MajorProject/delly/FibroALLDEL_somatic_filter.bcf  -s samplesFibroALL FibroALLDEL.bcf
delly filter -t DEL -f germline -o $SCRATCH/MajorProject/delly/FibroALLDEL_germline_filter.bcf  -s samplesFibroALL FibroALLDEL.bcf

delly call -t DUP -o $SCRATCH/MajorProject/delly/FibroALLDUP.bcf -g $SCRATCH/MajorProject/references/fasta/hg19.fa $SCRATCH/MajorProject/ALL_T2/bwa/ALL_sorted_dedup.bam $SCRATCH/MajorProject/fibro/bwa/S1Fibro_sorted_dedup.bam
delly filter -t DUP -f somatic -o $SCRATCH/MajorProject/delly/FibroALLDUP_somatic_filter.bcf  -s samplesFibroALL FibroALLDUP.bcf
delly filter -t DUP -f germline -o $SCRATCH/MajorProject/delly/FibroALLDUP_germline_filter.bcf  -s samplesFibroALL FibroALLDUP.bcf

delly call -t INV -o $SCRATCH/MajorProject/delly/FibroALLINV.bcf -g $SCRATCH/MajorProject/references/fasta/hg19.fa $SCRATCH/MajorProject/ALL_T2/bwa/ALL_sorted_dedup.bam $SCRATCH/MajorProject/fibro/bwa/S1Fibro_sorted_dedup.bam
delly filter -t INV -f somatic -o $SCRATCH/MajorProject/delly/FibroALLINV_somatic_filter.bcf  -s samplesFibroALL FibroALLINV.bcf
delly filter -t INV -f germline -o $SCRATCH/MajorProject/delly/FibroALLINV_germline_filter.bcf  -s samplesFibroALL FibroALLINV.bcf

delly call -t TRA -o $SCRATCH/MajorProject/delly/FibroALLTRA.bcf -g $SCRATCH/MajorProject/references/fasta/hg19.fa $SCRATCH/MajorProject/ALL_T2/bwa/ALL_sorted_dedup.bam $SCRATCH/MajorProject/fibro/bwa/S1Fibro_sorted_dedup.bam
delly filter -t TRA -f somatic -o $SCRATCH/MajorProject/delly/FibroALLTRA_somatic_filter.bcf  -s samplesFibroALL FibroALLTRA.bcf
delly filter -t TRA -f germline -o $SCRATCH/MajorProject/delly/FibroALLTRA_germline_filter.bcf  -s samplesFibroALL FibroALLTRA.bcf

delly call -t INS -o $SCRATCH/MajorProject/delly/FibroALLINS.bcf -g $SCRATCH/MajorProject/references/fasta/hg19.fa $SCRATCH/MajorProject/ALL_T2/bwa/ALL_sorted_dedup.bam $SCRATCH/MajorProject/fibro/bwa/S1Fibro_sorted_dedup.bam
delly filter -t INS -f somatic -o $SCRATCH/MajorProject/delly/FibroALLINS_somatic_filter.bcf  -s samplesFibroALL FibroALLINS.bcf
delly filter -t INS -f germline -o $SCRATCH/MajorProject/delly/FibroALLINS_germline_filter.bcf  -s samplesFibroALL FibroALLINS.bcf

echo $PBS_JOBID > $SCRATCH/MajorProject/DellyFibroALL.out.log
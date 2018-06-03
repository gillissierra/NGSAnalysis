#!/bin/bash
#PBS -l nodes=1:ppn=8
#PBS -l walltime=8:00:00
#PBS -N fibroHaplo
#PBS -o fibroAnalysisHaplo.out.log
#PBS -e fibroAnalysisHaplo.err.log
cd $PBS_O_WORKDIR

module load use.own 
cd $SCRATCH/MajorProject/fibro/bwa/gatk
module load java
module load GATK
### HAPLOTYPE CALLER ### The paths may need to be changed depending on your given path
# calls variants
java -jar $SCINET_GATK_JAR -T HaplotypeCaller -R $SCRATCH/MajorProject/references/gatk/hg19.fa -I $SCRATCH/MajorProject/fibro/bwa/gatk/S1Fibro_sorted_dedup_recal.bam --genotyping_mode DISCOVERY -et NO_ET -K $SCRATCH/MajorProject/gsamembers_broadinstitute.org.key -o $SCRATCH/MajorProject/fibro/bwa/gatk/S1Fibro_sorted_dedup_recal_raw_variants.vcf

### filters snvs
java -jar $SCINET_GATK_JAR -T SelectVariants -R $SCRATCH/MajorProject/references/gatk/hg19.fa -V $SCRATCH/MajorProject/fibro/bwa/gatk/S1Fibro_sorted_dedup_recal_raw_variants.vcf -selectType SNP -et NO_ET -K $SCRATCH/MajorProject/gsamembers_broadinstitute.org.key -o $SCRATCH/MajorProject/fibro/bwa/gatk/S1Fibro_sorted_dedup_recal_raw_variants_SNPs.vcf 

java -jar $SCINET_GATK_JAR -T VariantFiltration -R $SCRATCH/MajorProject/references/gatk/hg19.fa -V $SCRATCH/MajorProject/fibro/bwa/gatk/S1Fibro_sorted_dedup_recal_raw_variants_SNPs.vcf --filterExpression "QUAL < 100" --filterName "lowqual" -et NO_ET -K $SCRATCH/MajorProject/gsamembers_broadinstitute.org.key -o $SCRATCH/MajorProject/fibro/bwa/gatk/S1Fibro_sorted_dedup_recal_raw_variants_SNPs_filtered.vcf

# filters indels
java -jar $SCINET_GATK_JAR -T SelectVariants -R $SCRATCH/MajorProject/references/gatk/hg19.fa -V $SCRATCH/MajorProject/fibro/bwa/gatk/S1Fibro_sorted_dedup_recal_raw_variants.vcf -selectType INDEL -et NO_ET -K $SCRATCH/MajorProject/gsamembers_broadinstitute.org.key -o $SCRATCH/MajorProject/fibro/bwa/gatk/S1Fibro_sorted_dedup_recal_raw_variants_INDELs.vcf

java -jar $SCINET_GATK_JAR -T VariantFiltration -R $SCRATCH/MajorProject/references/gatk/hg19.fa -V $SCRATCH/MajorProject/fibro/bwa/gatk/S1Fibro_sorted_dedup_recal_raw_variants_INDELs.vcf --filterExpression "QUAL < 100" --filterName "lowqual" -et NO_ET -K $SCRATCH/MajorProject/gsamembers_broadinstitute.org.key -o $SCRATCH/MajorProject/fibro/bwa/gatk/S1Fibro_sorted_dedup_recal_raw_variants_INDELs_filtered.vcf

### ANNOVAR #########################################################################################
# make sure you have annovar module in private modules
module load annovar

table_annovar.pl $SCRATCH/MajorProject/fibro/bwa/gatk/S1Fibro_sorted_dedup_recal_raw_variants_SNPs_filtered.vcf /scratch/d/danfldhs/mchan/src/annovar/annovar/humandb/  -buildver hg19  -out $SCRATCH/MajorProject/fibro/bwa/gatk/S1Fibro_sorted_dedup_recal_raw_variants_SNPs_filtered_annotated.vcf -protocol refGene -operation g -vcfinput

echo $PBS_JOBID > $SCRATCH/MajorProject/fibroAnalysisHaplo.out.log
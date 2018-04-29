#!/usr/bin/perl

use strict;
use warnings;

# when running this perl program, run it in your "sample1" directory (MajorProject) that holds the normal and tumor directories

#expects input and output file to be passed as arguments when the program is run
my $name1 = shift(@ARGV); #File name base for normal eg S1N
my $name2 = shift(@ARGV); #File name base for tumor eg S1T
my $dirPath1 = shift(@ARGV); # directory path to your first bam file
my $dirPath2 = shift(@ARGV); # directory path to your second bam file
my $outFile = shift(@ARGV); # name of your bash script


#opens the file you created for your .sh
open(OUT, ">> $outFile") || die("Error: Could not open OUTPUT file... '$outFile' $!\n"); #append

#the following variable holds the bash script

print OUT "### MUTECT2 #######################################################################################
#copying the reference reference_files

cd \$SCRATCH/MajorProject/references
cp /scratch/d/danfldhs/mchan/public/reference_files/hg19_cosmic_v54_120711_2.vcf \$SCRATCH/MajorProject/references

cd \$SCRATCH/MajorProject
mkdir mutect2
cd \$SCRATCH/MajorProject/mutect2

java -jar \$SCINET_GATK_JAR -T MuTect2 -R \$SCRATCH/MajorProject/references/gatk/hg19.fa --cosmic \$SCRATCH/MajorProject/references/gatk/hg19_cosmic_v54_120711_2.vcf --dbsnp \$SCRATCH/MajorProject/references/gatk/dbsnp_138.hg19.vcf --input_file:normal \$SCRATCH/MajorProject/$dirPath1/bwa/$name1\_sorted_dedup.bam --input_file:tumor \$SCRATCH/MajorProject/$dirPath2/bwa/$name2\_sorted_dedup.bam -o \$SCRATCH/MajorProject/mutect2/$name1\_$name2\_somatic_variants.vcf

### ANNOVAR ########################################################################################
module load annovar 

table_annovar.pl \$SCRATCH/MajorProject/mutect2/$name1\_$name2\_somatic_variants.vcf /scratch/d/danfldhs/mchan/src/annovar/annovar/humandb/  -buildver hg19  -out \$SCRATCH/MajorProject/mutect2/$name1\_$name2\_somatic_variants_annotated.vcf -protocol refGene  -operation g -vcfinput";

#close script
close(OUT) || die "Could not close output file '$outFile' properly... $!";
#!/bin/bash
set -e; # this is for droping if any error. This should be run from runVEP_env

# this is to activate vep
source /gpfs/projects/bsc40/mschikora/anaconda3/etc/profile.d/conda.sh;
conda activate VarCall_CNV_env;
#conda activate runVEP_env;

# define the parent dir (this depends on where this is run)
mschikoraParentDir="/home/mschikora/samba"
if [ ! -d "$mschikoraParentDir" ]; then mschikoraParentDir="/gpfs/projects/bsc40/mschikora" ; fi

##############################################################################################################################
################## IF YOU WANT TO DOWNLOAD THE GITHUB PIPELINE YOU SHOULD CHANGE FunDir BY THE FULL PATH TO YOUR PIPELINE #####
FUNDIR="$mschikoraParentDir"/scripts/VarCall_CNV_ReadProcessing

## MAKE SURE THAT YOU ALSO CHANGE THIS PATH IN varcall_cnv_pipeline_miki.py, functions.py and run_vep.sh

###############################################################################################################################
##############################################################################################################################

# this script is useful to run vep
#VEP=/gpfs/projects/bsc40/mschikora/anaconda3/envs/runVEP_env/bin/vep
#BEDTOOLS=/gpfs/projects/bsc40/mschikora/anaconda3/envs/runVEP_env/bin/bedtools

VEP=/gpfs/projects/bsc40/mschikora/anaconda3/envs/VarCall_CNV_env/bin/vep
BEDTOOLS=/gpfs/projects/bsc40/mschikora/anaconda3/envs/VarCall_CNV_env/bin/bedtools
TABIX="$mschikoraParentDir"/software/tabix-0.2.6/tabix
CORRECT_BY_GENCODE="$FUNDIR"/vep_correct_ByGeneticCode.py
BGZIP=/gpfs/projects/bsc40/mschikora/anaconda3/envs/VarCall_CNV_env/bin/bgzip

# parse command line arguments
while [[ "$#" > 0 ]]; do case $1 in

	# arguments
	-i|--input_file) input_file="$2"; shift;; # a vcf
	-o|--output_file) output_file="$2"; shift;; # a modified vcf
	-f|--reference_fasta) reference_fasta="$2"; shift;; # the reference genome
	-gf|--gff) gff="$2"; shift;; # the gff
	-mch|--mitochondrial_chromosome) mitochondrial_chromosome="$2"; shift;; # the name of the mitochondrial chromosome
	-mcode|--mitochondrial_code) mitochondrial_code="$2"; shift;; # the code of the mitochondrial genetic code
	-gcode|--gDNA_code) gDNA_code="$2"; shift;; # the genetic code of the nuclear genes 

	# store true values as 1s:
	#-r|--repeat) repeat=1;;

	# debug
	*) echo "Unknown parameter passed: $1"; exit 1;;

esac; shift; done

gff_clean="$gff"_clean.gff; gff_clean_compressed="$gff_clean".gz

if [ ! -s "$gff_clean_compressed".tbi ]; then

	# eliminate strange lines,chromosomes and compress
	$BEDTOOLS sort -i $gff | egrep -v '^#' | egrep -v $'\tchromosome\t' > $gff_clean;
	$BGZIP -c $gff_clean > $gff_clean_compressed;

	# index with tabix
	$TABIX $gff_clean_compressed;

fi;

# run vep
echo "Running VEP";
$VEP --input_file $input_file --format "vcf" --output_file $output_file --fasta $reference_fasta --gff $gff_clean_compressed -v --force_overwrite --tab --fields "Uploaded_variation,Location,Allele,Gene,Feature,Feature_type,Consequence,cDNA_position,CDS_position,Protein_position,Amino_acids,Codons,Existing_variation,Extra";

# connect the protein changes taking into consideration the $mitochondrial_code
$CORRECT_BY_GENCODE -i $output_file --mito_chromosome $mitochondrial_chromosome --mito_code $mitochondrial_code --gDNA_code $gDNA_code;

echo "VEP finished correctly";



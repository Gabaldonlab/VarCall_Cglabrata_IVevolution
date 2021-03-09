#!/gpfs/projects/bsc40/mschikora/anaconda3/envs/VarCall_CNV_env/bin/python

# This is a SNP pipeline modified by Miki, intended to be run with python3. It has to be run with the VarCall_CNV_env (source activate VarCall_CNV_env)
import argparse, os
import pandas as pd
import numpy as np
from argparse import RawTextHelpFormatter
import copy as cp
import pickle
import string
import shutil 
from Bio import SeqIO
import random
import sys
from shutil import copyfile

# define the path of mschikora, were a lot of the dependecies of this pipeline are
mschikoraParentDir = "%s/samba"%(os.getenv("HOME")); # local
if os.path.exists(mschikoraParentDir): is_cluster = False    
else:
    is_cluster = True
    mschikoraParentDir = "/gpfs/projects/bsc40/mschikora"

##############################################################################################################################
################## IF YOU WANT TO DOWNLOAD THE GITHUB PIPELINE YOU SHOULD CHANGE FunDir BY THE FULL PATH TO YOUR PIPELINE #####

FunDir = "%s/scripts/VarCall_CNV_ReadProcessing"%mschikoraParentDir; sys.path.append(FunDir)
## MAKE SURE THAT YOU ALSO CHANGE THIS PATH IN varcall_cnv_pipeline_miki.py, functions.py and run_vep.sh

###############################################################################################################################
##############################################################################################################################

# import functions
print("Importing modules")
import functions as fun

# OTHER PATHS:

# define the path to the conda environmet where all packages are found
EnvDir = "/gpfs/projects/bsc40/mschikora/anaconda3/envs/VarCall_CNV_env"

# packages installed with BioConda into VarCall_CNV_env
bwa = "%s/bin/bwa"%EnvDir
picard = "%s/share/picard-2.18.26-0/picard.jar"%EnvDir
samtools = "%s/bin/samtools"%EnvDir
freebayes = "%s/bin/freebayes"%EnvDir
bcftools = "%s/bin/bcftools"%EnvDir
vcffilter = "%s/bin/vcffilter"%EnvDir
gatk = "%s/bin/gatk"%EnvDir # this is gatk4, and it has to be installed like that
vcfallelicprimitives = "%s/bin/vcfallelicprimitives"%EnvDir
qualimap = "%s/bin/qualimap"%EnvDir
sift4g = "%s/bin/sift4g"%EnvDir
java = "%s/bin/java"%EnvDir
gridss_run = "%s/software/gridss/downloaded_scripts/gridss.sh"%mschikoraParentDir
gridss_jar = "%s/software/gridss/downloaded_scripts/gridss-2.6.0-gridss-jar-with-dependencies.jar"%mschikoraParentDir
analyze_svVCF = "%s/generate_files_from_svVCF.R"%FunDir
ensembl_vep = "%s/run_vep.sh"%FunDir
clove = "%s/software/clove-0.17-jar-with-dependencies.jar"%mschikoraParentDir
#clove = "%s/software/clove-0.16-jar-with-dependencies.jar"%mschikoraParentDir
bcftools = "%s/bin/bcftools"%EnvDir
bgzip = "%s/bin/bgzip"%EnvDir
tabix = "%s/bin/tabix"%EnvDir
bedtools = "%s/bin/bedtools"%EnvDir

description = """
Run VariantCalling and/or CNV pipeline in cluster. It also runs indels and complex variants. It creates four folders under outdir (-o):

    - CNV_results: contains the results of the CNV analysis (gene_to_coverage.tab, gene_to_coverage_regions.tab). The regions file includes the regions +-10000  of the boundaries of the gene.

        Both files contain the following fields, for each of the features in the gff that have a "gene" tag (also includes tRNAs and rRNAs):
            1. median_reads_per_gene   
            2. chromosome      
            3. start   
            4. end     
            5. ID      
            6. fraction_covered_by_MoreThan1read: This is the fraction of the feature that is covered with more than one read. 
            7. relative_copy_number: Relative to the region. It can be 0 (if fraction_covered_by_MoreThan1read<0.1), 2 or 3 (duplication, if the median reads of the gene / median of the medians of the 30 (15 left and 15 right) surrounding genes is >1.9) or 1 otherwise. IMPORTANT: While losses (0) are absolute, the assignment of duplication (2 or 3) depends on the 30 surrounding genes' coverage. This is done to prevent smiley faces in telomeres. You won't see anything with whole-chromome duplications.

    - HaplotypeCaller_ploidy[1,2], freebayes_ploidy[1,2], bcftools_ploidy[1,2]: These contain the results of the three variant callers. The important output here is the *_annotated.tab vcf, which already contains the annotated variants (using ensembl Variant Effect Predictor) according to the given gff. It does not annotate non-gene RNAs.

This program should be run with a conda environment we call VarCall_CNV_env (everything in /gpfs/projects/bsc40/mschikora/anaconda3/envs/VarCall_CNV_env) to avoid problems with the versions. 

This program is tested for ploidies 1 and 2. We are not sure what would happen with higher ploidy.
Java version 8 needed. 
It will write all files under -o     

If the -sift_db option is indicated, it will require to have a sift database to be built for a given gff and genome

If a gff is provided, it will generate a table that integrates the variant annotation and the results of the three programs and, also their FILTERs, into outdir/integrated_variants_*_normalisation.tab, whre '*' is the program used for making the normalisation of the variants (vcflib or gatk)

Example command (runs for some reads against only chromosome A (which is fakely treated as mitochondrial for simplicity at the moment of testing the working of the mitocorrection)):

/gpfs/projects/bsc40/mschikora/scripts/VarCall_CNV_ReadProcessing/varcall_cnv_pipeline_miki.py -r /gpfs/projects/bsc40/mschikora/Cglabrata_antifungals/data/chromosomeA_Cglabrata/chromosomeA.fasta -thr 4 -caller all -o /gpfs/projects/bsc40/mschikora/Cglabrata_antifungals/VarCall/testingVarCallPipeline_chrA -p 1 -glm BOTH -gff /gpfs/projects/bsc40/mschikora/Cglabrata_antifungals/data/chromosomeA_Cglabrata/chromosomeA.gff  -cnv -f1 /gpfs/projects/bsc40/mschikora/Cglabrata_antifungals/data/trimmed_reads/RUN1_CST34_2G_FLZ_R1_trimmed.fq.gz -f2 /gpfs/projects/bsc40/mschikora/Cglabrata_antifungals/data/trimmed_reads/RUN1_CST34_2G_FLZ_R2_trimmed.fq.gz -mchr ChrA_C_glabrata_CBS138 -sv_gridss

Another way to test the pipeline is by running /gpfs/projects/bsc40/mschikora/scripts/VarCall_CNV_ReadProcessing/test_dir/test_pipeline.sh

IMPORTANT: You have to make sure that the gff, has only chromosomes that are in the reference
"""
              
parser = argparse.ArgumentParser(description=description, formatter_class=RawTextHelpFormatter)

# general args
parser.add_argument("-r", "--ref", dest="ref", required=True, help="Reference genome. Has to end with .fasta")
parser.add_argument("-thr", "--threads", dest="threads", default=2, type=int, help="Number of threads, Default: 16")
parser.add_argument("-o", "--outdir", dest="outdir", action="store", required=True, help="Directory where the data will be stored")
parser.add_argument("--replace", dest="replace", action="store_true", help="Replace existing files")
parser.add_argument("--replace_vep_integration", dest="replace_vep_integration", action="store_true", help="Replace existing files of the merging of the VEP output and the vcf files. This is for debugging purposes")
parser.add_argument("-p", "--ploidy", dest="ploidy", default=1, type=int, help="Ploidy, can be 1 or 2")

# alignment args
parser.add_argument("-f1", "--fastq1", dest="fastq1", default=None, help="fastq_1 file. Option required to obtain bam files")
parser.add_argument("-f2", "--fastq2", dest="fastq2", default=None, help="fastq_2 file. Option required to obtain bam files")
parser.add_argument("-sbam", "--sortedbam", dest="sortedbam", default=None, help="The path to the sorted bam file, which should have a bam.bai file in the same dir. This is mutually exclusive with providing reads")
parser.add_argument("--run_qualimap", dest="run_qualimap", action="store_true", help="Run qualimap for quality assessment of bam files")

# variant calling args
parser.add_argument("-caller", "--caller", dest="caller", required=False, default="no", help="SNP caller option to obtain vcf file. options: no/all/HaplotypeCaller/bcftools/freebayes.")
parser.add_argument("-glm", "--genotype_likelihood_model", dest="genotype_likelihood_model", default="BOTH", help="Genotype likelihoods calculation model. Can be SNP, INDEL or BOTH. Necessary for GATK")
parser.add_argument("-c", "--coverage", dest="coverage", default=20, type=int, help="minimum Coverage (int)")
parser.add_argument("-mchr", "--mitochondrial_chromosome", dest="mitochondrial_chromosome", default="mito_C_glabrata_CBS138", type=str, help="The name of the mitochondrial chromosome. This is important if you have mitochondrial proteins for which to annotate the impact of nonsynonymous variants, as the mitochondrial genetic code is different. This should be the same as in the gff. If there is no mitochondria just put no_mitochondria")
parser.add_argument("-mcode", "--mitochondrial_code", dest="mitochondrial_code", default=3, type=int, help="The code of the NCBI mitochondrial genetic code. For yeasts it is 3. You can find the numbers for your species here https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi")
parser.add_argument("-gcode", "--gDNA_code", dest="gDNA_code", default=1, type=int, help="The code of the NCBI gDNA genetic code. You can find the numbers for your species here https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi . For C. albicans it is 12. ")

# CNV args
parser.add_argument("-cnv", dest="run_cnv", action="store_true", default=False, help="Run the CNV pipeline, which outputs the number of copies that each gene has according to coverage. The gene ID's are based on the GFF3 files that have been provided in -gff")
parser.add_argument("--replace_cnv", dest="replace_cnv", action="store_true", help="Replace existing files of the cnv analysis")
parser.add_argument("-gff", "--gff-file", dest="gff", default=None, help="path to the GFF3 annotation of the reference genome. Make sure that the IDs are completely unique for each 'gene' tag. This is necessary for both the CNV analysis (it will look at genes there) and the annotation of the variants.")

# structural variation args
parser.add_argument("-sv_gridss", dest="run_gridss", action="store_true", default=False, help="Run gridss to call structural variation. In the c glabrata genome, it takes between 20 min and ")
parser.add_argument("-sv_gridss_purple_linx", dest="run_sv_gridss_purple_linx", action="store_true", default=False, help="Run the gridss_purple_linx pipeline to visualize and filter events. It requires at least one small variant caller. It requires specifying -f1_ref and f2_ref")
parser.add_argument("-f1_ref", "--fastq1_ref", dest="fastq1_ref", default=None, help="fastq_1 file of a reference sample. This is needed to run the gridss_purple_linx pipeline for visualization")
parser.add_argument("-f2_ref", "--fastq2_ref", dest="fastq2_ref", default=None, help="fastq_2 file of a reference sample. This is needed to run the gridss_purple_linx pipeline for visualization")
parser.add_argument("-kgSV_tab", "--known_genomes_withSV_and_shortReads_table", dest="known_genomes_withSV_and_shortReads_table", default=None, help="A path to a table that has ID, assembly, R1 and R2 fields. These are paths to an assembly and forward and reverse short reads, which are used for defining real structural variation, as oposed to --ref, as well as the performance of the gridss pipeline on these real variation.")

parser.add_argument("--replace_gridss", dest="replace_gridss", action="store_true", help="Replace existing files of the gridss analysis")

# an argument 

opt = parser.parse_args()

# debug commands
if not os.path.isdir(opt.outdir): os.mkdir(opt.outdir)
if not opt.gff is None and fun.file_is_empty(opt.gff): raise ValueError("%s is not a valid gff"%opt.gff)
if opt.run_gridss is True: opt.run_cnv = True

# check that the environment is correct
fun.run_cmd("echo 'This is a check of the environment in which the pipeline is running'; which bedtools")
if "/VarCall_CNV_env/" not in os.__file__: raise ValueError("The conda environment VarCall_CNV_env was not correctly activated")

# define the name that will be used as tag, it is the name of the outdir, without the full path
name_sample = opt.outdir.split("/")[-1]; print("working on %s"%name_sample)


# define files that may be used in many steps of the pipeline
if opt.sortedbam is None:

    bamfile = "%s/aligned_reads.bam"%opt.outdir
    sorted_bam = "%s.sorted"%bamfile
    index_bam = "%s.bai"%sorted_bam

else:

    # debug the fact that you prvided reads and bam. You should just provide one
    if any([not x is None for x in {opt.fastq1, opt.fastq2}]): raise ValueError("You have provided reads and a bam, you should only provide one")

    # get the files
    sorted_bam = opt.sortedbam
    index_bam = "%s.bai"%sorted_bam

# correct the gff file, so that it doesn't have lines starting with # and also add biotype (important for ensembl VEP)
if not opt.gff is None:
    correct_gff = "%s_corrected.gff"%(opt.gff); correct_gff_tmp = "%s_corrected_tmp.gff"%(opt.gff)

    if fun.file_is_empty(correct_gff) or opt.replace is True:
        print("correcting gff")
        correct_gff_cmd = "grep -v '^#' %s > %s"%(opt.gff, correct_gff_tmp); fun.run_cmd(correct_gff_cmd)
        os.rename(correct_gff_tmp, correct_gff)

    # modify gff to add biotype
    gff_with_biotype = "%s_with_biotype.gff"%correct_gff
    if fun.file_is_empty(gff_with_biotype) or opt.replace is True:
        print("adding biotype")

        starting_lines = [line for line in open(correct_gff, "r") if line.startswith("#")]
        df_gff3 = pd.read_csv(correct_gff, skiprows=list(range(len(starting_lines))), sep="\t", names=["chromosome", "source", "type_feature", "start", "end", "score", "strand", "phase", "attributes"])

        def add_biotype(row):
            if "biotype" not in row["attributes"] and "gene_biotype" not in row["attributes"]: row["attributes"] += ";biotype=%s"%row["type_feature"]
            return row["attributes"]

        # add biotype and save
        df_gff3["attributes"] = df_gff3.apply(lambda row: add_biotype(row), axis=1)
        df_gff3.to_csv(gff_with_biotype, sep="\t", index=False, header=False)

#####################################
############# BAM FILE ##############
#####################################

##### YOU NEED TO RUN THE BAM FILE #####

if all([not x is None for x in {opt.fastq1, opt.fastq2}]):

    print("WORKING ON ALIGNMENT")
    # delete bam file to debug

    ############# DEBUG #########
    #fun.remove_file(sorted_bam)
    #fun.remove_file(index_bam)
    ###########################

    fun.run_bwa_mem(opt.fastq1, opt.fastq2, opt.ref, opt.outdir, bamfile, sorted_bam, index_bam, name_sample, threads=opt.threads, replace=opt.replace)


else: print("Warning: No fastq file given, assuming that you provided a bam file")

# check that all the important files exist
if any([fun.file_is_empty(x) for x in {sorted_bam, index_bam}]): raise ValueError("You need the sorted and indexed bam files in ")

#### bamqc
if opt.run_qualimap is True:
    
    bamqc_outdir = "%s/bamqc_out"%opt.outdir
    if fun.file_is_empty("%s/qualimapReport.html"%bamqc_outdir) or opt.replace is True:
        print("Running bamqc to analyze the bam alignment")
        qualimap_std = "%s/std.txt"%bamqc_outdir
        try: bamqc_cmd = "%s bamqc -bam %s -outdir %s -nt %i > %s 2>&1"%(qualimap, sorted_bam, bamqc_outdir, opt.threads, qualimap_std); fun.run_cmd(bamqc_cmd)
        except: print("WARNING: qualimap failed likely due to memory errors, check %s"%qualimap_std)

# First create some files that are important for any program

# Create a reference dictionary
rstrip = opt.ref.split(".")[-1]
dictionary = "%sdict"%(opt.ref.rstrip(rstrip)); tmp_dictionary = "%s.tmp"%dictionary; 
if fun.file_is_empty(dictionary) or opt.replace is True:

    # remove any previously created tmp_file
    if not fun.file_is_empty(tmp_dictionary): os.unlink(tmp_dictionary)

    print("Creating picard dictionary")
    cmd_dict = "%s -jar %s CreateSequenceDictionary R=%s O=%s TRUNCATE_NAMES_AT_WHITESPACE=true"%(java, picard, opt.ref, tmp_dictionary); fun.run_cmd(cmd_dict)   
    os.rename(tmp_dictionary , dictionary)

# Index the reference
if fun.file_is_empty("%s.fai"%opt.ref) or opt.replace is True:
    print ("Indexing the reference...")
    cmd_indexRef = "%s faidx %s"%(samtools, opt.ref); fun.run_cmd(cmd_indexRef) # This creates a .bai file of the reference

#####################################
############### CNV #################
#####################################

if opt.run_cnv:

    print("Starting CNV analysis")

    # make a folder for the CNV anlysis
    cnv_outdir = "%s/CNV_results"%opt.outdir
    if not os.path.isdir(cnv_outdir): os.mkdir(cnv_outdir)

    # get the bed file, and also the one of the regions surrounding each gene
    print(opt.ref)
    bed_file = "%s.bed_index1"%correct_gff; bed_file_regions = fun.extract_BEDofGENES_of_gff3(correct_gff, bed_file, replace=opt.replace, reference=opt.ref)

    # define the interetsing beds
    gene_to_coverage_file = "%s/gene_to_coverage_genes.tab"%cnv_outdir
    gene_to_coverage_file_regions = "%s/gene_to_coverage_regions.tab"%cnv_outdir

    # go through each region of bed file
    for bed, final_coverge_file in [(bed_file, gene_to_coverage_file), (bed_file_regions, gene_to_coverage_file_regions)]: fun.write_coverage_per_gene_mosdepth_and_parallel(sorted_bam, opt.ref, cnv_outdir, bed, final_coverge_file, replace=(opt.replace or opt.replace_cnv))

 
    # write the integrated file
    integrated_coverage_file = "%s/genes_and_regions_coverage.tab"%cnv_outdir; integrated_coverage_file_tmp = "%s.tmp"%integrated_coverage_file
    if fun.file_is_empty(integrated_coverage_file) or opt.replace_cnv is True or opt.replace is True: 

       # integrate in one
        df_genes = pd.read_csv(gene_to_coverage_file, sep="\t")
        df_regions = pd.read_csv(gene_to_coverage_file_regions, sep="\t")
        df_integrated = df_genes.merge(df_regions, on="ID", validate="one_to_one", suffixes=("", "_+-10kb_region"))

        # write
        df_integrated.to_csv(integrated_coverage_file_tmp, sep="\t", header=True, index=False)
        os.rename(integrated_coverage_file_tmp, integrated_coverage_file)

    # remove everyhing that is not the coverage file
    for f in os.listdir(cnv_outdir): 
        if f not in {fun.get_file(gene_to_coverage_file), fun.get_file(gene_to_coverage_file_regions), fun.get_file(integrated_coverage_file)}: fun.delete_file_or_folder("%s/%s"%(cnv_outdir, f))
 
    # In Laia's script, she calculates coverage as the median reads per gene (cov per gene) / mean of the cov per gene across all genes

print("CNV analysis finished")


#####################################
##### STRUCTURAL VARIATION ##########
#####################################

if opt.run_gridss:

    print("Starting structural variation analysis with GRIDSS. This will call and annotate genomic rearrangements. This is not prepared for real work")

    # create the directories
    gridss_outdir = "%s/gridss_output"%opt.outdir
    #fun.delete_folder(gridss_outdir) # this is to debug

    # run pipeline, this has to be done with this if to run the pipeline
    if __name__ == '__main__':


        fun.run_GridssClove_optimising_parameters(sorted_bam, opt.ref, gridss_outdir, replace_covModelObtention=(opt.replace or opt.replace_gridss), threads=opt.threads, replace=(opt.replace or opt.replace_gridss), mitochondrial_chromosome=opt.mitochondrial_chromosome, simulation_types=["uniform"], n_simulated_genomes=3, target_ploidies=["haploid", "diploid_hetero"], range_filtering_benchmark="theoretically_meaningful", coverage=opt.coverage, known_genomes_withSV_and_shortReads_table=opt.known_genomes_withSV_and_shortReads_table, expected_ploidy=opt.ploidy)


    print("structural variation analysis with GRIDSS finished")


#####################################
######## VARIANTCALLING #############
#####################################

# initialize an array of files that have the VCF results filtered
filtered_vcf_results = []

# Go through the callers, creating in outdir a folder with the results of each
if opt.caller == "no": print("Stop. Doing the variant calling is not necessary.")
    
if opt.caller == "HaplotypeCaller" or opt.caller == "all":

    print("RUNNING GATK: HaplotypeCaller")

    # create a folder that will contain the output of VCF
    outdir_gatk = "%s/HaplotypeCaller_ploidy%i_out"%(opt.outdir, opt.ploidy)

    # run gatk and get the filtered filename
    gatk_out_filtered = fun.run_gatk_HaplotypeCaller(outdir_gatk, opt.ref, sorted_bam, opt.ploidy, opt.threads, opt.coverage, replace=opt.replace)

    # keep
    filtered_vcf_results.append(gatk_out_filtered)
    
    print("HaplotypeCaller is done")

if opt.caller == "bcftools" or opt.caller == "all":

    print("RUNNING bcftools")

    # create a folder that will contain the output of VCF
    outdir_bcftools = "%s/bcftools_ploidy%i_out"%(opt.outdir, opt.ploidy)
    if not os.path.isdir(outdir_bcftools): os.mkdir(outdir_bcftools)

    # look for the mpileup bcf in sister directories, as it is the same for any other ploidy
    mpileup_output = "%s/output.mpileup.bcf"%outdir_bcftools; mpileup_output_tmp = "%s.tmp"%mpileup_output
    for folder in os.listdir(opt.outdir):
        if folder.startswith("bcftools_ploidy") and folder.endswith("_out"):

            # look for the potential previously calculated mpielup outputs
            potential_previosuly_calculated_mpileup_output = "%s/%s/output.mpileup.bcf"%(opt.outdir, folder)
            if not fun.file_is_empty(potential_previosuly_calculated_mpileup_output): 
                print("taking %s from previous run"%potential_previosuly_calculated_mpileup_output)
                mpileup_output = potential_previosuly_calculated_mpileup_output; break

    # if there is no previous run
    if fun.file_is_empty(mpileup_output) or opt.replace is True:

        print("Running mpileup...")
        cmd_bcftools_mpileup = '%s mpileup -a "AD,DP" -O b -f %s -o %s --threads %i %s'%(bcftools, opt.ref, mpileup_output_tmp, opt.threads, sorted_bam); fun.run_cmd(cmd_bcftools_mpileup)
        os.rename(mpileup_output_tmp, mpileup_output)


    # run bcftools call
    call_output = "%s/output.raw.vcf"%outdir_bcftools; call_output_tmp = "%s.tmp"%call_output
    if fun.file_is_empty(call_output) or opt.replace is True:
        print("Running bcftools call ...")

        # define the ploidy specification
        if opt.ploidy==1: ploidy_cmd = "--ploidy %i"%opt.ploidy # This is all haploid
        else:
            # create a ploidy file if ploidy is 2. There's no way to simpli specify ploidy 2
            ploidy_file_bcftools = "%s/ploidy_file.tab"%outdir_bcftools
            open(ploidy_file_bcftools, "w").write("* * * * %i\n"%opt.ploidy) # CHROM, FROM, TO, SEX, PLOIDY

            ploidy_cmd = "--ploidy-file %s"%ploidy_file_bcftools

        cmd_bcftools_call = "%s call -m -f GQ,GP -v -O v --threads %i -o %s %s %s"%(bcftools, opt.threads, call_output_tmp, ploidy_cmd, mpileup_output); fun.run_cmd(cmd_bcftools_call)

        os.rename(call_output_tmp, call_output)
  
    #As there are no recommendations for bcftools, we decided to apply exclusively the filter for coverage. To apply harder filters please edit this command!
    
    # this generates a filtered, vcf, which only has the PASS ones.
    filtered_output = "%s/output.filt.vcf"%outdir_bcftools; filtered_output_tmp = "%s.tmp"%filtered_output
    if fun.file_is_empty(filtered_output) or opt.replace is True:
        print("Filtering bcftools ... ")
        cmd_filter = "%s filter -m x -e 'INFO/DP <= %i' -O v --threads %i -o %s %s"%(bcftools, opt.coverage, opt.threads, filtered_output_tmp, call_output); fun.run_cmd(cmd_filter)
        os.rename(filtered_output_tmp, filtered_output)

    # keep
    filtered_vcf_results.append(filtered_output)

    print("bcftools is done")

if opt.caller == "freebayes" or opt.caller == "all":

    print("RUNNING freebayes")

    # create a folder that will contain the output of VCF
    outdir_freebayes = "%s/freebayes_ploidy%i_out"%(opt.outdir, opt.ploidy)

    # run freebayes
    freebayes_filtered =  fun.run_freebayes_parallel(outdir_freebayes, opt.ref, sorted_bam, opt.ploidy, opt.coverage, replace=opt.replace) 

    # keep
    filtered_vcf_results.append(freebayes_filtered)
    
    print("freebayes is done")

# define the variants called by all the 
#variants_consensus_all_callers = 

#%s/variantsPASS.vcf"%cwd

##########

##############
# NORMALISATION OF THE VARIANTS. IT IS IMPORTANT TO REPRESENT THE VARIANTS IN THE PRIMITIVE FORM, AS IT ALLOWS TO COMPARE INDELS FROM SEVERAL PROGRAMS.
##############

print("Performing variant normalisation. For in/del variants, several programs may yield various variant representations. This is performed here with VCFLIB and GATK")

# first modify the vcfs so that the ref and alt alleles are uppercase. This is important for the vcflib normalisation
print("Puting all REF and ALT alleles to uppercase.")
for unnormalised_vcf in filtered_vcf_results:

    # load into df
    initial_lines_list = [line for line in open(unnormalised_vcf, "r", encoding='utf-8', errors='ignore') if line.startswith("##")]
    vcf_df = pd.read_csv(unnormalised_vcf, skiprows=list(range(len(initial_lines_list))), sep="\t", na_values=fun.vcf_strings_as_NaNs, keep_default_na=False)

    # put to uppercase
    vcf_df["REF"]  = vcf_df["REF"].apply(lambda x: x.upper())
    vcf_df["ALT"]  = vcf_df["ALT"].apply(lambda x: x.upper())

    # write to the same file, including the initial lines of the vcf for consistency
    open(unnormalised_vcf, "w").write("".join(initial_lines_list) + vcf_df.to_csv(sep="\t", index=False))

# initialize an array that will keep the path to the normalised VCFS
all_normalised_vcfs = set()

# normalise all the filtered_vcf_results with vcfallelicprimitives
for unnormalised_vcf in filtered_vcf_results:

    # define the normalised output
    folder = "/".join(unnormalised_vcf.split("/")[0:-1])
    normalised_vcf = "%s/output.filt.norm_vcflib.vcf"%folder; normalised_vcf_tmp = "%s.tmp"%normalised_vcf

    # generate an unifyed representation of the vcfs
    if fun.file_is_empty(normalised_vcf) or opt.replace is True:
        print("Running vcfallelicprimitives for vcf %s"%unnormalised_vcf)
        cmd_normalise = "%s --keep-geno %s > %s"%(vcfallelicprimitives, unnormalised_vcf, normalised_vcf_tmp); fun.run_cmd(cmd_normalise)
        os.rename(normalised_vcf_tmp, normalised_vcf)

    # keep
    all_normalised_vcfs.add(normalised_vcf)

print("VCFLIB Normalisation is done")



############################
# ANNOTATE VARIANTS WITH VEP
############################

if opt.gff is None: print("No gff provided. Skipping the annotation AND integration of the variants")

else:

    # create a dictionary with [typeNormalisation][sofware] = vep_df
    normalisation_to_software_to_vepDf = {}

    # go through each of the vcfs
    for normalised_vcf in all_normalised_vcfs:

        # get names
        software = normalised_vcf.split("/")[-2].split("_")[0]
        typeNormalisation = normalised_vcf.split("/")[-1].split(".")[2]
        #if software not in {"freebayes"}: continue # DEBUUUG

        # check if any of the integrated datasets are already created
        fileprefix = "%s/integrated_variants_%s_ploidy%i"%(opt.outdir, typeNormalisation, opt.ploidy)
        if any([fun.file_is_empty("%s.%s"%(fileprefix, x)) for x in ["py", "tab"]]) or opt.replace is True or opt.replace_vep_integration is True:

            print("Annotating variants with vep for %s"%normalised_vcf)

            # filter
            print("Loading vcf to calculate the allele frequencies")
            print(normalised_vcf)
            vcf_df , variant_to_frequency, variant_to_filter, var_to_GT, var_to_filters = fun.load_vcf_intoDF_GettingFreq_AndFilter(normalised_vcf)

            vcf_df.to_csv(normalised_vcf, sep="\t", index=False)

            # define an output file for VEP
            annotated_vcf = "%s_annotated.tab"%normalised_vcf; annotated_vcf_tmp = "%s.tmp"%annotated_vcf

            # run annotation by VEP
            if fun.file_is_empty(annotated_vcf) or opt.replace is True or opt.replace_vep_integration is True:
            #if True: # DEBUG
                print("Annotating with VEP %s"%normalised_vcf)
                fun.remove_file(annotated_vcf)
                fun.remove_file(annotated_vcf_tmp)

                vep_cmd = "%s -i %s -o %s -f %s -gf %s -mch %s -mcode %i -gcode %i"%(ensembl_vep, normalised_vcf, annotated_vcf_tmp, opt.ref, gff_with_biotype, opt.mitochondrial_chromosome, opt.mitochondrial_code, opt.gDNA_code); fun.run_cmd(vep_cmd)
                os.rename(annotated_vcf_tmp, annotated_vcf)

            # add the allele frequency of each variant, as calculated in alternative_allelle_frequencies

            print("getting vep table")
            vep_df = fun.load_vep_table_intoDF(annotated_vcf)
            vep_df["fraction_reads_coveringThisVariant"] = [variant_to_frequency[idx] for idx in vep_df.index]
            vep_df["FILTERtag"] = [variant_to_filter[idx] for idx in vep_df.index]
            vep_df["GT"] = [var_to_GT[idx] for idx in vep_df.index]
            vep_df["additional_filters"] = [var_to_filters[idx] for idx in vep_df.index]

            annotated_vcf_corrected = "%s.corrected"%annotated_vcf
            vep_df.to_csv(annotated_vcf_corrected, sep="\t", index=False)

            # add to the dictionary
            normalisation_to_software_to_vepDf.setdefault(typeNormalisation, {}).setdefault(software, vep_df)

    # generate integrated table with each software and also the FILTERtag
    for norm, software_to_vepDF in normalisation_to_software_to_vepDf.items():
        print("Working on, ", norm)

        # define the fileprefix and continur if it has not been generated
        fileprefix = "%s/integrated_variants_%s_ploidy%i"%(opt.outdir, norm, opt.ploidy)

        if any([fun.file_is_empty("%s.%s"%(fileprefix, x)) for x in ["py", "tab"]]) or opt.replace is True or opt.replace_vep_integration is True:
            print("stacking dfs...")

            # stack all the dfs with the common parts, indicating the software
            df = pd.DataFrame()
            for software, vepDF in software_to_vepDF.items(): df = df.append(vepDF[['#Uploaded_variation', 'Location', 'Allele', 'Gene', 'Feature', 'Feature_type', 'Consequence', 'cDNA_position', 'CDS_position', 'Protein_position', 'Amino_acids', 'Codons', 'Existing_variation', 'Extra', 'chromosome', 'position', 'ref', 'alt']])

            # remove the duplicated entries
            df["is_duplicate"] = df.duplicated(subset=['chromosome', 'position', 'ref', 'alt', 'Gene'], keep="first") # the first occurrence is False
            df = df[~(df.is_duplicate)] 
            df["var"] = df.index

            # check that the duplication also applies to the variant + gene combintaion
            if sum(df["is_duplicate"])!=sum(df.duplicated(subset=['#Uploaded_variation', 'Location', 'Allele', 'Gene', 'Feature', 'Feature_type', 'Consequence', 'cDNA_position', 'CDS_position', 'Protein_position', 'Amino_acids', 'Codons', 'Existing_variation', 'Extra', 'chromosome', 'position', 'ref', 'alt'], keep="first")): raise ValueError("The variants are not equaly annotated (with VEP) in all the callers")

            # add which programs called the variants and with the tags
            software_to_setVars = {software : set(vepDF.index) for software, vepDF in software_to_vepDF.items()}
            for software, setVars in software_to_setVars.items(): 
                print("Working on %s"%software)

                # add the sole calling
                df["%s_called"%software] = df["var"].apply(lambda x: x in setVars)

                # map each of the vars to a tag
                var_to_tag = {**dict(zip(software_to_vepDF[software].index, software_to_vepDF[software]["FILTERtag"])), 
                              **{var : "" for var in df[~df["%s_called"%software]].index}}

                # map the original uploaded variation and the GT filter
                var_to_GT_index = {**dict(zip(software_to_vepDF[software].index, software_to_vepDF[software]["GT_index"])), 
                    **{var : "" for var in df[~df["%s_called"%software]].index}}

                var_to_OriginalUploadedVar = {**dict(zip(software_to_vepDF[software].index, software_to_vepDF[software]["#Uploaded_variation_original"])), 
                    **{var : "" for var in df[~df["%s_called"%software]].index}}

                # map each of the vars to the frequency of calling
                var_to_freq = {**dict(zip(software_to_vepDF[software].index, software_to_vepDF[software]["fraction_reads_coveringThisVariant"])), 
                               **{var : 0.0 for var in df[~df["%s_called"%software]].index}}

                # map each of the vars to the GT
                var_to_GT = {**dict(zip(software_to_vepDF[software].index, software_to_vepDF[software]["GT"])), 
                               **{var : "" for var in df[~df["%s_called"%software]].index}}

                # map each of the vars to the additional_filters
                var_to_filters = {**dict(zip(software_to_vepDF[software].index, software_to_vepDF[software]["additional_filters"])), 
                               **{var : "" for var in df[~df["%s_called"%software]].index}}


                # add to df and define if True
                df["%s_FILTERtag"%software] = df["var"].apply(lambda x: var_to_tag[x])
                df["%s_fractionReadsCoveringThisVariant"%software] = df["var"].apply(lambda x: var_to_freq[x])
                df["%s_PASS"%software] = df["%s_FILTERtag"%software].apply(lambda x: x=="PASS")
                df["%s_GT"%software] = df["var"].apply(lambda x: var_to_GT[x])
                df["%s_additional_filters"%software] = df["var"].apply(lambda x: var_to_filters[x])
                df["%s_GT_index"%software] = df["var"].apply(lambda x: var_to_GT_index[x])
                df["%s_#Uploaded_variation_original"%software] = df["var"].apply(lambda x: var_to_OriginalUploadedVar[x])

                # check that all the variants that have been called by this software also have a fractionReadsCoveringThisVariant
                if sum(df["%s_fractionReadsCoveringThisVariant"%software]>0.0) <= sum(df["%s_called"%software])*0.95: 
                    raise ValueError("In %s, %s not almost all the variants that have been called also have a fraction of reads covering them. This may reflect an error in the way how read frequencies are calculated"%(name_sample, software))

            # write and save object
            df.to_csv("%s.tab.tmp"%fileprefix, sep="\t", index=False)
            fun.save_object(df, "%s.py.tmp"%fileprefix)
            os.rename("%s.tab.tmp"%fileprefix, "%s.tab"%fileprefix)
            os.rename("%s.py.tmp"%fileprefix, "%s.py"%fileprefix)

            print("files saved into %s"%("%s.tab"%fileprefix))

print("Variant annotation with VEP is done")


#### OUTPUT THE ALTERNATIVE GENOME IN FASTA FORMAT FOR SMALL VARIANTS (only use the vcflib) ####

alternative_genome = "%s/alternative_genome_vcflib_ploidy%i.fasta"%(opt.outdir, opt.ploidy); alternative_genome_tmp = "%s.tmp"%alternative_genome

# get the vcflib normailsed vcfs
vcflib_normailsed_vcfs = set([x for x in all_normalised_vcfs if "vcflib" in x])

if fun.file_is_empty(alternative_genome) or opt.replace is True:

    print("Getting alternative genome from consensus from vcflib")

    # get the alternative genome
    if len(vcflib_normailsed_vcfs)>0: fun.get_alternative_genome_from_intersetcion_VCFs(opt.outdir, vcflib_normailsed_vcfs, opt.ref, alternative_genome, threads=opt.threads, replace=opt.replace)

#########################################################################



####################################
# PREDICT IMPACT ON PROTEIN
###################################

# Use SIFT and polyphen to predict the impact of each of the vcfs

##### SIFT ######
if opt.gff is None: print("No gff provided. Skipping the analysis of SIFT")

else:

    # go through each of the vcfs and predict
    for normalised_vcf in all_normalised_vcfs:

        # get names
        software = normalised_vcf.split("/")[-2].split("_")[0]
        typeNormalisation = normalised_vcf.split("/")[-1].split(".")[2]

        
#################


print("VarCall Finished")

# create outfile
open("%s/finsihed_file_ploidy_good%i.txt"%(opt.outdir, opt.ploidy), "w").write("finsihed with pipeline\n")


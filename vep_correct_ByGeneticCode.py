#!/gpfs/projects/bsc40/mschikora/anaconda3/envs/VarCall_CNV_env/bin/python

# It used to be called from /gpfs/projects/bsc40/mschikora/anaconda3/envs/runVEP_env/bin/python

import argparse, os
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio import Seq 
import Bio.Seq as Seq
import Bio.SeqRecord as SeqRecord
import copy as cp
import itertools


### ARGUMENTS ###

parser = argparse.ArgumentParser(description="Takes an infile VEP file and writes the same one, but replacing the annotation of the mitochondrial genes, as indicated by those genes that are uin mito_chromosome and also the ones in the nuclear genome, as indicated by gDNA_code")

# general args
parser.add_argument("-i", "--infile", dest="infile", required=True, type=str, help="The input file that will be modified")
parser.add_argument("-mt", "--mito_chromosome", dest="mito_chromosome", required=False, default="mito_C_glabrata_CBS138", type=str, help="The chromosome")
parser.add_argument("-mc", "--mito_code", dest="mito_code", required=False, default=3, type=int, help="The code of translation of ncbi of mitochondrial proteins. Fungal mitochondrial by default.")
parser.add_argument("-gc", "--gDNA_code", dest="gDNA_code", required=False, default=1, type=int, help="The code of translation of ncbi of nuclear genes. Standard by default. C. albicans has 12")

opt = parser.parse_args()

# print the command line to run this
print("Correcting proteins for the genetic_code")

# define the strings that have to be considered as NaN in the VCF parsing
vcf_strings_as_NaNs = ['', '#N/A', '#N/A N/A', '#NA', '-1.#IND', '-1.#QNAN', '-NaN', '-nan', '1.#IND', '1.#QNAN', 'N/A', 'NULL', 'NaN', 'n/a', 'nan', 'null']

### FUNCTIONS ###

def chunks(l, n):
    
    """Yield successive n-sized chunks from a list l"""
    
    for i in range(0, len(l), n):
        yield l[i:i + n]

def get_aa(codons, genetic_code):

    """Takes a string of codons and returns the AA, according to the genetic_code"""

    # if there are no codons
    if codons=="-": return "-"

    # if it is a protein
    elif len(codons)%3==0: return str(Seq.Seq(codons).translate(table = genetic_code))

    # if not
    else: 
        if len(codons)<3: return "X"
        else: return (str(Seq.Seq("".join(list(chunks(codons, 3))[0:-1])).translate(table = genetic_code)) + "X")

def modify_DF_cols(row, genetic_code):

    """Takes a row of a VEP df according to the genetic_code and modifies the necessary rows"""

    if row["Codons"]=="-": codons_ref = "-"; codons_alt = "-"

    else:
        # get the codons
        codons_ref, codons_alt = row["Codons"].split("/")
        codons_ref = codons_ref.upper(); codons_alt = codons_alt.upper(); 

    # get refs
    aa_ref = get_aa(codons_ref, genetic_code); aa_alt = get_aa(codons_alt, genetic_code);

    # define the consequences
    consequences = set()

    # retainig stops
    if codons_ref in stop_codons and codons_alt in stop_codons: consequences.add('stop_retained_variant')

    # indels
    if (len(aa_ref.strip("-"))>len(aa_alt.strip("-")) and len(codons_ref)%3==0 and len(codons_alt)%3==0) or (aa_ref not in {"*", "X"} and aa_alt=="-"): consequences.add('inframe_deletion')
    if (len(aa_ref.strip("-"))<len(aa_alt.strip("-")) and len(codons_ref)%3==0 and len(codons_alt)%3==0) or (aa_alt not in {"*", "X"} and aa_ref=="-"): consequences.add('inframe_insertion')

    # frameshift
    if "X" in aa_alt: consequences.add('frameshift_variant')

    # SNPs
    if len(codons_ref)==3 and len(codons_alt)==3 and aa_ref!="*" and aa_alt!="*":
        if aa_ref==aa_alt: consequences.add('synonymous_variant')
        else: consequences.add('missense_variant')

    # stop codon info
    nStops_ref = aa_ref.count("*"); nStops_alt = aa_alt.count("*"); 
    if nStops_ref<nStops_alt: consequences.add('stop_gained')
    if nStops_ref>nStops_alt: consequences.add('stop_lost')

    # protein altering. This is equivalent to the thing of the row already. 
    if "protein_altering_variant" in row["Consequence"]: consequences.add('protein_altering_variant')

    # add any consequence that is not in any of the genCode_affected_vars
    consequences.update(set(row["Consequence"].split(",")).difference(genCode_affected_vars))

    # return the aa and the consequences
    return pd.Series({"Amino_acids": "%s/%s"%(aa_ref, aa_alt), "Consequence": ",".join(list(consequences))})

# get into df
vep_df = pd.read_csv(opt.infile, sep="\t", header=len([x for x in open(opt.infile, "r") if x.startswith("##")]), na_values=vcf_strings_as_NaNs, keep_default_na=False)

# define the types of variants that are affected by the genetic code
genCode_affected_vars = {'stop_retained_variant', 'inframe_deletion', 'inframe_insertion', 'frameshift_variant', 'synonymous_variant', 'missense_variant', 'stop_gained', 'stop_lost', 'protein_altering_variant'}

# define the idxs of each type of genes
typeGenes_to_idx = {"mito": vep_df.apply(lambda row: row["Location"].startswith(opt.mito_chromosome) and len(set(row["Consequence"].split(",")).intersection(genCode_affected_vars))>0, axis=1),
                    "nuclear": vep_df.apply(lambda row: row["Location"].startswith(opt.mito_chromosome) is False and len(set(row["Consequence"].split(",")).intersection(genCode_affected_vars))>0, axis=1)}

# define the code of each type of genes
typeGenes_to_code = {"mito":opt.mito_code, "nuclear":opt.gDNA_code}

# define a dataframe of non-affected rows in the VEP df
all_df = vep_df[~(typeGenes_to_idx["mito"]) & ~(typeGenes_to_idx["nuclear"])]

# go through each type of genes
for typeGenes, idx_affected_rows in typeGenes_to_idx.items():

    # define the affected df
    affected_df = vep_df[idx_affected_rows]    

    # define the genCode
    genCode = typeGenes_to_code[typeGenes]

    # define the stop codons
    stop_codons = set([codon for codon in ["".join(s) for s in itertools.product(["A", "C", "T", "G"], repeat=3)] if str(Seq.Seq(codon).translate(table = genCode))=="*"])

    # modify the rows if there's something to modify
    if len(affected_df)>0: affected_df[["Amino_acids","Consequence"]] = affected_df.apply(lambda r: modify_DF_cols(r, genCode), axis=1)

    # keep
    all_df = all_df.append(affected_df)

# write to the same as infile
all_df.to_csv(opt.infile, sep="\t", index=False, header=True)


#!/hpcf/apps/python/install/2.7.13/bin/python





"""
Program to get average editing frequency from CrispEsso. 
#Updated version 2.0
"""

import matplotlib
import pandas as pd
matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns
import string
import sys
import glob
import numpy as np
import os
import zipfile
import ast




sample_id = [sys.argv[1]]
#print(sample_id)
df = pd.read_csv(sys.argv[2],sep="\t",header=None) #amplicon- gRNA file
#print(df)
ref=sys.argv[3]
#print(ref)
alt=sys.argv[4]



# Reverse alignment & reference sequence if the strand is - or not found
def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    rev_comp_seq = []
    for base in reversed(seq):
        if base in complement:
            rev_comp_seq.append(complement[base])
        else:
            # other characters (e.g., '-'/delation)
            rev_comp_seq.append('*')  # Replace with *
    return ''.join(rev_comp_seq)



# Compare the gRNA with reference sequence, reverse if strand = "-"
def compare_gRNA_reference(strand, gRNA, aligned_sequence, reference_sequence):
    if len(gRNA) == 21:
        # If gRNA is exactly 21 long, remove the first base
        gRNA = gRNA[1:]

    
    if strand == '+' and gRNA in reference_sequence:
        return '+', aligned_sequence, reference_sequence
    
    
    elif strand == '-' and reverse_complement(gRNA) in reference_sequence:
        gRNA_reverse = reverse_complement(gRNA)
        if gRNA_reverse in reference_sequence:
            return '-', reverse_complement(aligned_sequence), reverse_complement(reference_sequence)
    

    return "gRNA NotFound", None, None # strand, alignment seq, reference seq, strand = 'gRNA NotFound' if no gRNA found


# Check if it's edited in gRNA_20: ref = C, alt = T
def is_edit(r,ref,alt,gRNA,strand):
    strand, read, ref_seq = compare_gRNA_reference(strand,gRNA, r.Aligned_Sequence, r.Reference_Sequence)
    if strand == "gRNA NotFound":
        return ["gRNA not found"], "gRNA NotFound"

    read = read[3:23]
    ref_seq = ref_seq[3:23]
    results = []
    edited = any(ref_seq[i] == ref and read[i] == alt for i in range(20)) 
    
    if edited:
        return ["EDITED", strand]
    else:
        return ["UNEDITED", strand]







# get #Reads from input frequency file, combine with strand
def get_info(f, gRNA, ref, alt,strand):
    df = pd.read_csv(f, sep="\t")
    #print(df.head())
    df[['edit_status', 'strand']] = df.apply(lambda r: is_edit(r, ref, alt, gRNA,strand), axis=1, result_type='expand')
    #print(df[['edit_status','strand']])

    N = df[df['edit_status'] == 'EDITED']['#Reads'].sum()
    denominator = df['#Reads'].sum()
    #strand = df[df['edit_status'] == 'EDITED']['strand'].iloc[0] if not df[df['edit_status'] == 'EDITED']['strand'].empty else 'None' 
    strand = df['strand'].iloc[0] 

    return N, denominator, strand





# Main function call
for s in sample_id:
	outfile = "%s.allele.edit_R37X.tsv"%(s)
	N_list = []
	#P_list = []
	gRNA_list = []
	name_list = []
	strand_list =[]
	denominator_list =[]

	for name, _, gRNA,strand in df.values:
		#print (name,gRNA)
		file = "/directory/CRISPResso_on_{1}/Alleles_frequency_table_around_sgRNA_{2}.txt".format(s,name,gRNA)
		#print (file)
		
		if os.path.isfile(file):
			#print (name,gRNA,file, ref,alt)
			N,denominator,strand = get_info(file,gRNA,ref,alt,strand)
												
			#print(N,strand)
			#print(denominator)

		else:
			N=0
			strand ='File not exist' # strand = File not exis if there's no output frequency file from CRISPREsso
			denominator =0
			#P=-1

		N_list.append(N)
		#P_list.append(P)
		gRNA_list.append(gRNA)
		name_list.append(name)
		strand_list.append(strand) # strand = 'gRNA NotFound' if gRNA doesn't match, strand still shown if only Unedited but still found
		denominator_list.append(denominator)




	out = pd.DataFrame()
	out['Sample'] = name_list
	out['gRNA'] = gRNA_list
	out['#Reads'] = N_list
	out['strand'] = strand_list
	out['denominator'] = denominator_list
	#out['%Reads'] = P_list
	output_dir = '/output_directory'	
	outfile = os.path.join(output_dir, outfile)
	out.to_csv(outfile,sep="\t",index=False)

	
	










#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 21 22:09:26 2021

Determines how c-values for targets in nucleotide sequence space convert to c-values in amino acid sequence space. 

@author: bmiller
"""
############################################################
# Import modules and functions

from random import (randint, sample, random)
import numpy as np
import statistics


############################################################
# Define Constants

TRIAL_RUNS = 10000
NUCSEQ_LEN = 789    
PROT_LEN = NUCSEQ_LEN/3
MUT_NUM = 263
BETA_TEST = 2
BETA_SEQ = 'GAAACGCTGGTGAAAGTAAAAGATGCTGAAGATCAGTTGGGTGCACGAGTGGGTTACATCGAACTGGATCTCAACAGCGGTAAGATCCTTGAGAGTTTTCGCCCCGAAGAACGTTTTCCAATGATGAGCACTTTTAAAGTTCTGCTATGTGGCGCGGTATTATCCCGTATTGACGCCGGGCAAGAGCAACTCGGTCGCCGCATACACTATTCTCAGAATGACTTGGTTGAGTACTCACCAGTCACAGAAAAGCATCTTACGGATGGCATGACAGTAAGAGAATTATGCAGTGCTGCCATAACCATGAGTGATAACACTGCGGCCAACTTACTTCTGACAACGATCGGAGGACCGAAGGAGCTAACCGCTTTTTTGCACAACATGGGGGATCATGTAACTCGCCTTGATCGTTGGGAACCGGAGCTGAATGAAGCCATACCAAACGACGAGCGTGACACCACGATGCCTGTAGCAATGGCAACAACGTTGCGCAAACTATTAACTGGCGAACTACTTACTCTAGCTTCCCGGCAACAATTAATAGACTGGATGGAGGCGGATAAAGTTGCAGGACCACTTCTGCGCTCGGCCCTTCCGGCTGGCTGGTTTATTGCTGATAAATCTGGAGCCGGTGAACGTGGGAGTCGCGGTATCATTGCAGCACTGGGGCCAGATGGTAAGCCCTCCCGTATCGTAGTTATCTACACGACGGGGAGTCAGGCAACTATGGATGAACGAAATAGACAGATCGCTGAGATAGGTGCCTCACTGATTAAGCATTGGTAATGA'
MUT_BIAS = [[.07, .16, .21], [.07, .21, .16], [.13, .35, .07], [.35, .13, .07]]
MUT_BIAS_REV = [[.07, .13, .35], [.07, .35, .13], [.16, .21, .07], [.21, .16, .07]]

TEST_SEQ = 'GAAACGCTG'

############################################################
# Functions for generating nucleotide and protein sequences

table = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                 
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
}

def DNA_prot(seq):
    protein = []
    if len(seq)%3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            codon_str = ''.join([str(elem) for elem in codon])
            protein.append(table[codon_str])
    return protein


def num_nuc(n):
    if n==1:
        return 'A'
    elif n==2:
        return 'T'
    elif n==3:
        return 'C'
    elif n==4:
        return 'G'
    else:
        print("num_nuc error")
        
def nuc_num(nuc):
    if nuc=='A':
        return 1
    elif nuc=='T':
        return 2
    elif nuc=='C':
        return 3
    elif nuc=='G':
        return 4
    else:
        print(nuc)
        print("nuc_num error")


def test_stop(codon):
    codon_str = ''.join([str(elem) for elem in codon])
    if codon_str in ['TAA', 'TAG', 'TGA']:
        return True
    else:
        return False
    

def generate_seq(n):
    DNA_seq = []
    codon = ['', '', '']
    for pos in range(0, n, 3):
        stop_codon = True
        while stop_codon:
            codon[0] = num_nuc(randint(1,4))
            codon[1] = num_nuc(randint(1,4))
            codon[2] = num_nuc(randint(1,4))
            stop_codon = test_stop(codon)
        DNA_seq.extend(codon)
    return DNA_seq


def mutbias_DNA(seq, m, mut_bias_mat, c):
    total = sum(item for subl in mut_bias_mat for item in subl)
    mut_bias = [[item / total for item in subl] for subl in mut_bias_mat]
    positions = list(range(len(seq)))
    mutated = seq.copy()
    for pos in positions:
        randnum1 = random()
        randnum2 = random()
        nucnum = nuc_num(seq[pos])
        brow = mut_bias[nucnum-1]
        rowsum = sum(brow)
        if randnum1 < rowsum * c/.25:
            if nucnum == 1: 
                if randnum2 <= brow[0]/rowsum:
                    mutated[pos] = 'T'
                elif randnum2 <= (brow[0] + brow[1])/rowsum:
                    mutated[pos] = 'C'
                else:
                    mutated[pos] = 'G'
            if nucnum == 2:  
                if randnum2 <= brow[0]/rowsum:
                    mutated[pos] = 'A'
                elif randnum2 <= (brow[0] + brow[1])/rowsum:
                    mutated[pos] = 'C'
                else:
                    mutated[pos] = 'G'
            if nucnum == 3:  
                if randnum2 <= brow[0]/rowsum:
                    mutated[pos] = 'A'
                elif randnum2 <= (brow[0] + brow[1])/rowsum:
                    mutated[pos] = 'T'
                else:
                    mutated[pos] = 'G'
            if nucnum == 4:  
                if randnum2 <= brow[0]/rowsum:
                    mutated[pos] = 'A'
                elif randnum2 <= (brow[0] + brow[1])/rowsum:
                    mutated[pos] = 'T'
                else:
                    mutated[pos] = 'C'
    return mutated


def mut_DNA(seq, m):
    positions = list(range(len(seq)))
    muts = sample(positions, m)
    mutated = seq.copy()
    for mut_pos in muts:
        new_nucn = nuc_num(mutated[mut_pos]) + randint(1,3)
        if new_nucn > 4:
            new_nucn -= 4
        mutated[mut_pos] = num_nuc(new_nucn)
    return mutated


def print_seq(seq):
    DNA_seq = ""
    prot_seq = ""
    for i in range(0, len(seq), 3):
        codon = seq[i:i + 3]
        codon_str = ''.join([str(elem) for elem in codon])
        amino = table[codon_str]
        DNA_seq += codon_str + " "
        prot_seq += " " + amino + "  "
    print(DNA_seq)
    print(prot_seq)

def change_stop(codon, mut_bias):
    print(codon)
    pos = randint(0,2)
    nucnum = nuc_num(codon[pos])
    brow = mut_bias[nucnum-1]
    rowsum = sum(brow)
    randnum = random()
    if nucnum == 1: 
        if randnum <= brow[0]/rowsum:
            codon[pos] = 'T'
        elif randnum <= (brow[0] + brow[1])/rowsum:
            codon[pos] = 'C'
        else:
            codon[pos] = 'G'
    if nucnum == 2:  
        if randnum <= brow[0]/rowsum:
            codon[pos] = 'A'
        elif randnum <= (brow[0] + brow[1])/rowsum:
            codon[pos] = 'C'
        else:
            codon[pos] = 'G'
    if nucnum == 3:  
        if randnum <= brow[0]/rowsum:
            codon[pos] = 'A'
        elif randnum <= (brow[0] + brow[1])/rowsum:
            codon[pos] = 'T'
        else:
            codon[pos] = 'G'
    if nucnum == 4:  
        if randnum <= brow[0]/rowsum:
            codon[pos] = 'A'
        elif randnum <= (brow[0] + brow[1])/rowsum:
            codon[pos] = 'T'
        else:
            codon[pos] = 'C'

    return test_stop(codon)



############################################################
# Run trials

amino_difs = []
for trial in range(TRIAL_RUNS):
    if BETA_TEST == 1:  
        nuc_seq = generate_seq(NUCSEQ_LEN)  
        mut_seq = mut_DNA(nuc_seq, MUT_NUM) 
    elif BETA_TEST == 2:
        nuc_seq = [char for char in BETA_SEQ]
        c = MUT_NUM/NUCSEQ_LEN
        mut_seq = mutbias_DNA(nuc_seq, MUT_NUM, MUT_BIAS, c)
    else:
        nuc_seq = [char for char in BETA_SEQ]
        c = MUT_NUM/NUCSEQ_LEN
        mut_seq = mutbias_DNA(nuc_seq, MUT_NUM, MUT_BIAS_REV, c)
    
    protein_wt = np.array(DNA_prot(nuc_seq))
    protein_mut = np.array(DNA_prot(mut_seq))
    amino_diff = PROT_LEN - np.sum(protein_wt == protein_mut)
    amino_difs.append(float(amino_diff))


print("Protein Differences: {:.2f} +/- {:.2f}".format(statistics.mean(amino_difs), 
      statistics.pstdev(amino_difs)))
print("Difference %: {:.0%} +/- {:.0%}".format(statistics.mean(amino_difs)/PROT_LEN, 
      statistics.pstdev(amino_difs)/PROT_LEN))


"""
# Percentage of mutations resulting in stop codons
stop_num = 0
for trial in range(20):
    nuc_seq = generate_seq(3)  
    if change_stop(nuc_seq, MUT_BIAS_REV):
        stop_num+=1

print("Stop Codon %: {:.2%}".format(stop_num/20))
"""





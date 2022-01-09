#
# Assignment 5/6 Biopython/Phylogenetics
# Bioinformatics BBE5, WS21/22
# Author: Dora Frauscher
# GitHub: https://github.com/dorafrau/BIN-Assignment-5-6.git
#

from Bio import Entrez, pairwise2
from Bio import SeqIO
import pandas as pd
from Bio.SeqUtils import GC

sequences = []
records = []

Entrez.email = 'A.N.Other@example.com'

accession_numbers = ['NM_001009392', 'NM_204628', 'NM_001009211', 'M26744', 'NM_001290980']
# creating the table that will be filled with the asked paramters after they were processed
table = pd.DataFrame(index=['1', '2', '3', '4', '5'],
                     columns=['title', 'acc_number', 'organism', 'length', 'GC_%', 'instability', 'aromaticity'])

for i in range(5):
    # retreiving information from genbenk about specific gene
    handle = Entrez.efetch(db="nucleotide", id=accession_numbers[i], rettype="gb", retmode="text")
    record = SeqIO.read(handle, 'genbank')
    records.append(record)

    ####################################################################
    # get nucleotide sequences
    ####################################################################
    sequence = record.seq
    sequences.append(sequence)
    # https://stackoverflow.com/questions/28354817/how-to-get-the-scientific-name-given-the-genbank-accession-code-to-biopython

    ####################################################################
    # Note down the basic info: accession number, title, organism, and length of the sequence
    ####################################################################
    table.loc[str(i + 1), 'organism'] = record.annotations['source']  # organism
    accessions = record.annotations['accessions']
    table.loc[str(i + 1), 'acc_number'] = accessions[0]  # accessionnumber
    table.loc[str(i + 1), 'title'] = record.description  # title
    table.loc[str(i + 1), 'length'] = len(str(sequence))

    ####################################################################
    # Determine the GC percentage
    ####################################################################
    table.loc[str(i + 1), 'GC_%'] = GC(str(sequence))  # https://biopython.org/docs/1.75/api/Bio.SeqUtils.html

####################################################################
# Perform a sequence alignment
####################################################################
seq1 = sequences[2]
seq2 = sequences[4]
# https://www.kaggle.com/mylesoneill/pairwise-alignment-using-biopython
a = pairwise2.align.localxx(seq1, seq2) # local alignment -  Identical characters are given 2 points, 1 point is deducted for each non-identical character.
print(pairwise2.format_alignment(*a[0]))
# I am aware that i should use "MultipleSeqAignment" but since the sequences are not the same length it doesnt work. Using pairwise, this problem is automaticalli taken care of.

####################################################################
# Create a phylogenetic tree with visualisation according to the scores of the sequence alignment
####################################################################
# Like this (using Pairwise alignmnet) though I cant get a score that i can use for the phylogenetic tree and I couldt find a solution to work alround

####################################################################
# Calculate: instability, aromaticity, isoelectric point
####################################################################
sequences_translated = [sec.transcribe().translate() for sec in sequences]
# the translated sequences contain of '*' characters and therefore cant be used to get insability, aromatzicity and isoelectric point SeqUtils.ProtParam module - I coult find a work-around
# furthermore i get warned about the lenght of the sequence not being a multiple of three

# print out whole table
pd.set_option('display.max_columns', None)
print(table)

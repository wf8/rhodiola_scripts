#! /usr/bin/env python
import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Align.Applications import ClustalwCommandline
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align import AlignInfo

"""
concatenate_matrices.py

Concatenate the chloroplast and nuclear matrices.

March 19, 2013
"""

# get the sequences out of each file
handle = open("chloroplast_total_alignment.fasta", "rU")
chloroplast_records = list(SeqIO.parse(handle, "fasta"))
handle.close()

handle = open("ITS_total_alignment_simple_ids.fas", "rU")
ITS_records = list(SeqIO.parse(handle, "fasta"))
handle.close()

'''
# final chloroplast matrix 1649 positions long
# final ITS matrix 681 long
i = 0
gap_chloroplast = ''
while i < 1649:
	gap_chloroplast = gap_chloroplast + '-'
	i += 1
i = 0
gap_ITS = ''
while i < 681:
	gap_ITS = gap_ITS + '-'
	i += 1
'''

combined_records = []

for ITS_record in ITS_records:
	found = False
	for chloroplast_record in chloroplast_records:
		if chloroplast_record.id.strip() == ITS_record.id.strip():
			combined_records.append(SeqRecord(ITS_record.seq + chloroplast_record.seq, id=ITS_record.id.strip(), description=''))
			found = True
			break

# write fasta file
SeqIO.write(combined_records, "combined_final_matrix.fasta", "fasta")
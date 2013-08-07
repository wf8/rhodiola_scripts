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
organize_chloroplast_sequences.py

Concatenate all the trnL, trnL-trnF, and trnS-trnG sequences from the Olfelt lab. Then
read in sequences from GenBank and align everthing. Uses ClustalW for alignments.

March 18, 2013
"""

# get the sequences out of each file

# get genbank sequences
handle = open("chloroplast_sequences_from_genbank_simple_ids.fasta", "rU")
final_records = list(SeqIO.parse(handle, "fasta"))
handle.close()

# get olfelt lab sequences 
handle = open("trnL_intron_final_matrix.fas", "rU")
trnL_records = list(SeqIO.parse(handle, "fasta"))
handle.close()

handle = open("trnL-trnF_spacer_final_matrix.fas", "rU")
trnLF_records = list(SeqIO.parse(handle, "fasta"))
handle.close()

handle = open("trnS-trnG_spacer__final_matrix.fas", "rU")
trnStrnG_records = list(SeqIO.parse(handle, "fasta"))
handle.close()

'''
# add the reverse complement of trnL records to genbank
for record in trnL_records:
	final_records.append(SeqRecord(record.seq.reverse_complement(), id=record.id, description=''))

# go through all trnL-trnF records, reverse complement and concatenate with final_records
for record in trnLF_records:
	found = False
	for final_record in final_records:
		if record.id.strip() == final_record.id.strip():
			final_record.seq = final_record.seq + record.seq.reverse_complement()
			found = True
			break
	if not found:
		final_records.append(SeqRecord(record.seq.reverse_complement(), id=record.id, description=''))

# write fasta file
SeqIO.write(final_records, "trnL-trnF.fasta", "fasta")
'''

# final trnL-trnL-trnF matrix 957 positions long
# final trnS-trnG 692 long
i = 0
gapLF = ''
while i<957:
	gapLF = gapLF + '-'
	i += 1
i = 0
gapSG = ''
while i<692:
	gapSG = gapSG + '-'
	i += 1

handle = open("trnL-trnL-trnF_final_matrix.fas", "rU")
trnLtrnLtrnF_records = list(SeqIO.parse(handle, "fasta"))
handle.close()

final_records = []

# go through all the trnLtrnLtrnF_records and concatenate trnS-trnG or a gap 
for LLF_record in trnLtrnLtrnF_records:
	found = False
	for S_record in trnStrnG_records:
		if LLF_record.id.strip() == S_record.id.strip():
			final_records.append(SeqRecord(LLF_record.seq + S_record.seq, id=LLF_record.id.strip(), description=''))
			found = True
			break
	if not found:
		final_records.append(SeqRecord(LLF_record.seq + gapSG, id=LLF_record.id.strip(), description=''))
		
for S_record in trnStrnG_records:
	found = False
	for LLF_record in trnLtrnLtrnF_records:
		if LLF_record.id.strip() == S_record.id.strip():
			found = True
			break
	if not found:
		final_records.append(SeqRecord(gapLF + S_record.seq, id=S_record.id.strip(), description=''))

# write fasta file
SeqIO.write(final_records, "chloroplast.fasta", "fasta")

# align
#clustalw_cline = ClustalwCommandline("clustalw", infile="chloroplast.fasta")
#clustalw_cline()
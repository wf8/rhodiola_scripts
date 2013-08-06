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
organize_sequences.py

Makes consensus sequences out of forward and reverse sequences, and then 
outputs aligned FASTA files for each DNA region. Uses ClustalW for alignments.

December 19, 2012
"""

# get a list of all files 
path = 'S1204413.NIU/SEQUENCES/'
filelist = os.listdir(path)
# alphabetize the list of filenames
filelist.sort()

def make_consensus( rev_string, for_string, seqfile):
    "function that accepts 2 sequence and returns the consensus sequence"
    # make fasta file for each paired sequence
    rev_sequence = Seq(rev_string.replace("\n", "").replace('\r', '').replace(' ', ''), IUPAC.ambiguous_dna)
    rev_sequence= rev_sequence.reverse_complement()
    for_sequence = Seq(for_string.replace("\n", "").replace('\r', '').replace(' ', ''), IUPAC.ambiguous_dna)
    paired_sequences = [SeqRecord(rev_sequence, id="rev"), SeqRecord(for_sequence, id="for")]
    if not os.path.exists("results/"):
        os.makedirs("results/")
    fasta_file = "results/" + seqfile + ".fasta"
    SeqIO.write(paired_sequences, fasta_file, "fasta")
    # align the paired sequences
    aln_file = "results/" + seqfile + ".aln"
    clustalw_cline = ClustalwCommandline("clustalw", infile=fasta_file, outfile=aln_file )
    clustalw_cline()
    # hack so that dumb_consensus will accept 1 base call against N
    f = open(aln_file, 'r+')
    contents = f.read()
    f.close()
    f = open(aln_file, 'w')
    f.write( contents.replace('N','.') )
    f.close()
    # read in alignment file and generate consensus
    alignment = AlignIO.read(aln_file, "clustal")
    summary_align = AlignInfo.SummaryInfo(alignment)
    return summary_align.dumb_consensus(ambiguous = "N", threshold=0.0, require_multiple=0)

ITS_records = []
ITS1_records = []
ITS2_records = []
trnS_trnG_records = []
trnL_trnF_records = []
trnL_records = []
i = 0
for seqfile in filelist:
    print "current file is: " + seqfile

    # read each seq file and turn it into a string
    seqstring = open(path + seqfile).read()    
    
    if i == 0:
        for_string = seqstring
        i += 1
    else:
        # send the for and rev strings to function that return seq object of consensus sequence
        consensusseq = make_consensus(seqstring, for_string, seqfile)
        # print "consensus sequence: " + consensusseq
        i = 0

		# get the sample id
        sample_id = seqfile[0:seqfile.find("_rev")]

        # make fasta files for each dna region of all the consensus sequences
        if "ITS1" in seqfile:
            ITS1_records.append( SeqRecord(consensusseq, id=sample_id) )
        elif "ITS2" in seqfile:
            ITS2_records.append( SeqRecord(consensusseq, id=sample_id) )
        elif "ITS" in seqfile:
            ITS_records.append( SeqRecord(consensusseq, id=sample_id) )
        elif "trnS_trnG" in seqfile:
            trnS_trnG_records.append( SeqRecord(consensusseq, id=sample_id) )
        elif "trnL_trnF" in seqfile:
            trnL_trnF_records.append( SeqRecord(consensusseq, id=sample_id) )
        elif "trnL" in seqfile:
            trnL_records.append( SeqRecord(consensusseq, id=sample_id) )

SeqIO.write(ITS_records, "ITS.fasta", "fasta")
SeqIO.write(ITS1_records, "ITS1.fasta", "fasta")
SeqIO.write(ITS2_records, "ITS2.fasta", "fasta")
SeqIO.write(trnS_trnG_records, "trnS_trnG.fasta", "fasta")
SeqIO.write(trnL_trnF_records, "trnL_trnF.fasta", "fasta")
SeqIO.write(trnL_records, "trnL.fasta", "fasta")

clustalw_cline = ClustalwCommandline("clustalw", infile="ITS.fasta")
clustalw_cline()
clustalw_cline = ClustalwCommandline("clustalw", infile="trnS_trnG.fasta")
clustalw_cline()
clustalw_cline = ClustalwCommandline("clustalw", infile="trnL_trnF.fasta")
clustalw_cline()
clustalw_cline = ClustalwCommandline("clustalw", infile="trnL.fasta")
clustalw_cline()
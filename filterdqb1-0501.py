#! /usr/bin/env python3

from Bio.Seq import Seq
from Bio import SeqIO
from HLA import HLA

DQB1_Seqs = list(SeqIO.parse("/scratch/okask/minion-hla/dbq1-05-01.fasta","fasta"))

HLASeqRecordList = []

for sequence in DQB1_Seqs:
    desc_list = sequence.description.split()
    gene_list = desc_list[1].split("*")
    hla_gene = gene_list[0]
    seq_str = str(sequence.seq)
    if ((hla_gene == "A" or hla_gene == "B" or hla_gene == "C" or hla_gene == "DQB1" or hla_gene == "DPB1" or hla_gene == "DRB1") and seq_str[:3].upper() == "ATG"):
        allele_list = gene_list[1].split(":")
        allele_len = len(allele_list)
        if allele_len < 3: 
            field_1 = allele_list[0]
            field_2 = allele_list[1]
            field_3 = "" 
            field_4 = ""
        elif allele_len < 4:
            field_1 = allele_list[0]                                                
            field_2 = allele_list[1]
            field_3 = allele_list[2]
            field_4 = ""
        else:
            field_1 = allele_list[0]                                                
            field_2 = allele_list[1]                                                
            field_3 = allele_list[2] 
            field_4 = allele_list[3]
        hla = HLA(hla_gene, field_1, field_2, field_3, field_4, seq_str, sequence.id)
        if hla.is01last():
            HLASeqRecordList.append(hla.toSeqRecord())

SeqIO.write(HLASeqRecordList, "DQB10501-filtered.fasta", "fasta")

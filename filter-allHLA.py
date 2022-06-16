#! venv/bin/python 
from Bio.Seq import Seq
from Bio import SeqIO
from HLA import HLA

HLA_Seqs = list(SeqIO.parse("hla_nuc.fasta","fasta"))

a_count = 0 
b_count = 0
c_count = 0 
dp_count = 0
dr_count = 0
dq_count = 0
other_count = 0 

HLASeqRecordList = []
HLA_A_List = []
HLA_B_List = []
HLA_C_List = []
HLA_DQB1_List = []
HLA_DRB1_List = []
HLA_DPB1_List = []

def bubblesortHLA(hla_list):
    swapped = False 
    for n in range(len(hla_list)-1, 0, -1):
        for i in range(n):
            if hla_list[i].compareHLAGene(hla_list[i+1]):
                swapped = True
                hla_list[i], hla_list[i+1] = hla_list[i+1], hla_list[i] 
        if not swapped:
            return 

    
for sequence in HLA_Seqs:
    desc_list = sequence.description.split()
    gene_list = desc_list[1].split("*")
    hla_gene = gene_list[0]
    if hla_gene == "A":
        a_count+=1 
    elif hla_gene == "B": 
        b_count+=1
    elif hla_gene == "C": 
        c_count+=1 
    elif hla_gene == "DQB1":
        dq_count+=1
    elif hla_gene == "DPB1":
        dp_count+=1
    elif hla_gene == "DRB1":
        dr_count+=1
    else: 
        other_count+=1
 
    seq_str = str(sequence.seq).upper()
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
        HLASeqRecordList.append(hla.toSeqRecord())
        if hla.getGene() == "A":
            HLA_A_List.append(hla)
        elif hla.getGene() == "B":
            HLA_B_List.append(hla)
        elif hla.getGene() == "C":
            HLA_C_List.append(hla)
        elif hla.getGene() == "DPB1":
            HLA_DPB1_List.append(hla)
        elif hla.getGene() == "DQB1":
            HLA_DQB1_List.append(hla)
        elif hla.getGene() == "DRB1":
            HLA_DRB1_List.append(hla)

print("Initial HLA A alleles: "+str(a_count)+", Filtered HLA A counts: "+str(len(HLA_A_List)))
print("Initial HLA B alleles: "+str(b_count)+", Filtered HLA B counts: "+str(len(HLA_B_List)))
print("Initial HLA C alleles: "+str(c_count)+", Filtered HLA C counts: "+str(len(HLA_C_List)))

print("Initial HLA DQB1 alleles: "+str(dq_count)+", Filtered HLA DQB1 counts: "+str(len(HLA_DQB1_List)))
print("Initial HLA DPB1 alleles: "+str(dp_count)+", Filtered HLA DPB1 counts: "+str(len(HLA_DPB1_List)))
print("Initial HLA DRB1 alleles: "+str(dr_count)+", Filtered HLA DRB1 counts: "+str(len(HLA_DRB1_List)))
print("Other count: "+ str(other_count))

#SeqIO.write(HLA_A_List, "HLA-A-filtered.fasta", "fasta")
#SeqIO.write(HLA_B_List, "HLA-B-filtered.fasta", "fasta")
#SeqIO.write(HLA_C_List, "HLA-C-filtered.fasta", "fasta")
#SeqIO.write(HLA_DRB1_List, "HLA-DRB1-filtered.fasta", "fasta")
#SeqIO.write(HLA_DQB1_List, "HLA-DQB1-filtered.fasta", "fasta")
#SeqIO.write(HLA_DPB1_List, "HLA-DPB1-filtered.fasta", "fasta")
SeqIO.write(HLASeqRecordList, "allHLA-filtered.fasta", "fasta")

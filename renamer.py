#!/usr/bin/env python3

from Bio import SeqIO
import sys

utr_records = []
protein_records = []
matched_records = [] 
output_utr_records = []
output_protein_records = []


utr_file = sys.argv[1]
prot_file = sys.argv[2]

for record in SeqIO.parse(utr_file, "fasta"):
    utr_records.append(record)

for record in SeqIO.parse(prot_file, "fasta"):
    protein_records.append(record)

def find_match(record):
    #print(record.id)
    for r in protein_records:
        if r.id == record.id:
            #print(record.id, r.id)
            return (record, r)

for record in utr_records:
    matched_records.append(find_match(record))

for r_pair in matched_records:
    desc = r_pair[0].description.strip()
    desc_words = desc.split()
    desc_words[0] = desc_words[0].replace("(2cells)","")
    desc_words[0] = desc_words[0].replace("Dcan.","")
    desc_words[0] = desc_words[0].replace(".","_")
    desc_genome = desc_words[-1].split("_")
    default = desc_genome[0] + "_" + desc_words[0]
    new_id = input(f'Renaming gene {default} ({r_pair[0].description[:45]}): ')
    r_pair[0].id = r_pair[0].description = new_id or default
    r_pair[1].id = r_pair[1].description = new_id or default
    output_utr_records.append(r_pair[0])
    output_protein_records.append(r_pair[1])

SeqIO.write(output_utr_records, "renamed_UTR_records.fasta", "fasta")
SeqIO.write(output_protein_records, "renamed_protein_records.fasta", "fasta")


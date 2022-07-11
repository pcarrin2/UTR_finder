#!/usr/bin/env python3

from Bio import SeqIO
import sys

records = []
clean_records = []

infile = sys.argv[1]
outfile = infile.rsplit(".", 1)[0] + "_renamed.fasta"

for record in SeqIO.parse(infile, "fasta"):
    records.append(record)

for r in records:
    desc = r.description.strip()
    desc_words = desc.split()
    desc_words[0] = desc_words[0].replace("(2cells)","")
    desc_words[0] = desc_words[0].replace("Dcan.","")
    desc_words[0] = desc_words[0].replace(".","_")
    desc_genome = desc_words[-1].split("_")
    default = desc_genome[0] + "_" + desc_words[0]
    new_id = input(f'Renaming gene {default} ({r.description[:45]}): ')
    r.id = r.description = new_id or default
    clean_records.append(r)

SeqIO.write(clean_records, outfile, "fasta")

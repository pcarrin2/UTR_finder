#!/usr/bin/env python3
help_string = """
This script finds 5' and 3' UTRs. Given a genome FASTA file and the locations of
genes in that file, it writes a series of multi-entry FASTA files containing UTRs
and CDS sequences of listed genes. Output is customizable with the -format option.

Usage: 
$ ./utr-fasta-creator.py -genome [genome file] -coords [coordinate file] 
    -outfile [output file prefix] -format [format string] -na_strategy [string]
    -rev_comp [reverse complement string]

Coordinate file should be in SnapGene "Export Features..." output format (TSV, one column
containing coordinates separated by "..").

Flags:

-genome       <file>      Genome file (FASTA format)

-coords       <file>      CDS coordinate file

-outfile      <file>      Output file prefix. If multiple outfiles are made,
                          they will be numbered sequentially.

-format      <string>     Format for output. Should be a three-letter string
                          specifying what actions to take for the 5' UTR, CDS,
                          and 3' UTR. Default is "NNN": leave all three parts
                          as nucleotide data, and concatenate them into a single
                          entry. "X" removes a part, "P" translates to protein,
                          "R" transcribes to RNA, and "|" between two letters
                          partitions the data into separate files. So "N|PX"
                          prints the 5' UTR as nucleotide data in one file, prints the
                          CDS as protein data in another file, and omits the 3' UTR.

-na_strategy <string>     A strategy for dealing with cases where TATA boxes or
                          poly-A sequences cannot be found in the genome file. 
                          Must be one of "o", "t", or "a". "o" omits the entry
                          entirely. "t" trims the entry to the start codon (if the
                          TATA box is absent) or stop codon (if the poly-A sequence
                          is absent), removing the affected UTR only. "a" approximates
                          the UTR with a 70-base-pair segment.

-rev_comp    <string>     A string in each gene's label that indicates that it is
                          transcribed in the reverse direction.
"""

import sys, re, argparse, warnings
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

parser = argparse.ArgumentParser(add_help=False)

parser.add_argument("-genome")
parser.add_argument("-coords")
parser.add_argument("-outfile")
parser.add_argument("-format", default="NNN")
parser.add_argument("-na_strategy", default="o")
parser.add_argument("-rev_comp", default="")
parser.add_argument("-h","--help", action='store_true', dest="help")

args = parser.parse_args()

if args.help:
    print(help_string)
    quit()

utr_len_guess = 70
genome_file = open(args.genome,"r")
genome_record = SeqIO.read(genome_file, "fasta")
genome_label = genome_record.id
genome = genome_record.seq
comp_genome = genome.complement()
cds_coord_file = open(args.coords, "r")

def na_parse(utr_coords, cds_coords, label):
    # takes UTR coordinates from find_utrs() as list, then CDS coordinates from input file.
    # if utr_coords contains None, parses according to na_strategy (o, t, or a).
    if None in utr_coords:
        print(f"gene {label.strip()}: couldn't find TATA box and/or polyA sequence")
        if args.na_strategy not in ["o","t","a"]:
            raise Exception(f"Script failed to find a UTR in gene {label.strip()}, and -na_strategy is improperly set. Please set to 'o', 't', or 'a'.")
        if args.na_strategy == "o":
            return None
        else:
            if args.na_strategy == "t":
                backup_coords = cds_coords
            else: #args.na_strategy = "a"
                backup_coords = [max(cds_coords[0]-70, 0), min(cds_coords[1]+70, len(genome.seq))]
            clean_coords = []
            for index, c in enumerate(utr_coords):
                if c == None:
                    clean_coords.append(backup_coords[index])
                else:
                    clean_coords.append(c)
            return clean_coords
    else:
        return utr_coords

def parse_coords(coord_string):
    # takes a line from a CDS coordinate file and parses it into a start coordinate (type int),
    # a stop coordinate (type int), a label (type str), and whether the gene is reverse-transcribed. 
    # returns (start, stop, label, rev) as a tuple of int, int, str, and bool.
    split_str = coord_string.split()
    label = ""
    coord_regex = r"^[0-9,]+\.\.[0-9,]+$" # example match: "12,345..67,890"
    for s in split_str:
        if re.match(coord_regex, s):
            coord_range = s
        else:
            label += s.strip() + " "
    label = label.strip()
    rev = False
    if args.rev_comp != "" and args.rev_comp in label:
        rev = True
    if not label:
        label = "(unlabeled)"
    coord_range = re.sub(",", "", coord_range) # take out the commas to make the coordinate numbers int-friendly
    coords = sorted([int(n) for n in coord_range.split("..")])
    return (coords[0], coords[1], label, rev)

def find_utrs(coords, rev):
    # scanning the complement instead of the reverse complement because I don't want to
    # deal with converting indices -- tbh it might be easier to use the reverse complement...
    if rev:
        gen = comp_genome
        feature_seqs = ["AAATAA","ATAT"]
    else:
        gen = genome
        feature_seqs = ["TATA","AATAAA"]
    
    results = []
    # finding first UTR (before the CDS)

    scan_iter = 1
    done = False
    while done == False:
        scan_start = max(0, coords[0] - scan_iter * utr_len_guess)
        scan_end = max(0, coords[0] - (scan_iter - 1) * utr_len_guess)
        seq_to_scan = gen[scan_start:scan_end]
        location = str(seq_to_scan).rfind(feature_seqs[0])
        # if we've gone too far and we haven't found anything, quit and print "n/a".
        # if we haven't found anything yet, but we aren't super far away, try again.
        # if we do find something, return its location.
        if scan_end == 0 or scan_iter > 20:
            done = True
            location_in_genome = None
        elif location == -1:
            scan_iter += 1
        else:
            done = True
            location_in_genome = scan_start + location
    results.append(location_in_genome)

    # finding second UTR (after the CDS)

    genome_end = len(gen) - 1
    scan_iter = 1
    done = False
    while done == False:
        scan_start = min(genome_end, coords[1] + (scan_iter - 1) * utr_len_guess)
        scan_end = min(genome_end, coords[1] + scan_iter * utr_len_guess)
        seq_to_scan = gen[scan_start:scan_end]
        location = str(seq_to_scan).find(feature_seqs[1])
        # if we've gone too far and we haven't found anything, quit and print "n/a".
        # if we haven't found anything yet, but we aren't super far away, try again.
        # if we do find something, return its location.
        if scan_start == genome_end or scan_iter > 20:
            done = True
            location_in_genome = None
        elif location == -1:
            scan_iter += 1
        else:
            done = True
            location_in_genome = scan_start + location + len(feature_seqs[1])
    results.append(location_in_genome)

    return results

genes_to_write = []

for i, line in enumerate(cds_coord_file):
    if not line.strip() or line[0]=="#": 
        # filters out empty lines, like trailing newline, and comments
        continue

    cds_coords = parse_coords(line)
    label = cds_coords[2]
    rev = cds_coords[3]
    # find UTRs
    coding_region = [cds_coords[0], cds_coords[1]]
    utr_coords = find_utrs(coding_region, rev)
    utrs = na_parse(utr_coords, cds_coords, label)

    if utrs == None: # --omit-nas has been set, and one or more UTRs couldn't be found
        continue

    # create FASTA entries based on those UTRs
    fiveprime = genome[utrs[0]:cds_coords[0]-1]
    cds = genome[cds_coords[0]-1:cds_coords[1]]
    threeprime = genome[cds_coords[1]:utrs[1]]

    if rev:
        orig_fiveprime = fiveprime
        fiveprime = threeprime.reverse_complement()
        cds = cds.reverse_complement()
        threeprime = orig_fiveprime.reverse_complement()

    format_str = args.format.upper()
    formatted_seqs = []
    for n, seq in enumerate([fiveprime, cds, threeprime]):
        if format_str[0] == "|":
            def add_seq(seq):
                formatted_seqs.append(seq)
            format_str = format_str[1:]
        else:
            def add_seq(seq):
                if formatted_seqs:
                    concat = Seq("").join([formatted_seqs[-1], seq])
                    del formatted_seqs[-1]
                    formatted_seqs.append(concat)
                else:
                    formatted_seqs.append(seq)

        if format_str[0] == "N":
            add_seq(seq)
        elif format_str[0] == "R":
            add_seq(seq.transcribe())
        elif format_str[0] == "P":
            mrna = seq.transcribe()
            start_codon_site = mrna.find("AUG")
            if start_codon_site > 2:
                warnings.warn(f"Gene {label}: Start codon site is farther away than expected ({start_codon_site} bp)")
            coding_mrna = mrna[start_codon_site:]
            add_seq(coding_mrna.translate())
        elif format_str[0] == "X":
            add_seq(Seq(""))
        else:
            raise Exception("Invalid format string.")
        
        format_str = format_str[1:]
        

    records = []
    for gene in formatted_seqs:
        records.append(SeqRecord(gene, id=label, description="from genome: " + genome_label))
    genes_to_write.append(records)

for n in range(len(genes_to_write[0])): # equal to number of outfiles needed
    outfile_name = args.outfile + "-" + str(n) + ".fa"
    SeqIO.write([seqlist[n] for seqlist in genes_to_write], outfile_name, "fasta")
    print(f"wrote partition to {outfile_name}")
print("done!")

genome_file.close()
cds_coord_file.close()


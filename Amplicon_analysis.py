import csv

"""
This is a python script that analyse Amplicon FASTQ file generated with Illumina Sequencing
The script classify which amino acids are present at a certain position
The purpose is to see which amino acid mutations have been generated
All the data are exported into a csv
"""

# Translate a codon of 3 bases into the appropriate amino acids
def amino_acids(aa):
    if aa == "TTT" or aa == "TTC":
        return "F"
    if aa == "TTA" or aa == "TTG" or aa == "CTT" or aa == "CTC" or aa == "CTA" or aa == "CTG":
        return "L"
    if aa == "ATT" or aa == "ATC" or aa == "ATA":
        return "I"
    if aa == "ATG":
        return "M"
    if aa == "GTT" or aa == "GTC" or aa == "GTA" or aa == "GTG":
        return "V"
    if aa == "TCT" or aa == "TCC" or aa == "TCA" or aa == "TCG" or aa == "AGT" or aa == "AGC":
        return "S"
    if aa == "CCT" or aa == "CCC" or aa == "CCA" or aa == "CCG":
        return "P"
    if aa == "ACT" or aa == "ACC" or aa == "ACA" or aa == "ACG":
        return "T"
    if aa == "GCT" or aa == "GCC" or aa == "GCA" or aa == "GCG":
        return "A"
    if aa == "TAT" or aa == "TAC":
        return "Y"
    if aa == "TAA" or aa == "TAG" or aa == "TGA":
        return "X"
    if aa == "CAT" or aa == "CAC":
        return "H"
    if aa == "CAA" or aa == "CAG":
        return "Q"
    if aa == "AAT" or aa == "AAC":
        return "N"
    if aa == "AAA" or aa == "AAG":
        return "K"
    if aa == "GAT" or aa == "GAC":
        return "D"
    if aa == "GAA" or aa == "GAG":
        return "E"
    if aa == "TGT" or aa == "TGC":
        return "C"
    if aa == "TGG":
        return "W"
    if aa == "CGT" or aa == "CGC" or aa == "CGA" or aa == "CGG" or aa == "AGA" or aa == "AGG":
        return "R"
    if aa == "GGT" or aa == "GGC" or aa == "GGA" or aa == "GGG":
        return "G"
    else:
        return "N/A"
    
# Start at amino acid position 55 because the amplicon is at that position
# Get the codon of of each position by getting 3 consecutive indices
def getCodon(position, seq):
    try:
        index = (position - 55)*3
        codon = str(seq[index]) + str(seq[index+1]) + str(seq[index+2])
        return codon
    except:
        return "N/A"

# Create csv fle to store all the amino acids present at each position
csv_file = open('Output_Amplicon.csv', 'w') 
writer = csv.writer(csv_file)
writer.writerow(['Position','AminoAcids', 'Count'])    

# Creat a list of positions to find the amino acids
list_positions = [90,111,115,122,102,131,142,146,166,168,175]

# Add row to the csv by iterating through the fastq file
for position in list_positions:
    list_AA = []
    with open("IDT_HM.assembled.fastq", 'r') as fh:
        lines = []
        for line in fh:
            lines.append(line.rstrip())
            # Only choose lines with length 4 which is how fastq file structures
            if len(lines) == 4:
                # Get the sequence
                seq = lines[1]
                # Get the codon
                codon = getCodon(position,seq)
                # Translate the codon
                amino_acid = amino_acids(codon)
                # If the amino acid not already present then add it to avoid duplicates
                if amino_acid not in list_AA:
                    list_AA.append(amino_acid)
                lines = []
    list_AA.remove('N/A')
    writer.writerow([position,list_AA,len(list_AA)])
                      
csv_file.close()
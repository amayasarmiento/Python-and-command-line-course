##Exam 



## Taking a file and reading the contents 
# Read the content of the DNA sequence file
#global variables 

#getting user input on the dna file they want to open 

dna_file = (input("Enter the name of the file you want to open 'no_frameDNA.txt' ; 'salmonella_sejplasmid.txt','salmonella_wesplasmid.txt':"))

with open('dna_file', 'r') as file:
    dna_strand = file.read()
    

#global variable as well - a translation dictionary 

codon_table = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
        'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
    }

##Getting a complementary strand for my second global variable, 
#first function definition - changes base for its complimentary

def complement_strand(dna_strand):
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    complemented_strand = ''
    for base in dna_strand:
        complemented_strand += complement.get(base, base)  # Get the complement of the base or leave it unchanged if not A, T, G, or C
    return complemented_strand

# Generate the complementary DNA strand by using function above
complementary_strand = complement_strand(dna_strand)

##Translate the DNA sequences into their potential protein sequences 
#This function uses the codon table above 

def dna_to_protein(dna_strand):
    
    protein_sequences = {}  # Start a dictionary to store protein sequences
    
    for frame_start in range(3):
        protein_sequence = ''  # Reset protein sequence for each reading frame
        for i in range(frame_start, len(dna_strand), 3):
            codon = dna_strand[i:i+3]
            amino_acid = codon_table.get(codon, '')    #This function uses the codon table above 
            if amino_acid:
                protein_sequence += amino_acid
            else:
                protein_sequence += 'X'
        protein_sequences[f'orf{frame_start + 1}'] = protein_sequence
    
    return protein_sequences

# Get protein sequences for the original DNA strand by using function with dna_strand input
protein_sequences = dna_to_protein(dna_strand)

# Get protein sequences for the complementary DNA strand by using output of first function 
complementary_protein_sequences = dna_to_protein(complementary_strand)

# Print protein sequences to see if they are different
for frame, sequence in protein_sequences.items():
    print(f"First 10 aminoacids for original DNA {frame}: {sequence[:10]}")

for frame, sequence in complementary_protein_sequences.items():
    print(f"First 10 aminoacids of sequence for complementary DNA {frame}: {sequence[:10]}")

#the dictionary above has made two different protein sequences because I had to call it twice for my two inputs
#so we make a single dictionary here to use it in the next function 

def merge_protein_sequences(protein_sequences, complementary_protein_sequences):
    merged_sequences = {}

    # Merge protein sequences for the original DNA strand
    for frame, sequence in protein_sequences.items():
        merged_sequences[f'original_{frame}'] = sequence

    # Merge protein sequences for the complementary DNA strand
    for frame, sequence in complementary_protein_sequences.items():
        merged_sequences[f'complementary_{frame}'] = sequence

    return merged_sequences
    
# Merge the protein sequences
merged_sequences = merge_protein_sequences(protein_sequences, complementary_protein_sequences)

##Identifying the most coding reading frame out of the six 
#first counting the aminoacids in between a START and STOP codon, then making a summary

def count_sequences_between_codons(merged_sequences):
    all_sequence_lengths = {}  # Dictionary to store sequence lengths for each ORF

    for frame, merged_sequence in merged_sequences.items():
        sequence_lengths = []
        in_sequence = False
        current_length = 0

        for amino_acid in merged_sequence:
            if amino_acid == 'M' and not in_sequence:
                # Start of a new sequence
                in_sequence = True
                current_length = 0
            elif amino_acid == '*' and in_sequence:
                # End of the current sequence
                in_sequence = False
                sequence_lengths.append(current_length)
            elif in_sequence:
                # Inside a sequence, increment length
                current_length += 1

        # Store sequence lengths for the current ORF
        all_sequence_lengths[frame] = sequence_lengths

    return all_sequence_lengths

sequence_lengths = count_sequences_between_codons(merged_sequences)

# Print summary outputs for each ORF
for frame, lengths in sequence_lengths.items():
    print(f"Summary for {frame}:")
    print(f"Number of sequences: {len(lengths)}")
    print(f"Max sequence length: {max(lengths) if lengths else 0}")
    print(f"Average sequence length: {sum(lengths) / len(lengths) if lengths else 0}\n")

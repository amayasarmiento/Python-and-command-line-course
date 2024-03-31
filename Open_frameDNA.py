{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": 1,
   "id": "f254da15-af25-4f0a-94aa-ab1f64c2cade",
   "metadata": {},
   "outputs": [],
   "source": [
    "#starting file with imports\n",
    "import numpy as np\n",
    "import pandas as pd\n",
    "import matplotlib.pyplot as plt"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 2,
   "id": "8ceb3140-985b-4fa7-be54-27253582fc0d",
   "metadata": {},
   "outputs": [],
   "source": [
    "## Taking a file and reading the contents \n",
    "# Read the content of the DNA sequence file\n",
    "#global variables \n",
    "\n",
    "with open('no_frameDNA.txt', 'r') as file:\n",
    "    dna_strand = file.read()\n",
    "    "
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 3,
   "id": "ad2ea756-9b07-43bc-8559-f461a1d6025a",
   "metadata": {},
   "outputs": [],
   "source": [
    "#global variable as well - a translation dictionary \n",
    "\n",
    "codon_table = {\n",
    "        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',\n",
    "        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',\n",
    "        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',\n",
    "        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',\n",
    "        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',\n",
    "        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',\n",
    "        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',\n",
    "        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',\n",
    "        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',\n",
    "        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',\n",
    "        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',\n",
    "        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',\n",
    "        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',\n",
    "        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',\n",
    "        'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',\n",
    "        'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',\n",
    "    }"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 4,
   "id": "cf8513d0-2df9-4e96-b808-b82b7adc275d",
   "metadata": {},
   "outputs": [],
   "source": [
    "##Getting a complementary strand for my second global variable, \n",
    "#first function definition - changes base for its complimentary\n",
    "\n",
    "def complement_strand(dna_strand):\n",
    "    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}\n",
    "    complemented_strand = ''\n",
    "    for base in dna_strand:\n",
    "        complemented_strand += complement.get(base, base)  # Get the complement of the base or leave it unchanged if not A, T, G, or C\n",
    "    return complemented_strand"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 5,
   "id": "407a1b59-7224-42aa-8d0e-8d8a79f6baae",
   "metadata": {},
   "outputs": [],
   "source": [
    "# Generate the complementary DNA strand by using function above\n",
    "complementary_strand = complement_strand(dna_strand)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 6,
   "id": "888d27c0-f618-4922-a024-52e2db283d3e",
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "First 10 aminoacids for original DNA orf1: SLNSVSSFSL\n",
      "First 10 aminoacids for original DNA orf2: LLTPCLVFR*\n",
      "First 10 aminoacids for original DNA orf3: S*LRV*FFVD\n",
      "First 10 aminoacids of sequence for complementary DNA orf1: RELRHRSKSN\n",
      "First 10 aminoacids of sequence for complementary DNA orf2: EN*GTDQKAT\n",
      "First 10 aminoacids of sequence for complementary DNA orf3: RIEAQIKKQL\n"
     ]
    }
   ],
   "source": [
    "##Translate the DNA sequences into their potential protein sequences \n",
    "#This function uses the codon table above \n",
    "\n",
    "def dna_to_protein(dna_strand):\n",
    "    \n",
    "    protein_sequences = {}  # Start a dictionary to store protein sequences\n",
    "    \n",
    "    for frame_start in range(3):\n",
    "        protein_sequence = ''  # Reset protein sequence for each reading frame\n",
    "        for i in range(frame_start, len(dna_strand), 3):\n",
    "            codon = dna_strand[i:i+3]\n",
    "            amino_acid = codon_table.get(codon, '')    #This function uses the codon table above \n",
    "            if amino_acid:\n",
    "                protein_sequence += amino_acid\n",
    "            else:\n",
    "                protein_sequence += 'X'\n",
    "        protein_sequences[f'orf{frame_start + 1}'] = protein_sequence\n",
    "    \n",
    "    return protein_sequences\n",
    "\n",
    "# Get protein sequences for the original DNA strand by using function with dna_strand input\n",
    "protein_sequences = dna_to_protein(dna_strand)\n",
    "\n",
    "# Get protein sequences for the complementary DNA strand by using output of first function \n",
    "complementary_protein_sequences = dna_to_protein(complementary_strand)\n",
    "\n",
    "# Print protein sequences to see if they are different\n",
    "for frame, sequence in protein_sequences.items():\n",
    "    print(f\"First 10 aminoacids for original DNA {frame}: {sequence[:10]}\")\n",
    "\n",
    "for frame, sequence in complementary_protein_sequences.items():\n",
    "    print(f\"First 10 aminoacids of sequence for complementary DNA {frame}: {sequence[:10]}\")"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 8,
   "id": "c80839ba-2e7d-49c6-8498-992f59dbac8f",
   "metadata": {},
   "outputs": [],
   "source": [
    "#the dictionary above has made two different protein sequences because I had to call it twice for my two inputs\n",
    "#so we make a single dictionary here to use it in the next function \n",
    "\n",
    "def merge_protein_sequences(protein_sequences, complementary_protein_sequences):\n",
    "    merged_sequences = {}\n",
    "\n",
    "    # Merge protein sequences for the original DNA strand\n",
    "    for frame, sequence in protein_sequences.items():\n",
    "        merged_sequences[f'original_{frame}'] = sequence\n",
    "\n",
    "    # Merge protein sequences for the complementary DNA strand\n",
    "    for frame, sequence in complementary_protein_sequences.items():\n",
    "        merged_sequences[f'complementary_{frame}'] = sequence\n",
    "\n",
    "    return merged_sequences\n",
    "    \n",
    "# Merge the protein sequences\n",
    "merged_sequences = merge_protein_sequences(protein_sequences, complementary_protein_sequences)\n",
    "\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 9,
   "id": "37067905-f101-4f7c-bcd3-5d3a5c30c820",
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Summary for original_orf1:\n",
      "Number of sequences: 69\n",
      "Max sequence length: 50\n",
      "Average sequence length: 19.782608695652176\n",
      "\n",
      "Summary for original_orf2:\n",
      "Number of sequences: 74\n",
      "Max sequence length: 69\n",
      "Average sequence length: 17.77027027027027\n",
      "\n",
      "Summary for original_orf3:\n",
      "Number of sequences: 64\n",
      "Max sequence length: 46\n",
      "Average sequence length: 19.078125\n",
      "\n",
      "Summary for complementary_orf1:\n",
      "Number of sequences: 52\n",
      "Max sequence length: 39\n",
      "Average sequence length: 12.615384615384615\n",
      "\n",
      "Summary for complementary_orf2:\n",
      "Number of sequences: 52\n",
      "Max sequence length: 58\n",
      "Average sequence length: 10.596153846153847\n",
      "\n",
      "Summary for complementary_orf3:\n",
      "Number of sequences: 53\n",
      "Max sequence length: 46\n",
      "Average sequence length: 11.037735849056604\n",
      "\n"
     ]
    }
   ],
   "source": [
    "##Identifying the most coding reading frame out of the six \n",
    "#first counting the aminoacids in between a START and STOP codon, then making a summary\n",
    "\n",
    "def count_sequences_between_codons(merged_sequences):\n",
    "    all_sequence_lengths = {}  # Dictionary to store sequence lengths for each ORF\n",
    "\n",
    "    for frame, merged_sequence in merged_sequences.items():\n",
    "        sequence_lengths = []\n",
    "        in_sequence = False\n",
    "        current_length = 0\n",
    "\n",
    "        for amino_acid in merged_sequence:\n",
    "            if amino_acid == 'M' and not in_sequence:\n",
    "                # Start of a new sequence\n",
    "                in_sequence = True\n",
    "                current_length = 0\n",
    "            elif amino_acid == '*' and in_sequence:\n",
    "                # End of the current sequence\n",
    "                in_sequence = False\n",
    "                sequence_lengths.append(current_length)\n",
    "            elif in_sequence:\n",
    "                # Inside a sequence, increment length\n",
    "                current_length += 1\n",
    "\n",
    "        # Store sequence lengths for the current ORF\n",
    "        all_sequence_lengths[frame] = sequence_lengths\n",
    "\n",
    "    return all_sequence_lengths\n",
    "\n",
    "sequence_lengths = count_sequences_between_codons(merged_sequences)\n",
    "\n",
    "# Print summary outputs for each ORF\n",
    "for frame, lengths in sequence_lengths.items():\n",
    "    print(f\"Summary for {frame}:\")\n",
    "    print(f\"Number of sequences: {len(lengths)}\")\n",
    "    print(f\"Max sequence length: {max(lengths) if lengths else 0}\")\n",
    "    print(f\"Average sequence length: {sum(lengths) / len(lengths) if lengths else 0}\\n\")\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "f5e33d10-c2e5-420b-b992-6c99fe47187f",
   "metadata": {},
   "outputs": [],
   "source": []
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3 (ipykernel)",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.10.12"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 5
}

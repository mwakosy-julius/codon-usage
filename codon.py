import streamlit as st
from collections import defaultdict

# Define the genetic code dictionary and ordered codons in table format
GENETIC_CODE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

# Define codon groups in table order for display
GENETIC_CODE_TABLE_ORDER = [
    ['TTT', 'TCT', 'TAT', 'TGT'],
    ['TTC', 'TCC', 'TAC', 'TGC'],
    ['TTA', 'TCA', 'TAA', 'TGA'],
    ['TTG', 'TCG', 'TAG', 'TGG'],
    ['CTT', 'CCT', 'CAT', 'CGT'],
    ['CTC', 'CCC', 'CAC', 'CGC'],
    ['CTA', 'CCA', 'CAA', 'CGA'],
    ['CTG', 'CCG', 'CAG', 'CGG'],
    ['ATT', 'ACT', 'AAT', 'AGT'],
    ['ATC', 'ACC', 'AAC', 'AGC'],
    ['ATA', 'ACA', 'AAA', 'AGA'],
    ['ATG', 'ACG', 'AAG', 'AGG'],
    ['GTT', 'GCT', 'GAT', 'GGT'],
    ['GTC', 'GCC', 'GAC', 'GGC'],
    ['GTA', 'GCA', 'GAA', 'GGA'],
    ['GTG', 'GCG', 'GAG', 'GGG']
]

# Function to calculate codon usage
def calculate_codon_usage(sequence):
    sequence = sequence.upper()
    codon_count = defaultdict(int)
    amino_acid_count = defaultdict(int)
    
    # Count occurrences of each codon
    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i+3]
        if codon in GENETIC_CODE:
            amino_acid = GENETIC_CODE[codon]
            codon_count[codon] += 1
            amino_acid_count[amino_acid] += 1

    # Calculate codon usage data
    codon_usage = {}
    for codon, amino_acid in GENETIC_CODE.items():
        count = codon_count[codon]
        total_usage = amino_acid_count[amino_acid]
        relative_usage = count / total_usage if total_usage > 0 else 0
        percentage = (count / sum(codon_count.values()) * 100) if sum(codon_count.values()) > 0 else 0.0
        codon_usage[codon] = {
            "amino_acid": amino_acid,
            "relative_usage": round(relative_usage, 2),
            "percentage": round(percentage, 1),
            "count": count
        }
    return codon_usage

# Streamlit UI
st.title("Codon Usage Calculator")
st.write("Enter a DNA sequence to calculate codon usage statistics.")

# Input DNA sequence
sequence_input = st.text_area("DNA Sequence", "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")

if st.button("Calculate Codon Usage"):
    result = calculate_codon_usage(sequence_input)

    # Display results in a genetic code table format
    st.write("## Codon Usage Results")
    st.write("Each entry shows codon, amino acid, relative usage, percentage, and count.")

    # Format table according to GENETIC_CODE_TABLE_ORDER
    table_data = []
    for row in GENETIC_CODE_TABLE_ORDER:
        row_data = []
        for codon in row:
            usage = result.get(codon, {"amino_acid": "", "relative_usage": 0.0, "percentage": 0.0, "count": 0})
            entry = f"{codon} {usage['amino_acid']: <2} {usage['relative_usage']:.2f} {usage['percentage']:>5.1f}% ({usage['count']:>5})"
            row_data.append(entry)
        table_data.append(row_data)

    # Display the table in Streamlit
    for row in table_data:
        st.write(" | ".join(row))
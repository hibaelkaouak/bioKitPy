# Dictionary of valid nucleotides for DNA and RNA
VALID_NUCLEOTIDES = {'DNA': ['A', 'T', 'C', 'G'],
                     'RNA': ['A', 'U', 'C', 'G']
}

# Dictionary for DNA complement mapping
COMPLEMENT_DNA = { 'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

CODON_TABLE = {'RNA': {
    'UCA': 'S', 'UCC': 'S', 'UCG': 'S', 'UCU': 'S',    # Serina
    'UUC': 'F', 'UUU': 'F',    # Fenilalanina
    'UUA': 'L', 'UUG': 'L',    # Leucina
    'UAC': 'Y', 'UAU': 'Y',    # Tirosina
    'UAA': '*', 'UAG': '*',    # Stop
    'UGC': 'C', 'UGU': 'C',    # Cisteina
    'UGA': '*',    # Stop
    'UGG': 'W',    # Triptofano
    'CUA': 'L', 'CUC': 'L', 'CUG': 'L', 'CUU': 'L',    # Leucina
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCU': 'P',    # Prolina
    'CAC': 'H', 'CAU': 'H',    # Histidina
    'CAA': 'Q', 'CAG': 'Q',    # Glutamina
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGU': 'R',    # Arginina
    'AUA': 'I', 'AUC': 'I', 'AUU': 'I',    # Isoleucina
    'AUG': 'M',    # Methionina
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACU': 'T',    # Treonina
    'AAC': 'N', 'AAU': 'N',    # Asparagina
    'AAA': 'K', 'AAG': 'K',    # Lisina
    'AGC': 'S', 'AGU': 'S',    # Serina
    'AGA': 'R', 'AGG': 'R',    # Arginina
    'GUA': 'V', 'GUC': 'V', 'GUG': 'V', 'GUU': 'V',    # Valina
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCU': 'A',    # Alanina
    'GAC': 'D', 'GAU': 'D',    # Acido Aspartico
    'GAA': 'E', 'GAG': 'E',    # Acido Glutamico
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGU': 'G'     # Glicina
}, 'DNA': {'TCA': 'S', 'TCC': 'S', 'TCT': 'S', 'TCG': 'S',    # Serina
    'TTC': 'F', 'TTT': 'F',    # Fenilalanina
    'TTA': 'L', 'TTG': 'L',    # Leucina
    'TAC': 'Y', 'TAT': 'Y',    # Tirosina
    'TAA': '*', 'TAG': '*',    # Stop
    'TGC': 'C', 'TGT': 'C',    # Cisteina
    'TGA': '*',    # Stop
    'TGG': 'W',    # Triptofano
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',    # Leucina
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',    # Prolina
    'CAC': 'H', 'CAT': 'H',    # Histidina
    'CAA': 'Q', 'CAG': 'Q',    # Glutamina
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',    # Arginina
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I',    # I  soleucina
    'ATG': 'M',    # Methionina
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',    # Treonina
    'AAC': 'N', 'AAT': 'N',    # Asparagina
    'AAA': 'K', 'AAG': 'K',    # Lisina
    'AGC': 'S', 'AGT': 'S',    # Serina
    'AGA': 'R', 'AGG': 'R',    # Arginina
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',    # Valina
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',    # Alanina
    'GAC': 'D', 'GAT': 'D',    # Acido Aspartico
    'GAA': 'E', 'GAG': 'E',    # Acido Glutamico
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G'     # Glicina
    }
}

VALID_AMINO_ACIDS = ['A', 'R', 'N', 'D', 'B', 'C', 'Q', 'E', 'Z', 'G',
                     'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

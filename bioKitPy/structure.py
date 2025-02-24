# Dictionary of valid nucleotides for DNA and RNA
VALID_NUCLEOTIDES = {'DNA': ['A', 'T', 'C', 'G'],
                     'RNA': ['A', 'U', 'C', 'G']
}

# Dictionary for DNA complement mapping
COMPLEMENT_DNA = { 'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

CODON_TABLE = {
    'UCA': 'S', 'UCC': 'S', 'UCU': 'S',    # Serina
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
}
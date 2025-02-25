import itertools

from .structure import *
import random
from collections import *

class Sequence:
    '''A class to rapresent a biological sequence (DNA or RNA).

    This class allows you to store, validate, and manipulate DNA or RNA sequences.
    '''
    def __init__(self, seq=None, seq_type='DNA', length=20, label='No Label'):
        """
        Initializes a Sequence object.
        If `seq` is None, a random sequence of length `length` is generated.

        Parameters:
        - __seq (str or None): The nucleotide sequence (default: None for random generation).
        - __seq_type (str): 'DNA' or 'RNA' (default: 'DNA').
        - label (str): A label for the sequence (default: 'No Label').
        - __length (int): Length of the sequence if generated randomly (default: 10).

        Raises:
        - ValueError: If seq_type is not 'DNA' or 'RNA'.
        - ValueError: If seq contains invalid characters for the given seq_type.
        - ValueError: If length is not a positive integer when generating a random sequence.
        """
        self.__seq_type = seq_type.upper()
        if self.__seq_type not in VALID_NUCLEOTIDES:
            raise ValueError("Invalid sequence type: must be 'DNA' or 'RNA'.")
        if seq is None or seq == "":
            if not isinstance(length, int) or length <= 0:
                raise ValueError("Length must be a positive integer greater than 0.")
            self.__seq = ''.join(random.choices(VALID_NUCLEOTIDES[self.__seq_type], k=length))
            print(self.__seq)
        else:
            self.__seq = seq.upper()
            if not self.__validate():
                raise ValueError(f"Invalid sequence: contains non-{self.__seq_type} nucleotides.")
        self.label = label
        self.__length = len(self.__seq)

    def __validate(self):
        """Checks if the sequence contains only valid nucleotides."""
        return set(VALID_NUCLEOTIDES[self.__seq_type]).issuperset(self.__seq)

    def __check_indices(self, start, end):
        """Validates the indices for substring operations"""
        if not isinstance(start, int) or (end is not None and not isinstance(end, int)):
            raise ValueError("Start and end must be integers or None.")
        if start < 0 or (end is not None and end > len(self.__seq)):
            raise ValueError(f"Start must be >= 0 and end must be <= {len(self.__seq)}.")
        if end is not None and start >= end:
            raise ValueError("Start must be less than end")

    def __str__(self):
        return f"{self.__seq}"

    def __repr__(self):
        return f"{self.__seq}"

    def get_info(self):
        """Returns information about the sequence."""
        return (f"Label: {self.label}\n"
                f"Type: {self.__seq_type}\n"
                f"Seqeunce: {self.__seq}\n"
                f"Length: {self.__length}"
                )

    def get_seq_type(self):
        """Returns the type of the sequence (DNA or RNA)."""
        return (f"Type: {self.__seq_type}")

    def get_seq(self):
        return self.__seq

    def count_nucleotides(self, start=0, end=None):
        """Counts the occurrences of each nucleotide in the sequence."""
        return {nuc: self.__seq[start:end].count(nuc) for nuc in VALID_NUCLEOTIDES[self.__seq_type]}

    def nucleotide_percentage(self, start=0, end=None):
        """Returns the percentage of each nucleotide in the sequence."""
        self.__check_indices(start, end)
        total_length = len(self.__seq[start:end])
        if total_length == 0:
            raise ValueError("Cannot calculate nucleotide percentage for an empty sequence.")
        return {nuc: round(self.__seq[start:end].count(nuc)/total_length, 3) for nuc in VALID_NUCLEOTIDES[self.__seq_type]}

    def gc_content(self, start=0, end=None):
        """Returns the GC content percentage of the sequence."""
        self.__check_indices(start, end)
        return round((self.__seq[start:end].count('G') + self.__seq[start:end].count('C')) / len(self.__seq[start:end]) * 100, 3)

    def get_complement(self, start=0, end=None):
        """Returns the complement of the DNA sequence (or a part of it)."""
        if self.__seq_type != 'DNA':
            raise ValueError("Complement function is only available for DNA sequences.")
        self.__check_indices(start, end)
        return ''.join(COMPLEMENT_DNA[nuc] for nuc in self.__seq[start:end])

    def get_reverse_complement(self, start=0, end=None):
        """Returns the reverse complement of the DNA sequence (or a part of it)."""
        if self.__seq_type != 'DNA':
            raise ValueError("Complement function is only available for DNA sequences.")
        self.__check_indices(start, end)
        return ''.join(COMPLEMENT_DNA[nuc] for nuc in self.__seq[start:end][::-1])

    def transcribe_dna(self, start=0, end=None):
        """Transcribes a DNA sequence into RNA (T → U)."""
        if self.__seq_type != 'DNA':
            raise ValueError("Transcribe function is only available for DNA sequences.")
        self.__check_indices(start, end)
        tmp_seq = self.__seq[start:end]
        return ''.join(tmp_seq.replace('T', 'U'))

    def translate_sequence(self, start=0, end=None):
        """Translates a DNA or RNA sequence into a protein sequence."""
        self.__check_indices(start, end)
        if end is None:
            end = self.__length
        return ''.join([CODON_TABLE[self.__seq_type][self.__seq[i:i+3]] for i in range(start,end-2,3)])

    def codon_usage(self, amino_acid, start=0, end=None):
        """
        Computes the relative frequency of codons that code for a given amino acid.

        Parameters:
        ----------
        amino_acid : str
            The amino acid for which to compute codon usage (one-letter code).
        start : int, optional
            The starting index of the sequence to analyze (default is 0).
        end : int, optional
            The ending index of the sequence to analyze (default is None, meaning full sequence).

        Returns:
        -------
        dict
            A dictionary where keys are codons and values are their relative frequency.

        Raises:
        ------
        ValueError
            If an invalid amino acid is provided.
        """
        self.__check_indices(start, end)
        amino_acid = amino_acid.upper()
        if not amino_acid in VALID_AMINO_ACIDS:
            raise ValueError("Invalid amino acid")
        if end is None:
            end = self.__length
        codon_counts = defaultdict(int)
        total_codons = 0
        for i in range(start, end-2, 3):
            codon = self.__seq[i:i+3]
            if CODON_TABLE[self.__seq_type][codon] == amino_acid:
                codon_counts[codon] += 1
                total_codons += 1
        codon_frequency = dict(codon_counts)
        if total_codons != 0:
            for codon in codon_frequency:
                codon_frequency[codon] /= total_codons
        return codon_frequency

    def locate_motif(self, motif, start=0, end=None):
        self.__check_indices(start, end)
        if end is None:
            end = self.__length
        motif_len = len(motif)
        if motif_len > len(self.__seq[start:end]):
            raise ValueError(f"Motif length ({motif_len}) exceeds selected sequence length ({len(self.__seq[start:end])}).")
        positions = []
        for i in range(start, end):
            if self.__seq[i] == motif[0]:
                if self.__seq[i:i+motif_len] == motif:
                    positions.append(i+1)
        return positions

    @staticmethod
    def find_consensus(sequences, type):
        matrix = []
        number_collums = sequences[0].__length
        for seq in sequences:
            if number_collums != seq.__length:
                raise ValueError()
            matrix.append(seq.__seq)
        number_rows = len(matrix)

        matrix_output = []
        for j in range(number_collums):
            current_collum = ''
            for i in range(number_rows):
                current_collum += matrix[i][j]
            matrix_output.append({nuc: current_collum.count(nuc) for nuc in VALID_NUCLEOTIDES[type]})
        print(''.join([max(dict,key=dict.get) for dict in matrix_output]))
        for nucl in VALID_NUCLEOTIDES[type]:
            string_output = nucl + ": "
            string_output += ' '.join([str(matrix_output[col][nucl]) for col in range(number_collums)])
            print(string_output)

    #@staticmethod
    #def rabbits(n, k):
        #f n == 1 or n == 2:
            #return 1
        #curr, prev = 1, 1
        #for _ in range(3, n + 1):
            #prev, curr = curr, curr + k * prev
        #return(curr)

    @staticmethod
    def build_overlap_graph(sequences, k=3):
        if k <= 0:
            raise ValueError("k must be a positive integer")
        edges = defaultdict(list)
        for s_id, t_id in itertools.combinations(sequences, 2):
            if sequences[s_id].__length < k or sequences[t_id].__length < k:
                raise ValueError(f"Sequences '{s_id}' or '{t_id}' are shorter than k={k}. Required: k ≤ sequence length.")
            s, t = sequences[s_id].__seq, sequences[t_id].__seq
            if s[-k:] == t[:k]:
                edges[s_id].append(t_id)
            if t[-k:] == s[:k]:
                edges[t_id].append(s_id)
        for start_node, adj_list in edges.items():
            for end_node in adj_list:
                print(f"{start_node} {end_node}")
        return edges














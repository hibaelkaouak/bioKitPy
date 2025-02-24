from .structure import *
import random
from collections import Counter

class Sequence:
    '''A class to rapresent a biological sequence (DNA or RNA).

    This class allows you to store, validate, and manipulate DNA or RNA sequences.
    '''
    def __init__(self, seq=None, seq_type='DNA', label='No Label', length=20):
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

    def






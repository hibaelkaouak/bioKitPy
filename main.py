# bioKitPy module
from bioKitPy import Sequence


strand = Sequence("", "DNA")
print(strand.get_info())
print(strand.count_nucleotides())
str2 = Sequence(strand.transcribe_dna(), "RNA")
print(strand.translate_sequence())
print(str2.translate_sequence())
print(strand.codon_usage('K'))

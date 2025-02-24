# bioKitPy module
from bioKitPy import Sequence

strand = Sequence("", "RNA")
print(strand.get_info())
print(strand.count_nucleotides())
print(strand.nucleotide_percentage())
print(strand.gc_content())
print(strand.nucleotide_percentage(7))
print(strand.gc_content(7))

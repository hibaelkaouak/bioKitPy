# bioKitPy module
from bioKitPy import Sequence, Utilities


#strand = Sequence("AATTTATTGTATATTGTATCGAATATGGCTGCATATTGTATTGCCGATTGTATGATATTGTATCAACTGCATTGTATACATTGTATTATTGTATAATTGTATATTGTATAACCATATTGTATACGGATTGTATGTTTATTGTATTATTGTATTCGAATTGTATTTAACAAGCGATTGTATATGGATTGTATATTGTATGCTGTATGCTTGATTGTATCGAATTGTATGATTGTATGTAGATTGTATATTGTATATTGTATGAATTGTATCACATTGTATGATTGTATGGATTGTATCGATTGTATTATTGTATCATAACAATTGTATATTGTATCATGGGCCCAATTGTATATTGTATATTGTATAATTGTATATTGTATATTGTATCATTGTATATTGTATATTGTATTCGTATTGTATTCATTGTATCTATTGTATATCTTTATTGTATAAATTGTATTCGTAATTGTATATATTGTATGATTGCAAAATTGTATATTGTATGAATCATTGTATCATTGTATCATTGTATTATTGTATACGATTGTATAATTGTATGATATTGTATCTATTGTATATATTGTATGAGAATTGTATACGCATTGTATATTGTATAAATTGTATCACCATTGTATTATTGTATGATTGTATATTGTATGATTGTATAAGATTGTATATATTGTATTGTATTGTATTAGATTGTATTTAAATTGTATTTGAAGATTGTATTCGGATTGTATTGCATCGTAATTGTATTTATTGTATTGACGATTGTATCTCGGATTGTATGGAATTGTATGACCATTGTATATTGTATATTGTATCTATAATTGTAT", "DNA")
#print(' '.join(str(item) for item in strand.locate_motif("ATTGTATAT")))
dictstrands = Utilities.parse_fasta("/Users/hibaelkaouak/bioKitPy/tests/rosalind_grph.txt")
Sequence.build_overlap_graph(dictstrands, 3)
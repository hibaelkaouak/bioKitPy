from .sequence import Sequence
class Utilities:

    @staticmethod
    def parse_fasta(file_path, type='DNA'):
        with open(file_path, 'r') as f:
            sequences = {}
            sequence_id = None
            current_seq = ""
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if sequence_id:
                        sequences[sequence_id] = Sequence(current_seq)
                    sequence_id = line.strip('>')
                    current_seq = ""
                else:
                    current_seq += line
        #for key, value in sequences.items():
            #print(value)
            sequences[sequence_id] = Sequence(current_seq)
            sequence_id = line.strip('>')
        return sequences
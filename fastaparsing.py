import fastaparser
import sys

class IdNotFoundException(Exception):

    def __init__(self, message="Id not found"):
        self.message = message
        super().__init__(self.message)


def get_ids():
    return [x.id for x in parser]

def print_length(seq):
    print(seq.id, len(seq.sequence_as_string()))
    

def print_all_lengths():
    for seq in parser:
        print_length(seq)

def print_desc_and_seq(seq):
    print('Description:', seq.description)
    print('Sequence:', seq.sequence_as_string())
    
    
def get_certain_sequence_by_index(index : int):
    return get_certain_sequence(ids[index])
    
    

def get_certain_sequence(name):
    first = next(filter(lambda x: x.id == name, parser), None)
    if first == None:
        raise IdNotFoundException()
    return first
         

def return_subseq(seq, lb:int, ub:int):
    return (seq.sequence_as_string()[3:10])


with open(sys.argv[1]) as fasta_file:
    parser = fastaparser.Reader(fasta_file)
    ids = get_ids()
    seq = get_certain_sequence_by_index(1)
    

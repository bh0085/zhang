nt_list = list('ATGC')
complement = {'G':'C',
              'A':'T',
              'T':'A',
              'C':'G',
              'N':'N'}

def reverse_complement(seq):
    return ''.join([complement[s] for s in seq[::-1]])

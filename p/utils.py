nt_list = list('ATGC')
complement = {'G':'C',
              'A':'T',
              'T':'A',
              'C':'G'}

def reverse_complement(seq):
    return [complement[s] for s in seq[::-1]]

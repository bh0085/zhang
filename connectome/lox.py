import Bio.Seq as seq
loxdm =  'taccgttcgtatagcatacattatacgaacggta'
loxp =   'ataacttcgtatagcatacattatacgaagttat'
lox66 =  'ataacttcgtatagcatacattatacgaacggta'
lox71 =  'taccgttcgtatagcatacattatacgaagttat'

def oligos_1():
    overhangs = ['CATG', 'ACAA']
    bc =  'GATGATTGA'
    kozak =  'gccacc' 
    start = 'atg'
    fwd =overhangs[0] + loxp + bc + kozak + start + lox71
    rev =  seq.reverse_complement(loxp+bc+kozak+start+lox71+overhangs[1])
    
    print fwd.upper()
    print rev.upper()


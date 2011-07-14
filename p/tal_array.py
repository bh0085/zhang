#!/usr/bin/env python
'''
Package of programs for my TAL characterization array and 
for the mismatch tolerance analysis that will lead up to it.

'''

import re

reporters = {'T8':'TACCANTNANTATA',
             'T4':'TAAGANTNANTATA',
             'T5':'TATGANTNANTATA',
             'T6':'TAGGANTNANTATA'}
ids = {'T8':'C','T4':'A','T5':'T','T6':'G'}
nt_list = list('ATGC')
complement = {'G':'C',
              'A':'T',
              'T':'A',
              'C':'G'}

def reverse_complement(seq):
    return [complement[s] for s in seq[::-1]]

def mismatch_reporters():
    enzymes = {'fwd': ['CTAGA', 'G'],
               'rev': ['GATCC', 'T']}

    primers, names = [], []
    for k, r in reporters.iteritems():
        for m_ct in 0, 2, 3:
            if m_ct > 0: nts = [n for n in nt_list if n!=ids[k]]
            else: nts = [ids[k]]
            for mm in nts:
                rep = re.sub('N', mm ,r, m_ct)
                rep = re.sub('N', ids[k], rep)
                
                fwd = list(enzymes['fwd'])
                fwd.insert(1,rep)
                rev = list(enzymes['rev'])
                rev.insert(1,''.join(reverse_complement(rep)))

                names.append('{0}_{1}mm={2}_rev'.format(k, m_ct, mm))
                primers.append(''.join(rev))
                names.append('{0}_{1}mm={2}_fwd'.format(k, m_ct, mm))
                primers.append(''.join(fwd))
              
    return dict([(n,p) for n,p in zip(names, primers)])

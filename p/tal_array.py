#!/usr/bin/env python
'''
Package of programs for my TAL characterization array and 
for the mismatch tolerance analysis that will lead up to it.

'''

import re, os
import utils as zutils
import compbio.config as cfg
import align as ali

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


def find_wells():
    mm = mismatch_reporters()
    ls = load_seqs()
    successes = {}
    for k,v in mm.iteritems():
        if not 'fwd' in k: continue
        successes[k] = [  (k1,v1) 
                      
                          for k1,v1 in ls.iteritems() 
                          if v[6:-1] in zutils.reverse_complement(v1)]
    print 
    print
    print
    print '              '+'\n              '.join(sorted([str((k,len([elt[0] for elt in v]),[elt[0] for elt in v])) for k,v in successes.iteritems()]))
    return successes
def load_seqs(seq_dir = cfg.dataPath('zhang/tal_array/seqs')):
    seqs = {}
    for f in os.listdir(seq_dir):
        fopen = open(os.path.join(seq_dir,f))
        lines = fopen.readlines()
        seqs[lines[0][1:].strip()] = ''.\
            join([l.strip() for l in lines[1:]]) 
        fopen.close()
    return seqs

def align_seqs():
    ls = load_seqs()
    return [ali.run_with_seqs(\
            ls.values()[0].replace('N',''),\
                ls.values()[i].replace('N',''))[1] 
            for i in range(5)]

def ali_stringarray(ls):
    return      [ali.run_with_seqs(\
            ls[0].replace('N',''),\
                ls[i].replace('N',''))[1] 
                 for i in range(len(ls))]

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
                rev.insert(1,zutils.reverse_complement(rep))

                names.append('{0}_{1}mm={2}_rev'.format(k, m_ct, mm))
                primers.append(''.join(rev))
                names.append('{0}_{1}mm={2}_fwd'.format(k, m_ct, mm))
                primers.append(''.join(fwd))
              
    return dict([(n,p) for n,p in zip(names, primers)])

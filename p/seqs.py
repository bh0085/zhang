#!/usr/bin/env python
'''
nt.py

Contains a few utilities for looking up nucleotide level info 
for zhang lab sequences of interest.

'''

import compbio.config as cfg
from Bio import SeqIO


ntfiles = {
    'nrx':cfg.dataPath('sequences/zhang/nt/nrx1_human_nt.gb'),
    'nlg':cfg.dataPath('sequences/zhang/nt/nlg1_human_nt.gb')                   
    }
aafiles = {
    'nrx':cfg.dataPath('sequences/zhang/aa/nrx1_human_aa.gb'),
    'nlg':cfg.dataPath('sequences/zhang/aa/nlg1_human_aa.gb')       
    }

def get_seq( name,  aa = True):
    seq = SeqIO.parse(open(aafiles[name]), 'genbank') if aa \
        else SeqIO.parse(open(ntfiles[name]),'genbank')

    return seq

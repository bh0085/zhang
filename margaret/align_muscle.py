#!/usr/bin/env python

import Bio.Align as align, Bio.SeqIO as sio, Bio, Bio.AlignIO as aio
import sys
import subprocess as spc

seqs = [align.SeqRecord(align.Seq(           'ATGATGGGGGATGATG')),\
              align.SeqRecord(align.Seq(           'ATGATGATGATG')),\
              ]





m_proc = spc.Popen('muscle -clw' ,
                   stdin = spc.PIPE,
                   stdout = spc.PIPE,
                   stderr=spc.PIPE,
                   shell = True)


  
sio.write(seqs, m_proc.stdin, "fasta")
m_proc.stdin.close()
align = aio.read(m_proc.stdout, "clustal")



print align

#!/usr/bin/env python
import re, sys, itertools as it
import zhang.p.utils as zutils

sw = sys.stdout.write
verbose = False
def readfiles_old(files):
    primers = {}


    for f in files:
        fopen = open(f)
        lines = fopen.read()
        pattern = re.compile('(?P<title>\S+.*)\n'+
                             '(?P<names>\S+.*)\n'+
                             '(?P<labels>\S+.*)\n'+
                             '(?P<seqs>(\S+.*\n)*)', re.M)
        hits = [(m.group('title'),
                 [elt.strip() for elt in m.group('names').split(',')],
                 [elt.strip() for elt in m.group('seqs').splitlines() ])
                for m in re.finditer(pattern, lines)]

        primers.update(dict(zip(*[[i.strip() 
                                   for i in it.chain(*zip(*hits)[j])] 
                                  for j in [0,1]])))
    return primers

def readfiles(files):
    primers = {}


    for f in files:
        fopen = open(f)
        lines = fopen.read()
        pattern = re.compile('(?P<title>\S+.*)\n'+
                             '(?P<names>\S+.*)\n'+
                             '(?P<labels>\S+.*)\n'+
                             '(?P<seqs>(\S+.*\n)*)', re.M)
        hits = [(m.group('title'),
                 [elt.strip()  for elt in m.group('names').split(',')],
                 [elt.strip() for elt in m.group('seqs').splitlines() ])
                for m in re.finditer(pattern, lines)]

        primers.update(dict([(e[0], 
                              {'fwd':e[2][0],'rev':e[2][1],
                               'fname':e[1][0], 'rname':e[1][1]},
                              )
                             for e in hits]))
        #primers.update(dict(zip(*[[i.strip() 
        #                           for i in it.chain(*zip(*hits)[j])] 
        #                          for j in [0,1]])))
    return primers

def bsai_overhangs(primers):
    return list(it.chain(*[((p['fname'],p['fwd'][7:11].upper()), 
                       (p['rname'],''.join(zutils.reverse_complement(\
                            p['rev'][7:11].upper()))))
                      for p in primers.values()]))
def bsai_print_final(primers):
    return list(it.chain(*[((p['fname'], 'GAGTAT' + p['fwd'].upper()),
                          (p['rname'], 'GAGTAT' + p['rev'].upper()))
                          for p in primers.values()]))

def primer_seqs(primers):
    return dict(list(it.chain(*[[(p['fname'],p['fwd']),
                                (p['rname'],p['rev'])]
                                for p in primers.values()])))

def bsai(files, verbose = verbose):

    for f in files:
      primers = readfiles([f])
      
      #CHECK BSAI SITES IN ONE FILE AT A TIME
      sw('\nFILE {0}:\n\n'.format(f))
      sw('checking for bsai sites:\n')

      overhangs = bsai_overhangs(primers)
      pseqs = primer_seqs(primers)
      ntot = len(primers)
      np = len(pseqs)

      #DO ANY PRIMERS LACK A BSAI SITE?
      fails = [(k,v) 
               for k,v in pseqs.iteritems() 
               if v.upper()[:6] != 'GGTCTC']
      if len(fails) == 0:
          sw('SUCCESS\n')
      else: 
          sw('FAILURE'+'{0}\n'.format(fails)); return(1)

      #CHECK FOR CORRECT LINKERS
      sw('\nchecking overhangs:\n')
      
      if verbose: sw('\n'.join([o[1] for o in overhangs])+'\n')
      
      #EVERY LINKER IS PAIRED ONCE (VECTOR IS PRIMED)
      if len(set([o[1] for o in overhangs])) == np / 2.:
          sw('SUCCESS\nevery primer appears paired:\n')
          sw('{0}\n'.format('\n'.join([' <--> '.join(['{0:5}'.format(elt[0]) for elt in list(g)] )
                                       for k , g in it.groupby(
                              sorted(overhangs, key = lambda x:x[1]),
                              key = lambda x: x[1])])))
      
      #TWO LINKERS ARE UNPAIRED (NO VECTOR PRIMERS)
      elif len(set([o[1] for o in overhangs])) ==np / 2. + 1:
          sw('SUCCESS\nevery primer appears paired:\n')
          sw('{0}\n'.format('\n'.join([' <--> '.join(['{0:5}'.format(elt[0]) for elt in list(g)] )
                                       for k , g in it.groupby(
                              sorted(overhangs, key = lambda x:x[1]),
                              key = lambda x: x[1])])))
      
      #ANOTHER CASE
      else:
          sw('FAILURE\n')
          grps =dict([(k,list(g))  for k , g in it.groupby(
                      sorted(overhangs, key = lambda x:x[1]),
                      key = lambda x: x[1])])
          sw('{0}\n'.format(
                  '\n'.join( '{0} share {1}'.format(
                          (', '.join([e[0] for e in elt[1]])),
                          elt[0]) 
                             for elt in grps.iteritems() 
                             if len(elt[1]) > 2)))
          if __name__ != '__main__': raise Exception()
          return(1)
      
      sw('\nEverything appears to be OK\n')

      
    primers = readfiles(files)  
    
    sw('printing primers:\n')
    sw('\n'.join(sorted(['>{0}\t{1}'.format(*p) for p in bsai_print_final(primers)])))
      
      

def usage(): 
    print '''usage: primers [run_method]
run_methods:
  print:   print the primers from the primer format files given to stdin
  bsai:    check for bsai compatibility, make ensure no repeated primers, 
           add 6 random nucleotides.

'''

def main():
    inp =  sys.stdin.read()
    files = re.findall(re.compile('\S+'),inp)

    args = sys.argv[1:]
    if args[0] == 'print':

        primers = readfiles(files)
        sys.stdout.write('\n'.join(['>{0:20}{1}'.format(k,v) 
                                    for k,v in primers.iteritems()]))
    elif args[0] == 'bsai':
        exit(bsai(files))
if __name__ == '__main__':
    main()

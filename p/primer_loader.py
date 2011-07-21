#!/usr/bin/env python
import re, sys, itertools as it
import zhang.p.utils as zutils

sw = sys.stdout.write
verbose = False
def readfiles_old(files):
    '''..deprecated. Around just in case.'''
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
    '''Read primer formatted files:
...{
                    (--empty line--)
product name        (arbitrary)
fwd_name, rev_name  (two entries, comma seperated)
[arbitrary notes    (exactly one line!)]
fwd_primer sequence (5' - 3')
rev_primer sequence (5' - 3')
                    (---empty line--)
}...
'''
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
    '''Given a dict of pcr products with primers, grab the 4 nt BsaI overhangs'''
    return list(it.chain(*[((p['fname'],p['fwd'][7:11].upper()), 
                       (p['rname'],''.join(zutils.reverse_complement(\
                            p['rev'][7:11].upper()))))
                      for p in primers.values()]))
def bsai_print_final(primers):
    '''Given a dict of pcr products with primers, print a list of primers with
extra nucleotides appended to allow cutting'''
    return list(it.chain(*[((p['fname'], 'GAGTAT' + p['fwd'].upper()),
                          (p['rname'], 'GAGTAT' + p['rev'].upper()))
                          for p in primers.values()]))

def primer_seqs(primers):
    '''Get a list of tuples containing all primer names/sequences'''
    return dict(list(it.chain(*[[(p['fname'],p['fwd']),
                                (p['rname'],p['rev'])]
                                for p in primers.values()])))

def sort_good(input_strings):
    return sorted(input_strings, key = lambda x: 
                  re.search(re.compile('[0-9]+'),x) and int(re.search(re.compile('[0-9]+'),x).group() ) +1 
                  or x)

def bsai(files, verbose = verbose):
    '''Check for bsaI compatibility, uniqueness etc and finally print a list of
primers with six nucleotides appended'''

    #CHECK EACH PRIMER FILE SEPERATELY
    for f in files:
      sw('\nFILE {0}:\n\n'.format(f))
      primers = readfiles([f])

      #CHECK BSAI SITES IN ONE FILE AT A TIME
      sw('checking for bsai sites:\n')
      overhangs = bsai_overhangs(primers)
      pseqs = primer_seqs(primers)
      ntot = len(primers);  np = len(pseqs)

      #DO ANY PRIMERS LACK A BSAI SITE?
      fails = [(k,v) 
               for k,v in pseqs.iteritems() 
               if v.upper()[:6] != 'GGTCTC']
      if len(fails) == 0: sw('SUCCESS\n')
      else:sw('FAILURE'+'{0}\n'.format(fails)); return(1)

      #CHECK LINKER PAIRING
      sw('\nchecking overhangs:\n')
      #EVERY LINKER IS PAIRED ONCE (VECTOR BACKBONE IS PRIMED)
      if verbose: sw('\n'.join([o[1] for o in overhangs])+'\n')     
      if len(set([o[1] for o in overhangs])) == np / 2.:
          sw('SUCCESS\nevery primer appears paired:\n')
          sw('{0}\n'.format('\n'.join([' <--> '.join(['{0:5}'.format(elt[0]) for elt in list(g)] )
                                       for k , g in it.groupby(
                              sorted(overhangs, key = lambda x:x[1]),
                              key = lambda x: x[1])])))
      
      #TWO LINKERS ARE UNPAIRED (NO VECTOR PRIMERS)
      elif len(set([o[1] for o in overhangs])) ==np / 2. + 1:
          sw('SUCCESS\nevery primer appears paired:\n')
          sw('{0}\n'.format('\n'.join(\
                      sort_good([

                                  ' <--> '.join(\
                                      sort_good(['{0:5}'.format(elt[0]) 
                                       for elt in list(g)]
                                      ))
                                  
                                  for k , g in it.groupby( sorted(overhangs, key = lambda x:x[1]),
                                                           key = lambda x: x[1])]))))

      #SOMETHING WEIRD IS GOING ON
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
          #IF THIS IS BEING RUN INTERACTIVELY, CRY FOUL!
          if __name__ != '__main__': raise Exception()
          return(1)
      
      sw('\nEverything appears to be OK\n')

    #NOW READ ALL OF THE FILES AT ONCE AND PRINT
    primers = readfiles(files)  
    sw('printing primers:\n')
    sw('\n'.join(sorted(['>{0}\t{1}'.format(*p) for p in bsai_print_final(primers)])))
      
      

def usage(): 
    print '''
usage: 

echo [filename1] ... [filenameN] | primers [run_method]

run_methods:
  print:   print the primers from the primer format files
  bsai:    check for bsai compatibility, make ensure no repeated primers, 
           add 6 random nucleotides. Print!

'''
    return 1

def main():
    '''Run from the command line. Takes one argument and stdin.'''

    args = sys.argv[1:]
    if len(args) != 1: exit(usage())
    if args[0] =='help':
        usage(); exit(0)
    elif args[0] == 'print':
        inp =  sys.stdin.read()
        files = re.findall(re.compile('\S+'),inp)
        primers = readfiles(files)
        sys.stdout.write('\n'.join(['>{0:20}{1}'.format(k,v) 
                                    for k,v in primers.iteritems()]))
    elif args[0] == 'bsai':
        inp =  sys.stdin.read()
        files = re.findall(re.compile('\S+'),inp)
        exit(bsai(files))

if __name__ == '__main__':
    main()

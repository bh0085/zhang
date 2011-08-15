import re, StringIO
import Bio.SeqIO as sio
def parse_gbfile(fpath = '/Users/bh0085/zhang/genbank/pCMV-Syn_N1_G1-10_NRX_MBD_(from ccdb).gbk'):

    fdata = open(fpath).readlines()

    print
    print 'subbing out any spaces in the locus line:'
    loc_re = re.compile('(LOCUS\s+)(.*[^ ])(\s+\d+ bp)')
    loc = re.search(loc_re, fdata[0]).groups()[1]
    new_loc = re.sub('[ \(\)]','_', loc)
    fdata[0] = re.sub(loc_re,'\g<1>{0}\g<3>'.format(new_loc), fdata[0])
                        
    print 'trying to load the file with no more mods'

    try:
        ioelt = StringIO.StringIO(
            '\n'.join(fdata)
            )
        seq = sio.parse(ioelt, 'genbank').next()
    except Exception, e:

        print 'initial parse failed, trying seqbuilder spacing fix'
        o_elt = [i for i,l in enumerate(fdata) if l[0:6]=='ORIGIN'][0]
        ioelt = StringIO.StringIO(\
            '\n'.join([l 
                       if ( not re.search(re.compile('^\s+\d+'),l)) or i < o_elt 
                       else l[2:] 
                       for i,l in enumerate(fdata)]))
        seq =  sio.parse(ioelt, 'genbank').next()
    return seq

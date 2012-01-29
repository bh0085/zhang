import align as ali
import compbio.config as cfg
import os 
import Bio.SeqIO as sio

default = cfg.dataPath('zhang/sequencing/piggybac')

def run_directory(directory = default):
    refs = [os.path.join(root, f) 
            for root, dirs, files in  os.walk(os.path.join(directory, 'refs')) 
            for f in files]
    all_refs = [sio.parse( f, format = 'fasta').next()
                for f in refs]
    

    results = [os.path.join(root, f)  
               for root, dirs, files in  os.walk(os.path.join(directory, 'results')) 
               for f in files]
    result_sequences = load_genewiz_seqs(results)

    for r in all_refs:
        print 'Aligning to reference: {0} '.format(r)
        out = align_seqs(r, result_sequences)
        for k,v in out.iteritems():
            print 'result {0}: {1}'.format(k,v)
    
def load_genewiz_seqs(filenames):
    seqs = {}
    for f in filenames:
        fopen = open(f)
        lines = fopen.readlines()
        seqs[lines[0][1:].strip()] = ''.\
            join([l.strip() for l in lines[1:]]) 
        fopen.close()
    return seqs

def align_seqs(ref,seqs):
    seqs = dict([(k,seqs[k]) for k in seqs.keys()[0:1]])
    return dict([(k, ali.run_with_seqs(\
            ref.seq.tostring().upper().replace('N',''),\
                seqs[k].replace('N',''))[1]) 
            for k,v in seqs.iteritems()])

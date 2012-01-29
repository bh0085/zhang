import cb.config as cfg
import os, re, itertools as it, json
import structs 
import cb.utils.plots as mp
import cb.utils.colors as mycolors
from numpy import *
import numpy as np

import Bio.PDB as PDB

figt = 'struct_{0}.pdf'

droot = cfg.dataPath('structs')
pdb_ann = os.path.join(droot, '3UGM.pdb')
pdb_plain = os.path.join(droot, '3UGM_.pdb')
pqr_file = os.path.join(droot, '3UGM.pqr')


def lookuptable():
    info = [['Alanine', 'Ala', 'A'],
 ['Arginine', 'Arg', 'R'],
 ['Asparagine', 'Asn', 'N'],
 ['Aspartic acid', 'Asp', 'D'],
 ['Cysteine', 'Cys', 'C'],
 ['Glutamic acid', 'Glu', 'E'],
 ['Glutamine', 'Gln', 'Q'],
 ['Glycine', 'Gly', 'G'],
 ['Histidine', 'His', 'H'],
 ['Isoleucine', 'Ile', 'I'],
 ['Leucine', 'Leu', 'L'],
 ['Lysine', 'Lys', 'K'],
 ['Methionine', 'Met', 'M'],
 ['Phenylalanine', 'Phe', 'F'],
 ['Proline', 'Pro', 'P'],
 ['Serine', 'Ser', 'S'],
 ['Threonine', 'Thr', 'T'],
 ['Tryptophan', 'Trp', 'W'],
 ['Tyrosine', 'Tyr', 'Y'],
 ['Valine', 'Val', 'V']]
    return dict([(e[1].upper(), e[2]) for e in info])



def parse_all():
    pdbr = structs.parse_pdb()
    pqrr = structs.parse_pqr()
    chains = []
    for i in range(len(pdbr)):
        b = pdbr[i]
        q = pqrr[i]
        chains.append({})
        for k in b.keys():
            chains[i][k] = {'pdb':b[k],
                         'pqr':q[k]}

    return chains
            
    
def compute_seq():
    fopen = open(pdb_ann)
    struct = structs.read_struct(fopen)
    fopen.close()
    residues = struct.get_residues()
    tal_len = 1048
    tstr = ['*'] * tal_len
    lookup = lookuptable()
    for r in residues:
        if r.parent.id == 'A':
            if not r.resname.upper() in lookup.keys(): continue
            tstr[r.id[1]-1] = lookup[r.resname.upper()]
    return ''.join(tstr)
              
    
def seq():
    return '''************************************************************************************************************************************************************************KPKVR*******E*LV******AHIVALSQHPAALGTVAVTYQHIITAL**ATHEDIVGVGKQWSGARALEALLTDAGELRGPPLQLDTGQLVKIAKRGGVTAMEAVHASRNALT**PLNLTPAQVVAIASNNGGKQALETVQRLLPVLCQAHGLTPAQVVAIASHDGGKQALETMQRLLPVLCQAHGLPPDQVVAIASNIGGKQALETVQRLLPVLCQAHGLTPDQVVAIASHGGGKQALETVQRLLPVLCQAHGLTPDQVVAIASHDGGKQALETVQRLLPVLCQAHGLTPDQVVAIASNGGGKQALETVQRLLPVLCQAHGLTPDQVVAIASNGGKQALETVQRLLPVLCQAHGLTPDQVVAIASHDGGKQALETVQRLLPVLCQTHGLTPAQVVAIASHDGGKQALETVQQLLPVLCQAHGLTPDQVVAIASNIGGKQALATVQRLLPVLCQAHGLTPDQVVAIASNGGGKQALETVQRLLPVLCQAHGLTPDQVVAIASNGGGKQALETVQRLLPVLCQAHGLTQVQVVAIASNIGGKQALETVQRLLPVLCQAHGLTPAQVVAIASHDGGKQALETVQRLLPVLCQAHGLTPDQVVAIASNGGGKQALETVQRLLP******GLTQEQVVAIASNNGGKQALETVQRLLPVLCQAHGLTPDQVVAIASNGGGKQALETVQRLLPVLCQAHGLTPAQVVAIASNIGGKQALETVQRLLPVLCQDHGLTLAQVVAIASNIGGKQALETVQRLLPVLCQ**GLTQDQVVAIASNIGGKQALETVQRLLPVLCQDHGLTPDQVVAIASNIGGKQALETVQRLLPVLCQDHGLTLDQVVAIASNGGKQALETVQRLLPVLCQDHGLTPDQVVAIASNSG'''

def rvds():
    tal_seq= seq()

    rep_seeds = [e
                 for e in re.compile('VVAIAS')
                 .finditer(tal_seq)]
    ends = [e.end() for e in rep_seeds]
    rvds = ends
    
    #rep1 = 197
    #rvds = []
    #for rep_count in range(22):
    #    rvds.append(rep1 + 34 * rep_count)
        
    return rvds[:-1]
    
def rvd_residues(chain, nbd_size = 0):
    rs = rvds()
   
    out = []
    for r in rs:
        keys = range( r - nbd_size+1, r+3 + nbd_size)
        out.append([chain[k] for k in keys])
        
    return out

def get_nn_locs():
    tags = [ seq()[r:r+5] for r in rvds() ]
    NN_locs = [rvds()[i]  
               for  i, e in enumerate(tags) 
               if e[:1] == 'N' and e[3] == 'G' and e[1] == 'N' ]
    return NN_locs

def run():
    
    chains = parse_all()
    
    c0 = chains[0]


    ksrt =sorted(c0.keys())
    charges = array([ sum([float(e['charge'])  for e in c0[k]['pqr']]) for k in ksrt])
    coords = array([ mean([e.get_coord() for e in c0[k]['pdb']], 0 ) for k in ksrt])
   
 
    c1,c2 = chains[1:]
    
    strand_charges = []
    strand_coords = []
    for c in [c1, c2]:
        ksrt =sorted(c.keys())
        strand_charges.append(array([ sum([float(e['charge'])  for e in c[k]['pqr']]) for k in ksrt]))
        strand_coords.append(array([ mean([e.get_coord() for e in c[k]['pdb']], 0 ) for k in ksrt]))
  
    k1 = sorted([nt for nt in c1])
    k2 = sorted([nt for nt in c2])
    s1_atoms = list(it.chain(*[ [e.get_coord() for e in c1[k]['pdb'] ] for k in k1]))
    s2_atoms = list(it.chain(*[ [e.get_coord() for e in c2[k]['pdb'] ] for k in k2]))
    dna_atoms = []
    dna_atoms.extend(s1_atoms)
    dna_atoms.extend(s2_atoms)
    dna_atoms = array(dna_atoms)
    
    #nearest neighbor params:
    kres =0 
    katoms_dna = 3
    kres_atoms = 3

    rvd_res = rvd_residues(c0,kres)

    xs = []
    ys = []
    cs = []
    ss = []
    ecs = []
    rdists = []

    rvd_groups = [];
    for i, rvd in enumerate(rvd_res):
        for r in rvd:
            atoms = array([e.get_coord() 
                           for e in r['pdb']])
            dists = sum(square(atoms),1)[:,newaxis] + \
                sum(square(dna_atoms),1)[newaxis,:] - \
                2 *sum(atoms[:,newaxis,:] * dna_atoms[newaxis,:,:],2)
            atom_srt_dists = np.sort(dists, 1)            
            atom_knn_avgdist = np.mean(atom_srt_dists[:,:katoms_dna],1)
            res_srt_dists = np.sort(atom_knn_avgdist)
            res_k_avgdist = res_srt_dists[:kres_atoms]

            xs.append(mean(atoms[:,0]))
            ys.append(mean(atoms[:,1]))
            #colors =  array([1,0,0]) * 1/atom_knn_avgdist[:,newaxis]
            cs.append(1/res_k_avgdist)
            rdists.append(res_k_avgdist)
            ss.append(50)
            ecs.append('none')
            rvd_groups.append(i)

    show_helix = False        
    if show_helix:
        cs = array(cs)
        cs /= np.max(cs)
        f = mp.fignum(1, (12,12))
        ax = f.add_subplot(111)
        ax.scatter(xs,ys,c = cs, s= ss, edgecolor = ecs)
        f.savefig(mp.figpath(figt.format('tal_rvd_neighborhoods')))
    
    rvd_dists =  [(k,[e[1] for e in list(g)]) 
                  for k,g in it.groupby(zip(rvd_groups,rdists), lambda x: x[0])]
    rs = rvds()
    tags = [ seq()[r:r+2] for r in rs ]
    nt = len(set(tags))
    
    tag_idx_in_ct = dict([(e,i) for i,e in enumerate(set(tags))])
    rvd_ct_map = dict([(i,tag_idx_in_ct[e]) for i,e in enumerate(tags)])
    ct = mycolors.getct(nt)

    f = mp.fignum(3, (12,12))
    ax= f.add_subplot(111)
    ax.set_xlabel('linear distance')
    ax.set_ylabel('nearest neighbor distance to DNA')
    labels_set = set([])
    for k, g in rvd_dists:
        if tags[k] in labels_set:
            ax.plot(g, color = ct[rvd_ct_map[k]])

        else:
            labels_set.add(tags[k])
            print 'labelling'
            ax.plot(g, color = ct[rvd_ct_map[k]],label = tags[k])
    ax.legend()
        
    f.savefig(mp.figpath(figt.format('tal_rvd_distances')))
        

    #plot_charges(coords, charges, strand_coords)
    
    return


def plot_charges(coords, charges, strand_coords):

    f0 = mp.fignum(1, (6,6))
    ax = f0.add_subplot(111)
    
    colors = [ 'red' if q > 0 else 'blue' for q in charges]
    ax.scatter(*coords[:,:2].T, c = colors, s = 50, zorder = 5)
    ax.scatter(*strand_coords[0][:,:2].T, c = 'black', alpha = 1, zorder = -1)
    ax.scatter(*strand_coords[1][:,:2].T, c = 'black', alpha = 1, zorder = -1)

    ax.scatter(*strand_coords[0][:,:2].T, c = 'black', alpha = .5, zorder = 6)
    ax.scatter(*strand_coords[1][:,:2].T, c = 'black', alpha = .5, zorder = 6)
    #ax.scatter(*txy.T, c = 'red')
    
    fp = mp.figpath(figt.format('tal_xy_charge_scatter'))
    f0.savefig(fp)
    
    return

def load_bonds(nbs = False):
    bfiles = dict([(rvd,cfg.dataPath('structs/bonds/3UGM/3UGM_{0}_{1}_NOFIX.bonds'\
                                         .format(rvd, 'NBS' if nbs else 'NONBS')))
              for rvd in ['NN', 'NH', 'NK']])
    all_bonds = {}
    for k,bf in bfiles.iteritems():
        fopen = open(bf)
        data = fopen.readlines()
        bonds = []
        for line in data:
            split = line.index('))') + 2
            bonds.append(dict(zip( ['atoms','dist'], 
                                   [json.loads(line[0:split].replace('(','[').replace(')',']').replace("'", '"')),json.loads(line[split:])]
                                   )))
        all_bonds[k.upper()] = bonds
        
    return all_bonds
        
def load_pdbfiles(nbs = False):
    bfiles = dict([(rvd, cfg.dataPath('structs/3ugm/3UGM_{0}_{1}_NOFIX.pdb'\
                                          .format(rvd,'NBS' if nbs else 'NONBS')))
                    for rvd in ['NN', 'NH', 'NK']])
    structs = {}
    for rvd, fname in bfiles.iteritems():
        fopen = open(fname)
        pparse = PDB.PDBParser(PERMISSIVE = 0)
        struct = pparse.get_structure('tmp', fopen)
        fopen.close()
        structs[rvd] = struct
        
    return structs

def get_bonded_structs():
    nbs = False
    bonds = load_bonds(nbs = nbs)
    structs = load_pdbfiles(nbs = nbs)

    struct_bonds = {}
    rvds = structs.keys()
    for rvd in rvds:
        bonded = bonds[rvd]
        struct = structs[rvd]
        
        atoms = [ a for a in struct.get_atoms() ]
        rvd_bonds = []
        for b in bonded: 
            ids_bonded =[ e[1] for e in b['atoms'] ]
            atoms_bonded = [
                [ a for a in atoms if a.serial_number == id] 
                for id in ids_bonded]
            dist = b['dist']
            assert len(atoms_bonded[0]) == len(atoms_bonded[1]) == 1
            atoms_bonded = [atoms_bonded[0][0], atoms_bonded[1][0]]
            ab = atoms_bonded 
            adist = sqrt(sum(square( ab[0].get_coord() - ab[1].get_coord())))

            
            rvd_bonds.append({'atoms':ab,
                              'dist':dist,
                              'adist':adist})
        struct_bonds[rvd] = rvd_bonds
    return struct_bonds


def analyze_bonded_structs(bonded_structs):
    tal_seq = seq()
    rvds = bonded_structs.keys()
    contexts = {}
    for rvd in rvds:
        contexts[rvd] = []
        sbond = bonded_structs[rvd]
        for bond in sbond:
            atoms = bond['atoms']
            res = atoms[0].parent
            pos = res.get_id()[1]
            ctx = tal_seq[pos-3:pos+3]

            if pos == 302:
                raise Exception()
            contexts[rvd].append('{0}: {1} in {2}'.format(pos, tal_seq[pos-1], ctx))
            

    return contexts
            

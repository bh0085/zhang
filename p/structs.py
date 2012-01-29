import Bio.PDB as PDB
import Bio.PDB.PDBIO as PDBIO
from Bio.PDB.PDBIO import Select

import cb.config as cfg
import os, itertools as it, re


droot = cfg.dataPath('structs')
pdb_ann = os.path.join(droot, '3UGM.pdb')
pdb_plain = os.path.join(droot, '3UGM_.pdb')
pqr_file = os.path.join(droot, '3UGM.pqr')



def test(filename = 'g_princeps.pdb'):
    st = read_struct(open(filename))
def read_struct(fopen):
    pparse = PDB.PDBParser(PERMISSIVE = 0)
    struct = pparse.get_structure('tmp', fopen)
    return struct
def res_pos(filename = 'g_princeps.pdb'):
    raise Exception()

def extract_chain(struct,fout):
    io=PDBIO()
    io.set_structure(struct)
    io.save(fout, select = ChainSelector('A'))
    return 

class ChainSelector(Select):
    id_accepted = None
    def __init__(self, id_accepted):
        self.id_accepted = id_accepted
    def accept_chain(self, chain):
        if chain.get_id()==self.id_accepted:
            return 1
        else:
            return 

def parse_pqr():
    fopen = open(pqr_file)
    pqlines = [l for l in fopen.readlines() if not l[0:6] == 'REMARK'][::-1]
    chains = []
    pqr_cols= ['atom_type', 'atom_number', 'atom_type', 'res_name', 'res_number', 'x', 'y', 'z', 'charge', 'radius']
    current_chain = []
    while pqlines:
        l = pqlines.pop()
        if l[0:3] == 'TER':
            chains.append(current_chain)
            current_chain = []
        elif l[0:3] == 'END':
            break
        else:
            current_chain.append(dict(zip(pqr_cols, re.split(re.compile('\s+'),l))))
        

    chain_residues = [];
    for i, c in enumerate(chains):
        chain_residues.append(dict([(k, list(g)) for k, g in it.groupby(c, lambda x: int(x['res_number']))]))
    
        

    return chain_residues
        
def parse_pdb():
    fopen = open(pdb_ann)
    struct = read_struct(fopen)
    fopen.close()

    chains = struct.get_chains()
    chain_atoms = []
    for c in chains:
        chain_atoms.append(c.get_atoms())
    
    chain_residues = [];
    for i, c in enumerate(chain_atoms):
        chain_residues.append(dict([(k, list(g)) for k, g in it.groupby(c, lambda x: int(x.parent.id[1]))]))
    
    return chain_residues

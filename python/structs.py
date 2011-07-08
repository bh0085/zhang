import Bio.PDB as PDB
import Bio.PDB.PDBIO as pio
from numpy import *
import numpy as np

def struct_center(filename = 'g_princeps.pdb'):
    pparse = PDB.PDBParser()
    struct == pparse.get_structure('tmp',open(filename))
    coords =  array([a.get_coord() for a in struct.get_atoms()])
    center = np.mean(coords, 0)
    print center
    
    

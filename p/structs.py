import Bio.PDB as PDB
import Bio.PDB.PDBIO as pio


def test(filename = 'g_princeps.pdb'):
    st = read_struct(open(filename))
def read_struct(fopen):
    pparse = PDB.PDBParser(PERMISSIVE = 0)
    struct = pparse.get_structure('tmp', fopen)
    return struct
def res_pos(filename = 'g_princeps.pdb'):
    raise Exception()


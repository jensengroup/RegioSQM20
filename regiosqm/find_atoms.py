import numpy as np
from rdkit import Chem

def find_identical_atoms(rdkit_mol, atom_list):
    
    len_list = len(atom_list)
    
    atom_rank = list(Chem.CanonicalRankAtoms(rdkit_mol, breakTies=False))
    for idx, atom in enumerate(rdkit_mol.GetAtoms()):
        if atom.GetIdx() in atom_list[:len_list]:
            sym_atoms = [int(atom_idx) for atom_idx, ranking in enumerate(atom_rank) if ranking == atom_rank[idx] and atom_idx not in atom_list] 
            atom_list.extend(sym_atoms)
    return atom_list


def remove_identical_atoms(rdkit_mol, name_list, smiles_list, atom_list):
    idx_list = []
    rank_kept = []
    atom_rank = list(Chem.CanonicalRankAtoms(rdkit_mol, breakTies=False))
    for idx, atom in enumerate(atom_list):
        if atom_rank[atom] not in rank_kept:
            rank_kept.append(atom_rank[atom])
            idx_list.append(idx)
    
    name_list = np.array(name_list)[idx_list].tolist()
    smiles_list = np.array(smiles_list)[idx_list].tolist()
    atom_list = np.array(atom_list)[idx_list].tolist()
    
    return name_list, smiles_list, atom_list



if __name__ == "__main__":
    
    smi = 'n1ccc(nc1Nc1ccc(cc1)C#N)Nc1c(cc(cc1C)/C=C/C#N)C'
    corr = [2,8,12,20] #18 is missing
    print(corr, find_identical_atoms(smi, corr))

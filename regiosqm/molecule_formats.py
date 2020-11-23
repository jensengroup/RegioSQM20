import subprocess

from rdkit.Chem import rdmolfiles, AllChem # TorsionFingerprints
from rdkit.ML.Cluster import Butina


def shell(cmd, shell=False):

    if shell:
        p = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    else:
        cmd = cmd.split()
        p = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    output, err = p.communicate()
    return output


def convert_xyz_to_sdf(xyzfile, sdffile):

    shell(f'obabel -ixyz {xyzfile} -osdf -xf > {sdffile}', shell=True)

    return


def convert_sdf_to_xyz(sdffile, xyzfile):
    
    shell(f'obabel -isdf {sdffile} -oxyz -xf > {xyzfile}', shell=True)

    return


def get_bonds(sdf_file):

    isav = 0
    atoms = 0
    bond_list = []

    searchlines = open(sdf_file, 'r').readlines()

    for i, line in enumerate(searchlines):
        words = line.split() #split line into words
        if len(words) < 1:
            continue
        if i == 3:
           atoms = int(words[0])
           bonds = int(words[1])
        if i > atoms+3 and i <= atoms+bonds+3:
           atom_1 = int(words[0])
           atom_2 = int(words[1])
           if atom_2 > atom_1:
              bond_list.append(tuple((atom_1,atom_2)))
           else:
              bond_list.append(tuple((atom_2,atom_1)))

    bond_list.sort()

    return bond_list


def compare_sdf_structure(start, end):
    """
    Returns True if structures are the same

    Return False if there has been a proton transfer
    """

    bond_start = get_bonds(start)
    bond_end = get_bonds(end)

    return bond_start == bond_end


def find_unique_confs(best_conformers, mol_files, threshold=0.5):
    """ Clustering conformers with RDKit's Butina algorithm
    to find unique conformer from a list of .sdf files
    using either heavy-atom root mean square deviation (RMSD) 
    or heavy-atom torsion fingerprint deviation (TFD) """

    rdkit_mol = next(rdmolfiles.ForwardSDMolSupplier(mol_files[0], sanitize=False, removeHs=True))
    for mol_file in mol_files[1:]:
        mol = next(rdmolfiles.ForwardSDMolSupplier(mol_file, sanitize=False, removeHs=True))
        rdkit_mol.AddConformer(mol.GetConformer(),assignId=True)

    # calculate difference matrix
    diffmat = AllChem.GetConformerRMSMatrix(rdkit_mol, prealigned=False) #threshold=0.5, sanitize=False, load AllChem
    # diffmat = TorsionFingerprints.GetTFDMatrix(rdkit_mol) #threshold=0.01, sanitize=True, load TorsionFingerprints

    # Cluster conformers
    num_confs = rdkit_mol.GetNumConformers()
    clt = Butina.ClusterData(diffmat, num_confs, threshold,
                             isDistData=True, reordering=True)

    # Get unique conformers
    centroid_idx = [c[0] for c in clt] # centroid indexes.
    unique_best_conformers = [best_conformers[i] for i in centroid_idx]
    
    return unique_best_conformers



if __name__ == "__main__":
    
    import sys

    print(compare_sdf_structure(sys.argv[1], sys.argv[2]))

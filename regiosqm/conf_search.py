import os
import numpy as np
import concurrent.futures
import subprocess
from operator import itemgetter

from ase.units import Hartree, mol, kcal

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

import molecule_formats as molfmt

# CPU usage
num_cpu = 4 #number of conformers being handled simultaneously (see also regiosqm.py)
OMP_NUM_THREADS = '1'
MKL_NUM_THREADS = '1'

# xTB path
XTBHOME = '/home/Ree/anaconda/envs/rdkit-env'
XTBPATH = '/home/Ree/anaconda/envs/rdkit-env/share/xtb'
MANPATH = '/home/Ree/anaconda/envs/rdkit-env/share/man'
LD_LIBRARY_PATH = '/home/Ree/anaconda/envs/rdkit-env/lib'


def run_xTB(xtb_args):

    global XTBHOME

    mol_obj, confname, chrg, method, solvent, original_path, dirprefix, use_sdf = xtb_args

    # Make a seperate folder for the calc and setup structure variables
    os.makedirs(f'{original_path}/{confname}', exist_ok=True)

    # Run calculation using previous optimized structure or write initial structure 
    if use_sdf:
        start_structure = f'{original_path}/../{dirprefix}_fast_conf_search/{confname}/xtbopt.xyz'
        final_structure = f'{original_path}/{confname}.xtbout.sdf'
    else:
        writer = Chem.rdmolfiles.SDWriter(f"{original_path}/{confname}/{confname}.sdf")
        writer.write(mol_obj)
        writer.close()

        molfmt.convert_sdf_to_xyz(f"{original_path}/{confname}/{confname}.sdf", f"{original_path}/{confname}/{confname}.xyz")

        start_structure = f'{original_path}/{confname}/{confname}.xyz'
        final_structure = f'{original_path}/{confname}.gfnff.sdf'

    # Run xTB calc
    cmd = f'{XTBHOME}/bin/xtb --gfn{method} {start_structure} --opt --gbsa {solvent} --chrg {chrg} --uhf 0'
    proc = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True, cwd=f'{original_path}/{confname}')
    output = proc.communicate()[0]
     
    # Save calc output
    with open(f'{original_path}/{confname}/{confname}.xtbout', 'w') as f:
        f.write(output)
    
    # Convert .xyz to .sdf using and check input/output connectivity
    if os.path.isfile(f'{original_path}/{confname}/xtbopt.xyz'):
        molfmt.convert_xyz_to_sdf(f'{original_path}/{confname}/xtbopt.xyz', final_structure) #convert optimized structure
    else:
        print(f'WARNING! xtbopt.xyz was not created => calc failed for {confname}')
        energy = 60000.0
        os.rename(start_structure, f"{original_path}/{confname}/xtbopt.xyz")
        return confname, energy

    same_structure = molfmt.compare_sdf_structure(f'{original_path}/../{dirprefix}_fast_conf_search/{confname}/{confname}.sdf', final_structure)
    if not same_structure:
        print(f'WARNING! Input/output mismatch for {confname}')
        energy = 60000.0
        os.remove(f"{original_path}/{confname}/xtbopt.xyz")
        os.rename(start_structure, f"{original_path}/{confname}/xtbopt.xyz")
        return confname, energy

    # Search for the molecular energy
    for i, line in enumerate(output.split('\n')):
        if 'TOTAL ENERGY' in line:
            energy = line.split()[3]

    try: #check if the molecular energy was found.
        energy = float(energy) * Hartree * mol/kcal #convert energy from Hartree to kcal/mol
    except Exception as e:
        print(e, confname)
        energy = 60000.0
        os.remove(f"{original_path}/{confname}/xtbopt.xyz")
        os.rename(start_structure, f"{original_path}/{confname}/xtbopt.xyz")

    return confname, energy
    


def fast_conf_search(rdkit_mol, name='mol', min_conf=1, rot_conf=3, max_conf=20, chrg=0, 
                    method='ff', solvent='Methanol', conf_cutoff=3, calc_dir=os.getcwd(), dirprefix=''): 
                    
    """ 
    Method options available: 0 (GFN-0), 1 (GFN-1), 2 (GFN-2), ff (GFN-FF)
    conf_cutoff has to be in units of kcal/mol. 
    """
    
    global XTBHOME
    global XTBPATH
    global MANPATH
    global LD_LIBRARY_PATH

    global OMP_NUM_THREADS
    global MKL_NUM_THREADS
    global num_cpu

    # Set env parameters for xTB
    os.environ['XTBHOME'] = XTBHOME
    os.environ['XTBPATH'] = XTBPATH
    os.environ['MANPATH'] = MANPATH
    os.environ['LD_LIBRARY_PATH'] = LD_LIBRARY_PATH

    os.environ["OMP_NUM_THREADS"] = OMP_NUM_THREADS
    os.environ['MKL_NUM_THREADS'] = MKL_NUM_THREADS
    
    # RDkit conf generator
    rdkit_mol = Chem.AddHs(rdkit_mol)

    rot_bond = rdMolDescriptors.CalcNumRotatableBonds(rdkit_mol)
    confs = min(min_conf + rot_conf*rot_bond, max_conf)

    # p = AllChem.ETKDGv3()
    # p.randomSeed = 90
    # p.useSmallRingTorsions=True
    # p.ETversion=2
    # p.useExpTorsionAnglePrefs=True
    # p.useBasicKnowledge=True
    # AllChem.EmbedMultipleConfs(rdkit_mol, numConfs=confs, params=p)
        
    AllChem.EmbedMultipleConfs(rdkit_mol, numConfs=confs,
            useExpTorsionAnglePrefs=True,
            useBasicKnowledge=True, ETversion=2) #, randomSeed=90)
    
    # Run calculations in parallel
    original_path = f'{calc_dir}/{dirprefix}_fast_conf_search'
    os.makedirs(original_path, exist_ok=True)
    os.chdir(original_path)
    xtb_args = [(Chem.Mol(rdkit_mol, False, idx), f"{name}-{idx}", chrg, method, solvent, original_path, dirprefix, False) for idx in range(rdkit_mol.GetNumConformers())]
    with concurrent.futures.ThreadPoolExecutor(max_workers=num_cpu) as executor:
        results = executor.map(run_xTB, xtb_args)
    os.chdir('..')

    # Save the results
    energy_list = list()
    conformer_list = list()
    for result in results:
        confname, energy = result
        conformer_list.append(confname)
        energy_list.append(energy)

    # Find the conformers below cutoff
    rel_energies = np.array(energy_list) - np.min(energy_list) #covert to relative energies
    below_cutoff = (rel_energies <= conf_cutoff).sum() #get number of conf below cutoff
    conf_tuble = list(zip(conformer_list,rel_energies)) #make a tuble
    conf_tuble = sorted(conf_tuble, key=itemgetter(1))[:below_cutoff] #get only the best conf below cutoff
    best_conformers = [item[0] for item in conf_tuble]

    # Find only unique conformers
    mol_files = [f'{original_path}/{confname}.gfnff.sdf' for confname in best_conformers]
    best_conformers = molfmt.find_unique_confs(best_conformers, mol_files, threshold=0.5)
    
    return best_conformers



def acc_conf_search(rdkit_mol, name='mol', min_conf=1, rot_conf=3, max_conf=20, chrg=0, 
                    method='1', solvent='Methanol', conf_cutoff=3, calc_dir=os.getcwd(), dirprefix=''):
    """  
    Run fast_conf_search and then higher level calculations.
    Method options available: 0 (GFN-0), 1 (GFN-1), 2 (GFN-2), ff (GFN-FF)
    """

    global num_cpu

    best_conformers = fast_conf_search(rdkit_mol, name=name, min_conf=min_conf, rot_conf=rot_conf, max_conf=max_conf, chrg=chrg, 
                                    method='ff', solvent=solvent, conf_cutoff=conf_cutoff, calc_dir=calc_dir, dirprefix=dirprefix)

    # Run calculations in parallel
    original_path = f'{calc_dir}/{dirprefix}_acc_conf_search'
    os.makedirs(original_path, exist_ok=True)
    os.chdir(original_path)
    xtb_args = [(None, confname, chrg, method, solvent, original_path, dirprefix, True) for confname in best_conformers]
    with concurrent.futures.ThreadPoolExecutor(max_workers=num_cpu) as executor:
        results = executor.map(run_xTB, xtb_args)
    os.chdir('..')

    # Get only the lowest energy conformer
    best_energy = 10e10
    for result in results:
        confname, energy = result

        if energy < best_energy:
            best_energy = energy
            best_conf = confname

    return best_conf, best_energy



if __name__ == "__main__":
    
    import time
    start = time.perf_counter()

    input_smiles = 'c1c(c2cc(sc2)C)n[nH]c1'
    
    print(acc_conf_search(Chem.MolFromSmiles(input_smiles), name='mol', min_conf=1, rot_conf=3, max_conf=20, chrg=0, 
                        method='1', solvent='Methanol', conf_cutoff=3, calc_dir=os.getcwd(), dirprefix='neutral'))

    finish = time.perf_counter()
    print(f'Finished in {round(finish-start, 2)} second(s)')

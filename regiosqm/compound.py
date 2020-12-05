import os
from rdkit import Chem

from protonate import generate_protonated_smiles
from tautomer_generator import rdkit_taut_gen
from reorder_atoms import get_atoms_in_order
from conf_search import acc_conf_search


class comp:

    def __init__(self, name=None, smiles=None, min_conf=1, rot_conf=3, max_conf=20, conf_cutoff=3, measured_atoms=None, calc_dir=os.getcwd()):
        """Compound object for collecting all data."""

        self.name = name
        self.smiles = smiles
        self.input_mol = Chem.MolFromSmiles(smiles)
        self.neutral_charge = Chem.GetFormalCharge(self.input_mol)
        self.measured_atoms = measured_atoms

        self.tauts = {}
        
        # Generate tautomers
        taut_names, taut_smiles, taut_energies = rdkit_taut_gen(self.input_mol, name=name, min_conf=min_conf, rot_conf=rot_conf, max_conf=max_conf, 
                                                            chrg=self.neutral_charge, method='1', solvent='Methanol', conf_cutoff=conf_cutoff, taut_cutoff=15, 
                                                            calc_dir=calc_dir, dirprefix='neutral')
        taut_mols = [Chem.MolFromSmiles(smi) for smi in taut_smiles]
        taut_mols = get_atoms_in_order(self.input_mol, taut_mols)
        
        # Generate protonated states for each tautomers
        for i, taut_name  in enumerate(taut_names):

            self.tauts[taut_name] = {}
            self.tauts[taut_name]['smiles'] = taut_smiles[i] #this is the smiles of the neutral tautomer
            self.tauts[taut_name]['energy'] = taut_energies[i] #this is the energy of the neutral tautomer
            self.tauts[taut_name]['rel_energy'] = taut_energies[i] - min(taut_energies) #this is the relative energy of the neutral tautomer
            self.tauts[taut_name]['prot'] = {}

            cnames, csmiles, catoms = generate_protonated_smiles(taut_mols[i], taut_name)

            prot_names_list = list()
            prot_energy_list = list()
            for cname, csmile, catom in zip(cnames, csmiles, catoms):

                best_conf, best_energy = acc_conf_search(Chem.MolFromSmiles(csmile), name=cname, min_conf=min_conf, rot_conf=rot_conf, max_conf=max_conf, 
                                                        chrg=int(self.neutral_charge+1), method='1', solvent='Methanol', conf_cutoff=conf_cutoff, 
                                                        calc_dir=calc_dir, dirprefix='prot')
                prot_names_list.append(best_conf) #this is the name of one of the protonated tautomers
                prot_energy_list.append(best_energy) #this is the energy of one of the protonated tautomers
                
            self.tauts[taut_name]['prot']['reac_sites'] = catoms
            self.tauts[taut_name]['prot']['smiles_list'] = csmiles
            self.tauts[taut_name]['prot']['names'] = prot_names_list
            self.tauts[taut_name]['prot']['energies'] = prot_energy_list



if __name__ == "__main__":
    import sys
    import time
    start = time.perf_counter()
    
    comp = comp(name="comp547", smiles='c1(cc(ccc1)S(F)(F)(F)(F)F)N')

    print(comp.name)
    print(comp.tauts)

    finish = time.perf_counter()
    print(f'Finished in {round(finish-start, 2)} second(s)')

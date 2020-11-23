import os
import time
import numpy as np
import concurrent.futures
import operator
from functools import reduce

from find_atoms import find_identical_atoms
import molecule_svg as molsvg
import compound

__version__ = "20"

# CPU usage
num_cpu = 2 #number of molecules being handled simultaneously (see also conf_search.py)


def analyse_results(comp, e_cut, e_cut2, calc_dir, exam):

    # Arrays for output
    output = []

    # Arrays for drawing
    smiles_list = []
    name_list = []
    pred_atoms_list = []
    pred_atoms2_list = []

    # Find the winners
    for taut in list(comp.tauts.keys()):
        
        # Read energies
        heats = comp.tauts[taut]['prot']['energies']
        heats = np.array(heats)
        
        # # Calculate absolute affinities
        e_neutral = comp.tauts[taut]['energy']
        heats = heats - e_neutral

        # Calculate relative energies
        min_heats = np.min(heats)
        rel_heats = heats - min_heats

        # Check for correct termination of calc
        if min_heats == 60000.0 - e_neutral or e_neutral == 60000.0:
            output.append(f'WARNING! No output conformers match input for {taut} \n') #print
            continue

        # Read reactive center
        atoms = comp.tauts[taut]['prot']['reac_sites']
        atoms = np.array(atoms)

        # Primary winners
        winners = np.where( rel_heats < e_cut )
        winners = winners[0]
        pred_atoms = np.unique(atoms[winners]).tolist()
        pred_atoms = find_identical_atoms(comp.smiles, pred_atoms)
        pred_atoms_list.append(pred_atoms)

        # Secondary winners
        winners2 = np.where( rel_heats < e_cut2 )
        winners2 = winners2[0]
        pred_atoms2 = np.unique(atoms[winners2]).tolist()
        pred_atoms2 = find_identical_atoms(comp.smiles, pred_atoms2)
        pred_atoms2_list.append(pred_atoms2)

        # Append smiles and name to lists, where the name is related to the reactivity
        smiles_list.append(comp.tauts[taut]['smiles'])
        if abs(min_heats) >= 100: #kcal/mol
            name_list.append(f"High | {'{:.2f}'.format(comp.tauts[taut]['rel_energy'])}")
        elif abs(min_heats) <= 70: #kcal/mol
            name_list.append(f"Low | {'{:.2f}'.format(comp.tauts[taut]['rel_energy'])}")
        else:
            name_list.append(f"Medium | {'{:.2f}'.format(comp.tauts[taut]['rel_energy'])}")

        ## Print results ##
        output.append(f'{taut} ') #print

        if exam:
            measure = comp.measured_atoms
            if set(measure).issubset(pred_atoms):
                output.append(f'corr {measure} == ') #print
            elif set(measure).issubset(pred_atoms2):
                output.append(f'semi {measure} == ') #print
            else:
                output.append(f'fail {measure} == ') #print

        output.append(','.join([str(x) for x in pred_atoms]) + ' ') #print
        output.append(','.join([str(x) for x in pred_atoms2]) + "\n") #print

        output.append(f"taut info: {taut} {comp.tauts[taut]['smiles']} {comp.tauts[taut]['energy']} \n") #print

        prot_names = comp.tauts[taut]['prot']['names']
        prot_names = np.array(prot_names)
        prot_smiles = comp.tauts[taut]['prot']['smiles_list']
        prot_smiles = np.array(prot_smiles)
        for winner in winners:
            output.append(f'1> {prot_names[winner]} {prot_smiles[winner]} {atoms[winner]} {heats[winner]} \n') #print prot_name, prot_smiles, prot_reaction_center, prot_energy
            
        for winner in winners2:
            if winner in winners: continue
            output.append(f'2> {prot_names[winner]} {prot_smiles[winner]} {atoms[winner]} {heats[winner]} \n') #print prot_name, prot_smiles, prot_reaction_center, prot_energy

    # Draw results - save SVG
    if exam:
        try:
            pred_atoms_list_concat = np.unique(reduce(operator.concat, pred_atoms_list))
            pred_atoms2_list_concat = np.unique(reduce(operator.concat, pred_atoms2_list))
        except:
            output.append(f'CRITICAL WARNING! No predictions for {comp.name} \n\n') #print
            return output
        
        output.append(f'FINAL SUMMARY: {comp.name} ') #print
        if set(measure).issubset(pred_atoms_list_concat):
            comp_status = "comp_corr"
        elif set(measure).issubset(pred_atoms2_list_concat):
            comp_status = "comp_semi"
        else:
            comp_status = "comp_fail"
        output.append(f'{comp_status} {measure} == ') #print
        output.append(','.join([str(x) for x in pred_atoms_list_concat]) + ' ') #print
        output.append(','.join([str(x) for x in pred_atoms2_list_concat]) + '\n') #print
        
        result_svg = molsvg.generate_structure(comp.smiles, smiles_list, name_list, [pred_atoms_list, pred_atoms2_list], highlight_measure=measure)
    else:
        result_svg = molsvg.generate_structure(comp.smiles, smiles_list, name_list, [pred_atoms_list, pred_atoms2_list])
    output.append("\n") #print

    f_draw = open(f'{calc_dir}/{comp.name}.svg','w')
    f_draw.write(result_svg)
    f_draw.close()
    return output


def calculate_and_analyse(args):

    # Time process - Start 
    start = time.perf_counter()

    # Initialize
    name, smiles, min_conf, rot_conf, max_conf, conf_cutoff, measured_atoms, calc_dir, e_cut, e_cut2, exam = args

    # Run calculations
    comp = compound.comp(name=name, smiles=smiles, min_conf=min_conf, rot_conf=rot_conf, max_conf=max_conf, conf_cutoff=conf_cutoff, 
                        measured_atoms=measured_atoms, calc_dir=calc_dir)

    # Analyse the results
    final_result = analyse_results(comp, e_cut, e_cut2, calc_dir, exam)

    # Time process - Finish 
    finish = time.perf_counter()
    timing_output = f'{name} finished in {round(finish-start, 2)} second(s) \n'

    return (final_result, timing_output)


def run_app(smiles_filename, e_cut=1.0, e_cut2=3.0, min_conf=5, rot_conf=5, max_conf=60, conf_cutoff=3, exam=False):
    """ read smiles smiles_filename in the format
    <compound name> <smiles>

    generate tautomers for each compound in the smiles_filename
    run conformational search on both neutral and protonated tautomers,
    and return a list consisting of a class object with data for each compound.
    Hereafter, run analyses on results.
    """

    global num_cpu

    # Read data from input file
    f = open(smiles_filename, "r")

    comps_name = []
    comps_smiles = []
    comps_measured_atoms = []
    for line in f:
        words = line.split()
        comps_name.append(words[0])
        comps_smiles.append(words[1])

        if exam:
            measured_atoms = [int(x) for x in words[2].split(",")]
            comps_measured_atoms.append(find_identical_atoms(words[1], measured_atoms))
        else:
            comps_measured_atoms.append(None)
    f.close() 

    # Run calculations for every compound
    calc_dir = os.getcwd()
    args = [(name, comps_smiles[i], min_conf, rot_conf, max_conf, conf_cutoff, comps_measured_atoms[i], calc_dir, e_cut, e_cut2, exam) for i, name in enumerate(comps_name)]
    with concurrent.futures.ThreadPoolExecutor(max_workers=num_cpu) as executor:
        results = executor.map(calculate_and_analyse, args)

    # Write the results to file
    f = open(f"{calc_dir}/{smiles_filename.split('.')[0]}.result",'w')
    g = open(f"{calc_dir}/{smiles_filename.split('.')[0]}.timing",'w')
    for result in results:
        final_result, timing_output = result
        f.write(''.join(final_result))
        g.write(''.join(timing_output))
    f.close()
    g.close()
    
    return



def main():

    import argparse
    import sys

    description = """  """
    epilog = """ """

    parser = argparse.ArgumentParser(
                    usage='%(prog)s [options]',
                    description=description,
                    formatter_class=argparse.RawDescriptionHelpFormatter,
                    epilog=epilog)

    parser.add_argument('-v', '--version', action='version', version='RegioSQM' + __version__ + "\nhttps://github.com/jensengroup/regiosqm20")
    parser.add_argument('-r', '--run', nargs=1, action='store', metavar='smiles_filename', help="")

    parser.add_argument('-c', '--cut', action='store', type=float, metavar='N', default=1.0, help="Energy cut-off for 'corr' classification (kcal/mol)")
    parser.add_argument('-c2', '--cut2', action='store', type=float, metavar='N', default=3.0, help="Energy cut-off for 'semi' classification (kcal/mol)")
    
    parser.add_argument('-mic', '--min_conf', action='store', type=int, metavar='N', default=5, help="Minimum number of conf")
    parser.add_argument('-rc', '--rot_conf', action='store', type=int, metavar='N', default=5, help="Number of conf pr. rotatable bonds")
    parser.add_argument('-mac', '--max_conf', action='store', type=int, metavar='N', default=60, help="Maximum number of conf")
    parser.add_argument('-cc', '--conf_cutoff', action='store', type=int, metavar='N', default=3, help="Energy cut-off for fast conf search (kcal/mol)")

    parser.add_argument('-e', '--exam', action='store_true', help='Check results vs reaction centers in smiles file')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    if args.run:
        run_app(*args.run, e_cut=args.cut, e_cut2=args.cut2, min_conf=args.min_conf, 
                rot_conf=args.rot_conf, max_conf=args.max_conf, conf_cutoff=args.conf_cutoff, exam=args.exam)
        return

    return 0



if __name__ == "__main__":

    import time

    start = time.perf_counter()
    main()
    finish = time.perf_counter()
    print(f'Finished in {round(finish-start, 2)} second(s)')

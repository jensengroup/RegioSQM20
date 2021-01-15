# RegioSQM20

RegioSQM20 predicts the regioselectivity of electrophilic aromatic substitution reactions in heteroaromatic systems.
The reactive sites are identified using semiempirical quantum mechanical calculations based on the open source software package xTB.

The program generates possible tautomers of the input molecule by employing RDKit (v. 2020.03.1 or newer), and highlights reaction sites with relative proton binding affinities below 1 kcal/mol (green circles) and 3 kcal/mol (red circles) for all tautomers with relative energies below 15 kcal/mol.
Furthermore, RegioSQM20 assings a qualitative prediction of the reactivity to each tautomer (low, medium, or high) based on the absolute proton affinity of the most stable protonated isomer.

More information is available at the [RegioSQM20 paper](https://doi.org/XX.XXXX/XXXXXXX).

# Installation

We recommend using anaconda to installed the following packages in a Python 3 environment:

    - xTB v. 6.3.2 or newer (https://github.com/grimme-lab/xtb)
    - RDKit v. 2020.03.1 or newer (http://www.rdkit.org/docs/Install.html)
    - obabel v. 2.4.1 or newer (https://openbabel.org/docs/dev/Installation/install.html)

RegioSQM20 depends on xTB for quantum calculations, OpenBabel for some
formation convertion and RDKit in the Python environment for everything else.

Note the following:

1. The path to your Python environment must be specified in the regiosqm/conf_search.py file.
2. The number of CPUs/molecules and CPUs/conformer has to be specified in the regiosqm/regiosqm.py and regiosqm/conf_search.py files, respectively. (Setting CPUs/molecules = 2 and CPUs/conformer = 4 will allow RegioSQM20 to use a maximum of 2*4=8 CPUs that you must have available.)

# Usage

An example of using RegioSQM20:

    # Create predictions for the molecules in example.smiles:
    
    cd example

    python ../regiosqm/regiosqm.py -r compounds.smiles


    # Create predictions for the molecules in example.smiles and highlight the experimentally observed reaction sites (black circles):

    cd example

    python ../regiosqm/regiosqm.py -e -r compounds.smiles


The results are now parseable from the results file, or viewable by 2D structures with regioselective indicators (in .svg format).



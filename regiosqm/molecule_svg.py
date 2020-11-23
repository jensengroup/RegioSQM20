
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
# from rdkit.Chem import rdDepictor
# rdDepictor.SetPreferCoordGen(True)

from reorder_atoms import get_atoms_in_order
from collections import defaultdict

# Drawing Options
color_predicted = (0.2, 1, 0.0) # Green
color_loseicted = (1.0, 0.1, 0.3) # Red
color_measured = (0.0, 0.0, 0.0) # Black
arad = 0.4 #0.25
#molsPerRow = 4 #change this in generate_structure()
subImgSize = (300,300)


def draw2d(mol, name, subImgSize, highlight_predicted, highlight_loseicted, measure=None):
    
    global color_predicted
    global color_loseicted
    global color_measured
    global arad
    
    d2d = rdMolDraw2D.MolDraw2DSVG(subImgSize[0], subImgSize[1])
    d2d.SetFontSize(1) #atom label font size
    dos = d2d.drawOptions()
    dos.legendFontSize=23 #legend font size
    dos.atomHighlightsAreCircles = False
    dos.fillHighlights = True

    atomHighlighs = defaultdict(list)
    highlightRads = {}
    for idx in highlight_predicted:
        atomHighlighs[idx].append(color_predicted)
        highlightRads[idx] = arad

    # did threshold find some predictions?
    # find ones not in predicted list
    highlight_loseicted = list(set(highlight_loseicted)-set(highlight_predicted))
    if len(highlight_loseicted):
        for idx in highlight_loseicted:
            atomHighlighs[idx].append(color_loseicted)
            highlightRads[idx] = arad

    if measure:
        for idx in measure:
            atomHighlighs[idx].append(color_measured)
            highlightRads[idx] = arad
    
    d2d.DrawMoleculeWithHighlights(mol, name, dict(atomHighlighs), {}, highlightRads, {})
    d2d.FinishDrawing()

    return d2d.GetDrawingText()


def generate_structure(ref_smi, smiles, names, predicted, highlight_measure=None):

    global subImgSize
    molsPerRow = 4

    highlight_predicted, highlight_loseicted = predicted

    if names == None:
        names = ['' for i in range(len(smiles))]

    nRows = len(smiles) // molsPerRow
    if len(smiles) % molsPerRow:
        nRows += 1
    if nRows == 1:
        molsPerRow = len(smiles)
    fullSize = (molsPerRow * subImgSize[0], nRows * subImgSize[1])

    header = """<svg version='1.1' baseProfile='full'
                xmlns='http://www.w3.org/2000/svg'
                        xmlns:rdkit='http://www.rdkit.org/xml'
                        xmlns:xlink='http://www.w3.org/1999/xlink'
                    xml:space='preserve'
    width='{0}px' height='{1}px' viewBox='0 0 {0} {1}'>
    <!-- END OF HEADER -->""".format(fullSize[0],fullSize[1])

    spacer = '<g transform="translate({0},{1})">\n{2}</g>'

    ### Make sure the atoms are in order ###
    mols = [Chem.MolFromSmiles(smi) for smi in smiles]
    mols = get_atoms_in_order(Chem.MolFromSmiles(ref_smi), mols)

    cwidth = 0
    cheight = 0
    drawed_mols = []
    for i in range(len(smiles)):
        res = draw2d(mols[i], names[i], subImgSize, highlight_predicted[i], highlight_loseicted[i], highlight_measure)
        res = res.split("\n")
        end_of_header = res.index("<!-- END OF HEADER -->") + 1 
        res = "\n".join(res[end_of_header:-2])
        
        res = "".join(spacer.format(int(cwidth*subImgSize[0]), int(cheight*subImgSize[1]), res))
        drawed_mols.append(res)
        
        if int(i+1) % molsPerRow == 0 and i != 0:
            cheight += 1
            cwidth = 0
        elif molsPerRow == 1:
            cheight += 1
            cwidth = 0
        else:
            cwidth += 1

    svg = header + "\n" + "\n".join(drawed_mols) + "\n</svg>"

    return svg


if __name__ == "__main__":

    ref_smi = 'c1c(c2cc(sc2)C)n[nH]c1'
    smiles = ['c1c(-c2cc(C)sc2)[nH]nc1', 'c1c(-c2cc(C)sc2)n[nH]c1']
    names = ['taut1', 'taut2']
    highlight_predicted = [[7,0], [10]]
    highlight_loseicted = [[7], [10]]
    highlight_measure = [0]
    
    result_svg = generate_structure(ref_smi, smiles, names, [highlight_predicted, highlight_loseicted], highlight_measure=highlight_measure)
    
    fd = open('test.svg','w')
    fd.write(result_svg)
    fd.close()



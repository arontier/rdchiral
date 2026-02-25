import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem
from rdkit.Chem.rdchem import ChiralType, BondStereo, BondDir
import rdchiral.utils as utils
utils.PLEVEL = 0

from rdchiral.initialization import rdchiralReaction, rdchiralReactants

rxn_smarts = "[CH3:1][CH2:2][Br:3].[OH-:4]>>[CH3:1][CH2:2][OH:4].[Br-:3]"
reactant_smiles = "CCBr.[OH-]"

rxn = rdchiralReaction(rxn_smarts)
reactants = rdchiralReactants(reactant_smiles)

# Run reactions manually to inspect
outcomes_raw = rxn.rxn.RunReactants(reactants.reactants_achiral_list)
print(f"RunReactants got {len(outcomes_raw)} outcomes")

for outcome in outcomes_raw:
    print("\n--- Outcome ---")
    for m in outcome:
        print(f"  Mol: {Chem.MolToSmiles(m)}")
        for a in m.GetAtoms():
            props = a.GetPropsAsDict()
            print(f"    atom idx={a.GetIdx()} sym={a.GetSymbol()} mapnum={a.GetAtomMapNum()} props={props}")

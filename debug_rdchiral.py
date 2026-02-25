import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem
import rdchiral.utils as utils
utils.PLEVEL = 5  # enable debug output

from rdchiral.main import rdchiralRunText, rdchiralRun
from rdchiral.initialization import rdchiralReaction, rdchiralReactants

rxn_smarts = "[CH3:1][CH2:2][Br:3].[OH-:4]>>[CH3:1][CH2:2][OH:4].[Br-:3]"
reactant_smiles = "CCBr.[OH-]"

print("=== Testing with rdchiralRun (PLEVEL=5) ===")
rxn = rdchiralReaction(rxn_smarts)
reactants = rdchiralReactants(reactant_smiles)

print("\nreactants_achiral_list has", len(reactants.reactants_achiral_list), "mols")
for i, m in enumerate(reactants.reactants_achiral_list):
    print(f"  mol {i}: {Chem.MolToSmiles(m)}")

print("\nRunning rdchiralRun ...")
results = rdchiralRun(rxn, reactants)
print("Results:", results)

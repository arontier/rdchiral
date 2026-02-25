import rdkit.Chem as Chem
import rdchiral.utils as utils
utils.PLEVEL = 5

from rdchiral.initialization import rdchiralReaction, rdchiralReactants
from rdchiral.main import rdchiralRun

rxn_smarts = "[CH3:1][CH2:2][Br:3].[OH-:4]>>[CH3:1][CH2:2][OH:4].[Br-:3]"
reactant_smiles = "CCBr.[OH-]"

rxn = rdchiralReaction(rxn_smarts)
reactants = rdchiralReactants(reactant_smiles)

print("atoms_rt_map (template mapnum -> template atom):")
for mapnum, atom in rxn.atoms_rt_map.items():
    print(f"  mapnum={mapnum} sym={atom.GetSymbol()} atom_mapnum={atom.GetAtomMapNum()}")

print("\natoms_r (reactant mapnum -> reactant atom):")
for mapnum, atom in reactants.atoms_r.items():
    print(f"  mapnum={mapnum} sym={atom.GetSymbol()}")

print("\nRunning rdchiralRun ...")
results = rdchiralRun(rxn, reactants)
print("Results:", results)

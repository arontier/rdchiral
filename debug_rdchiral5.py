import rdkit.Chem as Chem
import rdchiral.utils as utils
utils.PLEVEL = 5

from rdchiral.initialization import rdchiralReaction, rdchiralReactants
from rdchiral.main import rdchiralRun

rxn_smarts = "[CH3:1][CH2:2][Br:3].[OH-:4]>>[CH3:1][CH2:2][OH:4].[Br-:3]"

# Third test case: CCBr.[OH-]
reactant_smiles = "CCBr.[OH-]"

rxn = rdchiralReaction(rxn_smarts)
reactants = rdchiralReactants(reactant_smiles)

results = rdchiralRun(rxn, reactants)
print("Results:", results)

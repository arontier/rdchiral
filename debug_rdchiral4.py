import rdkit.Chem as Chem
import rdchiral.utils as utils
utils.PLEVEL = 0

from rdchiral.initialization import rdchiralReaction, rdchiralReactants

rxn_smarts = "[CH3:1][CH2:2][Br:3].[OH-:4]>>[CH3:1][CH2:2][OH:4].[Br-:3]"
reactant_smiles = "CCBr.[OH-]"

rxn = rdchiralReaction(rxn_smarts)
reactants = rdchiralReactants(reactant_smiles)

outcomes_raw = rxn.rxn.RunReactants(reactants.reactants_achiral_list)

for outcome in outcomes_raw:
    # Manually simulate the rdchiralRun flow
    combined = outcome[0]
    for j in range(1, len(outcome)):
        combined = Chem.AllChem.CombineMols(combined, outcome[j])
    
    print("After CombineMols SMILES:", Chem.MolToSmiles(combined, True))
    print("Atoms and their props:")
    for a in combined.GetAtoms():
        props = {k: v for k, v in a.GetPropsAsDict().items() if 'Query' in k}
        print(f"  idx={a.GetIdx()} sym={a.GetSymbol()} q_props={props} total_hs={a.GetTotalNumHs()} explicit_val={a.GetExplicitValence() if True else '?'}")
    
    # Try SanitizeMol before clearing _Query props
    print("\nTest 1: SanitizeMol WITHOUT clearing _Query props:")
    import copy
    m1 = copy.copy(combined)
    try:
        Chem.SanitizeMol(m1)
        print("  Success!")
    except Exception as e:
        print(f"  FAIL: {e}")
    
    # Clear _Query props then SanitizeMol
    print("\nTest 2: SanitizeMol AFTER clearing _Query props:")
    m2 = Chem.RWMol(combined)
    for a in m2.GetAtoms():
        for prop in list(a.GetPropsAsDict().keys()):
            if prop.startswith('_Query'):
                a.ClearProp(prop)
                print(f"  Cleared {prop} from atom {a.GetSymbol()}")
    m2 = m2.GetMol()
    try:
        Chem.SanitizeMol(m2)
        m2.UpdatePropertyCache()
        print("  Success! SMILES:", Chem.MolToSmiles(m2))
    except Exception as e:
        print(f"  FAIL: {e}")

    # Check what sanitize fails on specifically
    print("\nTest 3: Partial SanitizeMol to identify the issue:")
    m3 = Chem.RWMol(combined)
    try:
        Chem.SanitizeMol(m3, Chem.SanitizeFlags.SANITIZE_FINDRADICALS |
                         Chem.SanitizeFlags.SANITIZE_SETAROMATICITY |
                         Chem.SanitizeFlags.SANITIZE_SETCONJUGATION)
        print("  Partial sanitize OK")
    except Exception as e:
        print(f"  Partial sanitize failed: {e}")

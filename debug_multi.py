import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem

rxn_smarts = "[CH3:1][CH2:2][Br:3].[OH-:4]>>[CH3:1][CH2:2][OH:4].[Br-:3]"
rxn = AllChem.ReactionFromSmarts(rxn_smarts)
rxn.Initialize()
print("Num reactant templates:", rxn.GetNumReactantTemplates())
for i, rct in enumerate(rxn.GetReactants()):
    print(f"  Template {i}: {Chem.MolToSmiles(rct)}")

# Test with CCBr + [OH-]
reactant_smiles = "CCBr.[OH-]"
combined = Chem.MolFromSmiles(reactant_smiles)
Chem.AssignStereochemistry(combined, flagPossibleStereoCenters=True)
combined.UpdatePropertyCache(strict=False)
# Assign map numbers
for i, a in enumerate(combined.GetAtoms()):
    a.SetAtomMapNum(i+1)

print("\nCombined mol atoms:")
for a in combined.GetAtoms():
    print(f"  idx={a.GetIdx()}, symbol={a.GetSymbol()}, mapnum={a.GetAtomMapNum()}, Hs={a.GetTotalNumHs()}")

# Strip chirality like rdchiral does
from rdkit.Chem.rdchem import ChiralType, BondStereo, BondDir
for a in combined.GetAtoms():
    a.SetChiralTag(ChiralType.CHI_UNSPECIFIED)
for b in combined.GetBonds():
    b.SetStereo(BondStereo.STEREONONE)
    b.SetBondDir(BondDir.NONE)

print("\nAtoms after stripping chirality:")
for a in combined.GetAtoms():
    print(f"  idx={a.GetIdx()}, symbol={a.GetSymbol()}, mapnum={a.GetAtomMapNum()}, Hs={a.GetTotalNumHs()}")

# Try GetMolFrags
frags = Chem.GetMolFrags(combined, asMols=True)
print(f"\nGetMolFrags returned {len(frags)} frags:")
for i, frag in enumerate(frags):
    print(f"  Frag {i}: {Chem.MolToSmiles(frag)} - atoms:")
    for a in frag.GetAtoms():
        print(f"    local_idx={a.GetIdx()}, symbol={a.GetSymbol()}, mapnum={a.GetAtomMapNum()}, Hs={a.GetTotalNumHs()}, charge={a.GetFormalCharge()}")

# Try RunReactants
print("\nTrying RunReactants with frags tuple:")
try:
    outcomes = rxn.RunReactants(tuple(frags))
    print(f"  Got {len(outcomes)} outcomes")
    for o in outcomes:
        print("  ->", [Chem.MolToSmiles(m) for m in o])
except Exception as e:
    print(f"  Error: {e}")

print("\nTrying RunReactants with combined single mol:")
try:
    outcomes2 = rxn.RunReactants((combined,))
    print(f"  Got {len(outcomes2)} outcomes")
    for o in outcomes2:
        print("  ->", [Chem.MolToSmiles(m) for m in o])
except Exception as e:
    print(f"  Error: {e}")

# Try with plain mol (no map numbers, no rdchiral processing)
print("\nTrying RunReactants with plain mols (no map numbers):")
mol1 = Chem.MolFromSmiles("CCBr")
mol2 = Chem.MolFromSmiles("[OH-]")
try:
    outcomes3 = rxn.RunReactants((mol1, mol2))
    print(f"  Got {len(outcomes3)} outcomes")
    for o in outcomes3:
        print("  ->", [Chem.MolToSmiles(m) for m in o])
except Exception as e:
    print(f"  Error: {e}")

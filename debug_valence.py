import rdkit.Chem as Chem

# Test if SanitizeMol fails for a mapped [OH-]
# This simulates what rdchiral creates
m = Chem.RWMol()
# Add C-C-O([mapnum=4])-[H] and [Br-][mapnum=3]
m = Chem.MolFromSmiles('CC[OH-].[Br-]')
print("Basic mol:", Chem.MolToSmiles(m))
print("SanitizeMol on basic mol:")
try:
    Chem.SanitizeMol(m)
    print("  OK")
except Exception as e:
    print(f"  FAIL: {e}")

# Now set atom map numbers to simulate rdchiral outcome
m2 = Chem.RWMol(m)
atoms = list(m2.GetAtoms())
print(f"Atoms: {[(a.GetSymbol(), a.GetIdx()) for a in atoms]}")
for a in m2.GetAtoms():
    sym = a.GetSymbol()
    if sym == 'C' and a.GetIdx() == 0:
        a.SetAtomMapNum(1)
    elif sym == 'C' and a.GetIdx() == 1:
        a.SetAtomMapNum(2)
    elif sym == 'O':
        a.SetAtomMapNum(4)
    elif sym == 'Br':
        a.SetAtomMapNum(3)
m2_mol = m2.GetMol()
smiles_mapped = Chem.MolToSmiles(m2_mol, True)
print(f"\nMapped SMILES: {smiles_mapped}")

# Re-parse that SMILES and sanitize
print("SanitizeMol on re-parsed mapped mol:")
m3 = Chem.MolFromSmiles(smiles_mapped)
if m3 is None:
    print("  FAIL: Cannot parse")
else:
    try:
        Chem.SanitizeMol(m3)
        m3.UpdatePropertyCache()
        print("  OK:", Chem.MolToSmiles(m3))
    except Exception as e:
        print(f"  FAIL: {e}")

# What rdchiral does: take the mol directly (not re-parsed), sanitize it
print("\nSanitizeMol on mapped mol DIRECTLY (as rdchiral does):")
try:
    Chem.SanitizeMol(m2_mol)
    m2_mol.UpdatePropertyCache()
    print("  OK:", Chem.MolToSmiles(m2_mol))
except Exception as e:
    print(f"  FAIL: {e}")

# Check atoms before SanitizeMol
print("\nAtom details in m2_mol before SanitizeMol:")
m4 = Chem.RWMol(m2_mol)
for a in m4.GetAtoms():
    print(f"  {a.GetSymbol()} mapnum={a.GetAtomMapNum()} charge={a.GetFormalCharge()} Hs={a.GetNumExplicitHs()} noImplicit={a.GetNoImplicit()} props={a.GetPropsAsDict()}")
    
# Try without explicit valence validation
print("\nSanitize without SANITIZE_PROPERTIES:")
try:
    flags = (Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES)
    Chem.SanitizeMol(m2_mol, flags)
    print("  OK (without SANITIZE_PROPERTIES)")
except Exception as e:
    print(f"  FAIL: {e}")

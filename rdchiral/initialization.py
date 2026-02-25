import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem
from rdkit.Chem.rdchem import ChiralType, BondType, BondDir, BondStereo

from rdchiral.chiral import template_atom_could_have_been_tetra
from rdchiral.utils import vprint, PLEVEL
from rdchiral.bonds import enumerate_possible_cistrans_defs, bond_dirs_by_mapnum, \
    get_atoms_across_double_bonds

BondDirOpposite = {BondDir.ENDUPRIGHT: BondDir.ENDDOWNRIGHT,
                   BondDir.ENDDOWNRIGHT: BondDir.ENDUPRIGHT}

class rdchiralReaction(object):
    '''Class to store everything that should be pre-computed for a reaction. This
    makes library application much faster, since we can pre-do a lot of work
    instead of doing it for every mol-template pair

    Attributes:
        reaction_smarts (str): reaction SMARTS string
        rxn (rdkit.Chem.rdChemReactions.ChemicalReaction): RDKit reaction object.
            Generated from `reaction_smarts` using `initialize_rxn_from_smarts`
        template_r: Reaction reactant template fragments
        template_p: Reaction product template fragments
        atoms_rt_map (dict): Dictionary mapping from atom map number to RDKit Atom for reactants
        atoms_pt_map (dict): Dictionary mapping from atom map number to RDKit Atom for products
        atoms_rt_idx_to_map (dict): Dictionary mapping from atom idx to RDKit Atom for reactants
        atoms_pt_idx_to_map (dict): Dictionary mapping from atom idx to RDKit Atom for products

    Args:
        reaction_smarts (str): Reaction SMARTS string
    '''
    def __init__(self, reaction_smarts):
        # Keep smarts, useful for reporting
        self.reaction_smarts = reaction_smarts

        # Initialize - assigns stereochemistry and fills in missing rct map numbers
        self.rxn = initialize_rxn_from_smarts(reaction_smarts)

        # Combine template fragments so we can play around with mapnums
        self.template_r, self.template_p = get_template_frags_from_rxn(self.rxn)

        # Define molAtomMapNumber->atom dictionary for template rct and prd
        self.atoms_rt_map = {a.GetAtomMapNum(): a \
            for a in self.template_r.GetAtoms() if a.GetAtomMapNum()}
        self.atoms_pt_map = {a.GetAtomMapNum(): a \
            for a in self.template_p.GetAtoms() if a.GetAtomMapNum()}

        # Back-up the mapping for the reaction
        self.atoms_rt_idx_to_map = {a.GetIdx(): a.GetAtomMapNum()
            for a in self.template_r.GetAtoms()}
        self.atoms_pt_idx_to_map = {a.GetIdx(): a.GetAtomMapNum()
            for a in self.template_p.GetAtoms()}

        # Check consistency (this should not be necessary...)
        if any(self.atoms_rt_map[i].GetAtomicNum() != self.atoms_pt_map[i].GetAtomicNum() \
                for i in self.atoms_rt_map if i in self.atoms_pt_map):
            raise ValueError('Atomic identity should not change in a reaction!')

        # Call template_atom_could_have_been_tetra to pre-assign value to atom
        [template_atom_could_have_been_tetra(a) for a in self.template_r.GetAtoms()]
        [template_atom_could_have_been_tetra(a) for a in self.template_p.GetAtoms()]

        # Pre-list chiral double bonds (for copying back into outcomes/matching)
        self.rt_bond_dirs_by_mapnum = bond_dirs_by_mapnum(self.template_r)
        self.pt_bond_dirs_by_mapnum = bond_dirs_by_mapnum(self.template_p)

        # Enumerate possible cis/trans...
        self.required_rt_bond_defs, self.required_bond_defs_coreatoms = \
            enumerate_possible_cistrans_defs(self.template_r)

    def reset(self):
        '''Reset atom map numbers for template fragment atoms'''
        for (idx, mapnum) in self.atoms_rt_idx_to_map.items():
            self.template_r.GetAtomWithIdx(idx).SetAtomMapNum(mapnum)
        for (idx, mapnum) in self.atoms_pt_idx_to_map.items():
            self.template_p.GetAtomWithIdx(idx).SetAtomMapNum(mapnum)

class rdchiralReactants(object):
    '''Class to store everything that should be pre-computed for a reactant mol
    so that library application is faster.

    Supports both single-molecule and multi-molecule (dot-separated) SMILES.
    When the SMILES contains multiple fragments (e.g., "CCBr.[OH-]"), the
    reactants are stored as a combined molecule for atom-level operations,
    and also as a tuple of individual molecules for RunReactants.

    Attributes:
        reactant_smiles (str): Reactant SMILES string
        reactants (rdkit.Chem.rdchem.Mol): RDKit Molecule (combined)
        atoms_r (dict): Dictionary mapping from atom map number to atom in `reactants` Molecule
        idx_to_mapnum (callable): callable function that takes idx and returns atom map number
        reactants_achiral (rdkit.Chem.rdchem.Mol): achiral version of `reactants` (combined)
        reactants_achiral_list (tuple): tuple of individual achiral mols for RunReactants
        bonds_by_mapnum (list): List of reactant bonds
            (int, int, rdkit.Chem.rdchem.Bond)
        bond_dirs_by_mapnum (dict): Dictionary mapping from atom map number tuples to BondDir
        atoms_across_double_bonds (list): List of cis/trans specifications from `get_atoms_across_double_bonds`

    Args:
        reactant_smiles (str): Reactant SMILES string
        custom_reactant_mapping (bool): Whether custom atom mapping is provided
    '''
    def __init__(self, reactant_smiles, custom_reactant_mapping=False):
        # Keep original smiles, useful for reporting
        self.reactant_smiles = reactant_smiles
        self.custom_mapping = custom_reactant_mapping

        # Initialize into RDKit mol (combined)
        self.reactants = initialize_reactants_from_smiles(reactant_smiles, custom_reactant_mapping)

        # Set mapnum->atom dictionary
        # all reactant atoms must be mapped after initialization, so this is safe
        self.atoms_r = {a.GetAtomMapNum(): a for a in self.reactants.GetAtoms()}
        self.idx_to_mapnum = lambda idx: self.reactants.GetAtomWithIdx(idx).GetAtomMapNum()

        # Create combined achiral copy (used for atom-level operations)
        self.reactants_achiral = initialize_reactants_from_smiles(reactant_smiles, custom_reactant_mapping)
        [a.SetChiralTag(ChiralType.CHI_UNSPECIFIED) for a in self.reactants_achiral.GetAtoms()]
        [(b.SetStereo(BondStereo.STEREONONE), b.SetBondDir(BondDir.NONE))
            for b in self.reactants_achiral.GetBonds()]

        # Build per-fragment achiral mol list for multi-reactant RunReactants calls.
        # Use Chem.GetMolFrags to correctly split the combined achiral mol into
        # individual fragment Mols (handles charged species, brackets, etc. safely).
        #
        # Important: RDKit sets react_atom_idx = local atom index within each
        # individual reactant mol (NOT global index across all reactants).
        # So we need a mapping from (frag_position, local_idx) -> combined_mol_idx
        # to correctly look up atom map numbers via idx_to_mapnum.
        frag_atom_indices = Chem.GetMolFrags(self.reactants_achiral)  # tuple of tuples of combined-mol atom idx
        frags = Chem.GetMolFrags(self.reactants_achiral, asMols=True)
        if len(frags) > 1:
            self.reactants_achiral_list = tuple(frags)
            # For each fragment i, build a list where local_idx -> combined_mol_idx
            # frag_atom_indices[i][local_idx] = combined_mol_idx
            self._frag_local_to_global = frag_atom_indices  # tuple of tuples
            # Build a merged lookup: react_atom_idx (local per-frag) -> combined idx
            # When RunReactants((f0, f1, ...), react_atom_idx IS the local_idx within fn.
            # We need (frag_n, local_local_idx) -> combined_idx, but we don't know
            # which frag an outcome atom came from in react_atom_idx alone.
            # Instead build a single list: for frag i, local_idx j -> combined idx frag_atom_indices[i][j]
            # Store as a function that handles the lookup when we know the frag index.
            # However, since react_atom_idx is only local idx and we don't always know
            # which frag it belongs to, we build a combined react_atom_idx space by
            # using old_mapno to find the correct atom.
            self.reactant_frag_atom_indices = frag_atom_indices
        else:
            self.reactants_achiral_list = (self.reactants_achiral,)
            self.reactant_frag_atom_indices = None


        # Pre-list reactant bonds (for stitching broken products)
        self.bonds_by_mapnum = [
            (b.GetBeginAtom().GetAtomMapNum(), b.GetEndAtom().GetAtomMapNum(), b)
            for b in self.reactants.GetBonds()
        ]

        # Pre-list chiral double bonds (for copying back into outcomes/matching)
        self.bond_dirs_by_mapnum = {}
        for (i, j, b) in self.bonds_by_mapnum:
            if b.GetBondDir() != BondDir.NONE:
                self.bond_dirs_by_mapnum[(i, j)] = b.GetBondDir()
                self.bond_dirs_by_mapnum[(j, i)] = BondDirOpposite[b.GetBondDir()]

        # Get atoms across double bonds defined by mapnum
        self.atoms_across_double_bonds = get_atoms_across_double_bonds(self.reactants)


def initialize_rxn_from_smarts(reaction_smarts):
    '''Initialize RDKit reaction object from SMARTS string

    Args:
        reaction_smarts (str): Reaction SMARTS string

    Returns:
        rdkit.Chem.rdChemReactions.ChemicalReaction: RDKit reaction object
    '''
    # Initialize reaction
    rxn = AllChem.ReactionFromSmarts(reaction_smarts)
    rxn.Initialize()
    if rxn.Validate()[1] != 0:
        raise ValueError('validation failed')
    if PLEVEL >= 2: print('Validated rxn without errors')


    # Figure out if there are unnecessary atom map numbers (that are not balanced)
    # e.g., leaving groups for retrosynthetic templates. This is because additional
    # atom map numbers in the input SMARTS template may conflict with the atom map
    # numbers of the molecules themselves
    prd_maps = [a.GetAtomMapNum() for prd in rxn.GetProducts() for a in prd.GetAtoms() if a.GetAtomMapNum()]

    unmapped = 700
    for rct in rxn.GetReactants():
        rct.UpdatePropertyCache(strict=False)
        Chem.AssignStereochemistry(rct)
        # Fill in atom map numbers
        for a in rct.GetAtoms():
            if not a.GetAtomMapNum() or a.GetAtomMapNum() not in prd_maps:
                a.SetAtomMapNum(unmapped)
                unmapped += 1
    if PLEVEL >= 2: print('Added {} map nums to unmapped reactants'.format(unmapped-700))
    if unmapped > 800:
        raise ValueError('Why do you have so many unmapped atoms in the template reactants?')

    return rxn

def initialize_reactants_from_smiles(reactant_smiles, custom_reactant_mapping):
    '''Initialize RDKit molecule from SMILES string

    Args:
        reactant_smiles (str): Reactant SMILES string

    Returns:
        rdkit.Chem.rdchem.Mol: RDKit molecule
    '''
    # Initialize reactants
    reactants = Chem.MolFromSmiles(reactant_smiles)
    Chem.AssignStereochemistry(reactants, flagPossibleStereoCenters=True)
    reactants.UpdatePropertyCache(strict=False)
    # To have the product atoms match reactant atoms, we
    # need to populate the map number field, since this field
    # gets copied over during the reaction via reactant_atom_idx.
    if not custom_reactant_mapping:
        [a.SetAtomMapNum(i+1) for (i, a) in enumerate(reactants.GetAtoms())]
    if PLEVEL >= 2: print('Initialized reactants, assigned map numbers, stereochem, flagpossiblestereocenters')
    return reactants

def get_template_frags_from_rxn(rxn):
    '''Get template fragments from RDKit reaction object

    Args:
        rxn (rdkit.Chem.rdChemReactions.ChemicalReaction): RDKit reaction object

    Returns:
        (rdkit.Chem.rdchem.Mol, rdkit.Chem.rdchem.Mol): tuple of fragment molecules
    '''
    # Copy reaction template so we can play around with map numbers
    for i, rct in enumerate(rxn.GetReactants()):
        if i == 0:
            template_r = rct
        else:
            template_r = AllChem.CombineMols(template_r, rct)
    for i, prd in enumerate(rxn.GetProducts()):
        if i == 0:
            template_p = prd
        else:
            template_p = AllChem.CombineMols(template_p, prd)
    return template_r, template_p
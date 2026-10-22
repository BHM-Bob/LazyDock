
# which toolkit prepares the receptor/ligand pdbqt handed to the docking engine.
# Both are usable and the switch is NOT score-neutral: meeko and ADT classify a few
# nitrogen acceptors and polar hydrogens differently, so scores from the two backends
# are not comparable (see `prepare_receptor_with_meeko`).
PREPARE_BACKENDS = ('meeko', 'adfr')
PREPARE_BACKEND_DEFAULT = 'meeko'


# first meeko release where --flexres takes one comma-separated sequence
MEEKO_FLAT_FLEXRES_VERSION = (0, 7, 0)


# meeko breaks rings only when min_ring_size <= size <= max_ring_size; larger rings
# are SILENTLY ignored (no error, no warning) and stay rigid.  The default 33 misses
# most cyclic peptides: a measured scan gives the largest main ring at 56 atoms, and
# a per-molecule probe shows the threshold reacts sharply (a 48-mer ring is still
# rigid at max_ring_size=45 and only breaks at 48).  Hence 60, with headroom.
MACROCYCLE_MAX_RING_SIZE = 60

# meeko's ring breaking is ON BY DEFAULT for the whole 7..33 atom window, which is
# easy to miss: it does NOT only affect big macrocycles but also plain medium rings
# (a measured probe opens 7-, 10-, 12- and 20-membered rings with stock settings,
# while benzene stays closed only because the lower bound is 7).  That silently
# changes a ligand that ADT's `prepare_ligand` would deliver as a rigid ring, and
# the glue pseudo-atoms it adds (G0/CG0) have no grid in an AutoDock map.  Passing
# `rigid_macrocycles=True` makes every ring rigid, which is the ADT-equivalent
# behaviour; meeko documents it as "every ring is rigid without macrocycle option".
RIGID_MACROCYCLES = True

# disulfide bonds are covalent cross-links, so they must not become freely rotating
# torsions once the ring opens.  meeko matches the SMARTS and freezes the bond
# between its atoms: the indices are POSITIONS INSIDE THE SMARTS (0-based here,
# 1-based on the CLI), NEVER atom serial numbers.
# A SMARTS with no match is harmless, so this can be passed unconditionally.
SS_RIGIDIFY_SMARTS = ['[SX2]-[SX2]']
SS_RIGIDIFY_INDICES = [(0, 1)]

# macrocycle glue pseudo-atoms are named "G" with zero radius/charge (atom types
# G0..G3).  "G" is not an atom name of any standard residue, so dropping it needs no
# opt-in; it MUST happen before a pose reaches a pdb consumer: `ADModel(parse2std=True)`
# re-derives the element from the atom NAME and writes it into a 2-char column, so a
# glue atom would reach the output pdb as a fake element "G" and break
# PDBFixer/OpenMM and any re-preparation of the pose.
GLUE_ATOM_NAMES = ('G',)

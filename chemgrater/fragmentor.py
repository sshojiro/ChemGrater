# https://gist.github.com/sshojiro/251d60898844c99bba51d85162933cc3#file-extract-linear-fragments-ipynb
from rdkit.Chem import rdmolops
from rdkit import Chem
import networkx as nx


def bond_indices_to_nodes(m, bond_indices):
    """
    bond_indices_to_nodes
    ---
    Returns a sub-graph, indicating a Mol object and bond indices on RDKit.

    input: m, Mol object of RDKit
           bond_indices, indices of bonds, which is enumerated by RDKit
    return: g, graph object of NetworkX
    """
    G = nx.Graph()
    lst=[]
    for bond_ix in bond_indices:
        bond= m.GetBondWithIdx(bond_ix)
        lst+=[(bond.GetBeginAtomIdx(),bond.GetEndAtomIdx())]
    G.add_edges_from(lst)
    return G

def list_end_points_indices(g):
    """
    list_end_points_indices
    ---
    Detect end points in a sub-graph of a molecule
    expecting that the sub-graph has two end points.

    input: g, graph object of NetworkX
    return: indices. Returns blank array when the sub-graph is a ring.
    """
    return [k for k,v in g.adjacency()if len(v)==1]

def sort_atoms_in_indices(ix_start, g):
    indices_sorted = []
    atom_indices = set(g.nodes)
    atom_indices -= {ix_start}
    next_ix = ix_start
    indices_sorted += [next_ix]
    while any(k in atom_indices for k in list(g.adj[next_ix].keys())):
        if list(g.adj[next_ix].keys())[0] in atom_indices:
            next_ix = list(g.adj[next_ix].keys())[0]
        elif list(g.adj[next_ix].keys())[1] in atom_indices:
            next_ix = list(g.adj[next_ix].keys())[1]
        indices_sorted += [next_ix]
        atom_indices -= {next_ix}
    return indices_sorted

def generate_linear_fragments(mol, n_len=3):
    """
    generate_linear_fragments
    ---
    Extracting linear fragments from a Mol instance of RDKit, indicating the length.

    input: mol, Mol object of RDKit
           n_len, number of bonds
    return: smiles, pseudo-SMILES that expresses sub-graphs
            lst_indices, indexes of corresponding atoms
    """
    lst_paths = rdmolops.FindAllPathsOfLengthN(mol,n_len,useBonds=True,useHs=True)
    smiles = []
    lst_indices = []
    for p in lst_paths:
        lst=[]
        for bond_ix in p:
            bond = mol.GetBondWithIdx(bond_ix)
            lst+=[(bond.GetBeginAtomIdx(),bond.GetEndAtomIdx())]
        gg = bond_indices_to_nodes(mol, list(p))
        ind = list_end_points_indices(gg)
        if len(ind) == 0:
            # Note that the sub-graphs are ignored when the sub-graphs are rings.
            # For instance, when cyclo-ring happens to be fed as sub-graph,
            # the function `list_end_points_indices` returns a blank array [].
            continue
        atom_indices = list(sort_atoms_in_indices(ind[0], gg))
        lst_indices += [atom_indices]
        smi=''
        for ix in range(n_len):
            smi += mol.GetAtomWithIdx(atom_indices[ix]).GetSymbol()
            bond = mol.GetBondBetweenAtoms(atom_indices[ix],atom_indices[ix+1])
            if bond.GetBondType() is Chem.BondType.SINGLE:
                smi += '-'# explicitly shown, compared to aromatic bond.
            elif bond.GetBondType() is Chem.BondType.DOUBLE:
                smi += '='
            elif bond.GetBondType() is Chem.BondType.TRIPLE:
                smi += '#'
            elif bond.GetBondType() is Chem.BondType.QUADRUPLE:
                smi += '$'
        else:
            smi += mol.GetAtomWithIdx(atom_indices[ix+1]).GetSymbol()
        smiles += [smi]
    return smiles, lst_indices

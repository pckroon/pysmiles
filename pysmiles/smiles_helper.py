# -*- coding: utf-8 -*-
# Copyright 2018 Peter C Kroon

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import re
import operator

import networkx as nx

ISOTOPE_PATTERN = r'(?P<isotope>[\d]+)?'
ELEMENT_PATTERN = r'(?P<element>b|c|n|o|s|p|\*|[A-Z][a-z]{0,2})'
STEREO_PATTERN = r'(?P<stereo>@|@@|@TH[1-2]|@AL[1-2]|@SP[1-3]|@OH[\d]{1,2}|'\
                  r'@TB[\d]{1,2})?'
HCOUNT_PATTERN = r'(?P<hcount>H[\d]?)?'
CHARGE_PATTERN = r'(?P<charge>(-|\+)(\++|-+|[\d]{1,2})?)?'
CLASS_PATTERN = r'(?::(?P<class>[\d]+))?'
ATOM_PATTERN = re.compile('^' + ISOTOPE_PATTERN + ELEMENT_PATTERN +
                          STEREO_PATTERN + HCOUNT_PATTERN + CHARGE_PATTERN +
                          CLASS_PATTERN + '$')

VALENCES = {"B": (3,), "C": (4,), "N": (3, 5), "O": (2,), "P": (3, 5),
            "S": (2, 4, 6), "F": (1,), "Cl": (1,), "Br": (1,), "I": (1,)}

AROMATIC_ATOMS = "B C N O P S Se As".split()


def parse_atom(atom):
    """
    Parses a SMILES atom token, and returns a dict with the information.

    Note
    ----
    Can not deal with stereochemical information yet. This gets discarded.

    Parameters
    ----------
    atom : str
        The atom string to interpret. Looks something like one of the
        following: "C", "c", "[13CH3-1:2]"

    Returns
    -------
    dict
        A dictionary containing at least 'element' and 'charge'. If present,
        will also contain 'hcount', 'isotope', and 'class'.
    """
    defaults = {'charge': 0, 'hcount': 0}
    if not atom.startswith('[') and not atom.endswith(']'):
        if atom != '*':
            return {'element': atom, 'charge': 0}
        else:
            return defaults.copy()
    atom = atom.strip('[]')
    match = ATOM_PATTERN.match(atom)
    if match is None:
        raise ValueError('The atom {} is malformatted'.format(atom))
    out = defaults.copy()
    out.update({k: v for k, v in match.groupdict().items() if v is not None})

    parse_helpers = {'isotope': int,
                     'element': lambda x: x,
                     'stereo': lambda x: x,
                     'hcount': parse_hcount,
                     'charge': parse_charge,
                     'class': int}

    for attr, val_str in out.items():
        out[attr] = parse_helpers[attr](val_str)

    if out['element'] == '*':
        del out['element']

    if out.get('element') == 'H' and out.get('hcount', 0):
        raise ValueError("A hydrogen atom can't have hydrogens")

    if 'stereo' in out:
        print("I don't quite know how to handle stereo yet...")

    return out


def parse_hcount(hcount_str):
    """
    Parses a SMILES hydrogen count specifications.

    Parameters
    ----------
    hcount_str : str
        The hydrogen count specification to parse.

    Returns
    -------
    int
        The number of hydrogens specified.
    """
    if not hcount_str:
        return 0
    if hcount_str == 'H':
        return 1
    return int(hcount_str[1:])


def parse_charge(charge_str):
    """
    Parses a SMILES charge specification.

    Parameters
    ----------
    charge_str : str
        The charge specification to parse.

    Returns
    -------
    int
        The charge.
    """
    if not charge_str:
        return 0
    signs = {'-': -1, '+': 1}
    sign = signs[charge_str[0]]
    if len(charge_str) > 1 and charge_str[1].isdigit():
        charge = sign * int(charge_str[1:])
    else:
        charge = sign * charge_str.count(charge_str[0])
    return charge


def add_hydrogens(mol):
    """
    Adds explicit hydrogen nodes to `mol`, the amount is determined by the node
    attribute 'hcount'. Will remove the 'hcount' attribute.

    Parameters
    ----------
    mol : nx.Graph
        The molecule to which explicit hydrogens should be added. Is modified
        in-place.

    Returns
    -------
    None
        `mol` is modified in-place.
    """
    h_atom = parse_atom('[H]')
    if 'hcount' in h_atom:
        del h_atom['hcount']
    for n_idx in list(mol.nodes):
        hcount = mol.nodes[n_idx].get('hcount', 0)
        idxs = range(max(mol) + 1, max(mol) + hcount + 1)
        # Get the defaults from parse_atom.
        mol.add_nodes_from(idxs, **h_atom.copy())
        mol.add_edges_from([(n_idx, jdx) for jdx in idxs], order=1)
        if 'hcount' in mol.nodes[n_idx]:
            del mol.nodes[n_idx]['hcount']


def remove_hydrogens(mol):
    """
    Removes all explicit, simple hydrogens from `mol`. Simple means it is
    identical to the SMILES string "[H]", and has exactly one bond. Increments
    'hcount' where appropriate.

    Parameters
    ----------
    mol : nx.Graph
        The molecule whose explicit hydrogens should be removed. Is modified
        in-place.

    Returns
    -------
    None
        `mol` is modified in-place.
    """
    to_remove = set()
    defaults = parse_atom('[H]')
    for n_idx in mol.nodes:
        node = mol.nodes[n_idx]
        neighbors = mol[n_idx]
        if node == defaults and len(neighbors) == 1:
            to_remove.add(n_idx)
            neighbor = list(neighbors.keys())[0]
            mol.nodes[neighbor]['hcount'] = mol.nodes[neighbor].get('hcount', 0) + 1
    mol.remove_nodes_from(to_remove)


def fill_valence(mol, respect_hcount=True):
    """
    Sets the attribute 'hcount' on all nodes in `mol` that don't have it yet.
    The value to which it is set is based on the node's 'element', and the
    number of bonds it has. Default valences are as specified by the global
    variable VALENCES.

    Parameters
    ----------
    mol : nx.Graph
        The molecule whose nodes should get a 'hcount'. Is modified in-place.
    respect_hcount : bool
        If True, don't change the hcount on nodes that already have it set.

    Returns
    -------
    None
        `mol` is modified in-place.
    """
    for n_idx in mol:
        node = mol.nodes[n_idx]
        if 'hcount' in node and respect_hcount:
            continue
        missing = under_valence(mol, n_idx)
        node['hcount'] = missing


def under_valence(mol, node_idx, use_order=True):
    """
    Returns how much the specified node is under valence. If use_order is
    False, treat all bonds as if they are order 1. Returns `hcount` if it is
    set for the node.

    Parameters
    ----------
    mol : nx.Graph
        The molecule.
    node_idx : hashable
        The node to look at. Should be in mol.
    use_order : bool
        If False, treat all bonds as single.

    Returns
    -------
    int
        The number of missing bonds.
    """
    node = mol.nodes[node_idx]
    element = node.get('element')
    if element not in VALENCES:
        return 0
    val = VALENCES.get(element)
    if use_order:
        bond_orders = map(operator.itemgetter(2),
                          mol.edges(nbunch=node_idx, data='order', default=1))
        bonds = sum(bond_orders)
    else:
        bonds = len(mol[node_idx])
    bonds += node.get('hcount', 0)
    try:
        val = min(filter(lambda a: a >= bonds, val))
    except ValueError:  # More bonds than possible
        val = max(val)
    return max(int(val - bonds), 0)


def mark_aromatic_atoms(mol):
    """
    Sets the elements of all aromatic atoms in mol (according to Hückel's rule)
    to lowercase. Requires that the 'hcount' on atoms is correct.

    Parameters
    ----------
    mol : nx.Graph
        The molecule.

    Returns
    -------
    None
        `mol` is modified in-place.
    """
    # Only cycles can be aromiatic
    aromatic = set()
    for cycle in nx.cycle_basis(mol):
        electrons = 0
        maybe_aromatic = True
        for idx, jdx, order in mol.edges(nbunch=cycle, data='order', default=1):
            if idx in cycle and jdx in cycle:
                electrons += (order - 1) * 2
        for node_idx in cycle:
            node = mol.nodes[node_idx]
            element = node['element'].capitalize()
            hcount = node.get('hcount', 0)
            degree = mol.degree(node_idx) + hcount
            if element not in AROMATIC_ATOMS or degree not in (2, 3):
                maybe_aromatic = False
                break
            # missing cases:
            #   extracyclic sp2 heteroatom (e.g. =O)
            #   some charged cases
            if element in 'N P As'.split() and hcount == 1:
                electrons += 1
            elif element in 'O S Se'.split():
                electrons += 1
            if node.get('charge', 0) == +1 and not (element == 'C' and hcount == 0):
                electrons -= 1
        if maybe_aromatic and int(electrons) % 2 == 0:
            # definitely (anti) aromatic
            aromatic.update(cycle)
    for node_idx in mol:
        node = mol.nodes[node_idx]
        if node_idx not in aromatic:
            node['element'] = node['element'].capitalize()
        else:
            node['element'] = node['element'].lower()


def mark_aromatic_edges(mol):
    """
    Set all bonds between aromatic atoms (lowercase elements) to 1.5. Run
    `mark_aromatic_atoms` first.

    Parameters
    ----------
    mol : nx.Graph
        The molecule.

    Returns
    -------
    None
        `mol` is modified in-place.
    """
    for cycle in nx.cycle_basis(mol):
        for idx, jdx in mol.edges:
            if idx not in cycle or jdx not in cycle:
                continue
            if mol.nodes[idx]['element'].islower() and\
                  mol.nodes[jdx]['element'].islower():
                mol.edges[idx, jdx]['order'] = 1.5


def fix_aromatic(mol):
    """
    Sets hcount for all atoms, the element of all aromatic atoms to lowercase,
    and the order of all aromatic bonds to 1.5.

    Parameters
    ----------
    mol : nx.Graph
        The molecule.

    Returns
    -------
    None
        `mol` is modified in-place.
    """
    fill_valence(mol)
    mark_aromatic_atoms(mol)
    mark_aromatic_edges(mol)


def increment_bond_orders(molecule):
    """
    Increments bond orders up to what the atom's valence allows.

    Parameters
    ----------
    molecule : nx.Graph
        The molecule to process.

    Returns
    -------
    None
        molecule is modified in-place.
    """
    max_bond_order = 3
    # Gather the number of open spots for all atoms beforehand, since some
    # might have multiple oxidation states (e.g. S). We don't want to change
    # oxidation state halfway through for some funny reason. It shouln't be
    # nescessary, but it can't hurt.
    missing_bonds = {}
    for idx in molecule:
        missing_bonds[idx] = under_valence(molecule, idx)

    for idx, jdx in molecule.edges:
        missing_idx = missing_bonds[idx]
        missing_jdx = missing_bonds[jdx]
        edge_missing = min(missing_idx, missing_jdx)
        new_order = edge_missing + molecule.edges[idx, jdx].get("order", 1)
        new_order = min(new_order, max_bond_order)
        molecule.edges[idx, jdx]['order'] = new_order
        missing_bonds[idx] -= edge_missing
        missing_bonds[jdx] -= edge_missing

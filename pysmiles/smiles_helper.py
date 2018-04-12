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
CHARGE_PATTERN = r'(?P<charge>--|\+\+|-[\d]{0,2}|\+[\d]{0,2})?'
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
    if not atom.startswith('[') and not atom.endswith(']'):
        if atom != '*':
            return {'element': atom, 'charge': 0}
        else:
            return {'charge': 0, 'hcount': 0}
    atom = atom.strip('[]')
    match = ATOM_PATTERN.match(atom)
    if match is None:
        raise ValueError('The atom {} is malformatted'.format(atom))
    out = match.groupdict()
    out = {k: v for k, v in out.items() if v is not None}
    
    if 'isotope' in out:
        out['isotope'] = int(out['isotope'])
    if out['element'] == '*':
        del out['element']
    if 'hcount' in out:
        # Starts with H
        hcount = out['hcount']
        if hcount == 'H':
            out['hcount'] = 1
        else:
            out['hcount'] = int(hcount[1:])
        if out.get('element') == 'H':
            raise ValueError("a hydrogen atom can't have hydrogens")
    else:
        out['hcount'] = 0
    if 'stereo' in out:
        print("I don't quite know how to handle stereo yet...")
    if 'charge' in out:
        charge = out['charge']
        if charge == '--':
            charge = -2
        elif charge == '++':
            charge = +2
        elif charge == '-':
            charge = -1
        elif charge == '+':
            charge = +1
        else:
            charge = int(charge)
        out['charge'] = charge
    else:
        out['charge'] = 0
    if 'class' in out:
        out['class'] = int(out['class'])
    return out


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
    H_atom = parse_atom('[H]')
    if 'hcount' in H_atom:
        del H_atom['hcount']
    for n_idx in list(mol.nodes):
        hcount = mol.nodes[n_idx].get('hcount', 0)
        idxs = range(max(mol) + 1, max(mol) + hcount + 1)
        # Get the defaults from parse_atom.
        mol.add_nodes_from(idxs, **H_atom.copy())
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


def fill_valence(mol):
    """
    Sets the attribute 'hcount' on all nodes in `mol` that don't have it yet.
    The value to which it is set is based on the node's 'element', and the
    number of bonds it has. Default valences are as specified by the global
    variable VALENCES.

    Parameters
    ----------
    mol : nx.Graph
        The molecule whose nodes should get a 'hcount'. Is modified in-place.

    Returns
    -------
    None
        `mol` is modified in-place.
    """
    for n_idx in mol:
        node = mol.nodes[n_idx]
        element = node.get('element')
        if element not in VALENCES or 'hcount' in node:
            continue
        val = VALENCES.get(element)
        bond_orders = map(operator.itemgetter(2),
                          mol.edges(nbunch=n_idx, data='order', default=0))
        bonds = sum(bond_orders)
        try:
            val = min(filter(lambda a: a >= bonds, val))
        except ValueError:  # More bonds than possible
            val = max(val)

        if val - bonds > 0:
            node['hcount'] = int(val - bonds)
        else:
            node['hcount'] = 0


def under_valence(mol, node_idx, use_order=True):
    node = mol.nodes[node_idx]
    element = node.get('element')
    if element not in VALENCES:
        return 0
    val = VALENCES.get(element)
    bond_orders = map(operator.itemgetter(2),
                      mol.edges(nbunch=node_idx, data='order', default=0))
    bonds = sum(bond_orders)
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
    """
    # Only cycles can be aromiatic
    for cycle in nx.cycle_basis(mol):
        electrons = 0
        maybe_aromatic = True
        for idx, jdx, order in mol.edges(nbunch=cycle, data='order', default=1):
            if idx in cycle and jdx in cycle:
                electrons += (order - 1) * 2
        for node_idx in cycle:
            element = mol.nodes[node_idx]['element'].capitalize()
            hcount = mol.nodes[node_idx].get('hcount', 0)
            degree = mol.degree(node_idx) + hcount
            if element not in AROMATIC_ATOMS or degree not in (2, 3):
                maybe_aromatic = False
                break
            missing = under_valence(mol, node_idx)
            # missing cases:
            #   charged atoms
            #   extracyclic sp2 heteroatom (e.g. =O)
            if element in 'N P As'.split() and hcount == 1:
                electrons += 1
            elif element in 'O S Se'.split():
                electrons += 1
        if maybe_aromatic and int(electrons) % 2 == 0:
            print('Aromatic')
        else:
            print('Not aromatic')

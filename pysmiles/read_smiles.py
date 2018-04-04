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

import networkx as nx
import enum
import re
import operator

ISOTOPE_PATTERN = r'(?P<isotope>[\d]+)?'
ELEMENT_PATTERN = r'(?P<element>b|c|n|o|s|p|\*|[A-Z][a-z]{0,2})'
STEREO_PATTERN = r'(?P<stereo>@|@@|@TH[1-2]|@AL[1-2]|@SP[1-3]|@OH[\d]{1,2}|'\
                  '@TB[\d]{1,2})?'
HCOUNT_PATTERN = r'(?P<hcount>H[\d]?)?'
CHARGE_PATTERN = r'(?P<charge>--|\+\+|-[\d]{0,2}|\+[\d]{0,2})?'
CLASS_PATTERN = r'(?::(?P<class>[\d]+))?'
ATOM_PATTERN = re.compile(ISOTOPE_PATTERN + ELEMENT_PATTERN + STEREO_PATTERN +
                          HCOUNT_PATTERN + CHARGE_PATTERN + CLASS_PATTERN)

VALENCES = {"B": (3,), "C": (4,), "N": (3, 5), "O": (2,), "P": (3, 5),
            "S": (2, 4, 6), "F": (1,), "Cl": (1,), "Br": (1,), "I": (1,)}


class TokenType(enum.Enum):
    """Helper class defining the types of tokens possible in a SMILES string"""
    ATOM = enum.auto()
    BOND_TYPE = enum.auto()
    BRANCH_START = enum.auto()
    BRANCH_END = enum.auto()
    RING_NUM = enum.auto()


def tokenize(smiles):
    """
    Iterates over a SMILES string, yielding tokens.

    Parameters
    ----------
    smiles : iterable
        The SMILES string to iterate over

    Yields
    ------
    tuple(TokenType, str)
        A tuple describing the type of token and the associated data
    """
    organic_subset = 'B C N O P S F Cl Br I * b c n o s p'.split()
    smiles = iter(smiles)
    token = ''
    peek = None
    while True:
        try:
            char = peek if peek else next(smiles)
            peek = None
        except StopIteration:
            break
        if char == '[':
            token = char
            for char in smiles:
                token += char
                if char == ']':
                    break
            yield TokenType.ATOM, token
        elif char in organic_subset:
            peek = next(smiles, '')
            if char + peek in organic_subset:
                yield TokenType.ATOM, char + peek
                peek = None
            else:
                yield TokenType.ATOM, char
        elif char in '-=#$:.':
            yield TokenType.BOND_TYPE, char
        elif char == '(':
            yield TokenType.BRANCH_START, '('
        elif char == ')':
            yield TokenType.BRANCH_END, ')'
        elif char == '%':
            yield TokenType.RING_NUM, int(next(smiles) + next(smiles))
        elif char.isdigit():
            yield TokenType.RING_NUM, int(char)


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
            return {}
    atom = atom.strip('[]')
    match = ATOM_PATTERN.fullmatch(atom)
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
    for n_idx in list(mol.nodes):
        hcount = mol.nodes[n_idx].get('hcount', 0)
        idxs = range(max(mol) + 1, max(mol) + hcount + 1)
        # Get the defaults from parse_atom.
        mol.add_nodes_from(idxs, **parse_atom('[H]'))
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
        val = min(filter(lambda a: a >= bonds, val))

        if val - bonds > 0:
            node['hcount'] = int(val - bonds)
        else:
            node['hcount'] = 0


def aromatize_bonds(mol):
    """
    Sets bond orders between atoms which are specified to be aromatic to 1.5.
    Atoms are aromatic if their element is lowercase.

    Paratmeters
    -----------
    mol : nx.Graph
        The molecule to be made aromatic.

    Returns
    -------
    None
        `mol` is modified in-place.

    Raises
    ------
    ValueError
        If there are atoms which are specified to be aromatic that are not in a
        ring.
    """
    elements = nx.get_node_attributes(mol, 'element')
    for cycle in nx.cycle_basis(mol):
        # If all elements in cycle are lowercase (or missing, *) it's aromatic
        if all(elements.get(n_idx, 'x').islower() for n_idx in cycle):
            for u, v in mol.edges(nbunch=cycle):
                if not (u in cycle and v in cycle):
                    continue
                mol.edges[u, v]['order'] = 1.5
            for n_idx in cycle:
                mol.nodes[n_idx]['element'] = mol.nodes[n_idx]['element'].upper()

    for n_idx in mol:
        if mol.nodes[n_idx].get('element', 'X').islower():
            raise ValueError("You specified an aromatic atom outside of a"
                             " ring. This is impossible")


def read_smiles(smiles, explicit_H=True):
    """
    Parses a SMILES string.

    Parameters
    ----------
    smiles : iterable
        The SMILES string to parse. Should conform to the OpenSMILES
        specification.
    explicit_H : bool
        Whether hydrogens should be explicit nodes in the outout graph, or be
        implicit in 'hcount' attributes.

    Returns
    -------
    nx.Graph
        A graph describing a molecule. Nodes will have an 'element' and a
        'charge', and if `explicit_H` is False a 'hcount'. Depending on the
        input, they will also have 'isotope' and 'class' information.
        Edges will have an 'order'.
    """
    bond_to_order = {'-': 1, '=': 2, '#': 3, '$': 4, ':': 1.5, '.': 0}
    mol = nx.Graph()
    anchor = None
    idx = 0
    default_bond = 1
    next_bond = None
    branches = []
    ring_nums = {}
    for tokentype, token in tokenize(smiles):
        if tokentype == TokenType.ATOM:
            mol.add_node(idx, **parse_atom(token))
            if anchor is not None:
                if next_bond is None:
                    next_bond = default_bond
                mol.add_edge(anchor, idx, order=next_bond)
                next_bond = None
            anchor = idx
            idx += 1
        elif tokentype == TokenType.BRANCH_START:
            branches.append(anchor)
        elif tokentype == TokenType.BRANCH_END:
            anchor = branches.pop()
        elif tokentype == TokenType.BOND_TYPE:
            next_bond = bond_to_order[token]
        elif tokentype == TokenType.RING_NUM:
            if token in ring_nums:
                jdx, order = ring_nums[token]
                if next_bond is None and order is None:
                    next_bond = default_bond
                elif order is None:  # Note that the check is needed,
                    next_bond = next_bond  # But this could be pass.
                elif next_bond is None:
                    next_bond = order
                elif next_bond != order:  # Both are not None
                    raise ValueError('Conflicting bond orders for ring '
                                     'between indices {}'.format(token))
                # idx is the index of the *next* atom we're adding. So: -1.
                mol.add_edge(idx - 1, jdx, order=next_bond)
                next_bond = None
                del ring_nums[token]
            else:
                # idx is the index of the *next* atom we're adding. So: -1.
                ring_nums[token] = (idx - 1, next_bond)
                next_bond = None

    # Time to deal with aromaticity
    aromatize_bonds(mol)

    # Add Hydrogens
    fill_valence(mol)

    if explicit_H:
        add_hydrogens(mol)
    else:
        remove_hydrogens(mol)
    return mol


if __name__ == '__main__':
    example = 'OCC(CCC)C(C(C)C)CCC'
    mol = read_smiles(example)
    mol = read_smiles('c1ccccc1')
    mol = read_smiles('C(C(C(C(C(C(C(C(C(C(C(C(C(C(C(C(C(C(C(C(C))))))))))))))))))))C')
    mol = read_smiles('N1CC2CCCC2CC1')
    mol = read_smiles('C%25CCCCC%25')
    mol3 = read_smiles('C1CCCCC1C1CCCCC1')
    mol4 = read_smiles('C1CC[13CH2]CC1C1CCCCC1')
    print(mol4.nodes(data=True))
    mol = read_smiles('[Rh-](Cl)(Cl)(Cl)(Cl)$[Rh-](Cl)(Cl)(Cl)Cl')
    mol2 = read_smiles('[15OH1-:4][HoH3]')
    molx = read_smiles('[O--][13CH2+3][14CH3+2][C@H4][Rh@OH19]')
    spiro = read_smiles('C12(CCCCC1)CC([H])CCC2*', False)

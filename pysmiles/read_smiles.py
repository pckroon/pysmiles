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

import enum

import networkx as nx

from .smiles_helper import (add_hydrogens, remove_hydrogens, parse_atom,
                            fill_valence)


@enum.unique
class TokenType(enum.Enum):
    ATOM = 1
    BOND_TYPE = 2
    BRANCH_START = 3
    BRANCH_END = 4
    RING_NUM = 5


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
        char = peek if peek else next(smiles, '')
        peek = None
        if not char:
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
            # If smiles is too short this will raise a ValueError, which is
            # (slightly) prettier than a StopIteration.
            yield TokenType.RING_NUM, int(next(smiles, '') + next(smiles, ''))
        elif char.isdigit():
            yield TokenType.RING_NUM, int(char)


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
            for idx, jdx in mol.edges(nbunch=cycle):
                if not (idx in cycle and jdx in cycle):
                    continue
                mol.edges[idx, jdx]['order'] = 1.5
            for n_idx in cycle:
                mol.nodes[n_idx]['element'] = mol.nodes[n_idx]['element'].upper()

    for n_idx in mol:
        if mol.nodes[n_idx].get('element', 'X').islower():
            raise ValueError("You specified an aromatic atom outside of a"
                             " ring. This is impossible")


def read_smiles(smiles, explicit_H=True, zero_order_bonds=True):
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
                if next_bond or zero_order_bonds:
                    mol.add_edge(anchor, idx, order=next_bond)
                next_bond = None
            anchor = idx
            idx += 1
        elif tokentype == TokenType.BRANCH_START:
            branches.append(anchor)
        elif tokentype == TokenType.BRANCH_END:
            anchor = branches.pop()
        elif tokentype == TokenType.BOND_TYPE:
            if next_bond is not None:
                raise ValueError('Previous bond (order {}) not used. '
                                 'Overwritten by "{}"'.format(next_bond, token))
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
                if mol.has_edge(idx-1, jdx):
                    raise ValueError('Edge specified by marker {} already '
                                     'exists'.format(token))
                elif idx-1 == jdx:
                    raise ValueError('Marker {} specifies a bond between an '
                                     'atom and itself'.format(token))
                if next_bond or zero_order_bonds:
                    mol.add_edge(idx - 1, jdx, order=next_bond)
                next_bond = None
                del ring_nums[token]
            else:
                if idx == 0:
                    raise ValueError("Can't have a marker ({}) before an atom"
                                     "".format(token))
                # idx is the index of the *next* atom we're adding. So: -1.
                ring_nums[token] = (idx - 1, next_bond)
                next_bond = None
    if ring_nums:
        raise KeyError('Unmatched ring indices {}'.format(list(ring_nums.keys())))
    # Time to deal with aromaticity
    aromatize_bonds(mol)

    # Add Hydrogens
    fill_valence(mol)

    if explicit_H:
        add_hydrogens(mol)
    else:
        remove_hydrogens(mol)
    return mol

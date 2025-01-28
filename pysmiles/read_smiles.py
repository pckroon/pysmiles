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

"""
Exposes functionality needed for parsing SMILES strings.
"""

import enum
import logging

import networkx as nx

from . import PTE
from .smiles_helper import (add_explicit_hydrogens, remove_explicit_hydrogens,
                            parse_atom, fill_valence, bonds_missing, format_atom,
                            correct_aromatic_rings,
                            _mark_chiral_atoms, _annotate_ez_isomers)
from .write_smiles import write_smiles_component

LOGGER = logging.getLogger(__name__)


@enum.unique
class TokenType(enum.Enum):
    """Possible SMILES token types"""
    ATOM = 1
    BOND_TYPE = 2
    BRANCH_START = 3
    BRANCH_END = 4
    RING_NUM = 5
    EZSTEREO = 6
    CHIRAL = 7


def _tokenize(smiles):
    """
    Iterates over a SMILES string, yielding tokens.

    Parameters
    ----------
    smiles : iterable
        The SMILES string to iterate over

    Yields
    ------
    tuple(TokenType, int, str)
        A tuple describing the type of token, its position, and the associated
        data.
    """
    organic_subset = 'B C N O P S F Cl Br I * b c n o s p'.split()
    smiles = iter(smiles)
    token = ''
    idx = -1
    peek = None
    while True:
        idx += 1
        char = peek if peek else next(smiles, '')
        peek = None
        if not char:
            break
        if char == '[':
            token = char
            for char in smiles:  # pragma: no branch
                token += char
                if char == ']':
                    break
            yield TokenType.ATOM, idx, token
        elif char in organic_subset:
            peek = next(smiles, '')
            if char + peek in organic_subset:
                yield TokenType.ATOM, idx, char + peek
                peek = None
            else:
                yield TokenType.ATOM, idx, char
        elif char in '-=#$:.':
            yield TokenType.BOND_TYPE, idx, char
        elif char == '(':
            yield TokenType.BRANCH_START, idx, '('
        elif char == ')':
            yield TokenType.BRANCH_END, idx, ')'
        elif char == '%':
            # If smiles is too short this will raise a ValueError, which is
            # (slightly) prettier than a StopIteration.
            yield TokenType.RING_NUM, idx, int(next(smiles, '') + next(smiles, ''))
        elif char in '/\\':
            yield TokenType.EZSTEREO, idx, char
        elif char.isdigit():  # pragma: no branch
            yield TokenType.RING_NUM, idx, int(char)


def base_smiles_parser(smiles, strict=True, node_attr='desc', edge_attr='desc'):
    """
    Parse but do not interpret a SMILES string to a graph. Nodes and edges will
    be annotated with their constructing string description as ``node_attr`` and
    ``edge_attr``. The indices of these strings are annotated in the ``'_pos'``
    attribute. Finally, the smiles string itself is annotated as a graph
    attribute under ``'smiles'``.

    This functions also returns stereochemical information that cannot be
    directly annotated in the graph, since it depends on interpreting nodes and
    edges.

    Parameters
    ----------
    smiles : iterable
    strict : bool
    node_attr : collections.abc.Hashable
    edge_attr : collections.abc.Hashable

    Returns
    -------
    nx.Graph
        The created graph. Nodes are annotated with the string/token that define
        their attributes in the ``node_attr`` key, edges in the ``edge_attr``
        key. The indices at which these strings can be found is annotated as
        ``_pos``.
    List[Tuple[Tuple[int, int, str], Tuple[int, int, str]]]
        A list of E/Z isomer pairs, where the integers are node indices and last
        element is the directional token from the smiles string (*i.e.* / or \).
    List[Tuple[int, int, {edge_attr: str, '_pos': int}]
        A list of created ring bonds, to be used for R/S isomer annotation. It
        contains tuples consisting of 2 node indices, and a dictionary
        containing the edge attributes.
    """
    mol = nx.Graph(smiles=smiles)
    anchor = None
    idx = 0
    next_bond = None
    branches = []
    ring_nums = {}
    current_ez = None
    ez_isomer_atoms = {}
    created_ring_bonds = []
    prev_token = None
    prev_type = None
    for tokentype, token_idx, token in _tokenize(smiles):
        if tokentype == TokenType.ATOM:
            mol.add_node(idx, **{node_attr: token, '_pos': token_idx})
            if anchor is not None:
                if next_bond is None:
                    next_bond = ""
                mol.add_edge(anchor, idx, **{edge_attr: next_bond, '_pos': token_idx})
                next_bond = None
            anchor = idx
            idx += 1
        elif tokentype == TokenType.BRANCH_START:
            if anchor is None:
                raise SyntaxError('You cannot start a branch before an anchor.')
            branches.append(anchor)
        elif tokentype == TokenType.BRANCH_END:
            anchor = branches.pop()
        elif tokentype == TokenType.BOND_TYPE:
            if next_bond is not None:
                raise ValueError('Previous bond (order {}) not used. '
                                 'Overwritten by "{}"'.format(next_bond, token))
            next_bond = token
        elif tokentype == TokenType.RING_NUM:
            if token in ring_nums:
                jdx, order = ring_nums[token]
                if next_bond is None and order is None:
                    next_bond = ""
                elif order is None:  # Note that the check is needed,
                    next_bond = next_bond  # But this could be pass.
                elif next_bond is None:
                    next_bond = order
                elif next_bond != order:  # Both are not None
                    raise ValueError('Conflicting bond orders for ring '
                                     'between indices {}'.format(token))
                # idx is the index of the *next* atom we're adding. So: -1.
                if mol.has_edge(idx - 1, jdx):
                    raise ValueError('Edge specified by marker {} already '
                                     'exists'.format(token))
                if idx - 1 == jdx:
                    raise ValueError('Marker {} specifies a bond between an '
                                     'atom and itself'.format(token))
                mol.add_edge(idx - 1, jdx, **{edge_attr: next_bond, '_pos': token_idx})
                next_bond = None
                del ring_nums[token]
                # we need to keep track of ring bonds here for the
                # chirality assignment
                created_ring_bonds.append((idx-1, jdx, {edge_attr: next_bond, '_pos': token_idx}))
            else:
                if idx == 0:
                    raise ValueError("Can't have a marker ({}) before an atom"
                                     "".format(token))
                # idx is the index of the *next* atom we're adding. So: -1.
                ring_nums[token] = (idx - 1, next_bond)
                next_bond = None
        elif tokentype == TokenType.EZSTEREO:  # pragma: no branch
            ez_isomer_atoms[anchor] = token
            ez_isomer_atoms[idx] = token
        prev_token = token
        prev_type = tokentype

    if ring_nums and strict:
        raise KeyError('Unmatched ring indices {}'.format(list(ring_nums.keys())))

    if current_ez and strict:
        raise ValueError('There is an unmatched stereochemical token.')

    return mol, ez_isomer_atoms, created_ring_bonds

def read_smiles(smiles, explicit_hydrogen=False, zero_order_bonds=True, 
                reinterpret_aromatic=True, strict=True):
    """
    Parses a SMILES string.

    Parameters
    ----------
    smiles : iterable
        The SMILES string to parse. Should conform to the OpenSMILES
        specification.
    explicit_hydrogen : bool
        Whether hydrogens should be explicit nodes in the output graph, or be
        implicit in 'hcount' attributes.
    zero_order_bonds : bool
        Whether zero-order bonds (".") should be added as edges with an order of
        0.
    reinterpret_aromatic : bool
        Whether aromaticity should be determined from the created molecule,
        instead of taken from the SMILES string.
    strict : bool
        Whether to be more strict in accepting what is a valid SMILES string and
        what is not.

    Returns
    -------
    nx.Graph
        A graph describing a molecule. Nodes will have an 'element', 'aromatic'
        and a 'charge', and if `explicit_hydrogen` is False a 'hcount'.
        Depending on the input, they will also have 'isotope' and 'class'
        information.
        Edges will have an 'order'.
    """
    bond_to_order = {'-': 1, '=': 2, '#': 3, '$': 4, ':': 1.5, '.': 0}
    default_bond = 1
    default_aromatic_bond = 1.5
    mol, ez_isomer_atoms, ring_bonds = base_smiles_parser(smiles, strict=strict,
                                                          node_attr='_atom_str', edge_attr='_bond_str')
    for node in mol:
        mol.nodes[node].update(parse_atom(mol.nodes[node]['_atom_str']))
    for idx, jdx, attrs in ring_bonds:
        if mol.nodes[idx].get('rs_isomer'):
            mol.nodes[idx]['rs_isomer'][1].append(jdx)
        if mol.nodes[jdx].get('rs_isomer'):
            mol.nodes[jdx]['rs_isomer'][1].append(idx)
    for edge in mol.edges:
        if mol.edges[edge]['_bond_str']:
            order = bond_to_order[mol.edges[edge]['_bond_str']]
            if not order and not zero_order_bonds:
                mol.remove_edge(*edge)
                continue
            mol.edges[edge]['order'] = order
        elif mol.nodes[edge[0]].get('aromatic') and mol.nodes[edge[1]].get('aromatic'):
            mol.edges[edge]['order'] = default_aromatic_bond
        else:
            mol.edges[edge]['order'] = default_bond

    if reinterpret_aromatic:
        correct_aromatic_rings(mol, strict=strict)

    # This is a bit of an overreach, we should only add implicit hydrogens to
    # atoms in the organic subset. However, all non-organic atoms already have
    # a hcount, so all is well.
    fill_valence(mol)

    for node in mol:
        element = mol.nodes[node].get('element', '*')
        if element != '*' and element not in PTE:
            msg = f'Unknown element {element}'
            if strict:
                raise KeyError(msg)
            else:
                LOGGER.warning(msg)
        elif element != '*' and bonds_missing(mol, node):
            attrs = mol.nodes[node]
            # Grab some of the smiles string for debugging purposes
            surrounding_characters = 5
            node_position = (attrs['_pos'] - surrounding_characters,
                             attrs['_pos'] + len(attrs['_atom_str']) + surrounding_characters)
            if node_position[0] < 0:
                left_idx = 0
                prefix = ''
            else:
                left_idx = node_position[0]
                prefix = '...'
            if node_position[1] >= len(smiles):
                right_idx = len(smiles) - 1
                suffix = ''
            else:
                right_idx = node_position[1]
                suffix = '...'
            debug_smiles = prefix + smiles[left_idx:right_idx] + suffix
            msg = (f'Node {node} ({format_atom(mol, node)}) has non-standard'
                   f' valence: {debug_smiles}')
            if strict:
                raise KeyError(msg)
            else:
                LOGGER.warning(msg)

    # post-processing of E/Z isomerism
    _annotate_ez_isomers(mol, ez_isomer_atoms)
    # post-processing of chiral atoms
    _mark_chiral_atoms(mol)

    if explicit_hydrogen:
        add_explicit_hydrogens(mol)
    else:
        remove_explicit_hydrogens(mol)

    return mol

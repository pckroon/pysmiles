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
                            parse_atom, fill_valence, mark_aromatic_edges,
                            mark_aromatic_atoms,  bonds_missing, format_atom,
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
            for char in smiles:  # pragma: no branch
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
        elif char in '/\\':
            yield TokenType.EZSTEREO, char
        elif char.isdigit():  # pragma: no branch
            yield TokenType.RING_NUM, int(char)


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
    mol = nx.Graph()
    anchor = None
    idx = 0
    default_bond = 1
    next_bond = None
    branches = []
    ring_nums = {}
    current_ez = None
    prev_token = None
    prev_type = None
    ez_isomer_pairs = []
    for tokentype, token in _tokenize(smiles):
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
                if idx-1 == jdx:
                    raise ValueError('Marker {} specifies a bond between an '
                                     'atom and itself'.format(token))
                if next_bond or zero_order_bonds:
                    mol.add_edge(idx - 1, jdx, order=next_bond)
                    # we need to keep track of ring bonds here for the
                    # chirality assignment
                    if mol.nodes[idx-1].get('stereo', False):
                        mol.nodes[idx-1]['stereo'][1].append(jdx)
                    if mol.nodes[jdx].get('stereo', False):
                        mol.nodes[jdx]['stereo'][1].append(idx-1)

                next_bond = None
                del ring_nums[token]
            else:
                if idx == 0:
                    raise ValueError("Can't have a marker ({}) before an atom"
                                     "".format(token))
                # idx is the index of the *next* atom we're adding. So: -1.
                ring_nums[token] = (idx - 1, next_bond)
                next_bond = None
        elif tokentype == TokenType.EZSTEREO:  # pragma: no branch
            # FIXME "It is permissible, but not required, that every atom attached to a double bond be marked."
            # we found the second ez reference and
            # annotate the molecule
            if current_ez:
                ez_isomer_pairs.append((current_ez, (idx, anchor, token)))
                current_ez = None
            # current_ez is formatted as:
            # ligand, anchor, token where ligand is the atom defining CIS/TRANS
            # we found a token belonging to a branch (e.g. C(\F))
            elif prev_type == TokenType.BRANCH_START:
                current_ez = (idx, anchor, token)
#            elif prev_type == TokenType.BRANCH_STOP:
#                raise ValueError(f"The E/Z token {token} may not follow a closing bracket.")
            else:
                current_ez = (anchor, idx, token)

        prev_token = token
        prev_type = tokentype

    if ring_nums and strict:
        raise KeyError('Unmatched ring indices {}'.format(list(ring_nums.keys())))

    if current_ez and strict:
        raise ValueError('There is an unmatched stereochemical token.')

    if reinterpret_aromatic:
        mark_aromatic_atoms(mol, strict=strict)
        mark_aromatic_edges(mol)
    else:
        mark_aromatic_edges(mol)

    # This is a bit of an overreach, we should only add implicit hydrogens to
    # atoms in the organic subset. However, all non-organic atoms already have
    # a hcount, so all is well.
    fill_valence(mol)

    if strict:
        for node in mol:
            element = mol.nodes[node].get('element', '*')
            if element != '*' and element not in PTE:
                raise KeyError(f'Unknown element {element}')
            elif element != '*' and bonds_missing(mol, node):
                debug_smiles = write_smiles_component(nx.ego_graph(mol, node))
                raise KeyError(f'Node {node} ({format_atom(mol, node)}) has'
                               f' non-standard valence: ...{debug_smiles}...')

    # post-processing of E/Z isomerism
    _annotate_ez_isomers(mol, ez_isomer_pairs)

    # post-processing of chiral atoms
    # potentially we need all hydrogen in place
    _mark_chiral_atoms(mol)

    if explicit_hydrogen:
        add_explicit_hydrogens(mol)
    else:
        remove_explicit_hydrogens(mol)

    return mol

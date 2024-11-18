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


def base_smiles_parser(smiles, strict=True, node_attr='desc', edge_attr='desc'):
    """
    Parse but do not interpret a SMILES string to a graph. Nodes and edges will
    be annotated with their constructing string description.

    Parameters
    ----------
    smiles : str

    Returns
    -------
    nx.Graph
    """
    mol = nx.Graph(smiles=smiles)
    anchor = None
    idx = 0
    next_bond = None
    branches = []
    ring_nums = {}
    current_ez = None
    ez_isomer_pairs = []
    prev_token = None
    prev_type = None
    for tokentype, token in _tokenize(smiles):
        if tokentype == TokenType.ATOM:
            mol.add_node(idx, **{node_attr: token, 'rs_bond_order': []})
            if anchor is not None:
                if next_bond is None:
                    next_bond = ""
                mol.add_edge(anchor, idx, **{edge_attr: next_bond})
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
                mol.add_edge(idx - 1, jdx, **{edge_attr: next_bond})
                next_bond = None
                del ring_nums[token]
                # we need to keep track of ring bonds here for the
                # chirality assignment
                mol.nodes[idx - 1]['rs_bond_order'].append(jdx)
                mol.nodes[jdx]['rs_bond_order'].append(idx - 1)
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
            else:
                current_ez = (anchor, idx, token)

        prev_token = token
        prev_type = tokentype

    if ring_nums and strict:
        raise KeyError('Unmatched ring indices {}'.format(list(ring_nums.keys())))

    if current_ez and strict:
        raise ValueError('There is an unmatched stereochemical token.')

    return mol, ez_isomer_pairs

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
    mol, ez_isomer_pairs = base_smiles_parser(smiles, strict=strict, node_attr='atom_str', edge_attr='bond_str')
    for node in mol:
        mol.nodes[node].update(parse_atom(mol.nodes[node]['atom_str']))
        if mol.nodes[node].get('rs_isomer') and mol.nodes[node]['rs_bond_order']:
            mol.nodes[node]['rs_isomer'][1].extend(mol.nodes[node]['rs_bond_order'])
        del mol.nodes[node]['atom_str']
        if 'rs_bond_order' in mol.nodes[node]:
            del mol.nodes[node]['rs_bond_order']
    for edge in mol.edges:
        if mol.edges[edge]['bond_str']:
            order = bond_to_order[mol.edges[edge]['bond_str']]
            if not order and not zero_order_bonds:
                mol.remove_edge(*edge)
                continue
            mol.edges[edge]['order'] = order
        elif mol.nodes[edge[0]].get('aromatic') and mol.nodes[edge[1]].get('aromatic'):
            mol.edges[edge]['order'] = default_aromatic_bond
        else:
            mol.edges[edge]['order'] = default_bond
        del mol.edges[edge]['bond_str']

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
            debug_smiles = write_smiles_component(nx.ego_graph(mol, node))
            msg = (f'Node {node} ({format_atom(mol, node)}) has non-standard'
                   f' valence: ...{debug_smiles}...')
            if strict:
                raise KeyError(msg)
            else:
                LOGGER.warning(msg)

    # post-processing of E/Z isomerism
    _annotate_ez_isomers(mol, ez_isomer_pairs)
    # post-processing of chiral atoms
    _mark_chiral_atoms(mol)

    if explicit_hydrogen:
        add_explicit_hydrogens(mol)
    else:
        remove_explicit_hydrogens(mol)

    return mol

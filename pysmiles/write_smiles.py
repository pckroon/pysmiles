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
Exposes functionality for writing SMILES strings
"""

from collections import defaultdict

import networkx as nx

from .smiles_helper import remove_explicit_hydrogens, format_atom


def _get_ring_marker(used_markers):
    """
    Returns the lowest number larger than 0 that is not in `used_markers`.

    Parameters
    ----------
    used_markers : Container
        The numbers that can't be used.

    Returns
    -------
    int
        The lowest number larger than 0 that's not in `used_markers`.
    """
    new_marker = 1
    while new_marker in used_markers:
        new_marker += 1
    return new_marker


def _write_edge_symbol(molecule, n_idx, n_jdx):
    """
    Determines whether a symbol should be written for the edge between `n_idx`
    and `n_jdx` in `molecule`. It should not be written if it's a bond of order
    1 or an aromatic bond between two aromatic atoms; unless it's a single bond
    between two aromatic atoms.

    Parameters
    ----------
    molecule : nx.Graph
        The molecule.
    n_idx : Hashable
        The first node key describing the edge.
    n_jdx : Hashable
        The second node key describing the edge.

    Returns
    -------
    bool
        Whether an explicit symbol is needed for this edge.
    """
    order = molecule.edges[n_idx, n_jdx].get('order', 1)
    aromatic_atoms = molecule.nodes[n_idx].get('element', '*').islower() and\
                     molecule.nodes[n_jdx].get('element', '*').islower()
    aromatic_bond = aromatic_atoms and order == 1.5
    cross_aromatic = aromatic_atoms and order == 1
    single_bond = order == 1
    return cross_aromatic or not (aromatic_bond or single_bond)


def write_smiles(molecule, default_element='*', start=None):
    """
    Creates a SMILES string describing `molecule` according to the OpenSMILES
    standard.

    Parameters
    ----------
    molecule : nx.Graph
        The molecule for which a SMILES string should be generated.
    default_element : str
        The element to write if the attribute is missing for a node.
    start : Hashable
        The atom at which the depth first traversal of the molecule should
        start. A sensible one is chosen: preferably a terminal heteroatom.

    Returns
    -------
    str
        The SMILES string describing `molecule`.
    """
    molecule = molecule.copy()
    remove_explicit_hydrogens(molecule)

    if start is None:
        # Start at a terminal atom, and if possible, a heteroatom.
        def keyfunc(idx):
            """Key function for finding the node at which to start."""
            return (molecule.degree(idx),
                    # True > False
                    molecule.nodes[idx].get('element', default_element) == 'C',
                    idx)
        start = min(molecule.nodes, key=keyfunc)


    order_to_symbol = {0: '.', 1: '-', 1.5: ':', 2: '=', 3: '#', 4: '$'}

    dfs_successors = nx.dfs_successors(molecule, source=start)

    predecessors = defaultdict(list)
    for node_key, successors in dfs_successors.items():
        for successor in successors:
            predecessors[successor].append(node_key)
    predecessors = dict(predecessors)
    # We need to figure out which edges we won't cross when doing the dfs.
    # These are the edges we'll need to add to the smiles using ring markers.
    edges = set()
    for n_idx, n_jdxs in dfs_successors.items():
        for n_jdx in n_jdxs:
            edges.add(frozenset((n_idx, n_jdx)))
    total_edges = set(map(frozenset, molecule.edges))
    ring_edges = total_edges - edges

    atom_to_ring_idx = defaultdict(list)
    ring_idx_to_bond = {}
    ring_idx_to_marker = {}
    for ring_idx, (n_idx, n_jdx) in enumerate(ring_edges, 1):
        atom_to_ring_idx[n_idx].append(ring_idx)
        atom_to_ring_idx[n_jdx].append(ring_idx)
        ring_idx_to_bond[ring_idx] = (n_idx, n_jdx)

    branch_depth = 0
    branches = set()
    to_visit = [start]
    smiles = ''

    while to_visit:
        current = to_visit.pop()
        if current in branches:
            branch_depth += 1
            smiles += '('
            branches.remove(current)

        if current in predecessors:
            # It's not the first atom we're visiting, so we want to see if the
            # edge we last crossed to get here is interesting.
            previous = predecessors[current]
            assert len(previous) == 1
            previous = previous[0]
            if _write_edge_symbol(molecule, previous, current):
                order = molecule.edges[previous, current].get('order', 1)
                smiles += order_to_symbol[order]
        smiles += format_atom(molecule, current, default_element)
        if current in atom_to_ring_idx:
            # We're going to need to write a ring number
            ring_idxs = atom_to_ring_idx[current]
            for ring_idx in ring_idxs:
                ring_bond = ring_idx_to_bond[ring_idx]
                if ring_idx not in ring_idx_to_marker:
                    marker = _get_ring_marker(ring_idx_to_marker.values())
                    ring_idx_to_marker[ring_idx] = marker
                    new_marker = True
                else:
                    marker = ring_idx_to_marker.pop(ring_idx)
                    new_marker = False

                if _write_edge_symbol(molecule, *ring_bond) and new_marker:
                    order = molecule.edges[ring_bond].get('order', 1)
                    smiles += order_to_symbol[order]
                smiles += str(marker) if marker < 10 else '%{}'.format(marker)

        if current in dfs_successors:
            # Proceed to the next node in this branch
            next_nodes = dfs_successors[current]
            # ... and if needed, remember to return here later
            branches.update(next_nodes[1:])
            to_visit.extend(next_nodes)
        elif branch_depth:
            # We're finished with this branch.
            smiles += ')'
            branch_depth -= 1

    smiles += ')' * branch_depth
    return smiles

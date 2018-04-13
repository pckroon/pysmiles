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

from collections import defaultdict

import networkx as nx

from .smiles_helper import remove_hydrogens, fix_aromatic


def write_smiles(molecule, default_element='*', start=None):
    if start is None:
        start = min(molecule.nodes)

    molecule = molecule.copy()
    remove_hydrogens(molecule)
    fix_aromatic(molecule)

    def get_symbol(node_key):
        name = molecule.nodes[node_key].get('element', default_element)
        charge = molecule.nodes[node_key].get('charge', 0)
        stereo = molecule.nodes[node_key].get('stereo', None)
        isotope = molecule.nodes[node_key].get('isotope', '')
        class_ = molecule.nodes[node_key].get('class', '')
        if stereo is not None:
            raise NotImplementedError
        if stereo is None and isotope == '' and charge == 0 and\
                name in 'B C N O P S F Cl Br I b c n o p s se as'.split():
            return name
        if charge > 0:
            chargestr = '+'
            if charge > 1:
                chargestr += str(charge)
        elif charge < 0:
            chargestr = '-'
            if charge < -1:
                chargestr += str(-charge)
        else:
            chargestr = ''
        if class_ != '':
            class_ = ':{}'.format(class_)
        return '[{isotope}{name}{charge}{class_}]'.format(isotope=isotope,
                                                          name=name,
                                                          charge=chargestr,
                                                          class_=class_)

    order_to_symbol = {0: '.', 1: '', 1.5: '', 2: '=', 3: '#', 4: '$'}

    dfs_successors = nx.dfs_successors(molecule)
    predecessors = defaultdict(list)
    for node_key, successors in dfs_successors.items():
        for successor in successors:
            predecessors[successor].append(node_key)
    predecessors = dict(predecessors)

    edges = set()
    for u, vs in dfs_successors.items():
        for v in vs:
            edges.add(frozenset((u, v)))
    total_edges = set(map(frozenset, molecule.edges))
    ring_edges = total_edges - edges

    ring_markers = {}
    marker_to_bond = {}
    for marker, (u, v) in enumerate(ring_edges, 1):
        marker = str(marker) if marker < 10 else '%{}'.format(marker)
        ring_markers[u] = marker
        ring_markers[v] = marker
        marker_to_bond[marker] = (u, v)
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
            previous = predecessors[current]
            assert len(previous) == 1
            previous = previous[0]
            order = molecule.edges[previous, current].get('order', 1)
            smiles += order_to_symbol[order]
        smiles += get_symbol(current)
        if current in ring_markers:
            marker = ring_markers[current]
            ring_bond = marker_to_bond[marker]
            order = molecule.edges[ring_bond].get('order', 1)
            smiles += order_to_symbol[order]
            smiles += marker

        if current in dfs_successors:
            next_nodes = dfs_successors[current]
            branches.update(next_nodes[1:])
            to_visit.extend(next_nodes)
        elif branch_depth:
            smiles += ')'
            branch_depth -= 1

    smiles += ')' * branch_depth
    return smiles

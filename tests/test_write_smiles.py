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

from pysmiles import write_smiles
from pysmiles.smiles_helper import mark_aromatic_atoms

mol = nx.Graph()
mol.add_edges_from([(0, 1), (1, 2), (1, 3), (3, 4)], order=1)
for idx, (ele, count) in enumerate(zip('CCOCC', (3, 0, 0, 2, 3))):
    mol.nodes[idx]['element'] = ele
#    mol.nodes[idx]['hcount'] = count
mol.edges[1, 2]['order'] = 2

print(write_smiles(mol))

mol = nx.cycle_graph(6)
mol.add_edge(3, 6, order=1)
for node_key in mol:
    mol.nodes[node_key]['element'] = 'C'
mol.nodes[6]['element'] = 'O'
for edge in [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)]:
    mol.edges[edge]['order'] = 1.5
print(write_smiles(mol))

mol = nx.Graph()
mol.add_edges_from([(0, 1), (1, 2), (1, 3), (3, 4), (1, 5)], order=1)
for idx, ele in enumerate('ABCDEF'):
    mol.nodes[idx]['element'] = ele
print(write_smiles(mol))

mol = nx.Graph()
mol.add_edges_from([(0, 1), (1, 2), (1, 3), (3, 4), (1, 5), (3, 6)], order=1)
for idx, ele in enumerate('CCCCOCO'):
    mol.nodes[idx]['element'] = ele
mol.nodes[4]['charge'] = -1
mol.edges[3, 6]['order'] = 2
print(write_smiles(mol))

mol = nx.Graph()
mol.add_cycle(range(6), order=1.5)
mol.add_cycle(range(6, 12), order=1.5)
mol.add_edge(5, 7, order=1)
for n_idx in mol:
    mol.nodes[n_idx]['element'] = 'C'

print(mol.edges(data=True))
print(write_smiles(mol))

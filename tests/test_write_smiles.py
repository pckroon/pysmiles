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

import unittest

import networkx as nx

from pysmiles import write_smiles, read_smiles
from pysmiles.smiles_helper import correct_aromatic_rings, fill_valence
from pysmiles.testhelper import GraphTest


class Tests(GraphTest):
    def test_simple(self):
        mol = nx.Graph()
        nx.add_path(mol, range(4), order=1)
        for n_idx in mol:
            mol.nodes[n_idx].update(element='C', charge=0)
        fill_valence(mol)

        smiles = write_smiles(mol)
        found = read_smiles(smiles)
        self.assertEqualGraphs(mol, found)

    def test_hetero(self):
        mol2 = nx.Graph()
        mol2.add_edges_from([(0, 1), (1, 2), (1, 3), (3, 4)], order=1)
        for idx, (ele, count) in enumerate(zip('CCOCC', (3, 0, 0, 2, 3))):
            mol2.nodes[idx]['element'] = ele
            mol2.nodes[idx]['hcount'] = count
            mol2.nodes[idx]['charge'] = 0
        mol2.edges[1, 2]['order'] = 2

        smiles = write_smiles(mol2)
        found = read_smiles(smiles)
        self.assertEqualGraphs(found, mol2)

    def test_cycle(self):
        mol = nx.cycle_graph(6)
        mol.add_edge(3, 6, order=1)
        for node_key in mol:
            mol.nodes[node_key]['element'] = 'C'
            mol.nodes[node_key]['charge'] = 0
        mol.nodes[6]['element'] = 'O'
        for edge in [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)]:
            mol.edges[edge]['order'] = 1.5
        fill_valence(mol)
        correct_aromatic_rings(mol)
        smiles = write_smiles(mol)
        found = read_smiles(smiles)
        self.assertEqualGraphs(mol, found)

    def test_aromatic(self):
        mol = nx.Graph()
        mol.add_edges_from([(0, 1), (1, 2), (1, 3), (3, 4), (1, 5), (3, 6)], order=1)
        for idx, ele in enumerate('CCCCOCO'):
            mol.nodes[idx]['element'] = ele
            mol.nodes[idx]['charge'] = 0
        mol.nodes[4]['charge'] = -1
        mol.nodes[4]['hcount'] = 0
        mol.edges[3, 6]['order'] = 2
        fill_valence(mol)
        correct_aromatic_rings(mol)
        smiles = write_smiles(mol)
        found = read_smiles(smiles)
        self.assertEqualGraphs(found, mol)

    def test_multiring(self):
        mol = nx.Graph()
        mol.add_cycle(range(6), order=1.5)
        mol.add_cycle(range(6, 12), order=1.5)
        mol.add_edge(5, 7, order=1)
        for n_idx in mol:
            mol.nodes[n_idx]['element'] = 'C'
            mol.nodes[n_idx]['charge'] = 0
        fill_valence(mol)
        correct_aromatic_rings(mol)
        smiles = write_smiles(mol)
        found = read_smiles(smiles)
        self.assertEqualGraphs(found, mol)


if __name__ == '__main__':
    unittest.main()

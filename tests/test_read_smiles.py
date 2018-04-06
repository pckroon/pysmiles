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

import operator
import unittest

import networkx as nx

from pysmiles import read_smiles


class Tests(unittest.TestCase):
    def assertEqualGraphs(self, graph1, graph2):
        out = nx.is_isomorphic(graph1, graph2,
                               node_match=operator.eq, edge_match=operator.eq)
        if not out:
            simple_out = nx.is_isomorphic(graph1, graph2)
            if not simple_out:
                raise AssertionError("Graphs not isomorphic")
            else:
                raise AssertionError('Node or edge attributes do no match')

    def test_single_chain(self):
        smiles = 'CCCC'
        found = read_smiles(smiles, explicit_H=False)

        expected = nx.Graph()
        nx.add_path(expected, (0, 1, 2, 3))
        data = [(0, {'element': 'C', 'charge': 0, 'hcount': 3}),
                (1, {'element': 'C', 'charge': 0, 'hcount': 2}),
                (2, {'element': 'C', 'charge': 0, 'hcount': 2}),
                (3, {'element': 'C', 'charge': 0, 'hcount': 3})]
        data = dict(data)
        nx.set_node_attributes(expected, data)
        data = [((0, 1), {'order': 1}),
                ((1, 2), {'order': 1}),
                ((2, 3), {'order': 1})]
        data = dict(data)
        nx.set_edge_attributes(expected, data)

        self.assertEqualGraphs(found, expected)

    def test_branched(self):
        smiles = 'CCC(CC)CC'
        found = read_smiles(smiles, explicit_H=False)

        expected = nx.Graph()
        nx.add_path(expected, (0, 1, 2, 3, 4))
        expected.add_edges_from([(2, 5), (5, 6)])
        data = {0: {'element': 'C', 'charge': 0, 'hcount': 3},
                1: {'element': 'C', 'charge': 0, 'hcount': 2},
                2: {'element': 'C', 'charge': 0, 'hcount': 1},
                3: {'element': 'C', 'charge': 0, 'hcount': 2},
                4: {'element': 'C', 'charge': 0, 'hcount': 3},
                5: {'element': 'C', 'charge': 0, 'hcount': 2},
                6: {'element': 'C', 'charge': 0, 'hcount': 3}}
        nx.set_node_attributes(expected, data)
        data = [((0, 1), {'order': 1}),
                ((1, 2), {'order': 1}),
                ((2, 3), {'order': 1}),
                ((2, 5), {'order': 1}),
                ((3, 4), {'order': 1}),
                ((5, 6), {'order': 1})]
        data = dict(data)
        nx.set_edge_attributes(expected, data)
        self.assertEqualGraphs(found, expected)

    def test_double_bond(self):
        smiles = 'C=C'
        found = read_smiles(smiles, explicit_H=False)

        expected = nx.Graph()
        data = [(0, 1, {'order': 2})]
        expected.add_edges_from(data)
        data = {0: {'element': 'C', 'charge': 0, 'hcount': 2},
                1: {'element': 'C', 'charge': 0, 'hcount': 2}}
        nx.set_node_attributes(expected, data)

        self.assertEqualGraphs(found, expected)

    def test_ring(self):
        smiles = 'C1CC1'
        found = read_smiles(smiles, explicit_H=False)

        expected = nx.Graph()
        data = [(0, 1, {'order': 1}),
                (0, 2, {'order': 1}),
                (1, 2, {'order': 1})]
        expected.add_edges_from(data)
        data = {0: {'charge': 0, 'element': 'C', 'hcount': 2},
                1: {'charge': 0, 'element': 'C', 'hcount': 2},
                2: {'charge': 0, 'element': 'C', 'hcount': 2}}
        nx.set_node_attributes(expected, data)

        self.assertEqualGraphs(found, expected)

        smiles = "C%11CC%11"
        found = read_smiles(smiles, explicit_H=False)
        self.assertEqualGraphs(found, expected)

    def test_aromatic(self):
        smiles = 'c1ccccc1'
        found = read_smiles(smiles, explicit_H=False)

        expected = nx.Graph()
        data = [(0, 1, {'order': 1.5}),
                (0, 5, {'order': 1.5}),
                (1, 2, {'order': 1.5}),
                (2, 3, {'order': 1.5}),
                (3, 4, {'order': 1.5}),
                (4, 5, {'order': 1.5})]
        expected.add_edges_from(data)
        data = {0: {'charge': 0, 'element': 'C', 'hcount': 1},
                1: {'charge': 0, 'element': 'C', 'hcount': 1},
                2: {'charge': 0, 'element': 'C', 'hcount': 1},
                3: {'charge': 0, 'element': 'C', 'hcount': 1},
                4: {'charge': 0, 'element': 'C', 'hcount': 1},
                5: {'charge': 0, 'element': 'C', 'hcount': 1}}
        nx.set_node_attributes(expected, data)

        self.assertEqualGraphs(found, expected)

    def test_ring_bond(self):
        smiles = 'C1CC=1'
        found = read_smiles(smiles, explicit_H=False)

        expected = nx.Graph()
        data = [(0, 1, {'order': 1}),
                (0, 2, {'order': 2}),
                (1, 2, {'order': 1})]
        expected.add_edges_from(data)
        data = {0: {'charge': 0, 'element': 'C', 'hcount': 1},
                1: {'charge': 0, 'element': 'C', 'hcount': 2},
                2: {'charge': 0, 'element': 'C', 'hcount': 1}}
        nx.set_node_attributes(expected, data)

        self.assertEqualGraphs(found, expected)

        smiles = 'C=1CC=1'
        found = read_smiles(smiles, explicit_H=False)
        self.assertEqualGraphs(found, expected)

        smiles = 'C=1CC1'
        found = read_smiles(smiles, explicit_H=False)
        self.assertEqualGraphs(found, expected)

        smiles = 'C%11CC=%11'
        found = read_smiles(smiles, explicit_H=False)
        self.assertEqualGraphs(found, expected)

    def test_branch_bond(self):
        smiles = 'OS(=O)(=S)O'
        found = read_smiles(smiles, explicit_H=False)

        expected = nx.Graph()
        data = [(0, 1, {'order': 1}),
                (1, 2, {'order': 2}),
                (1, 3, {'order': 2}),
                (1, 4, {'order': 1})]
        expected.add_edges_from(data)
        data = {0: {'charge': 0, 'element': 'O', 'hcount': 1},
                1: {'charge': 0, 'element': 'S', 'hcount': 0},
                2: {'charge': 0, 'element': 'O', 'hcount': 0},
                3: {'charge': 0, 'element': 'S', 'hcount': 0},
                4: {'charge': 0, 'element': 'O', 'hcount': 1}}
        nx.set_node_attributes(expected, data)

        self.assertEqualGraphs(found, expected)

    def test_deep_branch(self):
        smiles = 'C(C(C(C(C(C(C(C(C(C(C(C(C(C(C(C(C(C(C(C(C))))))))))))))))))))C'
        found = read_smiles(smiles, explicit_H=False)

        expected = nx.Graph()
        data = [(0, 1, {'order': 1}),
                (0, 21, {'order': 1}),
                (1, 2, {'order': 1}),
                (2, 3, {'order': 1}),
                (3, 4, {'order': 1}),
                (4, 5, {'order': 1}),
                (5, 6, {'order': 1}),
                (6, 7, {'order': 1}),
                (7, 8, {'order': 1}),
                (8, 9, {'order': 1}),
                (9, 10, {'order': 1}),
                (10, 11, {'order': 1}),
                (11, 12, {'order': 1}),
                (12, 13, {'order': 1}),
                (13, 14, {'order': 1}),
                (14, 15, {'order': 1}),
                (15, 16, {'order': 1}),
                (16, 17, {'order': 1}),
                (17, 18, {'order': 1}),
                (18, 19, {'order': 1}),
                (19, 20, {'order': 1})]
        expected.add_edges_from(data)
        data = {0: {'charge': 0, 'element': 'C', 'hcount': 2},
                1: {'charge': 0, 'element': 'C', 'hcount': 2},
                2: {'charge': 0, 'element': 'C', 'hcount': 2},
                3: {'charge': 0, 'element': 'C', 'hcount': 2},
                4: {'charge': 0, 'element': 'C', 'hcount': 2},
                5: {'charge': 0, 'element': 'C', 'hcount': 2},
                6: {'charge': 0, 'element': 'C', 'hcount': 2},
                7: {'charge': 0, 'element': 'C', 'hcount': 2},
                8: {'charge': 0, 'element': 'C', 'hcount': 2},
                9: {'charge': 0, 'element': 'C', 'hcount': 2},
                10: {'charge': 0, 'element': 'C', 'hcount': 2},
                11: {'charge': 0, 'element': 'C', 'hcount': 2},
                12: {'charge': 0, 'element': 'C', 'hcount': 2},
                13: {'charge': 0, 'element': 'C', 'hcount': 2},
                14: {'charge': 0, 'element': 'C', 'hcount': 2},
                15: {'charge': 0, 'element': 'C', 'hcount': 2},
                16: {'charge': 0, 'element': 'C', 'hcount': 2},
                17: {'charge': 0, 'element': 'C', 'hcount': 2},
                18: {'charge': 0, 'element': 'C', 'hcount': 2},
                19: {'charge': 0, 'element': 'C', 'hcount': 2},
                20: {'charge': 0, 'element': 'C', 'hcount': 3},
                21: {'charge': 0, 'element': 'C', 'hcount': 3}}
        nx.set_node_attributes(expected, data)

        self.assertEqualGraphs(found, expected)

    def test_multiring(self):
        smiles = 'N1CC2CCCC2CC1'
        found = read_smiles(smiles, explicit_H=False)

        expected = nx.Graph()
        data = [(0, 1, {'order': 1}),
                (0, 8, {'order': 1}),
                (1, 2, {'order': 1}),
                (2, 3, {'order': 1}),
                (2, 6, {'order': 1}),
                (3, 4, {'order': 1}),
                (4, 5, {'order': 1}),
                (5, 6, {'order': 1}),
                (6, 7, {'order': 1}),
                (7, 8, {'order': 1})]
        expected.add_edges_from(data)
        data = {0: {'charge': 0, 'element': 'N', 'hcount': 1},
                1: {'charge': 0, 'element': 'C', 'hcount': 2},
                2: {'charge': 0, 'element': 'C', 'hcount': 1},
                3: {'charge': 0, 'element': 'C', 'hcount': 2},
                4: {'charge': 0, 'element': 'C', 'hcount': 2},
                5: {'charge': 0, 'element': 'C', 'hcount': 2},
                6: {'charge': 0, 'element': 'C', 'hcount': 1},
                7: {'charge': 0, 'element': 'C', 'hcount': 2},
                8: {'charge': 0, 'element': 'C', 'hcount': 2}}
        nx.set_node_attributes(expected, data)

        self.assertEqualGraphs(found, expected)

    def test_spiro(self):
        smiles = 'C12(CCCCC1)CCCCC2'
        found = read_smiles(smiles, explicit_H=False)

        expected = nx.Graph()
        data = [(0, 1, {'order': 1}),
                (0, 5, {'order': 1}),
                (0, 6, {'order': 1}),
                (0, 10, {'order': 1}),
                (1, 2, {'order': 1}),
                (2, 3, {'order': 1}),
                (3, 4, {'order': 1}),
                (4, 5, {'order': 1}),
                (6, 7, {'order': 1}),
                (7, 8, {'order': 1}),
                (8, 9, {'order': 1}),
                (9, 10, {'order': 1})]
        expected.add_edges_from(data)
        data = {0: {'charge': 0, 'element': 'C', 'hcount': 0},
                1: {'charge': 0, 'element': 'C', 'hcount': 2},
                2: {'charge': 0, 'element': 'C', 'hcount': 2},
                3: {'charge': 0, 'element': 'C', 'hcount': 2},
                4: {'charge': 0, 'element': 'C', 'hcount': 2},
                5: {'charge': 0, 'element': 'C', 'hcount': 2},
                6: {'charge': 0, 'element': 'C', 'hcount': 2},
                7: {'charge': 0, 'element': 'C', 'hcount': 2},
                8: {'charge': 0, 'element': 'C', 'hcount': 2},
                9: {'charge': 0, 'element': 'C', 'hcount': 2},
                10: {'charge': 0, 'element': 'C', 'hcount': 2}}
        nx.set_node_attributes(expected, data)

        self.assertEqualGraphs(found, expected)

    def test_muffle_H(self):
        smiles = '[H]C([H])([H])[H]'
        found = read_smiles(smiles, explicit_H=False)

        expected = nx.Graph()
        data = []
        expected.add_edges_from(data)
        data = {1: {'charge': 0, 'element': 'C', 'hcount': 4}}
        expected.add_nodes_from(data)
        nx.set_node_attributes(expected, data)
        self.assertEqualGraphs(found, expected)

        smiles = '[H]O([H])[H]O([H])[H]'
        found = read_smiles(smiles, explicit_H=False)
        expected = nx.Graph()
        data = [(1, 3, {'order': 1}), (3, 4, {'order': 1})]
        expected.add_edges_from(data)
        data = {1: {'charge': 0, 'element': 'O', 'hcount': 2},
                3: {'charge': 0, 'element': 'H', 'hcount': 0},
                4: {'charge': 0, 'element': 'O', 'hcount': 2}}
        nx.set_node_attributes(expected, data)
        self.assertEqualGraphs(found, expected)

    def test_atom_spec(self):
        # Little bit of everything
        smiles = '[013CH3-1:1]'
        found = read_smiles(smiles, explicit_H=False)

        expected = nx.Graph()
        data = {'charge': -1, 'class': 1, 'element': 'C', 'hcount': 3, 'isotope': 13}
        expected.add_node(0, **data)
        self.assertEqualGraphs(found, expected)

        # Leading 0 isotope; just - for charge
        smiles = '[013CH3-:1]'
        found = read_smiles(smiles, explicit_H=False)
        self.assertEqualGraphs(found, expected)

        # (Deprecated) ++ syntax; 2 letter element
        smiles = '[Cu++]'
        found = read_smiles(smiles, explicit_H=False)

        expected = nx.Graph()
        data = {'charge': 2, 'element': 'Cu', 'hcount': 0}
        expected.add_node(0, **data)
        self.assertEqualGraphs(found, expected)

        smiles = '[Cu+2]'
        found = read_smiles(smiles, explicit_H=False)
        self.assertEqualGraphs(found, expected)

        # 3 letter element
        smiles = '[Uuo+4]'
        found = read_smiles(smiles, explicit_H=False)
        expected = nx.Graph()
        data = {'charge': 4, 'element': 'Uuo', 'hcount': 0}
        expected.add_node(0, **data)
        self.assertEqualGraphs(found, expected)

        smiles = '[2H][CH2]'
        found = read_smiles(smiles, explicit_H=False)

        expected = nx.Graph()
        data = [(0, 1, {'order': 1})]
        expected.add_edges_from(data)
        data = {0: {'charge': 0, 'element': 'H', 'hcount': 0, 'isotope': 2},
                1: {'charge': 0, 'element': 'C', 'hcount': 2}}
        nx.set_node_attributes(expected, data)
        self.assertEqualGraphs(found, expected)

        smiles = '[HH]'
        with self.assertRaises(ValueError):
            read_smiles(smiles, explicit_H=False)

    def test_rhodium(self):
        smiles = '[Rh-](Cl)(Cl)(Cl)(Cl)$[Rh-](Cl)(Cl)(Cl)Cl'
        found = read_smiles(smiles, explicit_H=False)

        expected = nx.Graph()
        data = [(0, 1, {'order': 1}),
                (0, 2, {'order': 1}),
                (0, 3, {'order': 1}),
                (0, 4, {'order': 1}),
                (0, 5, {'order': 4}),
                (5, 6, {'order': 1}),
                (5, 7, {'order': 1}),
                (5, 8, {'order': 1}),
                (5, 9, {'order': 1})]
        expected.add_edges_from(data)
        data = {0: {'charge': -1, 'element': 'Rh', 'hcount': 0},
                1: {'charge': 0, 'element': 'Cl', 'hcount': 0},
                2: {'charge': 0, 'element': 'Cl', 'hcount': 0},
                3: {'charge': 0, 'element': 'Cl', 'hcount': 0},
                4: {'charge': 0, 'element': 'Cl', 'hcount': 0},
                5: {'charge': -1, 'element': 'Rh', 'hcount': 0},
                6: {'charge': 0, 'element': 'Cl', 'hcount': 0},
                7: {'charge': 0, 'element': 'Cl', 'hcount': 0},
                8: {'charge': 0, 'element': 'Cl', 'hcount': 0},
                9: {'charge': 0, 'element': 'Cl', 'hcount': 0}}
        nx.set_node_attributes(expected, data)
        self.assertEqualGraphs(found, expected)

    def test_explicit_H(self):
        smiles = 'c1occc1'
        found = read_smiles(smiles, explicit_H=True)
        expected = nx.Graph()
        data = [(0, 1, {'order': 1.5}),
                (0, 4, {'order': 1.5}),
                (0, 5, {'order': 1}),
                (1, 2, {'order': 1.5}),
                (2, 3, {'order': 1.5}),
                (2, 6, {'order': 1}),
                (3, 4, {'order': 1.5}),
                (3, 7, {'order': 1}),
                (4, 8, {'order': 1})]
        expected.add_edges_from(data)
        data = {0: {'charge': 0, 'element': 'C'},
                1: {'charge': 0, 'element': 'O'},
                2: {'charge': 0, 'element': 'C'},
                3: {'charge': 0, 'element': 'C'},
                4: {'charge': 0, 'element': 'C'},
                5: {'charge': 0, 'element': 'H'},
                6: {'charge': 0, 'element': 'H'},
                7: {'charge': 0, 'element': 'H'},
                8: {'charge': 0, 'element': 'H'}}
        nx.set_node_attributes(expected, data)
        self.assertEqualGraphs(found, expected)

        smiles = '[O-]C(O)CN'
        found = read_smiles(smiles, explicit_H=True)
        expected = nx.Graph()
        data = [(0, 1, {'order': 1}),
                (1, 2, {'order': 1}),
                (1, 3, {'order': 1}),
                (1, 5, {'order': 1}),
                (2, 6, {'order': 1}),
                (3, 4, {'order': 1}),
                (3, 7, {'order': 1}),
                (3, 8, {'order': 1}),
                (4, 9, {'order': 1}),
                (4, 10, {'order': 1})]
        expected.add_edges_from(data)
        data = {0: {'charge': -1, 'element': 'O'},
                1: {'charge': 0, 'element': 'C'},
                2: {'charge': 0, 'element': 'O'},
                3: {'charge': 0, 'element': 'C'},
                4: {'charge': 0, 'element': 'N'},
                5: {'charge': 0, 'element': 'H'},
                6: {'charge': 0, 'element': 'H'},
                7: {'charge': 0, 'element': 'H'},
                8: {'charge': 0, 'element': 'H'},
                9: {'charge': 0, 'element': 'H'},
                10: {'charge': 0, 'element': 'H'}}
        nx.set_node_attributes(expected, data)
        self.assertEqualGraphs(found, expected)

    def test_wrong_cycle(self):
        smiles = 'c1ccccc2'
        with self.assertRaises(KeyError):
            read_smiles(smiles, explicit_H=False)

    def test_wrong_marker(self):
        smiles = "C%1CC1"
        with self.assertRaises(ValueError):
            read_smiles(smiles, explicit_H=False)

    def test_wrong_noncycle(self):
        smiles = 'c1c1CC'
        with self.assertRaises(ValueError):
            read_smiles(smiles, explicit_H=False)

    def test_wrong_noncycle2(self):
        smiles = 'CC11C'
        with self.assertRaises(ValueError):
            read_smiles(smiles, explicit_H=False)

    def test_marker_start(self):
        smiles = '1CCC1'
        with self.assertRaises(ValueError):
            read_smiles(smiles, explicit_H=False)

    def test_wrong_aromatic(self):
        smiles = 'cccccc'
        with self.assertRaises(ValueError):
            read_smiles(smiles, explicit_H=False)

    def test_wrong_ring_bond(self):
        smiles = 'C=1CC-1'
        with self.assertRaises(ValueError):
            read_smiles(smiles, explicit_H=False)


if __name__ == '__main__':
    unittest.main()

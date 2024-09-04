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
from hypothesis import strategies as st
from hypothesis.stateful import (RuleBasedStateMachine, rule, invariant,
                                 initialize, precondition, Bundle)
from hypothesis import note, settings, assume, HealthCheck
from hypothesis_networkx import graph_builder

import networkx as nx

from pysmiles import read_smiles
from pysmiles import write_smiles
from pysmiles.smiles_helper import (
    add_explicit_hydrogens, remove_explicit_hydrogens, correct_aromatic_rings,
    fill_valence,
    increment_bond_orders, kekulize, dekekulize)
from pysmiles.testhelper import assertEqualGraphs

isotope = st.integers(min_value=1)
element = st.sampled_from('C N O S P'.split())
charge = st.integers(min_value=-2, max_value=2)
class_ = st.integers(min_value=1)

node_data = st.fixed_dictionaries({
    'isotope': st.one_of(st.none(), isotope),
    'element': st.one_of(element),
    'charge': charge,  # Charge can not be none for comparison reasons
    'aromatic': st.just(False),
    'class': st.one_of(st.none(), class_)
}).map(lambda d: {k: v for k, v in d.items() if v is not None})
FRAGMENTS = st.sampled_from('C O N P S c1ccccc1 C(=O)[O-]'.split())

@settings(max_examples=500, stateful_step_count=100, deadline=None)
class SMILESTest(RuleBasedStateMachine):
    molecule = Bundle("molecule")

    @initialize(fragment=FRAGMENTS)
    def setup(self, fragment):
        self.mol = read_smiles(fragment)
        note(self.mol.nodes(data=True))

    @rule(data=st.data())
    def add_edge(self, data):
        nodes_with_hydrogen = set(n for n in self.mol if self.mol.nodes[n].get('hcount'))
        assume(nodes_with_hydrogen)
        first_node = data.draw(st.sampled_from(sorted(nodes_with_hydrogen)), label='first_node')
        possibles = sorted(nodes_with_hydrogen - {first_node} - set(self.mol[first_node]))
        if possibles:
            second_node = data.draw(st.one_of(st.sampled_from(possibles), st.none()), label='second_node')
        else:
            second_node = None
        if second_node is None:
            # Make a new fragment, and merge the two.
            new_fragment = data.draw(FRAGMENTS, label='new_fragment')
            new_fragment = read_smiles(new_fragment)
            mapping = dict(zip(range(len(new_fragment)), range(len(self.mol), len(self.mol)+len(new_fragment))))
            nx.relabel_nodes(new_fragment, mapping, copy=False)
            second_node = data.draw(st.sampled_from(sorted(n for n in new_fragment if new_fragment.nodes[n].get('hcount'))), label='second_node')
            self.mol.update(new_fragment)
        max_order = min(self.mol.nodes[first_node].get('hcount', 0), self.mol.nodes[second_node].get('hcount', 0), 4)
        order = data.draw(st.one_of(st.integers(min_value=0, max_value=max_order)), label='order')
        self.mol.add_edge(first_node, second_node, order=order)
        self.mol.nodes[first_node]['hcount'] -= order
        self.mol.nodes[second_node]['hcount'] -= order

    @rule()
    def add_explicit_hydrogens(self):
        add_explicit_hydrogens(self.mol)
        remove_explicit_hydrogens(self.mol)

    @rule(function=st.sampled_from([
        (kekulize, {}),
        (dekekulize, {}),
        (correct_aromatic_rings, {'strict': False}),
        (fill_valence, {}),
        (increment_bond_orders, {}),
    ]))
    def helper(self, function):
        function[0](self.mol, **function[1])

    @invariant()
    def write_read_cycle(self):
        smiles = write_smiles(self.mol)
        note(self.mol.nodes(data=True))
        note(self.mol.edges(data=True))
        note(smiles)
        found = read_smiles(smiles, reinterpret_aromatic=False)
        note(found.nodes(data=True))
        note(found.edges(data=True))
        assertEqualGraphs(self.mol, found)

Tester = SMILESTest.TestCase

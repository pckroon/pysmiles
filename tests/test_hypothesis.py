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
                                 initialize)
from hypothesis import note, settings, assume

import networkx as nx

from pysmiles import read_smiles
from pysmiles import write_smiles
from pysmiles.smiles_helper import (
    add_explicit_hydrogens, remove_explicit_hydrogens, fill_valence,
    correct_aromatic_rings, increment_bond_orders, kekulize, dekekulize)
from pysmiles.testhelper import assertEqualGraphs

isotope = st.integers(min_value=1)
class_ = st.integers(min_value=1)

NODE_DATA = st.fixed_dictionaries({
    'isotope': st.one_of(st.none(), isotope),
    'class': st.one_of(st.none(), class_)
}).map(lambda d: {k: v for k, v in d.items() if v is not None})
FRAGMENTS = st.sampled_from('C O N P S c1ccccc1 C(=O)[O-] C(=O)O *'.split())

@settings(max_examples=500, stateful_step_count=50, deadline=None)
class SMILESTest(RuleBasedStateMachine):
    @initialize(fragment=FRAGMENTS)
    def setup(self, fragment):
        self.mol = read_smiles(fragment)
        note(self.mol.nodes(data=True))

    @staticmethod
    def _get_bonding_nodes(mol):
        return set(n for n in mol if mol.nodes[n].get('hcount') or mol.nodes[n].get('element', '*') == '*')

    @rule(data=st.data())
    def add_edge(self, data):
        nodes_with_hydrogen = self._get_bonding_nodes(self.mol)
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
            # need to have at least one hydrogen or * to add an edge.
            candidates = self._get_bonding_nodes(new_fragment)
            assume(candidates)
            second_node = data.draw(st.sampled_from(sorted(candidates)), label='second_node')
            self.mol.update(new_fragment)
        orders = [self.mol.nodes[first_node].get('hcount', 0), self.mol.nodes[second_node].get('hcount', 0), 4]
        if self.mol.nodes[first_node].get('element', '*') == '*':
            orders[0] = 4
        if self.mol.nodes[second_node].get('element', '*') == '*':
            orders[1] = 4
        max_order = min(orders)
        order = data.draw(st.one_of(st.integers(min_value=0, max_value=max_order)), label='order')
        self.mol.add_edge(first_node, second_node, order=order)
        self.mol.nodes[first_node]['hcount'] -= order
        self.mol.nodes[second_node]['hcount'] -= order

        # *-atoms could get negative hcounts.
        for node in (first_node, second_node):
            if self.mol.nodes[node]['hcount'] < 0:
                self.mol.nodes[node]['hcount'] = 0

    @rule(data=st.data())
    def increment_edge_order(self, data):
        candidates = self._get_bonding_nodes(self.mol)
        options = []
        for n1, n2 in self.mol.edges:
            order = self.mol.edges[(n1, n2)].get('order', 1)
            if order < 4 and order != 1.5 and n1 in candidates and n2 in candidates:
                options.append((n1, n2))
        if not options:
            return
        edge = data.draw(st.sampled_from(options), label='edge')
        self.mol.edges[edge]['order'] = self.mol.edges[edge].get('order', 1) + 1
        self.mol.nodes[edge[0]]['hcount'] -= 1
        self.mol.nodes[edge[1]]['hcount'] -= 1
        for node in edge:
            if self.mol.nodes[node]['hcount'] < 0:
                self.mol.nodes[node]['hcount'] = 0

    @rule(data=st.data())
    def decrement_edge_order(self, data):
        options = []
        for edge in self.mol.edges:
            order = self.mol.edges[edge].get('order', 1)
            if order >= 1 and order != 1.5:  # 0-order is allowed, because why not.
                options.append(edge)
        if not options:
            return
        edge = data.draw(st.sampled_from(options), label='edge')
        self.mol.edges[edge]['order'] = self.mol.edges[edge].get('order', 1) - 1
        self.mol.nodes[edge[0]]['hcount'] += 1
        self.mol.nodes[edge[1]]['hcount'] += 1

    @rule(data=st.data())
    def decorate_node(self, data):
        node_data = st.dictionaries(
            keys=st.sampled_from(list(self.mol.nodes)),
            values=NODE_DATA
        )
        node_data = data.draw(node_data, label='node_data')
        for idx, attrs in node_data.items():
            self.mol.nodes[idx].update(attrs)

    @rule(data=st.data())
    def change_charge(self, data):
        candidates = [(idx, self.mol.nodes[idx]['hcount']) for idx in self._get_bonding_nodes(self.mol)]
        delta = st.fixed_dictionaries({}, optional={idx: st.integers(min_value=-hcount, max_value=hcount) for (idx, hcount) in candidates})
        delta = data.draw(delta, label='charge change')
        for idx, charge in delta.items():
            self.mol.nodes[idx]['charge'] += charge
            del self.mol.nodes[idx]['hcount']
        fill_valence(self.mol)

    @rule()
    def cycle_explicit_hydrogens(self):
        add_explicit_hydrogens(self.mol)
        remove_explicit_hydrogens(self.mol)

    @rule(function=st.sampled_from([
        (kekulize, {}),
        (dekekulize, {}),
        (correct_aromatic_rings, {'strict': True}),
        (fill_valence, {}),
        (increment_bond_orders, {}),
    ]))
    def helper(self, function):
        function[0](self.mol, **function[1])

    @invariant()
    def write_read_cycle(self):
        smiles = write_smiles(self.mol)
        # note(self.mol.nodes(data=True))
        # note(self.mol.edges(data=True))
        note(smiles)
        found = read_smiles(smiles, reinterpret_aromatic=True)
        reference = self.mol.copy()
        correct_aromatic_rings(reference)
        # note(found.nodes(data=True))
        # note(found.edges(data=True))
        assertEqualGraphs(reference, found, excluded_node_attrs=['_pos', '_atom_str'], excluded_edge_attrs=['_pos', '_bond_str'])

Tester = SMILESTest.TestCase

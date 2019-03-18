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
                                 initialize, precondition)
from hypothesis import note, settings
from hypothesis_networkx import graph_builder

from pysmiles import read_smiles
from pysmiles import write_smiles
from pysmiles.smiles_helper import (
    add_explicit_hydrogens, remove_explicit_hydrogens, mark_aromatic_atoms,
    mark_aromatic_edges, fill_valence, correct_aromatic_rings,
    increment_bond_orders, )
from pysmiles.testhelper import assertEqualGraphs


def no_none_value(mapping):
    return {k: v for k, v in mapping.items() if v is not None}


def no_hydrogen_hcount(mapping):
    if mapping.get('element') == 'H' and 'hcount' in mapping:
        mapping['hcount'] = 0
    return mapping


def valid_hydrogen_count(mol):
    for node in mol.nodes:
        H_neighbors = sum(1 for neighbor in mol[node] if mol.nodes[neighbor].get('element') == 'H')
        hcount = mol.nodes[node].get('hcount', 0)
        if H_neighbors + hcount > 9:
            mol.nodes[node]['hcount'] = hcount + H_neighbors - 9


isotope = st.integers(min_value=1)
element = st.sampled_from('C N O S P H'.split())
hcount = st.integers(min_value=0, max_value=9)
charge = st.integers(min_value=-99, max_value=99)
class_ = st.integers(min_value=1)

node_data = st.fixed_dictionaries({
    'isotope': st.one_of(st.none(), isotope),
    'hcount': st.one_of(st.none(), hcount),
    'element': st.one_of(st.none(), element),
    'charge': charge,  # Charge can not be none for comparison reasons
    'aromatic': st.just(False),
    'class': st.one_of(st.none(), class_)
}).map(no_none_value).map(no_hydrogen_hcount)
order = st.sampled_from([1, 2, 3, 0, 4, 1.5])
edge_data = st.fixed_dictionaries({'order': order})


arom_triangle = read_smiles('[*]1[*][*]1')
for edge in arom_triangle.edges:
    arom_triangle.edges[edge]['order'] = 1


class SMILESTest(RuleBasedStateMachine):
    @initialize(mol=graph_builder(node_data=node_data, edge_data=edge_data,
                                  min_nodes=1, max_nodes=7))
    def setup(self, mol):
        valid_hydrogen_count(mol)
        self.mol = mol
        note(self.mol.nodes(data=True))
        note(self.mol.edges(data=True))

    @rule()
    def add_explicit_hydrogens(self):
        add_explicit_hydrogens(self.mol)

    @rule()
    def remove_explicit_hydrogens(self):
        remove_explicit_hydrogens(self.mol)

    # We can't run this halfway through, because aromaticity does not always get
    # encoded in the SMILES string. In particular when * elements are involved.
    # @rule()
    # def correct_aromatic_rings(self):
    #     correct_aromatic_rings(self.mol)

    @rule()
    def fill_valence(self):
        fill_valence(self.mol)

    @rule()
    def increment_bond_orders(self):
        increment_bond_orders(self.mol)

    @precondition(lambda self: hasattr(self, 'mol'))
    @invariant()
    def write_read_cycle(self):
        smiles = write_smiles(self.mol)
        note(self.mol.nodes(data=True))
        note(self.mol.edges(data=True))
        note(smiles)

        # self.mol can exist in a mixed implicit/explicit H style. The reference
        # must be one or the other, since we can't read in mixed mode. We want
        # to be sure we produce the correct answer in both cases though.
        for expl_H in (False, True):
            ref_mol = self.mol.copy()
            defaults = {'charge': 0, 'hcount': 0}
            for node in ref_mol:
                for key, val in defaults.items():
                    if key not in ref_mol.nodes[node]:
                        ref_mol.nodes[node][key] = val
            if expl_H:
                add_explicit_hydrogens(ref_mol)
            else:
                remove_explicit_hydrogens(ref_mol)
            found = read_smiles(smiles, explicit_hydrogen=expl_H,
                                reinterpret_aromatic=False)
            note(found.nodes(data=True))
            note(found.edges(data=True))
            assertEqualGraphs(ref_mol, found)


SMILESTest.TestCase.settings = settings(max_examples=100, stateful_step_count=10)
Tester = SMILESTest.TestCase

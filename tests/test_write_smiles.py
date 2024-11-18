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

import pytest

from pysmiles import write_smiles, read_smiles
from pysmiles.testhelper import assertEqualGraphs, make_mol


@pytest.mark.parametrize('node_data, edge_data, kwargs', (
    (
        [(0, {'element': 'C', 'charge': 0, 'aromatic': False, 'hcount': 3}),
         (1, {'element': 'C', 'charge': 0, 'aromatic': False, 'hcount': 2}),
         (2, {'element': 'C', 'charge': 0, 'aromatic': False, 'hcount': 2}),
         (3, {'element': 'C', 'charge': 0, 'aromatic': False, 'hcount': 3})],
        [(0, 1, {'order': 1}),
         (1, 2, {'order': 1}),
         (2, 3, {'order': 1})],
        dict(explicit_hydrogen=False, zero_order_bonds=True, reinterpret_aromatic=True),
    ),
    (
        [(0, {'element': 'C', 'charge': 0, 'aromatic': False, 'hcount': 3}),
         (1, {'element': 'C', 'charge': 0, 'aromatic': False, 'hcount': 1}),
         (2, {'element': 'C', 'charge': 0, 'aromatic': False, 'hcount': 3}),
         (3, {'element': 'C', 'charge': 0, 'aromatic': False, 'hcount': 3})],
        [(0, 1, {'order': 1}),
         (1, 2, {'order': 1}),
         (1, 3, {'order': 1})],
        dict(explicit_hydrogen=False, zero_order_bonds=True, reinterpret_aromatic=True),
    ),
    (
        [(0, {'element': 'C', 'charge': 0, 'aromatic': False, 'hcount': 2}),
         (1, {'element': 'C', 'charge': 0, 'aromatic': False, 'hcount': 1}),
         (2, {'element': 'O', 'charge': 0, 'aromatic': False, 'hcount': 0}),
         (3, {'element': 'C', 'charge': 0, 'aromatic': False, 'hcount': 1}),
         (4, {'element': 'C', 'charge': 0, 'aromatic': False, 'hcount': 2})],
        [(0, 1, {'order': 2}),
         (1, 2, {'order': 1}),
         (2, 3, {'order': 1}),
         (3, 4, {'order': 2})],
        dict(explicit_hydrogen=False, zero_order_bonds=True, reinterpret_aromatic=True),
    ),
    (
        [(0, {'element': 'C', 'charge': 0, 'aromatic': True, 'hcount': 1}),
         (1, {'element': 'C', 'charge': 0, 'aromatic': True, 'hcount': 1}),
         (2, {'element': 'C', 'charge': 0, 'aromatic': True, 'hcount': 1}),
         (3, {'element': 'C', 'charge': 0, 'aromatic': True, 'hcount': 0}),
         (4, {'element': 'C', 'charge': 0, 'aromatic': True, 'hcount': 1}),
         (5, {'element': 'C', 'charge': 0, 'aromatic': True, 'hcount': 1}),
         (6, {'element': 'O', 'charge': 0, 'aromatic': False, 'hcount': 1})],
        [(0, 1, {'order': 1.5}),
         (1, 2, {'order': 1.5}),
         (2, 3, {'order': 1.5}),
         (3, 4, {'order': 1.5}),
         (4, 5, {'order': 1.5}),
         (5, 0, {'order': 1.5}),
         (3, 6, {'order': 1})],
        dict(explicit_hydrogen=False, zero_order_bonds=True, reinterpret_aromatic=False),
    ),
    (
        [(0, {'element': 'C', 'charge': 0, 'aromatic': False, 'hcount': 1}),
         (1, {'element': 'C', 'charge': 0, 'aromatic': False, 'hcount': 1}),
         (2, {'element': 'C', 'charge': 0, 'aromatic': False, 'hcount': 1}),
         (3, {'element': 'C', 'charge': 0, 'aromatic': False, 'hcount': 0}),
         (4, {'element': 'O', 'charge': 0, 'aromatic': False, 'hcount': 0}),
         (6, {'element': 'O', 'charge': 0, 'aromatic': False, 'hcount': 1})],
        [(0, 1, {'order': 2}),
         (1, 2, {'order': 1}),
         (2, 3, {'order': 2}),
         (3, 4, {'order': 1}),
         (4, 0, {'order': 1}),
         (3, 6, {'order': 1})],
        dict(explicit_hydrogen=False, zero_order_bonds=True, reinterpret_aromatic=True),
    ),
    (
        [(0, {'element': 'C', 'charge': 0, 'aromatic': True, 'hcount': 1}),
         (1, {'element': 'C', 'charge': 0, 'aromatic': True, 'hcount': 1}),
         (2, {'element': 'C', 'charge': 0, 'aromatic': True, 'hcount': 0}),
         (3, {'element': 'C', 'charge': 0, 'aromatic': True, 'hcount': 1}),
         (4, {'element': 'C', 'charge': 0, 'aromatic': True, 'hcount': 1}),
         (5, {'element': 'C', 'charge': 0, 'aromatic': True, 'hcount': 1}),
         (6, {'element': 'C', 'charge': 0, 'aromatic': True, 'hcount': 1}),
         (7, {'element': 'C', 'charge': 0, 'aromatic': True, 'hcount': 0}),
         (8, {'element': 'C', 'charge': 0, 'aromatic': True, 'hcount': 1}),
         (9, {'element': 'C', 'charge': 0, 'aromatic': True, 'hcount': 1}),],
        [(0, 1, {'order': 1.5}),
         (1, 2, {'order': 1.5}),
         (2, 3, {'order': 1.5}),
         (3, 4, {'order': 1.5}),
         (4, 5, {'order': 1.5}),
         (5, 6, {'order': 1.5}),
         (6, 7, {'order': 1.5}),
         (7, 8, {'order': 1.5}),
         (8, 9, {'order': 1.5}),
         (9, 0, {'order': 1.5}),
         (2, 7, {'order': 1.5})],
        dict(explicit_hydrogen=False, zero_order_bonds=True, reinterpret_aromatic=True),
    ),
    (
        [(0, {'element': 'C', 'charge': -1, 'aromatic': False, 'hcount': 3}),
         (1, {'element': 'Cl', 'charge': 0, 'aromatic': False, 'hcount': 2}),
         (2, {'element': 'Tc', 'charge': 0, 'aromatic': False, 'hcount': 5}),
         (3, {'element': 'C', 'charge': 1, 'aromatic': False, 'hcount': 2})],
        [(0, 1, {'order': 1}),
         (0, 2, {'order': 1}),
         (2, 3, {'order': 1})],
        dict(explicit_hydrogen=False, zero_order_bonds=True, reinterpret_aromatic=True),
    ),
    (
        [(0, {'element': 'C', 'charge': -1, 'aromatic': False, 'hcount': 4}),
         (1, {'element': 'C', 'charge': -2, 'aromatic': False, 'hcount': 2}),
         (2, {'element': 'C', 'charge': -2, 'aromatic': False, 'hcount': 4, 'class': 1}),
         (3, {'element': 'C', 'charge': 1, 'aromatic': False, 'hcount': 2})],
        [(0, 1, {'order': 1}),
         (1, 2, {'order': 1}),
         (2, 3, {'order': 1})],
        dict(explicit_hydrogen=False, zero_order_bonds=True, reinterpret_aromatic=True),
    ),
    (
        [(0, {'element': 'C', 'charge': 0, 'aromatic': False, 'hcount': 4}),
         (1, {'element': 'C', 'charge': 0, 'aromatic': False, 'hcount': 4})],
        [],
        dict(explicit_hydrogen=False, zero_order_bonds=False, reinterpret_aromatic=True),
    ),
    (
        [('a', {'element': 'C', 'charge': 0, 'aromatic': False, 'hcount': 4}),
         (1, {'element': 'C', 'charge': 0, 'aromatic': False, 'hcount': 4})],
        [],
        dict(explicit_hydrogen=False, zero_order_bonds=False, reinterpret_aromatic=True),
    ),
    (
      [(0, {'element': 'Se', 'charge': 0, 'aromatic': False, 'hcount': 2})],
      [],
      dict(explicit_hydrogen=False, zero_order_bonds=False, reinterpret_aromatic=True),
    ),
    (
      [(0, {'element': 'As', 'charge': 0, 'aromatic': False, 'hcount': 3})],
      [],
      dict(explicit_hydrogen=False, zero_order_bonds=False, reinterpret_aromatic=True),
    ),
))
def test_write_smiles(node_data, edge_data, kwargs):
    mol = make_mol(node_data, edge_data)
    smiles = write_smiles(mol)
    print(smiles)
    print(kwargs)
    found = read_smiles(smiles, **kwargs)
    assertEqualGraphs(mol, found, excluded_node_attrs=['_pos', '_atom_str'], excluded_edge_attrs=['_pos', '_bond_str'])

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
import networkx as nx

from pysmiles.smiles_helper import (
    fill_valence, remove_explicit_hydrogens,
    add_explicit_hydrogens, correct_aromatic_rings,
    valence, kekulize, dekekulize, increment_bond_orders,
    _bonds, _reorder_cycle, _check_for_ez_conflicts,
    _interpret_cis_trans_tokens,
)
from pysmiles.testhelper import assertEqualGraphs, make_mol


@pytest.mark.parametrize('helper, kwargs, n_data_in, e_data_in, n_data_out, e_data_out', (
    # 1
    (
        add_explicit_hydrogens, {},
        [(0, {'element': 'C'})],
        [],
        [(0, {'element': 'C'})],
        [],
    ),
    # 2
    (
        add_explicit_hydrogens, {},
        [(0, {'element': 'C', 'hcount': 2})],
        [],
        [(0, {'element': 'C'}),
         (1, {'element': 'H', 'aromatic': False, 'charge': 0}),
         (2, {'element': 'H', 'aromatic': False, 'charge': 0})],
        [(0, 1, {'order': 1}),
         (0, 2, {'order': 1})],
    ),
    # 3
    (
        add_explicit_hydrogens, {},
        [(0, {'element': 'C', 'hcount': 2}),
         (1, {'element': 'C', 'hcount': 2})],
        [(0, 1, {'order': 2})],
        [(0, {'element': 'C'}),
         (1, {'element': 'H', 'aromatic': False, 'charge': 0}),
         (2, {'element': 'H', 'aromatic': False, 'charge': 0}),
         (3, {'element': 'C'}),
         (4, {'element': 'H', 'aromatic': False, 'charge': 0}),
         (5, {'element': 'H', 'aromatic': False, 'charge': 0}),],
        [(0, 1, {'order': 1}),
         (0, 2, {'order': 1}),
         (3, 4, {'order': 1}),
         (3, 5, {'order': 1}),
         (0, 3, {'order': 2})],
    ),
    # 4
    (
        remove_explicit_hydrogens, {},
        [(0, {'element': 'C'}),
         (1, {'element': 'H', 'aromatic': False, 'charge': 0}),
         (2, {'element': 'H', 'aromatic': False, 'charge': 0}),
         (3, {'element': 'C'}),
         (4, {'element': 'H', 'aromatic': False, 'charge': 0}),
         (5, {'element': 'H', 'aromatic': False, 'charge': 0}),
         (6, {'element': 'C'}),],
        [(0, 1, {'order': 1}),
         (0, 2, {'order': 1}),
         (3, 4, {'order': 1}),
         (3, 5, {'order': 1}),
         (0, 3, {'order': 2}),
         (3, 6, {'order': 1})],
        [(0, {'element': 'C', 'hcount': 2}),
         (1, {'element': 'C', 'hcount': 2}),
         (2, {'element': 'C', 'hcount': 0})],
        [(0, 1, {'order': 2}),
         (1, 2, {'order': 1})],
    ),
    # 5
    (
        remove_explicit_hydrogens, {},
        [(0, {'element': 'H'}),
         (1, {'element': 'H'}),],
        [(0, 1, {'order': 1})],
        [(0, {'element': 'H', 'hcount': 0}),
         (1, {'element': 'H', 'hcount': 0}),],
        [(0, 1, {'order': 1})],
    ),
    # 6
    (
        remove_explicit_hydrogens, {},
        [(0, {'element': 'C'}),
         (1, {'element': 'H'}),],
        [(0, 1, {'order': 2})],
        [(0, {'element': 'C', 'hcount': 0}),
         (1, {'element': 'H', 'hcount': 0}),],
        [(0, 1, {'order': 2})],
    ),
    (
        remove_explicit_hydrogens, {},
        [(0, {'element': 'H', 'hcount': 0}),
         (1, {'element': 'C', 'hcount': 0, 'rs_isomer': (0, 2, 3, 4)}),
         (2, {'element': 'Br', 'hcount': 0}),
         (3, {'element': 'F', 'hcount': 0}),
         (4, {'element': 'C', 'hcount': 3})],
        [(0, 1, {'order': 1}),
         (1, 2, {'order': 1}),
         (1, 3, {'order': 1}),
         (1, 4, {'order': 1}),],
        [(1, {'element': 'C', 'hcount': 1, 'rs_isomer': (1, 2, 3, 4)}),
         (2, {'element': 'Br', 'hcount': 0}),
         (3, {'element': 'F', 'hcount': 0}),
         (4, {'element': 'C', 'hcount': 3})],
        [(1, 2, {'order': 1}),
         (1, 3, {'order': 1}),
         (1, 4, {'order': 1}),],
    ),
    (
        remove_explicit_hydrogens, {},
        [(0, {'element': 'H', 'ez_isomer': [0, 2, 3, 4, 'trans']}),
         (1, {'element': 'F', 'hcount': 0}),
         (2, {'element': 'C'}),
         (3, {'element': 'C', 'hcount': 1}),
         (4, {'element': 'F', 'ez_isomer': [4, 3, 2, 0, 'trans'], 'hcount': 0}),],
        [(0, 2, {'order': 1}),
         (1, 2, {'order': 1}),
         (2, 3, {'order': 2}),
         (3, 4, {'order': 1}),],
        [(1, {'element': 'F', 'hcount': 0, 'ez_isomer': [1, 2, 3, 4, 'cis']}),
         (2, {'element': 'C', 'hcount': 1}),
         (3, {'element': 'C', 'hcount': 1}),
         (4, {'element': 'F', 'ez_isomer': [4, 3, 2, 1, 'cis'], 'hcount': 0}),],
        [(1, 2, {'order': 1}),
         (2, 3, {'order': 2}),
         (3, 4, {'order': 1}),],
    ),
    (
        remove_explicit_hydrogens, {},
        [(0, {'element': 'H', 'ez_isomer': [0, 2, 3, 4, 'trans']}),
         (1, {'element': 'H', 'hcount': 0}),
         (2, {'element': 'C'}),
         (3, {'element': 'C', 'hcount': 1}),
         (4, {'element': 'F', 'ez_isomer': [4, 3, 2, 0, 'trans'], 'hcount': 0}),],
        [(0, 2, {'order': 1}),
         (1, 2, {'order': 1}),
         (2, 3, {'order': 2}),
         (3, 4, {'order': 1}),],
        [(1, {'element': 'H', 'hcount': 0, 'ez_isomer': [1, 2, 3, 4, 'cis']}),
         (2, {'element': 'C', 'hcount': 1}),
         (3, {'element': 'C', 'hcount': 1}),
         (4, {'element': 'F', 'ez_isomer': [4, 3, 2, 1, 'cis'], 'hcount': 0}),],
        [(1, 2, {'order': 1}),
         (2, 3, {'order': 2}),
         (3, 4, {'order': 1}),],
    ),
    # 6
    (
        fill_valence,
        {'respect_hcount': True, 'respect_bond_order': True, 'max_bond_order': 3},
        [(0, {'element': 'C'}),
         (1, {'element': 'C'})],
        [(0, 1, {'order': 1})],
        [(0, {'element': 'C', 'hcount': 3}),
         (1, {'element': 'C', 'hcount': 3})],
        [(0, 1, {'order': 1})],
    ),
    # 6
    (
        fill_valence,
        {'respect_hcount': True, 'respect_bond_order': True, 'max_bond_order': 3},
        [(0, {'element': 'C', 'hcount': 1}),
         (1, {'element': 'C'})],
        [(0, 1, {'order': 1})],
        [(0, {'element': 'C', 'hcount': 1}),
         (1, {'element': 'C', 'hcount': 3})],
        [(0, 1, {'order': 1})],
    ),
    # 7
    (
        fill_valence,
        {'respect_hcount': False, 'respect_bond_order': True, 'max_bond_order': 3},
        [(0, {'element': 'C', 'hcount': 1}),
         (1, {'element': 'C'})],
        [(0, 1, {'order': 1})],
        [(0, {'element': 'C', 'hcount': 3}),
         (1, {'element': 'C', 'hcount': 3})],
        [(0, 1, {'order': 1})],
    ),
    # 8
    (
        fill_valence,
        {'respect_hcount': True, 'respect_bond_order': False, 'max_bond_order': 3},
        [(0, {'element': 'C'}),
         (1, {'element': 'C'})],
        [(0, 1, {'order': 1})],
        [(0, {'element': 'C', 'hcount': 1}),
         (1, {'element': 'C', 'hcount': 1})],
        [(0, 1, {'order': 3})],
    ),
    # 9
    (
        # This case sort of stinks, since there's a single aromatic bond not in
        # a cycle.
        fill_valence,
        {'respect_hcount': True, 'respect_bond_order': False, 'max_bond_order': 3},
        [(0, {'element': 'C'}),
         (1, {'element': 'C'})],
        [(0, 1, {'order': 1.5})],
        [(0, {'element': 'C', 'hcount': 2}),
         (1, {'element': 'C', 'hcount': 2})],
        [(0, 1, {'order': 1.5})],
    ),
    # 10
    (
        fill_valence,
        {'respect_hcount': False, 'respect_bond_order': True, 'max_bond_order': 3},
        [(0, {'element': 'Tc'}),
         (1, {'element': 'C', 'hcount': 5})],
        [(0, 1, {'order': 1})],
        [(0, {'element': 'Tc', 'hcount': 0}),
         (1, {'element': 'C', 'hcount': 5})],
        [(0, 1, {'order': 1})],
    ),
    # 11
    (
        correct_aromatic_rings, {},
        [(0, {'element': 'C', 'hcount': 1, 'charge': 0}),
         (1, {'element': 'C', 'hcount': 1, 'charge': 0}),
         (2, {'element': 'C', 'hcount': 1, 'charge': 0}),
         (3, {'element': 'C', 'hcount': 0, 'charge': 0}),
         (4, {'element': 'C', 'hcount': 3, 'charge': 0})],
        [(0, 1, {'order': 2}),
         (1, 2, {'order': 1}),
         (2, 3, {'order': 2}),
         (3, 0, {'order': 1}),
         (3, 4, {'order': 1})],
        [(0, {'element': 'C', 'hcount': 1, 'charge': 0, 'aromatic': True}),
         (1, {'element': 'C', 'hcount': 1, 'charge': 0, 'aromatic': True}),
         (2, {'element': 'C', 'hcount': 1, 'charge': 0, 'aromatic': True}),
         (3, {'element': 'C', 'hcount': 0, 'charge': 0, 'aromatic': True}),
         (4, {'element': 'C', 'hcount': 3, 'charge': 0, 'aromatic': False})],
        [(0, 1, {'order': 1.5}),
         (1, 2, {'order': 1.5}),
         (2, 3, {'order': 1.5}),
         (3, 0, {'order': 1.5}),
         (3, 4, {'order': 1})],
    ),
    # 12
    (
        correct_aromatic_rings, {},
        [(0, {'element': 'C', 'hcount': 2, 'charge': 0}),
         (1, {'element': 'C', 'hcount': 2, 'charge': 0}),
         (2, {'element': 'C', 'hcount': 2, 'charge': 0}),
         (3, {'element': 'C', 'hcount': 1, 'charge': 0}),
         (4, {'element': 'C', 'hcount': 3, 'charge': 0})],
        [(0, 1, {'order': 1}),
         (1, 2, {'order': 1}),
         (2, 3, {'order': 1}),
         (3, 0, {'order': 1}),
         (3, 4, {'order': 1})],
        [(0, {'element': 'C', 'hcount': 2, 'charge': 0, 'aromatic': False}),
         (1, {'element': 'C', 'hcount': 2, 'charge': 0, 'aromatic': False}),
         (2, {'element': 'C', 'hcount': 2, 'charge': 0, 'aromatic': False}),
         (3, {'element': 'C', 'hcount': 1, 'charge': 0, 'aromatic': False}),
         (4, {'element': 'C', 'hcount': 3, 'charge': 0, 'aromatic': False})],
        [(0, 1, {'order': 1}),
         (1, 2, {'order': 1}),
         (2, 3, {'order': 1}),
         (3, 0, {'order': 1}),
         (3, 4, {'order': 1})],
    ),
    # 13
    (
        correct_aromatic_rings, {},
        [(0, {'element': 'C', 'hcount': 1, 'charge': 0}),
         (1, {'element': 'C', 'hcount': 1, 'charge': 0}),
         (2, {'element': 'C', 'hcount': 1, 'charge': 0}),
         (3, {'hcount': 1, 'charge': 0}),],
        [(0, 1, {'order': 2}),
         (1, 2, {'order': 1}),
         (2, 3, {'order': 2}),
         (3, 0, {'order': 1}),],
        [(0, {'element': 'C', 'hcount': 1, 'charge': 0, 'aromatic': True}),
         (1, {'element': 'C', 'hcount': 1, 'charge': 0, 'aromatic': True}),
         (2, {'element': 'C', 'hcount': 1, 'charge': 0, 'aromatic': True}),
         (3, {'hcount': 1, 'charge': 0, 'aromatic': True}),],
        [(0, 1, {'order': 1.5}),
         (1, 2, {'order': 1.5}),
         (2, 3, {'order': 1.5}),
         (3, 0, {'order': 1.5}),],
    ),
    # 16
    (
        correct_aromatic_rings, {},
        [(0, {'charge': 1, 'aromatic': False}),
         (1, {'charge': 0, 'aromatic': True}),
         (2, {'charge': 0, 'aromatic': True}),],
        [(0, 1, {'order': 1}),
         (1, 2, {'order': 1}),
         (2, 0, {'order': 1}),],
        [(0, {'charge': 1, 'aromatic': False}),
         (1, {'charge': 0, 'aromatic': False}),
         (2, {'charge': 0, 'aromatic': False}),],
        [(0, 1, {'order': 1}),
         (1, 2, {'order': 1}),
         (2, 0, {'order': 1}),],
    ),
    # 18
    (
        correct_aromatic_rings, {},
        [(0, {'element': 'C'}),
         (1, {'element': 'C'}),
         (2, {'element': 'C'}),
         (3, {'element': 'C'}),],
        [(0, 1, {}),
         (1, 2, {}),
         (2, 3, {}),
         (3, 0, {})],
        [(0, {'element': 'C', 'aromatic': False}),
         (1, {'element': 'C', 'aromatic': False}),
         (2, {'element': 'C', 'aromatic': False}),
         (3, {'element': 'C', 'aromatic': False}),],
        [(0, 1, {}),
         (1, 2, {}),
         (2, 3, {}),
         (3, 0, {})],
    ),
    # 19
    (
        correct_aromatic_rings, {},
        [(0, {'element': 'C', 'hcount': 1}),
         (1, {'element': 'C', 'hcount': 1}),
         (2, {'element': 'C', 'hcount': 1}),
         (3, {'element': 'C', 'hcount': 1}),],
        [(0, 1, {}),
         (1, 2, {'order': 2}),
         (2, 3, {}),
         (3, 0, {'order': 2}),],
        [(0, {'element': 'C', 'hcount': 1, 'aromatic': True}),
         (1, {'element': 'C', 'hcount': 1, 'aromatic': True}),
         (2, {'element': 'C', 'hcount': 1, 'aromatic': True}),
         (3, {'element': 'C', 'hcount': 1, 'aromatic': True}),],
        [(0, 1, {'order': 1.5}),
         (1, 2, {'order': 1.5}),
         (2, 3, {'order': 1.5}),
         (3, 0, {'order': 1.5})],
    ),
(
        correct_aromatic_rings, {},
        [(0, {'element': 'C', 'hcount': 1, 'aromatic': True}),
         (1, {'element': 'C', 'hcount': 1, 'aromatic': True}),
         (2, {'element': 'C', 'hcount': 1, 'aromatic': True}),
         (3, {'element': 'C', 'hcount': 1, 'aromatic': True}),],
        [(0, 1, {}),
         (1, 2, {}),
         (2, 3, {}),
         (3, 0, {}),],
        [(0, {'element': 'C', 'hcount': 1, 'aromatic': True}),
         (1, {'element': 'C', 'hcount': 1, 'aromatic': True}),
         (2, {'element': 'C', 'hcount': 1, 'aromatic': True}),
         (3, {'element': 'C', 'hcount': 1, 'aromatic': True}),],
        [(0, 1, {'order': 1.5}),
         (1, 2, {'order': 1.5}),
         (2, 3, {'order': 1.5}),
         (3, 0, {'order': 1.5})],
    ),
    # 20
    (
        correct_aromatic_rings, {},
        [(0, {'element': 'C', 'hcount': 1}),
         (1, {'element': 'C', 'hcount': 1}),
         (2, {'element': 'C', 'hcount': 1}),
         (3, {'element': 'C', 'hcount': 1}),],
        [(0, 1, {}),
         (1, 2, {}),
         (2, 3, {}),],
        [(0, {'element': 'C', 'hcount': 1, 'aromatic': False}),
         (1, {'element': 'C', 'hcount': 1, 'aromatic': False}),
         (2, {'element': 'C', 'hcount': 1, 'aromatic': False}),
         (3, {'element': 'C', 'hcount': 1, 'aromatic': False}),],
        [(0, 1, {}),
         (1, 2, {}),
         (2, 3, {}),],
    ),
    # 21
    (
        correct_aromatic_rings, {},
        [(0, {'element': 'C', 'hcount': 1, 'aromatic': True}),
         (1, {'element': 'C', 'hcount': 1, 'aromatic': True}),
         (2, {'element': 'C', 'hcount': 1, 'aromatic': True}),
         (3, {'element': 'C', 'hcount': 1, 'aromatic': True}),
         (4, {'element': 'O', 'hcount': 0, 'aromatic': True}),],
        [(0, 1, {}),
         (1, 2, {}),
         (2, 3, {}),
         (3, 4, {}),
         (4, 0, {})],
        [(0, {'element': 'C', 'hcount': 1, 'aromatic': False}),
         (1, {'element': 'C', 'hcount': 1, 'aromatic': False}),
         (2, {'element': 'C', 'hcount': 1, 'aromatic': False}),
         (3, {'element': 'C', 'hcount': 1, 'aromatic': False}),
         (4, {'element': 'O', 'hcount': 0, 'aromatic': False}),],
        [(0, 1, {'order': 2}),
         (1, 2, {}),
         (2, 3, {'order': 2}),
         (3, 4, {}),
         (4, 0, {})],
    ),
    (
        correct_aromatic_rings, {},
        [(0, {'element': 'C', 'hcount': 1, 'aromatic': True}),
         (1, {'element': 'C', 'hcount': 1, 'aromatic': True}),
         (2, {'element': 'C', 'hcount': 1, 'aromatic': True}),
         (3, {'element': 'C', 'hcount': 1, 'aromatic': True}),
         (4, {'element': 'N', 'hcount': 1, 'aromatic': True}),],
        [(0, 1, {}),
         (1, 2, {}),
         (2, 3, {}),
         (3, 4, {}),
         (4, 0, {}),],
        [(0, {'element': 'C', 'hcount': 1, 'aromatic': False}),
         (1, {'element': 'C', 'hcount': 1, 'aromatic': False}),
         (2, {'element': 'C', 'hcount': 1, 'aromatic': False}),
         (3, {'element': 'C', 'hcount': 1, 'aromatic': False}),
         (4, {'element': 'N', 'hcount': 1, 'aromatic': False}),],
        [(0, 1, {'order': 2}),
         (1, 2, {}),
         (2, 3, {'order': 2}),
         (3, 4, {}),
         (4, 0, {}),],
    ),
    (
        correct_aromatic_rings, {},
        [(0, {'element': 'C', 'hcount': 1, 'aromatic': True}),
         (1, {'element': 'C', 'hcount': 1, 'aromatic': True}),
         (2, {'element': 'C', 'hcount': 1, 'aromatic': True}),
         (3, {'element': 'C', 'hcount': 1, 'aromatic': True}),
         (4, {'element': 'N', 'aromatic': True}),
         (5, {'element': 'H', 'aromatic': False})],
        [(0, 1, {}),
         (1, 2, {}),
         (2, 3, {}),
         (3, 4, {}),
         (4, 0, {}),
         (4, 5, {})],
        [(0, {'element': 'C', 'hcount': 1, 'aromatic': False}),
         (1, {'element': 'C', 'hcount': 1, 'aromatic': False}),
         (2, {'element': 'C', 'hcount': 1, 'aromatic': False}),
         (3, {'element': 'C', 'hcount': 1, 'aromatic': False}),
         (4, {'element': 'N', 'aromatic': False}),
         (5, {'element': 'H', 'aromatic': False})],
        [(0, 1, {'order': 2}),
         (1, 2, {}),
         (2, 3, {'order': 2}),
         (3, 4, {}),
         (4, 0, {}),
         (4, 5, {})],
    ),
    (
        correct_aromatic_rings, {},
        [(0, {'element': 'C', 'hcount': 1, 'aromatic': True}),
         (1, {'element': 'C', 'hcount': 1, 'aromatic': True}),
         (2, {'element': 'C', 'hcount': 1, 'aromatic': True}),
         (3, {'element': '*', 'hcount': 1, 'aromatic': True}),],
        [(0, 1, {'order': 1}),
         (1, 2, {'order': 1}),
         (2, 3, {'order': 1}),
         (0, 3, {'order': 1}),],
        [(0, {'element': 'C', 'hcount': 1, 'aromatic': True}),
         (1, {'element': 'C', 'hcount': 1, 'aromatic': True}),
         (2, {'element': 'C', 'hcount': 1, 'aromatic': True}),
         (3, {'element': '*', 'hcount': 1, 'aromatic': True}),],
        [(0, 1, {'order': 1.5}),
         (1, 2, {'order': 1.5}),
         (2, 3, {'order': 1.5}),
         (0, 3, {'order': 1.5}),],
    ),
    (
        kekulize, {},
        [(0, {'element': 'C', 'hcount': 1, 'aromatic': True}),
         (1, {'element': 'C', 'hcount': 1, 'aromatic': True}),
         (2, {'element': 'C', 'hcount': 1, 'aromatic': True}),
         (3, {'element': 'C', 'hcount': 1, 'aromatic': True}),],
        [(0, 1, {'order': 1.5}),
         (1, 2, {'order': 1.5}),
         (2, 3, {'order': 1.5}),
         (0, 3, {'order': 1.5}),],
        [(0, {'element': 'C', 'hcount': 1, 'aromatic': False}),
         (1, {'element': 'C', 'hcount': 1, 'aromatic': False}),
         (2, {'element': 'C', 'hcount': 1, 'aromatic': False}),
         (3, {'element': 'C', 'hcount': 1, 'aromatic': False}), ],
        [(0, 1, {'order': 1}),
         (1, 2, {'order': 2}),
         (2, 3, {'order': 1}),
         (0, 3, {'order': 2}), ],
    ),
    (
        dekekulize, {},
        [(0, {'element': 'C', 'hcount': 1, 'aromatic': False}),
         (1, {'element': 'C', 'hcount': 1, 'aromatic': False}),
         (2, {'element': 'C', 'hcount': 1, 'aromatic': False}),
         (3, {'element': 'C', 'hcount': 1, 'aromatic': False}), ],
        [(0, 1, {'order': 1}),
         (1, 2, {'order': 2}),
         (2, 3, {'order': 1}),
         (0, 3, {'order': 2}), ],
        [(0, {'element': 'C', 'hcount': 1, 'aromatic': True}),
         (1, {'element': 'C', 'hcount': 1, 'aromatic': True}),
         (2, {'element': 'C', 'hcount': 1, 'aromatic': True}),
         (3, {'element': 'C', 'hcount': 1, 'aromatic': True}),],
        [(0, 1, {'order': 1.5}),
         (1, 2, {'order': 1.5}),
         (2, 3, {'order': 1.5}),
         (0, 3, {'order': 1.5}),],
    ),
))
def test_helper(helper, kwargs, n_data_in, e_data_in, n_data_out, e_data_out):
    print(f'{helper=} {kwargs=}')
    print(f'{n_data_in=} {e_data_in=}')
    print(f'{n_data_out=} {e_data_out=}')
    mol = make_mol(n_data_in, e_data_in)
    helper(mol, **kwargs)
    ref_mol = make_mol(n_data_out, e_data_out)
    print(mol.nodes(data=True))
    print(mol.edges(data=True))
    assertEqualGraphs(mol, ref_mol)


@pytest.mark.parametrize('atom, expected', [
    ({'element': 'C'}, [4]),
    ({'element': 'N'}, [3, 5]),
    ({'element': 'S'}, [2, 4, 6]),
    ({'element': 'P'}, [3, 5]),
    ({'element': 'N', 'charge': 1}, [4]),
    ({'element': 'B', 'charge': 1}, [2]),
])
def test_valence(atom, expected):
    found = valence(atom)
    assert found == expected


@pytest.mark.parametrize('atom', [
    {"element": "H", "charge": 5},
    {"element": "H", "charge": -10000}
])
def test_valence_error(atom):
    with pytest.raises(ValueError):
        valence(atom)


def test_bonds_without_order_counts_neighbors():
    mol = make_mol(
        [(0, {'element': 'C'}), (1, {'element': 'N'})],
        [(0, 1, {'order': 3})],
    )

    assert _bonds(mol, 0, use_order=False) == 1


def test_correct_aromatic_rings_warns_for_unmatched_aromatic_region(caplog):
    mol = make_mol(
        [(0, {'element': 'C', 'hcount': 1, 'aromatic': True}),
         (1, {'element': 'C', 'hcount': 1, 'aromatic': True}),
         (2, {'element': 'C', 'hcount': 1, 'aromatic': True})],
        [(0, 1, {'order': 1.5}),
         (1, 2, {'order': 1.5}),
         (2, 0, {'order': 1.5})],
    )

    with caplog.at_level('WARNING'):
        correct_aromatic_rings(mol, strict=False)

    assert 'Your molecule is invalid and cannot be kekulized.' in caplog.text
    assert sorted(edge['order'] for *_, edge in mol.edges(data=True)) == [1, 1, 2]


def test_kekulize_rejects_non_kekulizable_aromatic_region():
    mol = make_mol(
        [(0, {'element': 'C', 'hcount': 1, 'aromatic': True}),
         (1, {'element': 'C', 'hcount': 1, 'aromatic': True}),
         (2, {'element': 'C', 'hcount': 1, 'aromatic': True})],
        [(0, 1, {'order': 1.5}),
         (1, 2, {'order': 1.5}),
         (2, 0, {'order': 1.5})],
    )

    with pytest.raises(ValueError, match='Aromatic region cannot be kekulized.'):
        kekulize(mol)


def test_reorder_cycle_returns_empty_for_no_nodes():
    assert _reorder_cycle(nx.Graph(), []) == []


def test_reorder_cycle_raises_for_non_cycle():
    graph = nx.path_graph([0, 1, 2])
    graph.add_node(3)

    with pytest.raises(nx.NetworkXNoCycle, match='Not a cycle'):
        _reorder_cycle(graph, [0, 1, 2, 3])


def test_increment_bond_orders_skips_edges_once_a_node_is_satisfied(caplog):
    mol = make_mol(
        [(0, {'element': 'C'}),
         (1, {'element': 'C'}),
         (2, {'element': 'C'})],
        [(0, 1, {'order': 1}),
         (1, 2, {'order': 1})],
    )

    with caplog.at_level('WARNING'):
        increment_bond_orders(mol)

    assert sorted(edge['order'] for *_, edge in mol.edges(data=True)) == [1, 3]
    assert 'Unable to completely assign bond orders from valence.' in caplog.text


@pytest.mark.parametrize('anchor, tagged_nodes, ez_isomer_class', [
    (5, (1, 2), {1: '/', 2: '/'}),
    (2, (1, 3), {1: '/', 3: '\\'}),
])
def test_check_for_ez_conflicts_raises(anchor, tagged_nodes, ez_isomer_class):
    with pytest.raises(ValueError, match='Conflicting cis/trans assignment'):
        _check_for_ez_conflicts(anchor, tagged_nodes, ez_isomer_class)


def test_check_for_ez_conflicts_allows_distinct_same_side_tokens():
    _check_for_ez_conflicts(5, (1, 2), {1: '/', 2: '\\'})


def test_interpret_cis_trans_tokens_marks_case_three_as_cis():
    mol = make_mol([(0, {}), (1, {}), (2, {}), (3, {})], [])

    _interpret_cis_trans_tokens(mol, [[[0, 1, '/'], [3, 2, '\\']]])

    assert mol.nodes[0]['ez_isomer'] == [(0, 1, 2, 3, 'cis')]
    assert mol.nodes[3]['ez_isomer'] == [(3, 2, 1, 0, 'cis')]

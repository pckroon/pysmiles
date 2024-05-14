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

from pysmiles.smiles_helper import (
    correct_aromatic_rings, fill_valence, remove_explicit_hydrogens,
    add_explicit_hydrogens, mark_aromatic_atoms, mark_aromatic_edges,
    valence
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
        [(0, {'element': 'Tc', 'hcount': 6}),
         (1, {'element': 'C', 'hcount': 5})],
        [(0, 1, {'order': 1})],
    ),
    # 11
    (
        mark_aromatic_atoms, {},
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
        [(0, 1, {'order': 2}),
         (1, 2, {'order': 1}),
         (2, 3, {'order': 2}),
         (3, 0, {'order': 1}),
         (3, 4, {'order': 1})],
    ),
    # 12
    (
        mark_aromatic_atoms, {},
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
        mark_aromatic_atoms, {},
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
        [(0, 1, {'order': 2}),
         (1, 2, {'order': 1}),
         (2, 3, {'order': 2}),
         (3, 0, {'order': 1}),],
    ),
    # 14
    (
        mark_aromatic_edges, {},
        [(0, {'charge': 1, 'aromatic': True}),
         (1, {'charge': 0, 'aromatic': True}),
         (2, {'charge': 0, 'aromatic': True}),],
        [(0, 1, {'order': 1}),
         (1, 2, {'order': 1}),
         (2, 0, {'order': 1}),],
        [(0, {'charge': 1, 'aromatic': True}),
         (1, {'charge': 0, 'aromatic': True}),
         (2, {'charge': 0, 'aromatic': True}),],
        [(0, 1, {'order': 1.5}),
         (1, 2, {'order': 1.5}),
         (2, 0, {'order': 1.5}),],
    ),
    # 15
    (
        mark_aromatic_edges, {},
        [(0, {'charge': 1, 'aromatic': True}),
         (1, {'charge': 0, 'aromatic': True}),
         (2, {'charge': 0, 'aromatic': True}),],
        [(0, 1, {'order': 1}),
         (1, 2, {'order': 1}),],
        [(0, {'charge': 1, 'aromatic': True}),
         (1, {'charge': 0, 'aromatic': True}),
         (2, {'charge': 0, 'aromatic': True}),],
        [(0, 1, {'order': 1.5}),
         (1, 2, {'order': 1.5}),],
    ),
    # 16
    (
        # This case smells a bit. Not all atoms in a cycle are aromatic, so only
        # some of the bonds become aromatic.
        mark_aromatic_edges, {},
        [(0, {'charge': 1, 'aromatic': False}),
         (1, {'charge': 0, 'aromatic': True}),
         (2, {'charge': 0, 'aromatic': True}),],
        [(0, 1, {'order': 1}),
         (1, 2, {'order': 1}),
         (2, 0, {'order': 1}),],
        [(0, {'charge': 1, 'aromatic': False}),
         (1, {'charge': 0, 'aromatic': True}),
         (2, {'charge': 0, 'aromatic': True}),],
        [(0, 1, {'order': 1}),
         (1, 2, {'order': 1.5}),
         (2, 0, {'order': 1}),],
    ),
    # 17
    (
        mark_aromatic_edges, {},
        [(0, {'charge': 1, 'aromatic': True}),
         (1, {'charge': 0, 'aromatic': True}),
         (2, {'charge': 0, 'aromatic': True}),
         (3, {'aromatic': False})],
        [(0, 1, {'order': 1}),
         (1, 2, {'order': 1}),
         (2, 0, {'order': 1}),
         (2, 3, {'order': 1})],
        [(0, {'charge': 1, 'aromatic': True}),
         (1, {'charge': 0, 'aromatic': True}),
         (2, {'charge': 0, 'aromatic': True}),
         (3, {'aromatic': False})],
        [(0, 1, {'order': 1.5}),
         (1, 2, {'order': 1.5}),
         (2, 0, {'order': 1.5}),
         (2, 3, {'order': 1})],
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
        [(0, {'element': 'C', 'hcount': 2, 'aromatic': False}),
         (1, {'element': 'C', 'hcount': 2, 'aromatic': False}),
         (2, {'element': 'C', 'hcount': 2, 'aromatic': False}),
         (3, {'element': 'C', 'hcount': 2, 'aromatic': False}),],
        [(0, 1, {'order': 1}),
         (1, 2, {'order': 1}),
         (2, 3, {'order': 1}),
         (3, 0, {'order': 1})],
    ),
    # 19
    (
        correct_aromatic_rings, {},
        [(0, {'element': 'C', 'hcount': 1}),
         (1, {'element': 'C', 'hcount': 1}),
         (2, {'element': 'C', 'hcount': 1}),
         (3, {'element': 'C', 'hcount': 1}),],
        [(0, 1, {}),
         (1, 2, {}),
         (2, 3, {}),
         (3, 0, {})],
        [(0, {'element': 'C', 'hcount': 1, 'aromatic': True}),
         (1, {'element': 'C', 'hcount': 1, 'aromatic': True}),
         (2, {'element': 'C', 'hcount': 1, 'aromatic': True}),
         (3, {'element': 'C', 'hcount': 1, 'aromatic': True}),],
        [(0, 1, {'order': 1.5}),
         (1, 2, {'order': 1.5}),
         (2, 3, {'order': 1.5}),
         (3, 0, {'order': 1.5})],
    ),
    # 20 - this should lead to bond-orders of three ...
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
        [(0, 1, {'order': 2}),
         (1, 2, {'order': 1}),
         (2, 3, {'order': 2}),],
    ),
    # 21
    (
        correct_aromatic_rings, {},
        [(0, {'element': 'C', 'hcount': 1}),
         (1, {'element': 'C', 'hcount': 1}),
         (2, {'element': 'C', 'hcount': 1}),
         (3, {'element': 'C', 'hcount': 1}),
         (4, {'element': 'O', 'hcount': 0}),],
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
         (1, 2, {'order': 1}),
         (2, 3, {'order': 2}),
         (3, 4, {'order': 1}),
         (4, 0, {'order': 1})],
    ),
    (
        correct_aromatic_rings, {},
        [(0, {'element': 'C', 'hcount': 1}),
         (1, {'element': 'C', 'hcount': 1}),
         (2, {'element': 'C', 'hcount': 1}),
         (3, {'element': 'C', 'hcount': 1}),
         (4, {'element': 'N', 'hcount': 1}),],
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
         (1, 2, {'order': 1}),
         (2, 3, {'order': 2}),
         (3, 4, {'order': 1}),
         (4, 0, {'order': 1}),],
    ),
    (
        correct_aromatic_rings, {},
        [(0, {'element': 'C', 'hcount': 1}),
         (1, {'element': 'C', 'hcount': 1}),
         (2, {'element': 'C', 'hcount': 1}),
         (3, {'element': 'C', 'hcount': 1}),
         (4, {'element': 'N'}),
         (5, {'element': 'H'})],
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
         (4, {'element': 'N', 'hcount': 0, 'aromatic': False}),
         (5, {'element': 'H', 'aromatic': False})],
        [(0, 1, {'order': 2}),
         (1, 2, {'order': 1}),
         (2, 3, {'order': 2}),
         (3, 4, {'order': 1}),
         (4, 0, {'order': 1}),
         (4, 5, {'order': 1})],
    ),
))
def test_helper(helper, kwargs, n_data_in, e_data_in, n_data_out, e_data_out):
    mol = make_mol(n_data_in, e_data_in)
    helper(mol, **kwargs)
    ref_mol = make_mol(n_data_out, e_data_out)
    assertEqualGraphs(mol, ref_mol)


@pytest.mark.parametrize('atom, expected', [
    ({'element': 'C'}, [4]),
    ({'element': 'N'}, [3]),
    ({'element': 'S'}, [2, 4, 6]),
    ({'element': 'P'}, [3, 5]),
    ({'element': 'N', 'charge': 1}, [4]),
    ({'element': 'B', 'charge': 1}, [2]),
])
def test_valence(atom, expected):
    found = valence(atom)
    assert found == expected


@pytest.mark.parametrize('atom', [
    {"charge": 1},
    {"charge": -10000}
])
def test_valence_error(atom):
    with pytest.raises(ValueError):
        valence(atom)

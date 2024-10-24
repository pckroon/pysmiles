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
from pysmiles import read_smiles
from pysmiles.testhelper import assertEqualGraphs, make_mol


@pytest.mark.parametrize('smiles, node_data, edge_data, explicit_h', (
    (
        'CCCC',
        [(0, {'element': 'C', 'charge': 0, 'hcount': 3, 'aromatic': False}),
         (1, {'element': 'C', 'charge': 0, 'hcount': 2, 'aromatic': False}),
         (2, {'element': 'C', 'charge': 0, 'hcount': 2, 'aromatic': False}),
         (3, {'element': 'C', 'charge': 0, 'hcount': 3, 'aromatic': False})],
        [(0, 1, {'order': 1}),
         (1, 2, {'order': 1}),
         (2, 3, {'order': 1})],
        False
    ),
    (
        'CCC(CC)CC',
        [(0, {'element': 'C', 'charge': 0, 'hcount': 3, 'aromatic': False}),
         (1, {'element': 'C', 'charge': 0, 'hcount': 2, 'aromatic': False}),
         (2, {'element': 'C', 'charge': 0, 'hcount': 1, 'aromatic': False}),
         (3, {'element': 'C', 'charge': 0, 'hcount': 2, 'aromatic': False}),
         (4, {'element': 'C', 'charge': 0, 'hcount': 3, 'aromatic': False}),
         (5, {'element': 'C', 'charge': 0, 'hcount': 2, 'aromatic': False}),
         (6, {'element': 'C', 'charge': 0, 'hcount': 3, 'aromatic': False})],
        [(0, 1, {'order': 1}),
         (1, 2, {'order': 1}),
         (2, 3, {'order': 1}),
         (2, 5, {'order': 1}),
         (3, 4, {'order': 1}),
         (5, 6, {'order': 1})],
        False
    ),
    (
        'C=C',
        [(0, {'element': 'C', 'charge': 0, 'hcount': 2, 'aromatic': False}),
         (1, {'element': 'C', 'charge': 0, 'hcount': 2, 'aromatic': False})],
        [(0, 1, {'order': 2})],
        False
    ),
    (
        'C1CC1',
        [(0, {'charge': 0, 'element': 'C', 'hcount': 2, 'aromatic': False}),
         (1, {'charge': 0, 'element': 'C', 'hcount': 2, 'aromatic': False}),
         (2, {'charge': 0, 'element': 'C', 'hcount': 2, 'aromatic': False})],
        [(0, 1, {'order': 1}),
         (0, 2, {'order': 1}),
         (1, 2, {'order': 1})],
        False
    ),
    (
        'C%11CC%11',
        [(0, {'charge': 0, 'element': 'C', 'hcount': 2, 'aromatic': False}),
         (1, {'charge': 0, 'element': 'C', 'hcount': 2, 'aromatic': False}),
         (2, {'charge': 0, 'element': 'C', 'hcount': 2, 'aromatic': False})],
        [(0, 1, {'order': 1}),
         (0, 2, {'order': 1}),
         (1, 2, {'order': 1})],
        False
    ),
    (
        'c1ccccc1',
        [(0, {'charge': 0, 'element': 'C', 'hcount': 1, 'aromatic': True}),
         (1, {'charge': 0, 'element': 'C', 'hcount': 1, 'aromatic': True}),
         (2, {'charge': 0, 'element': 'C', 'hcount': 1, 'aromatic': True}),
         (3, {'charge': 0, 'element': 'C', 'hcount': 1, 'aromatic': True}),
         (4, {'charge': 0, 'element': 'C', 'hcount': 1, 'aromatic': True}),
         (5, {'charge': 0, 'element': 'C', 'hcount': 1, 'aromatic': True})],
        [(0, 1, {'order': 1.5}),
         (0, 5, {'order': 1.5}),
         (1, 2, {'order': 1.5}),
         (2, 3, {'order': 1.5}),
         (3, 4, {'order': 1.5}),
         (4, 5, {'order': 1.5})],
        False
    ),
    (
        'C1CC=1',
        [(0, {'charge': 0, 'element': 'C', 'hcount': 1, 'aromatic': False}),
         (1, {'charge': 0, 'element': 'C', 'hcount': 2, 'aromatic': False}),
         (2, {'charge': 0, 'element': 'C', 'hcount': 1, 'aromatic': False})],
        [(0, 1, {'order': 1}),
         (0, 2, {'order': 2}),
         (1, 2, {'order': 1})],
        False
    ),
    (
        'C=1CC=1',
        [(0, {'charge': 0, 'element': 'C', 'hcount': 1, 'aromatic': False}),
         (1, {'charge': 0, 'element': 'C', 'hcount': 2, 'aromatic': False}),
         (2, {'charge': 0, 'element': 'C', 'hcount': 1, 'aromatic': False})],
        [(0, 1, {'order': 1}),
         (0, 2, {'order': 2}),
         (1, 2, {'order': 1})],
        False
    ),
    (
        'C=1CC1',
        [(0, {'charge': 0, 'element': 'C', 'hcount': 1, 'aromatic': False}),
         (1, {'charge': 0, 'element': 'C', 'hcount': 2, 'aromatic': False}),
         (2, {'charge': 0, 'element': 'C', 'hcount': 1, 'aromatic': False})],
        [(0, 1, {'order': 1}),
         (0, 2, {'order': 2}),
         (1, 2, {'order': 1})],
        False
    ),
    (
        'C%11CC=%11',
        [(0, {'charge': 0, 'element': 'C', 'hcount': 1, 'aromatic': False}),
         (1, {'charge': 0, 'element': 'C', 'hcount': 2, 'aromatic': False}),
         (2, {'charge': 0, 'element': 'C', 'hcount': 1, 'aromatic': False})],
        [(0, 1, {'order': 1}),
         (0, 2, {'order': 2}),
         (1, 2, {'order': 1})],
        False
    ),
    (
        'OS(=O)(=S)O',
        [(0, {'charge': 0, 'element': 'O', 'hcount': 1, 'aromatic': False}),
         (1, {'charge': 0, 'element': 'S', 'hcount': 0, 'aromatic': False}),
         (2, {'charge': 0, 'element': 'O', 'hcount': 0, 'aromatic': False}),
         (3, {'charge': 0, 'element': 'S', 'hcount': 0, 'aromatic': False}),
         (4, {'charge': 0, 'element': 'O', 'hcount': 1, 'aromatic': False})],
        [(0, 1, {'order': 1}),
         (1, 2, {'order': 2}),
         (1, 3, {'order': 2}),
         (1, 4, {'order': 1})],
        False
    ),
    (
        'C(C(C(C(C(C(C(C(C(C(C(C(C(C(C(C(C(C(C(C(C))))))))))))))))))))C',
        [(0, {'charge': 0, 'element': 'C', 'hcount': 2, 'aromatic': False}),
         (1, {'charge': 0, 'element': 'C', 'hcount': 2, 'aromatic': False}),
         (2, {'charge': 0, 'element': 'C', 'hcount': 2, 'aromatic': False}),
         (3, {'charge': 0, 'element': 'C', 'hcount': 2, 'aromatic': False}),
         (4, {'charge': 0, 'element': 'C', 'hcount': 2, 'aromatic': False}),
         (5, {'charge': 0, 'element': 'C', 'hcount': 2, 'aromatic': False}),
         (6, {'charge': 0, 'element': 'C', 'hcount': 2, 'aromatic': False}),
         (7, {'charge': 0, 'element': 'C', 'hcount': 2, 'aromatic': False}),
         (8, {'charge': 0, 'element': 'C', 'hcount': 2, 'aromatic': False}),
         (9, {'charge': 0, 'element': 'C', 'hcount': 2, 'aromatic': False}),
         (10, {'charge': 0, 'element': 'C', 'hcount': 2, 'aromatic': False}),
         (11, {'charge': 0, 'element': 'C', 'hcount': 2, 'aromatic': False}),
         (12, {'charge': 0, 'element': 'C', 'hcount': 2, 'aromatic': False}),
         (13, {'charge': 0, 'element': 'C', 'hcount': 2, 'aromatic': False}),
         (14, {'charge': 0, 'element': 'C', 'hcount': 2, 'aromatic': False}),
         (15, {'charge': 0, 'element': 'C', 'hcount': 2, 'aromatic': False}),
         (16, {'charge': 0, 'element': 'C', 'hcount': 2, 'aromatic': False}),
         (17, {'charge': 0, 'element': 'C', 'hcount': 2, 'aromatic': False}),
         (18, {'charge': 0, 'element': 'C', 'hcount': 2, 'aromatic': False}),
         (19, {'charge': 0, 'element': 'C', 'hcount': 2, 'aromatic': False}),
         (20, {'charge': 0, 'element': 'C', 'hcount': 3, 'aromatic': False}),
         (21, {'charge': 0, 'element': 'C', 'hcount': 3, 'aromatic': False})],
        [(0, 1, {'order': 1}),
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
         (19, 20, {'order': 1})],
        False
    ),
    (
        'N1CC2CCCC2CC1',
        [(0, {'charge': 0, 'element': 'N', 'hcount': 1, 'aromatic': False}),
         (1, {'charge': 0, 'element': 'C', 'hcount': 2, 'aromatic': False}),
         (2, {'charge': 0, 'element': 'C', 'hcount': 1, 'aromatic': False}),
         (3, {'charge': 0, 'element': 'C', 'hcount': 2, 'aromatic': False}),
         (4, {'charge': 0, 'element': 'C', 'hcount': 2, 'aromatic': False}),
         (5, {'charge': 0, 'element': 'C', 'hcount': 2, 'aromatic': False}),
         (6, {'charge': 0, 'element': 'C', 'hcount': 1, 'aromatic': False}),
         (7, {'charge': 0, 'element': 'C', 'hcount': 2, 'aromatic': False}),
         (8, {'charge': 0, 'element': 'C', 'hcount': 2, 'aromatic': False})],
        [(0, 1, {'order': 1}),
         (0, 8, {'order': 1}),
         (1, 2, {'order': 1}),
         (2, 3, {'order': 1}),
         (2, 6, {'order': 1}),
         (3, 4, {'order': 1}),
         (4, 5, {'order': 1}),
         (5, 6, {'order': 1}),
         (6, 7, {'order': 1}),
         (7, 8, {'order': 1})],
        False
    ),
    (
        'c1ccc2ccccc2c1',
        [(0, {'charge': 0, 'element': 'C', 'aromatic': True}),
         (1, {'charge': 0, 'element': 'C', 'aromatic': True}),
         (2, {'charge': 0, 'element': 'C', 'aromatic': True}),
         (3, {'charge': 0, 'element': 'C', 'aromatic': True}),
         (4, {'charge': 0, 'element': 'C', 'aromatic': True}),
         (5, {'charge': 0, 'element': 'C', 'aromatic': True}),
         (6, {'charge': 0, 'element': 'C', 'aromatic': True}),
         (7, {'charge': 0, 'element': 'C', 'aromatic': True}),
         (8, {'charge': 0, 'element': 'C', 'aromatic': True}),
         (9, {'charge': 0, 'element': 'C', 'aromatic': True}),
         (10, {'charge': 0, 'element': 'H', 'aromatic': False}),
         (11, {'charge': 0, 'element': 'H', 'aromatic': False}),
         (12, {'charge': 0, 'element': 'H', 'aromatic': False}),
         (13, {'charge': 0, 'element': 'H', 'aromatic': False}),
         (14, {'charge': 0, 'element': 'H', 'aromatic': False}),
         (15, {'charge': 0, 'element': 'H', 'aromatic': False}),
         (16, {'charge': 0, 'element': 'H', 'aromatic': False}),
         (17, {'charge': 0, 'element': 'H', 'aromatic': False})],
        [(0, 1, {'order': 1.5}),
         (0, 9, {'order': 1.5}),
         (0, 10, {'order': 1}),
         (1, 2, {'order': 1.5}),
         (1, 11, {'order': 1}),
         (2, 3, {'order': 1.5}),
         (2, 12, {'order': 1}),
         (3, 4, {'order': 1.5}),
         (3, 8, {'order': 1.5}),
         (4, 5, {'order': 1.5}),
         (4, 13, {'order': 1}),
         (5, 6, {'order': 1.5}),
         (5, 14, {'order': 1}),
         (6, 7, {'order': 1.5}),
         (6, 15, {'order': 1}),
         (7, 8, {'order': 1.5}),
         (7, 16, {'order': 1}),
         (8, 9, {'order': 1.5}),
         (9, 17, {'order': 1})],
        True
    ),
    (
        'C12(CCCCC1)CCCCC2',
        [(0, {'charge': 0, 'element': 'C', 'hcount': 0, 'aromatic': False}),
         (1, {'charge': 0, 'element': 'C', 'hcount': 2, 'aromatic': False}),
         (2, {'charge': 0, 'element': 'C', 'hcount': 2, 'aromatic': False}),
         (3, {'charge': 0, 'element': 'C', 'hcount': 2, 'aromatic': False}),
         (4, {'charge': 0, 'element': 'C', 'hcount': 2, 'aromatic': False}),
         (5, {'charge': 0, 'element': 'C', 'hcount': 2, 'aromatic': False}),
         (6, {'charge': 0, 'element': 'C', 'hcount': 2, 'aromatic': False}),
         (7, {'charge': 0, 'element': 'C', 'hcount': 2, 'aromatic': False}),
         (8, {'charge': 0, 'element': 'C', 'hcount': 2, 'aromatic': False}),
         (9, {'charge': 0, 'element': 'C', 'hcount': 2, 'aromatic': False}),
         (10, {'charge': 0, 'element': 'C', 'hcount': 2, 'aromatic': False})],
        [(0, 1, {'order': 1}),
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
         (9, 10, {'order': 1})],
        False
    ),
    (
        '[H]C([H])([H])[H]',
        [(0, {'charge': 0, 'element': 'C', 'hcount': 4, 'aromatic': False})],
        [],
        False
    ),
    (
        '[H][O+]([H])[H].O([H])[H]',
        [(1, {'charge': 1, 'element': 'O', 'hcount': 2, 'aromatic': False}),
         (3, {'charge': 0, 'element': 'H', 'hcount': 0, 'aromatic': False}),
         (4, {'charge': 0, 'element': 'O', 'hcount': 2, 'aromatic': False})],
        [(1, 3, {'order': 1}),
         (3, 4, {'order': 0})],
        False
    ),
    (
        '[O]=C',
        [(0, {'charge': 0, 'element': 'O', 'hcount': 0, 'aromatic': False}),
         (1, {'charge': 0, 'element': 'C', 'hcount': 2, 'aromatic': False}),],
        [(0, 1, {'order': 2})],
        False
    ),
    (
        '[H][H]',
        [(0, {'charge': 0, 'element': 'H', 'hcount': 0, 'aromatic': False}),
         (1, {'charge': 0, 'element': 'H', 'hcount': 0, 'aromatic': False}),],
        [(0, 1, {'order': 1})],
        False
    ),
    (
        '[13CH3-1:1]',
        [(0, {'charge': -1, 'class': 1, 'element': 'C', 'hcount': 3, 'isotope': 13, 'aromatic': False})],
        [],
        False
    ),
    (
        '[013CH3-:1]',
        [(0, {'charge': -1, 'class': 1, 'element': 'C', 'hcount': 3, 'isotope': 13, 'aromatic': False})],
        [],
        False
    ),
    (
        'Cl[Cu++]Cl',
        [
            (0, {'charge': 0, 'element': 'Cl', 'hcount': 0, 'aromatic': False}),
            (1, {'charge': 2, 'element': 'Cu', 'hcount': 0, 'aromatic': False}),
            (2, {'charge': 0, 'element': 'Cl', 'hcount': 0, 'aromatic': False}),
        ],
        [(0, 1, {'order': 1}), (1, 2, {'order': 1})],
        False
    ),
    (
        '[Cu+2](Cl)Cl',
        [
            (0, {'charge': 2, 'element': 'Cu', 'hcount': 0, 'aromatic': False}),
            (1, {'charge': 0, 'element': 'Cl', 'hcount': 0, 'aromatic': False}),
            (2, {'charge': 0, 'element': 'Cl', 'hcount': 0, 'aromatic': False}),
        ],
        [(0, 1, {'order': 1}), (0, 2, {'order': 1})],
        False
    ),
    (
        '[OgH4+4]',
        [(0, {'charge': 4, 'element': 'Og', 'hcount': 4, 'aromatic': False})],
        [],
        False
    ),
    (
        '[2H][CH2+]',
        [(0, {'charge': 0, 'element': 'H', 'hcount': 0, 'isotope': 2, 'aromatic': False}),
         (1, {'charge': 1, 'element': 'C', 'hcount': 2, 'aromatic': False})],
        [(0, 1, {'order': 1})],
        False
    ),
    (
        '*',
        [(0, {'charge': 0, 'aromatic': False, 'hcount': 0})],
        [],
        False
    ),
    (
        '[*--]',
        [(0, {'charge': -2, 'aromatic': False, 'hcount': 0})],
        [],
        False
    ),
    (
        '[*-]',
        [(0, {'charge': -1, 'aromatic': False, 'hcount': 0})],
        [],
        False
    ),
    (
        '[H+]',
        [(0, {'element': 'H', 'charge': 1, 'hcount': 0, 'aromatic': False})],
        [],
        False
    ),
    (
        '[cH-1]1[cH-1][cH-1]1',
        [(0, {'element': 'C', 'aromatic': False, 'charge': -1, 'hcount': 1}),
         (1, {'element': 'C', 'aromatic': False, 'charge': -1, 'hcount': 1}),
         (2, {'element': 'C', 'aromatic': False, 'charge': -1, 'hcount': 1})],
        [(0, 1, {'order': 1}),
         (1, 2, {'order': 1}),
         (2, 0, {'order': 1})],
        False
    ),
    (
        'c1cscc1',
        [(0, {'charge': 0, 'element': 'C', 'aromatic': False, 'hcount': 1}),
         (1, {'charge': 0, 'element': 'C', 'aromatic': False, 'hcount': 1}),
         (2, {'charge': 0, 'element': 'S', 'aromatic': False, 'hcount': 0}),
         (3, {'charge': 0, 'element': 'C', 'aromatic': False, 'hcount': 1}),
         (4, {'charge': 0, 'element': 'C', 'aromatic': False, 'hcount': 1}),],
        [(0, 1, {'order': 2}),
         (0, 4, {'order': 1}),
         (1, 2, {'order': 1}),
         (2, 3, {'order': 1}),
         (3, 4, {'order': 2}),],
        False
    ),
    (
        'c12ccccc1[nH]cc2',
        [(0, {'charge': 0, 'element': 'C', 'aromatic': True, 'hcount': 0}),
         (1, {'charge': 0, 'element': 'C', 'aromatic': True, 'hcount': 1}),
         (2, {'charge': 0, 'element': 'C', 'aromatic': True, 'hcount': 1}),
         (3, {'charge': 0, 'element': 'C', 'aromatic': True, 'hcount': 1}),
         (4, {'charge': 0, 'element': 'C', 'aromatic': True, 'hcount': 1}),
         (5, {'charge': 0, 'element': 'C', 'aromatic': True, 'hcount': 0}),
         (6, {'charge': 0, 'element': 'N', 'aromatic': False, 'hcount': 1}),
         (7, {'charge': 0, 'element': 'C', 'aromatic': False, 'hcount': 1}),
         (8, {'charge': 0, 'element': 'C', 'aromatic': False, 'hcount': 1}),],
        [(0, 1, {'order': 1.5}),
         (1, 2, {'order': 1.5}),
         (2, 3, {'order': 1.5}),
         (3, 4, {'order': 1.5}),
         (4, 5, {'order': 1.5}),
         (5, 0, {'order': 1.5}),
         (5, 6, {'order': 1}),
         (6, 7, {'order': 1}),
         (7, 8, {'order': 2}),
         (8, 0, {'order': 1}),],
        False
    ),
    (
        'c1cc2ccc3c2c1cc3',
        [(0, {'charge': 0, 'element': 'C', 'aromatic': True, 'hcount': 1}),
         (1, {'charge': 0, 'element': 'C', 'aromatic': True, 'hcount': 1}),
         (2, {'charge': 0, 'element': 'C', 'aromatic': True, 'hcount': 0}),
         (3, {'charge': 0, 'element': 'C', 'aromatic': True, 'hcount': 1}),
         (4, {'charge': 0, 'element': 'C', 'aromatic': True, 'hcount': 1}),
         (5, {'charge': 0, 'element': 'C', 'aromatic': True, 'hcount': 0}),
         (6, {'charge': 0, 'element': 'C', 'aromatic': True, 'hcount': 0}),
         (7, {'charge': 0, 'element': 'C', 'aromatic': True, 'hcount': 0}),
         (8, {'charge': 0, 'element': 'C', 'aromatic': True, 'hcount': 1}),
         (9, {'charge': 0, 'element': 'C', 'aromatic': True, 'hcount': 1}),],
        [(0, 1, {'order': 1.5}),
         (0, 7, {'order': 1.5}),
         (1, 2, {'order': 1.5}),
         (2, 3, {'order': 1.5}),
         (3, 4, {'order': 1.5}),
         (4, 5, {'order': 1.5}),
         (6, 2, {'order': 1.5}),
         (5, 6, {'order': 1.5}),
         (5, 9, {'order': 1.5}),
         (6, 7, {'order': 1.5}),
         (6, 2, {'order': 1.5}),
         (9, 8, {'order': 1.5}),
         (7, 8, {'order': 1.5}),],
        False
    ),
    (
        '[Rh-](Cl)(Cl)(Cl)(Cl)$[Rh-](Cl)(Cl)(Cl)Cl',
        [(0, {'charge': -1, 'element': 'Rh', 'hcount': 0, 'aromatic': False}),
         (1, {'charge': 0, 'element': 'Cl', 'hcount': 0, 'aromatic': False}),
         (2, {'charge': 0, 'element': 'Cl', 'hcount': 0, 'aromatic': False}),
         (3, {'charge': 0, 'element': 'Cl', 'hcount': 0, 'aromatic': False}),
         (4, {'charge': 0, 'element': 'Cl', 'hcount': 0, 'aromatic': False}),
         (5, {'charge': -1, 'element': 'Rh', 'hcount': 0, 'aromatic': False}),
         (6, {'charge': 0, 'element': 'Cl', 'hcount': 0, 'aromatic': False}),
         (7, {'charge': 0, 'element': 'Cl', 'hcount': 0, 'aromatic': False}),
         (8, {'charge': 0, 'element': 'Cl', 'hcount': 0, 'aromatic': False}),
         (9, {'charge': 0, 'element': 'Cl', 'hcount': 0, 'aromatic': False})],
        [(0, 1, {'order': 1}),
         (0, 2, {'order': 1}),
         (0, 3, {'order': 1}),
         (0, 4, {'order': 1}),
         (0, 5, {'order': 4}),
         (5, 6, {'order': 1}),
         (5, 7, {'order': 1}),
         (5, 8, {'order': 1}),
         (5, 9, {'order': 1})],
        False
    ),
    (
        'c1occc1',
        [(0, {'charge': 0, 'element': 'C', 'aromatic': False}),
         (1, {'charge': 0, 'element': 'O', 'aromatic': False}),
         (2, {'charge': 0, 'element': 'C', 'aromatic': False}),
         (3, {'charge': 0, 'element': 'C', 'aromatic': False}),
         (4, {'charge': 0, 'element': 'C', 'aromatic': False}),
         (5, {'charge': 0, 'element': 'H', 'aromatic': False}),
         (6, {'charge': 0, 'element': 'H', 'aromatic': False}),
         (7, {'charge': 0, 'element': 'H', 'aromatic': False}),
         (8, {'charge': 0, 'element': 'H', 'aromatic': False})],
        [(0, 1, {'order': 1}),
         (0, 4, {'order': 2}),
         (0, 5, {'order': 1}),
         (1, 2, {'order': 1}),
         (2, 3, {'order': 2}),
         (2, 6, {'order': 1}),
         (3, 4, {'order': 1}),
         (3, 7, {'order': 1}),
         (4, 8, {'order': 1})],
        True
    ),
    (
        'c1[asH]ccc1',
        [(0, {'charge': 0, 'element': 'C', 'aromatic': False}),
         (1, {'charge': 0, 'element': 'As', 'aromatic': False}),
         (2, {'charge': 0, 'element': 'C', 'aromatic': False}),
         (3, {'charge': 0, 'element': 'C', 'aromatic': False}),
         (4, {'charge': 0, 'element': 'C', 'aromatic': False}),
         (5, {'charge': 0, 'element': 'H', 'aromatic': False}),
         (6, {'charge': 0, 'element': 'H', 'aromatic': False}),
         (7, {'charge': 0, 'element': 'H', 'aromatic': False}),
         (8, {'charge': 0, 'element': 'H', 'aromatic': False}),
         (9, {'charge': 0, 'element': 'H', 'aromatic': False})],
        [(0, 1, {'order': 1}),
         (0, 4, {'order': 2}),
         (0, 5, {'order': 1}),
         (1, 2, {'order': 1}),
         (2, 3, {'order': 2}),
         (2, 6, {'order': 1}),
         (3, 4, {'order': 1}),
         (3, 7, {'order': 1}),
         (4, 8, {'order': 1}),
         (1, 9, {'order': 1}),],
        True
    ),
    (
        'c1[se]ccc1',
        [(0, {'charge': 0, 'element': 'C', 'aromatic': False}),
         (1, {'charge': 0, 'element': 'Se', 'aromatic': False}),
         (2, {'charge': 0, 'element': 'C', 'aromatic': False}),
         (3, {'charge': 0, 'element': 'C', 'aromatic': False}),
         (4, {'charge': 0, 'element': 'C', 'aromatic': False}),
         (5, {'charge': 0, 'element': 'H', 'aromatic': False}),
         (6, {'charge': 0, 'element': 'H', 'aromatic': False}),
         (7, {'charge': 0, 'element': 'H', 'aromatic': False}),
         (8, {'charge': 0, 'element': 'H', 'aromatic': False})],
        [(0, 1, {'order': 1}),
         (0, 4, {'order': 2}),
         (0, 5, {'order': 1}),
         (1, 2, {'order': 1}),
         (2, 3, {'order': 2}),
         (2, 6, {'order': 1}),
         (3, 4, {'order': 1}),
         (3, 7, {'order': 1}),
         (4, 8, {'order': 1})],
        True
    ),
    (
        '[O-]C(O)CN',
        [(0, {'charge': -1, 'element': 'O', 'aromatic': False}),
         (1, {'charge': 0, 'element': 'C', 'aromatic': False}),
         (2, {'charge': 0, 'element': 'O', 'aromatic': False}),
         (3, {'charge': 0, 'element': 'C', 'aromatic': False}),
         (4, {'charge': 0, 'element': 'N', 'aromatic': False}),
         (5, {'charge': 0, 'element': 'H', 'aromatic': False}),
         (6, {'charge': 0, 'element': 'H', 'aromatic': False}),
         (7, {'charge': 0, 'element': 'H', 'aromatic': False}),
         (8, {'charge': 0, 'element': 'H', 'aromatic': False}),
         (9, {'charge': 0, 'element': 'H', 'aromatic': False}),
         (10, {'charge': 0, 'element': 'H', 'aromatic': False})],
        [(0, 1, {'order': 1}),
         (1, 2, {'order': 1}),
         (1, 3, {'order': 1}),
         (1, 5, {'order': 1}),
         (2, 6, {'order': 1}),
         (3, 4, {'order': 1}),
         (3, 7, {'order': 1}),
         (3, 8, {'order': 1}),
         (4, 9, {'order': 1}),
         (4, 10, {'order': 1})],
        True
    ),
    (
        '[*+]1[*][*]1',
        [(0, {'charge': 1, 'aromatic': False, 'hcount': 0}),
         (1, {'charge': 0, 'aromatic': False, 'hcount': 0}),
         (2, {'charge': 0, 'aromatic': False, 'hcount': 0})],
        [(0, 1, {'order': 1}),
         (1, 2, {'order': 1}),
         (2, 0, {'order': 1}),],
        False
    ),
    (
        'N1[*][*]1',
        [(0, {'element': 'N', 'charge': 0, 'aromatic': False, 'hcount': 1}),
         (1, {'charge': 0, 'aromatic': False, 'hcount': 0}),
         (2, {'charge': 0, 'aromatic': False, 'hcount': 0})],
        [(0, 1, {'order': 1}),
         (1, 2, {'order': 1}),
         (2, 0, {'order': 1}),],
        False
    ),
    (
        'c12ccccc1=CC=CC=2',
        [(0, {'element': 'C', 'charge': 0, 'aromatic': True, 'hcount': 0}),
         (1, {'element': 'C', 'charge': 0, 'aromatic': True, 'hcount': 1}),
         (2, {'element': 'C', 'charge': 0, 'aromatic': True, 'hcount': 1}),
         (3, {'element': 'C', 'charge': 0, 'aromatic': True, 'hcount': 1}),
         (4, {'element': 'C', 'charge': 0, 'aromatic': True, 'hcount': 1}),
         (5, {'element': 'C', 'charge': 0, 'aromatic': True, 'hcount': 0}),
         (6, {'element': 'C', 'charge': 0, 'aromatic': True, 'hcount': 1}),
         (7, {'element': 'C', 'charge': 0, 'aromatic': True, 'hcount': 1}),
         (8, {'element': 'C', 'charge': 0, 'aromatic': True, 'hcount': 1}),
         (9, {'element': 'C', 'charge': 0, 'aromatic': True, 'hcount': 1})],
        [(0, 1, {'order': 1.5}),
         (0, 5, {'order': 1.5}),
         (0, 9, {'order': 1.5}),
         (1, 2, {'order': 1.5}),
         (2, 3, {'order': 1.5}),
         (3, 4, {'order': 1.5}),
         (4, 5, {'order': 1.5}),
         (5, 6, {'order': 1.5}),
         (6, 7, {'order': 1.5}),
         (7, 8, {'order': 1.5}),
         (8, 9, {'order': 1.5})],
        False
    ),
    # chiral center S/L alanine
    (
        'C[C@@H](C(=O)O)N',
        [(0, {'element': 'C', 'charge': 0, 'aromatic': False, 'hcount': 3}),
         (1, {'charge': 0, 'aromatic': False, 'element': 'C', 'rs_isomer': (0, 1, 5, 2), 'hcount': 1}),
         (2, {'element': 'C', 'charge': 0, 'aromatic': False, 'hcount': 0}),
         (3, {'element': 'O', 'charge': 0, 'aromatic': False, 'hcount': 0}),
         (4, {'element': 'O', 'charge': 0, 'aromatic': False, 'hcount': 1}),
         (5, {'element': 'N', 'charge': 0, 'aromatic': False, 'hcount': 2}),],
        [(0, 1, {'order': 1}),
         (1, 2, {'order': 1}),
         (1, 5, {'order': 1}),
         (2, 3, {'order': 2}),
         (2, 4, {'order': 1}),],
        False
    ),
    # chiral center R/D alanine
    (
        'C[C@H](C(=O)O)N',
        [(0, {'element': 'C', 'charge': 0, 'aromatic': False, 'hcount': 3}),
         (1, {'charge': 0, 'aromatic': False, 'element': 'C', 'rs_isomer': (0, 1, 2, 5), 'hcount': 1}),
         (2, {'element': 'C', 'charge': 0, 'aromatic': False, 'hcount': 0}),
         (3, {'element': 'O', 'charge': 0, 'aromatic': False, 'hcount': 0}),
         (4, {'element': 'O', 'charge': 0, 'aromatic': False, 'hcount': 1}),
         (5, {'element': 'N', 'charge': 0, 'aromatic': False, 'hcount': 2}), ],
        [(0, 1, {'order': 1}),
         (1, 2, {'order': 1}),
         (1, 5, {'order': 1}),
         (2, 3, {'order': 2}),
         (2, 4, {'order': 1}), ],
        False
    ),
    (  # R/D ALA, to see if explicit hydrogens work
        'C[C@H](C(=O)O)N',
        [(0, {'element': 'C', 'charge': 0, 'aromatic': False}),
         (1, {'charge': 0, 'aromatic': False, 'element': 'C', 'rs_isomer': (0, 9, 2, 5)}),
         (2, {'element': 'C', 'charge': 0, 'aromatic': False}),
         (3, {'element': 'O', 'charge': 0, 'aromatic': False}),
         (4, {'element': 'O', 'charge': 0, 'aromatic': False}),
         (5, {'element': 'N', 'charge': 0, 'aromatic': False}),
         (6, {'charge': 0, 'aromatic': False, 'element': 'H'}),
         (7, {'charge': 0, 'aromatic': False, 'element': 'H'}),
         (8, {'charge': 0, 'aromatic': False, 'element': 'H'}),
         (9, {'charge': 0, 'aromatic': False, 'element': 'H'}),
         (10, {'charge': 0, 'aromatic': False, 'element': 'H'}),
         (11, {'charge': 0, 'aromatic': False, 'element': 'H'}),
         (12, {'charge': 0, 'aromatic': False, 'element': 'H'})],
        [(0, 1, {'order': 1}),
         (0, 6, {'order': 1}),
         (0, 7, {'order': 1}),
         (0, 8, {'order': 1}),
         (1, 2, {'order': 1}),
         (1, 5, {'order': 1}),
         (1, 9, {'order': 1}),
         (2, 3, {'order': 2}),
         (2, 4, {'order': 1}),
         (4, 10, {'order': 1}),
         (5, 11, {'order': 1}),
         (5, 12, {'order': 1})],
        True
    ),
    # test with smiles and ring bond
    ('[C@]1(Br)(Cl)CCCC(F)C1',
        [(0, {'element': 'C', 'charge': 0, 'aromatic': False, 'rs_isomer': (8, 1, 2, 3), 'hcount': 0}),
         (1, {'element': 'Br', 'charge': 0, 'aromatic': False, 'hcount': 0}),
         (2, {'element': 'Cl', 'charge': 0, 'aromatic': False, 'hcount': 0}),
         (3, {'element': 'C', 'charge': 0, 'aromatic': False, 'hcount': 2}),
         (4, {'element': 'C', 'charge': 0, 'aromatic': False, 'hcount': 2}),
         (5, {'element': 'C', 'charge': 0, 'aromatic': False, 'hcount': 2}),
         (6, {'element': 'C', 'charge': 0, 'aromatic': False, 'hcount': 1}),
         (7, {'element': 'F', 'charge': 0, 'aromatic': False, 'hcount': 0}),
         (8, {'element': 'C', 'charge': 0, 'aromatic': False, 'hcount': 2}),],
        [(0, 1, {'order': 1}),
         (0, 2, {'order': 1}),
         (0, 3, {'order': 1}),
         (0, 8, {'order': 1}),
         (3, 4, {'order': 1}),
         (4, 5, {'order': 1}),
         (5, 6, {'order': 1}),
         (6, 7, {'order': 1}),
         (6, 8, {'order': 1}),],
        False,
    ),
    ('[C@H]1(Br)CC1',
        [(0, {'element': 'C', 'charge': 0, 'aromatic': False, 'hcount': 1, 'rs_isomer': (3, 0, 1, 2)}),
         (1, {'element': 'Br', 'charge': 0, 'aromatic': False, 'hcount': 0}),
         (2, {'element': 'C', 'charge': 0, 'aromatic': False, 'hcount': 2}),
         (3, {'element': 'C', 'charge': 0, 'aromatic': False, 'hcount': 2}),],
        [(0, 1, {'order': 1}),
         (0, 2, {'order': 1}),
         (0, 3, {'order': 1}),
         (2, 3, {'order': 1}),],
        False,
    ),
    ('C1C[C@H]1Br',
        [(0, {'element': 'C', 'charge': 0, 'aromatic': False, 'hcount': 2}),
         (1, {'element': 'C', 'charge': 0, 'aromatic': False, 'hcount': 2}),
         (2, {'element': 'C', 'charge': 0, 'aromatic': False, 'hcount': 1, 'rs_isomer': (0, 2, 1, 3)}),
         (3, {'element': 'Br', 'charge': 0, 'aromatic': False, 'hcount': 0}),],
        [(0, 1, {'order': 1}),
         (1, 2, {'order': 1}),
         (0, 2, {'order': 1}),
         (2, 3, {'order': 1}),],
        False,
    ),
    (r'F\C=C\F',
        [(0, {'element': 'F', 'charge': 0, 'aromatic': False, 'hcount': 0, 'ez_isomer': (0, 1, 2, 3, 'trans')}),
         (1, {'element': 'C', 'charge': 0, 'aromatic': False, 'hcount': 1}),
         (2, {'element': 'C', 'charge': 0, 'aromatic': False, 'hcount': 1}),
         (3, {'element': 'F', 'charge': 0, 'aromatic': False, 'hcount': 0, 'ez_isomer': (3, 2, 1, 0, 'trans')}),],
        [(0, 1, {'order': 1}),
         (1, 2, {'order': 2}),
         (2, 3, {'order': 1}),],
        False,
    ),
    (r'C(/F)=C\F',
        [(0, {'element': 'C', 'charge': 0, 'aromatic': False, 'hcount': 1}),
         (1, {'element': 'F', 'charge': 0, 'aromatic': False, 'hcount': 0, 'ez_isomer': (1, 0, 2, 3, 'trans')}),
         (2, {'element': 'C', 'charge': 0, 'aromatic': False, 'hcount': 1}),
         (3, {'element': 'F', 'charge': 0, 'aromatic': False, 'hcount': 0, 'ez_isomer': (3, 2, 0, 1, 'trans')}),],
        [(0, 1, {'order': 1}),
         (0, 2, {'order': 2}),
         (2, 3, {'order': 1}),],
        False,
    ),
    (r'C(\F)=C\F',
        [(0, {'element': 'C', 'charge': 0, 'aromatic': False, 'hcount': 1}),
         (1, {'element': 'F', 'charge': 0, 'aromatic': False, 'hcount': 0, 'ez_isomer': (1, 0, 2, 3, 'cis')}),
         (2, {'element': 'C', 'charge': 0, 'aromatic': False, 'hcount': 1}),
         (3, {'element': 'F', 'charge': 0, 'aromatic': False, 'hcount': 0, 'ez_isomer': (3, 2, 0, 1, 'cis')}),],
        [(0, 1, {'order': 1}),
         (0, 2, {'order': 2}),
         (2, 3, {'order': 1}),],
        False,
    ),
))
def test_read_smiles(smiles, node_data, edge_data, explicit_h):
    found = read_smiles(smiles, explicit_hydrogen=explicit_h)
    print(found.nodes(data=True))
    print(found.edges(data=True))
    print(smiles)
    expected = make_mol(node_data, edge_data)
    assertEqualGraphs(found, expected)

@pytest.mark.parametrize('smiles, error_type', (
    ('[CL-]', ValueError),
    ('[HH]', ValueError),
    ('c1ccccc2', KeyError),
    ('C%1CC1', ValueError),
    ('c1c1CC', ValueError),
    ('CC11C', ValueError),
    ('1CCC1', ValueError),
    ('ccc', SyntaxError),
    ('C=1CC-1', ValueError),
))
def test_invalid_smiles(smiles, error_type):
    with pytest.raises(error_type):
        read_smiles(smiles)



@pytest.mark.parametrize('smiles,expected', [
    ('N[C@](Br)(O)C', 'N Br O C'),
    ('Br[C@](O)(N)C', 'Br O N C'),
    ('O[C@](Br)(C)N', 'O Br C N'),
    ('Br[C@](C)(O)N', 'Br C O N'),
    ('C[C@](Br)(N)O', 'C Br N O'),
    ('Br[C@](N)(C)O', 'Br N C O'),
    ('C[C@@](Br)(O)N', 'C Br N O'),
    ('Br[C@@](N)(O)C', 'Br N C O'),
    ('[C@@](C)(Br)(O)N', 'C Br N O'),
    ('[C@@](Br)(N)(O)C', 'Br N C O'),
    ('N[C@H](O)C', 'N C O C'),
    ('FC1C[C@](Br)(Cl)CCC1', 'C Br Cl C'),
    ('[C@]1(Br)(Cl)CCCC(F)C1', 'C Br Cl C'),
])
def test_chiral(smiles, expected):
    molecule = read_smiles(smiles)
    found = None
    for n_idx in molecule:
        if 'rs_isomer' in molecule.nodes[n_idx]:
            found = molecule.nodes[n_idx]['rs_isomer']
            break
    assert [molecule.nodes[n]['element'] for n in found] == expected.split()


@pytest.mark.parametrize('smiles, records',[
    (r'F/C=C/F', [(0, 1, 2, 3, 'trans'), (3, 2, 1, 0, 'trans')]),
    (r'C(\F)=C/F', [(1, 0, 2, 3, 'trans'), (3, 2, 0, 1, 'trans')]),
    (r'F\C=C/F', [(0, 1, 2, 3, 'cis'), (3, 2, 1, 0, 'cis')]),
    (r'C(/F)=C/F', [(1, 0, 2, 3, 'cis'), (3, 2, 0, 1, 'cis')]),
    (r'F/C(CC)=C/F', [(0, 1, 4, 5, 'trans'), (5, 4, 1, 0, 'trans')]),
    (r'F/C=C=C=C/F', [(0, 1, 4, 5, 'trans'), (5, 4, 1, 0, 'trans')]),
    (r'F/C=C=C=C\F', [(0, 1, 4, 5, 'cis'), (5, 4, 1, 0, 'cis')]),
    ('c1ccccc1', None)
])
def test_cis_trans_reading(smiles, records):
    mol = read_smiles(smiles, explicit_hydrogen=False)
    if records:
        for record in records:
            assert mol.nodes[record[0]]['ez_isomer'] == record
    else:
        assert len(nx.get_node_attributes(mol, 'ez_isomer')) == 0


@pytest.mark.parametrize('smiles', (
    'c1c[nH]cc1',
    'c1cNcc1',
    'c1cScc1',
    'c1cnc[nH]1',
    'c1cncN1',
    'c1cscc1',))
def test_kekulize(smiles):
    g = read_smiles(smiles)
    assert len(g) > 0


@pytest.mark.parametrize('smiles', (
    'cc',
    'cn',))
def test_skip_kekulize(smiles):
    g = read_smiles(smiles, reinterpret_aromatic=False)
    for node in g.nodes:
        assert g.nodes[node]['aromatic']
    assert g.edges[(0, 1)]['order'] == 1.5

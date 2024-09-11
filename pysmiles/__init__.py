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

"""
PySMILES: The lightweight python module for reading and writing SMILES strings.
"""

import builtins
from json import load


try:
    from importlib.resources import files, as_file
    import atexit
    from contextlib import ExitStack
except ImportError:  # pragma: no cover
    from pathlib import Path
    PTE_FILE_NAME = Path(__file__).parent / "PTE.json"
    del Path
else:
    ref = files('pysmiles') / "PTE.json"
    file_manager = ExitStack()
    atexit.register(file_manager.close)
    PTE_FILE_NAME = file_manager.enter_context(as_file(ref))
    del files, as_file, atexit, ExitStack


def _read_pte(pte_file_name):
    pte = {}
    with open(pte_file_name) as pte_file:
        data = load(pte_file)

    for row in data['Rows']:
        element = dict(zip(data['Columns'], row))
        pte[element['Symbol']] = element
        pte[element['AtomicNumber']] = element
    return pte


PTE = _read_pte(PTE_FILE_NAME)

from .read_smiles import read_smiles
from .write_smiles import write_smiles
from .smiles_helper import (fill_valence, add_explicit_hydrogens,
                            remove_explicit_hydrogens, correct_aromatic_rings)

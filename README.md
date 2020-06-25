[![Build Status](https://travis-ci.org/pckroon/pysmiles.svg?branch=master)](https://travis-ci.org/pckroon/pysmiles)
[![Coverage Status](https://coveralls.io/repos/github/pckroon/pysmiles/badge.svg?branch=master)](https://coveralls.io/github/pckroon/pysmiles?branch=master)

# pysmiles: The lightweight and pure-python SMILES reader and writer

This is a small project I started because I couldn't find any SMILES reader or
writer that was easy to install (read: Python only). Currently, the writer is
extremely basic, and although it should produce valid SMILES they won't be
pretty. The reader is in a better state, and should be usable.

SMILES strings are assumed to be as specified by the
[OpenSmiles standard][opensmiles].

## Molecules
Molecules are depicted as [Networkx][networkx] graphs. Atoms are the nodes of
the graph, and bonds are the edges. Nodes can have the following attributes:
- element: str. This describes the element of the atom. Defaults to '\*'
    meaning unknown.
- aromatic: bool. Whether the atom is part of an (anti)-aromatic system. 
    Defaults to False.
- isotope: float. The mass of the atom. Defaults to unknown.
- hcount: int. The number of implicit hydrogens attached to this atom.
    Defaults to 0.
- charge: int. The charge of this atom. Defaults to 0.
- class: int. The "class" of this atom. Defaults to 0.

Edges have the following attributes:
- order: Number. The bond order. 1.5 is used for aromatic bonds. Defaults to 1.

There is currently no way of specifying stereo chemical information, and this
is discarded upon reading. Somewhere in the future this will probably be
stored in the "stereo" attribute of nodes.

## Reading SMILES
The function `read_smiles(smiles, explicit_hydrogen=False,
zero_order_bonds=True, reinterpret_aromatic=True)` can be used to parse a
SMILES string. It should not be used to validate whether a string is a valid
SMILES string --- the function does very little validation whether your SMILES string makes chemical sense.
Edges in the created molecule will always have an 'order'
attribute. Nodes will have the relevant attributes in so far they are
specified. Atoms for which the element is not known (\*) will not have an
element attribute.
- `explicit_hydrogen` determines whether hydrogen atoms should be
    represented as explicit nodes in the created molecule, or implicit in the
    'hcount' attribute.
- `zero_order_bonds` determines whether zero-order bonds (.) in the SMILES
    string should result in edges in the produced molecule.
- `reinterpret_aromatic` determines whether aromaticity should be 
    reinterpreted, and determined from the constructed molecule, or whether
    the aromaticity specifications from the SMILES string (lower case 
    elements) should be taken as leading. If `True`, will also set bond orders 
    to 1 for bonds that are not part of an aromatic ring and have a bond order 
    of 1.5. If `False`, will create a molecule using *only* the information in 
    the SMILES string.

### Stereochemical information
Currently the library cannot handle stereochemical information, neither E/Z nor
R/S. Any stereochemical information that was in the SMILES string will be
*discarded* upon parsing. This means there will be no difference between
parsing *e.g.* `N[C@](Br)(O)C`, `N[C@@](Br)(O)C` and `NC(Br)(O)C`. Parsing
these *will result in the same molecule*. The same holds for *e.g.* `F/C=C/F`
and `FC=CF`. These will result in the same molecule.

Whenever stereochemical information is being discarded a warning will be
logged using the built-in `logging` module. If you want to disable all the
messages logged by `pysmiles` you can add the following snippet to your code,
without interfering with any logging by your own code:

```python
import logging
logging.getLogger('pysmiles').setLevel(logging.CRITICAL)  # Anything higher than warning
```


## Writing SMILES
The function `write_smiles(molecule, default_element='*', start=None)` can be
used to write SMILES strings from a molecule. The function does *not* check 
whether your molecule makes chemical sense. Instead, it writes a SMILES 
representation of the molecule you provided, and nothing else.
- `default_element` is the element to use for nodes that do not have an 
    'element' attribute.
- `start` is the key of the node where the depth first traversal should be 
    started. Something clever is done if not specified.

## Additional functions
In addition to these two core functions, four more functions are exposed that
can help in creating chemically relevant molecules with minimal work.

- `fill_valence(mol, respect_hcount=True, respect_bond_order=True,
                 max_bond_order=3)`
    This function will fill the valence of all atoms in your molecule by 
    incrementing the 'hcount' and, if specified, bond orders. Note that it does
    not use 'charge' attribute to find the correct valence.
    - `repect_hcount`: bool. Whether existing hcounts can be overwritten.
    - `respect_bond_order`: bool. Whether bond orders can be changed
    - `max_bond_order`: int. The maximum bond order that will be set.
- `add_explicit_hydrogens(mol)`
    This function transforms implicit hydrogens, specified by 'hcount' 
    attributes, to explicit nodes.
- `remove_explicit_hydrogens(mol)`
    This function does the inverse of `add_explicit_hydrogens`: it will remove
    explicit hydrogen nodes and add them to the relevant 'hcount' attributes.
- `correct_aromatic_rings(mol)`
    This function marks all (anti)-aromatic atoms in your molecule, and sets 
    all bonds between (anti)-aromatic atoms to order 1.5.
    It fills the valence of all atoms (see also `fill_valence`) before trying
    to figure our which atoms are aromatic. It works by first finding all 
    atoms that are in a ring. Next, for every atom in every ring it is checked
    whether the atoms are sp2 hybridized (note that this is a vague term. 
    Strictly speaking we check whether their element is something that *could*
    be aromatic, and whether they have 2 or 3 bonds.). Finally, the number of 
    electrons per ring is counted, and if this is even, the atoms in the ring
    are said to be aromatic.
    This function is the most fragile in the whole library, and I expect it to
    produce wrong answers in some cases. In particular for fused (aromatic)
    ring systems (such as indole) and rings with extracyclic heteroatoms
    (O=C1C=CC=C1). Buyer beware.

## Examples
### Reading
```python
from pysmiles import read_smiles

smiles = 'C1CC[13CH2]CC1C1CCCCC1'
mol = read_smiles(smiles)

print(mol.nodes(data='element'))
# [(0, 'C'),
#  (1, 'C'),
#  (2, 'C'),
#  (3, 'C'),
#  (4, 'C'),
#  (5, 'C'),
#  (6, 'C'),
#  (7, 'C'),
#  (8, 'C'),
#  (9, 'C'),
#  (10, 'C'),
#  (11, 'C')]
print(mol.nodes(data='hcount'))
# [(0, 2),
#  (1, 2),
#  (2, 2),
#  (3, 2),
#  (4, 2),
#  (5, 1),
#  (6, 1),
#  (7, 2),
#  (8, 2),
#  (9, 2),
#  (10, 2),
#  (11, 2)]

mol_with_H = read_smiles(smiles, explicit_hydrogen=True)
print(mol_with_H.nodes(data='element'))
# [(0, 'C'),
#  (1, 'C'),
#  (2, 'C'),
#  (3, 'C'),
#  (4, 'C'),
#  (5, 'C'),
#  (6, 'C'),
#  (7, 'C'),
#  (8, 'C'),
#  (9, 'C'),
#  (10, 'C'),
#  (11, 'C'),
#  (12, 'H'),
#  (13, 'H'),
#  (14, 'H'),
#  (15, 'H'),
#  (16, 'H'),
#  (17, 'H'),
#  (18, 'H'),
#  (19, 'H'),
#  (20, 'H'),
#  (21, 'H'),
#  (22, 'H'),
#  (23, 'H'),
#  (24, 'H'),
#  (25, 'H'),
#  (26, 'H'),
#  (27, 'H'),
#  (28, 'H'),
#  (29, 'H'),
#  (30, 'H'),
#  (31, 'H'),
#  (32, 'H'),
# (33, 'H')]
```

### Writing
```python
import networkx as nx
from pysmiles import write_smiles, fill_valence

mol = nx.Graph()
mol.add_edges_from([(0, 1), (1, 2), (1, 3), (3, 4), (1, 5), (3, 6)])
for idx, ele in enumerate('CCCCOCO'):
    mol.nodes[idx]['element'] = ele
mol.nodes[4]['charge'] = -1
mol.nodes[4]['hcount'] = 0
mol.edges[3, 6]['order'] = 2

print(write_smiles(mol))
# [O-]C(=O)C([C])([C])[C]
fill_valence(mol, respect_hcount=True)
print(write_smiles(mol))
# [O-]C(=O)C(C)(C)C
```

## Limitations
- The writer produces non-recommended SMILES strings (as per OpenSmiles).
- `fill_valence` does not use 'charge' to find the correct valence.
- `correct_aromatic_rings` is fragile.
- There is currently no way of specifying stereo chemical information. The 
    parser can deal with it, but it will be discarded.
- It only processes SMILES. This might later be extended to e.g. InChi, SLN,
    SMARTS, etc.

## Requirements
- [networkx][networkx]

## Similar projects
There are more python projects that deal with SMILES, and I try to list at 
least some of them here. If yours is missing, feel free to open up a PR.
- [PySMILE](https://github.com/jhosmer/PySmile): A similar named project, 
    capable of encoding/decoding SMILE format objects. Doesn't deal with 
    SMILES.
- [RDKit](https://github.com/rdkit/rdkit): A collection of cheminformatics and 
    machine-learning software, capable of reading and writing SMILES, InChi, 
    and others.
- [OpenEye Chem toolkit](https://www.eyesopen.com/oechem-tk): The OpenEye 
    chemistry toolkit is a programming library for chemistry and 
    cheminformatics. It is capable of dealing with (canonical) SMILES and 
    InChi.

## License
PySmiles is distributed under the Apache 2.0 license.
    Copyright 2018 Peter C Kroon

    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

        http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.


[opensmiles]: http://opensmiles.org/
[networkx]: https://networkx.github.io/

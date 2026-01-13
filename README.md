[![Build Status](https://travis-ci.org/pckroon/pysmiles.svg?branch=master)](https://travis-ci.org/pckroon/pysmiles)
[![Coverage Status](https://coveralls.io/repos/github/pckroon/pysmiles/badge.svg?branch=master)](https://coveralls.io/github/pckroon/pysmiles?branch=master)

# pysmiles: The lightweight and pure-python SMILES reader and writer

This is a small project I started because I couldn't find any SMILES reader or
writer that was easy to install (read: Python only). Currently, the writer is
extremely basic, and although it should produce valid SMILES they won't be
pretty, but see also issue #17. The reader is in a better state, and should be usable.

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
- rs_isomer and ez_isomer: See below.
- stereo: Not yet used, but reserved for more interpretable stereochemical 
  information

Edges have the following attributes:
- order: Number. The bond order. 1.5 is used for aromatic bonds. Defaults to 1.

## Reading SMILES
The function `read_smiles(smiles, explicit_hydrogen=False,
zero_order_bonds=True, reinterpret_aromatic=True, strict=True)` can be used to 
parse a SMILES string. It should not be used to validate whether a string is 
a valid SMILES string --- the function does very little validation whether 
your SMILES string makes chemical sense.
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
    of 1.5. See also below. If `False`, will create a molecule using *only* the 
    information in the SMILES string.
- `strict` determines whether the function be a bit more critical about what 
    constitutes a valid SMILES string. In particular, checks whether the 
    specified aromatic fragments can be kekulized, and whether atom valency 
    is sensible. 

### Stereochemical information
Tetrahedral chirality is stored on nodes in the 'rs_isomer' attribute. It 
consists of 4 node indices. The first, together with the central atom, 
indicates the axis of rotation/observation. The last three indices indicate 
the direction of rotation in counter-clockwise order. For example, when 
parsing the SMILES `N[C@](Br)(O)C` node 1 (the central carbon) will have 
'rs_isomer' attribute `(0, 2, 3, 4)`. This means when looking along the axis 
(0, 1) (nitrogen to carbon), the nodes 2 (bromine), 3 (oxygen), 4 (methyl) 
will be in counter-clockwise order. If the central atom has an implicit hydrogen
(which of course does not have its own node index) we use the index of the 
central atom instead.
There is currently no way to easily convert this to R/S labels.

E/Z chirality around double bonds is encoded on nodes in the 'ez_isomer' 
attribute. This attribute consists of 4 node indices and 1 string indication 
either 'cis' or 'trans'. The node indices indicate the dihedral angle. For 
example, when parsing the SMILES `Br/C=C\F` node 0 (the bromine) will have
'ez_isomer' attribute `[0, 1, 2, 3, 'cis']`. This means the dihedral angle of
the axes (0, 1) and (2, 3) around the central axis (1, 2) is 0 degrees. 
There is currently no way to easily convert this to E/Z labels.

### Aromaticity
Aromaticity in SMILES has little to do with how an organic chemist might 
interpret the term. Here we use "aromatic" to mean "delocalization induced 
molecular equivalence" (DIME). DIME means that molecular fragments that are 
equivalent under delocalization should get the same representation. As an 
example, `C1=CC=CC=C1` and `C=1C=CC=CC1` are equivalent, and should both be 
represented as `c1ccccc1`. On the other hand, pyrolle (`c1c[nH]cc1`) can only
be represented as `C1=CNC=C1`, and is not considered aromatic. A series of 
blogposts with a decent writeup on the subject can be found [here][depth_first].

Despite how simple it seems at a glance, here be dragons lurking in the details.
The main problem is that in order to determine whether an edge is aromatic, 
you need to generate *all* simple cycles. To illustrate why this is a problem,
let's look at buckminsterfullerene: it contains a staggering *374,237,206* 
cycles! It's really not feasible to investigate all of those separately.

Therefor we approximate the aromaticity for ring systems that are too large. By
default we consider 30 atoms large. For the exact algorithm, please see the
docstring of `dekekulize`. The approximation usually works fine, but may 
erroneously mark some triangle edges as aromatic. If you have a SMILES 
string with a large ring system for which the approximation does not work, 
please open an issue. As a workaround, you can read the SMILES with 
`reinterpret_aromatic=False`, and manually call `correct_aromatic_rings` 
with more appropriate thresholds.

A second complication comes from wildcard (\*) atoms. If possible we try to 
kekulize the provided molecule without using them, but there may be edge 
cases where this does not work as intended.

## Writing SMILES
The function `write_smiles(molecule, default_element='*', start=None)` can be
used to write SMILES strings from a molecule. The function does *not* check 
whether your molecule makes chemical sense. Instead, it writes a SMILES 
representation of the molecule you provided, and nothing else. It does not 
yet know how to write stereochemical information that is present in your 
molecule, and will log a warning if this happens. See below if you need to 
silence these.
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
    incrementing the 'hcount' and, if specified, bond orders.
    - `repect_hcount`: bool. Whether existing hcounts can be overwritten.
    - `respect_bond_order`: bool. Whether bond orders can be changed.
    - `max_bond_order`: int. The maximum bond order that will be set.
- `add_explicit_hydrogens(mol)`
    This function transforms implicit hydrogens, specified by 'hcount' 
    attributes, to explicit nodes.
- `remove_explicit_hydrogens(mol)`
    This function does the inverse of `add_explicit_hydrogens`: it will remove
    explicit hydrogen nodes and add them to the relevant 'hcount' attributes.
- `correct_aromatic_rings(mol, strict=True, estimation_threshold=None,
                          max_ring_size=None)`
    This function marks all (anti)-aromatic atoms in your molecule, and sets 
    all bonds between (anti)-aromatic atoms to order 1.5.

    It works by first finding all atoms that could participate in an aromatic
    system. These are atoms that can make a double bond without introducing 
    a charge. Then, it tries to assign alternating double and single bonds  
    between the atoms found. If it cannot and `strict` is `True` an error 
    will be raised. Otherwise, it starts to look for cycles that consist of 
    such alternating double and single bonds. Atoms in these cycles will be 
    marked as aromatic, and these edges will get a bond order of 1.5. The 
    relevant atoms that are not in any cycle will be marked as not-aromatic, 
    and these bonds will be assigned alternating bond orders of 1 and 2.
- `kekulize(mol)`
    Assigns alternating single and double bonds to all aromatic regions in
    ``mol``. Will raise an error if it cannot.
- `dekekulize(mol, estimation_threshold=None, max_ring_size=None)`
    Finds all cycles in `mol` that consist of alternating single and double
    bonds, and marks them as aromatic.

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

## Logging
Pysmiles uses the python logging module to log warnings. If you need to silence
these, you can use the following snippet in your code:
```python
import logging
logging.getLogger('pysmiles').setLevel(logging.CRITICAL)  # Anything higher than warning
```

## Limitations
- The writer produces non-recommended SMILES strings (as per OpenSmiles).
- The writer is better described as a "serializer": if the graph provided
    doesn't make chemical sense the produced "SMILES" string will be an
    exact representation of that graph. Because of this, the SMILES string
    will be invalid though.
- The writer cannot deal with stereochemistry, and ignores it.
- There is no way to easily transform stereochemical information to 
  human-readable labels like R/S/E/Z.
- Extended stereochemical centers, such as described by `NC(Br)=[C@]=C(O)C` 
    are not supported.
- Although the SMILES `F/C(/Br)=C\F` is valid, pysmiles does not understand it.
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
[depth_first]: https://depth-first.com/articles/2020/02/10/a-comprehensive-treatment-of-aromaticity-in-the-smiles-language/

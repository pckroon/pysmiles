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
"element", "isotope", "hcount", "charge", and "class". There is currently no
way of specifying stereochemical information, and this is discarded upon
reading. Somewhere in the future this will probably be stored in the "stereo"
attribute of nodes. Edges will have an "order" attribute.

## Examples
### Reading
```python
from pysmiles import read_smiles

smiles = 'C1CC[13CH2]CC1C1CCCCC1'
mol = read_smiles(smiles)

print(mol.nodes(data='element'))
#  [(0, 'C'), (1, 'C'), (2, 'C'), (3, 'C'), (4, 'C'), (5, 'C'), (6, 'C'),
#   (7, 'C'), (8, 'C'), (9, 'C'), (10, 'C'), (11, 'C'), (12, 'H'), (13, 'H'),
#   (14, 'H'), (15, 'H'), (16, 'H'), (17, 'H'), (18, 'H'), (19, 'H'),
#   (20, 'H'), (21, 'H'), (22, 'H'), (23, 'H'), (24, 'H'), (25, 'H'),
#   (26, 'H'), (27, 'H'), (28, 'H'), (29, 'H'), (30, 'H'), (31, 'H'),
#   (32, 'H'), (33, 'H')]

mol_without_H = read_smiles(smiles, explicit_H=False)
print(mol_without_H.nodes(data='element'))
# [(0, 'C'), (1, 'C'), (2, 'C'), (3, 'C'), (4, 'C'), (5, 'C'), (6, 'C'),
#  (7, 'C'), (8, 'C'), (9, 'C'), (10, 'C'), (11, 'C')]
print(mol_without_H.nodes(data='hcount'))
# [(0, 2), (1, 2), (2, 2), (3, 2), (4, 2), (5, 1), (6, 1), (7, 2), (8, 2),
#  (9, 2), (10, 2), (11, 2)]
```

### Writing
```python
import networkx as nx
from pysmiles import write_smiles

mol = nx.Graph()
mol.add_edges_from([(0, 1), (1, 2), (1, 3), (3, 4), (1, 5), (3, 6)])
for idx, ele in enumerate('CCCCOCO'):
    mol.nodes[idx]['element'] = ele
mol.nodes[4]['charge'] = -1
mol.edges[3, 6]['order'] = 2

print(write_smiles(mol))
# CC(C)(C(=O)[O-])C
```

## Limitations
- The writer produces non-recommended SMILES strings (as per OpenSmiles). In
    addition, it doesn't use the 'hcount' and valence to determine bond orders.
- There is currently no way of specifying stereochemical information. The parser
can deal with it, but it will be discarded.
- There are no tests.
- It is not pip-installable
- It only reads SMILES. This might later be extended to e.g. InChi, SLN, etc.

## Requirements
- [networkx][networkx]

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
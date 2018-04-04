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

import networkx as nx
import enum
import re


ISOTOPE_PATTERN = r'(?P<isotope>[\d]+)?'
ELEMENT_PATTERN = r'(?P<element>b|c|n|o|s|p|\*|[A-Z][a-z]{0,2})'
STEREO_PATTERN = r'(?P<stereo>@|@@|@TH[1-2]|@AL[1-2]|@SP[1-3]|@OH[\d]{1,2}|'\
                  '@TB[\d]{1,2})?'
HCOUNT_PATTERN = r'(?P<hcount>H[\d]?)?'
CHARGE_PATTERN = r'(?P<charge>--|\+\+|-[\d]{0,2}|\+[\d]{0,2})?'
CLASS_PATTERN = r'(?::(?P<class>[\d]+))?'
ATOM_PATTERN = re.compile(ISOTOPE_PATTERN + ELEMENT_PATTERN + STEREO_PATTERN +
                          HCOUNT_PATTERN + CHARGE_PATTERN + CLASS_PATTERN)


class TokenType(enum.Enum):
    ATOM = enum.auto()
    BOND_TYPE = enum.auto()
    BRANCH_START = enum.auto()
    BRANCH_END = enum.auto()
    RING_NUM = enum.auto()


def tokenize(smiles):
    organic_subset = 'B C N O P S F Cl Br I * b c n o s p'.split()
    smiles = iter(smiles)
    token = ''
    peek = None
    while True:
        try:
            char = peek if peek else next(smiles)
            peek = None
        except StopIteration:
            break
        if char == '[':
            token = char
            for char in smiles:
                token += char
                if char == ']':
                    break
            yield TokenType.ATOM, token
        elif char in organic_subset:
            peek = next(smiles, '')
            if char + peek in organic_subset:
                yield TokenType.ATOM, char + peek
                peek = None
            else:
                yield TokenType.ATOM, char
        elif char in '-=#$:.':
            yield TokenType.BOND_TYPE, char
        elif char == '(':
            yield TokenType.BRANCH_START, '('
        elif char == ')':
            yield TokenType.BRANCH_END, ')'
        elif char == '%':
            yield TokenType.RING_NUM, int(next(smiles) + next(smiles))
        elif char.isdigit():
            yield TokenType.RING_NUM, int(char)


def parse_atom(atom):
    atom = atom.strip('[]')
    match = ATOM_PATTERN.fullmatch(atom)
    if match is None:
        raise ValueError('The atom {} is malformatted'.format(atom))
    out = match.groupdict()
    out = {k: v for k, v in out.items() if v is not None}
    if 'isotope' in out:
        out['isotope'] = int(out['isotope'])
    if out['element'] == '*':
        del out['element']
    if 'hcount' in out:
        # Starts with H
        hcount = out['hcount']
        if hcount == 'H':
            out['hcount'] = 1
        else:
            out['hcount'] = int(hcount[1:])
    if 'stereo' in out:
        print("I don't quite know how to handle stereo yet...")
    if 'charge' in out:
        charge = out['charge']
        if charge == '--':
            charge = -2
        elif charge == '++':
            charge = +2
        elif charge == '-':
            charge = -1
        elif charge == '+':
            charge = +1
        else:
            charge = int(charge)
        out['charge'] = charge
    if 'class' in out:
        out['class'] = int(out['class'])
    return out


def read_smiles(smiles):
    bond_to_order = {'-': 1, '=': 2, '#': 3, '$': 4, ':': 1.5, '.': 0}
    mol = nx.Graph()
    anchor = None
    idx = 0
    default_bond = 1
    next_bond = None
    branches = []
    ring_nums = {}
    for tokentype, token in tokenize(smiles):
        if tokentype == TokenType.ATOM:
            mol.add_node(idx, **parse_atom(token))
            if anchor is not None:
                if next_bond is None:
                    next_bond = default_bond
                mol.add_edge(anchor, idx, order=next_bond)
                next_bond = None
            anchor = idx
            idx += 1
        elif tokentype == TokenType.BRANCH_START:
            branches.append(anchor)
        elif tokentype == TokenType.BRANCH_END:
            anchor = branches.pop()
        elif tokentype == TokenType.BOND_TYPE:
            next_bond = bond_to_order[token]
        elif tokentype == TokenType.RING_NUM:
            if token in ring_nums:
                jdx, order = ring_nums[token]
                if next_bond is None and order is None:
                    next_bond = default_bond
                elif order is None:  # Note that the check is needed,
                    next_bond = next_bond  # But this could be pass.
                elif next_bond is None:
                    next_bond = order
                elif next_bond != order:  # Both are not None
                    raise ValueError('Conflicting bond orders for ring '
                                     'between indices {}'.format(token))
                # idx is the index of the *next* atom we're adding. So: -1.
                mol.add_edge(idx - 1, jdx, order=next_bond)
                next_bond = None
                del ring_nums[token]
            else:
                ring_nums[token] = (idx - 1, next_bond)
                next_bond = None
    return mol


if __name__ == '__main__':
    example = 'OCC(CCC)C(C(C)C)CCC'
    mol = read_smiles(example)
    mol = read_smiles('c1ccccc1')
    mol = read_smiles('C(C(C(C(C(C(C(C(C(C(C(C(C(C(C(C(C(C(C(C(C))))))))))))))))))))C')
    mol = read_smiles('N1CC2CCCC2CC1')
    mol = read_smiles('C%25CCCCC%25')
    mol3 = read_smiles('C1CCCCC1C1CCCCC1')
    mol4 = read_smiles('C1CC[13CH2]CC1C1CCCCC1')
    print(mol4.nodes(data=True))
    mol = read_smiles('[Rh-](Cl)(Cl)(Cl)(Cl)$[Rh-](Cl)(Cl)(Cl)Cl')
    mol2 = read_smiles('[15OH1-:4][HoH3]')
    molx = read_smiles('[O--][13CH2+3][14CH3+2][C@H4][Rh@OH19]')
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 10 17:08:10 2018

@author: Peter Kroon
"""

import networkx as nx
import itertools
import enum

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
            try:
                peek = next(smiles)
            except StopIteration:
                peek = ''
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
    if not atom.startswith('['):
        return {'element': atom} if atom != '*' else {}
    # '[' isotope? symbol chiral? hcount? charge? class? ']'
    # isotope ::= NUMBER
    # chiral ::= '@' | '@@'
    # hcount ::= 'H' | 'H' DIGIT
    # charge ::= '-' | '-' DIGIT? DIGIT | '+' | '+' DIGIT? DIGIT | '--' deprecated | '++' deprecated
    # class ::= ':' NUMBER
    out = {}
    atom = atom.strip('[]')
    atom, *class_ = atom.split(':')

    # Isotope?
    idx = 0
    while idx < len(atom):
        if not atom[idx].isdigit():
            break
        idx += 1
    isotope = atom[:idx]
    atom = atom[idx:]
    if isotope:
        out['isotope'] = int(isotope)

    # Symbol
    idx = 0
    while idx < len(atom):
        if (not atom[idx].isalpha() and atom[idx] != '*') or (idx != 0 and atom[idx].isupper()):
            break
        idx += 1
    element = atom[:idx]
    atom = atom[idx:]
    if element != '*':
        out['element'] = element

    # Stereo?
    idx = 0
    while idx < len(atom):
        if atom[idx] in 'H+-':
            break
        idx += 1
    stereo = atom[:idx]
    atom = atom[idx:]
    if stereo:
        print("I don't know how to deal with stereochemistry yet...")
        out['stereo'] = stereo

    # Hcount?
    idx = 0
    while idx < len(atom):
        if atom[idx] in '+-':
            break
        idx += 1
    hcount = atom[:idx]
    atom = atom[idx:]
    if len(hcount) == 2:
        hcount = int(hcount[1])
    elif hcount:
        hcount = 1
    if hcount:
        out['hcount'] = hcount

    # Charge?
    charge = atom
    if charge == '--':
        charge = -2
    elif charge == '++':
        charge = +2
    elif charge == '-':
        charge = -1
    elif charge == '+':
        charge = +1
    elif charge:
        charge = int(charge)
    else:
        charge = 0

    if charge:
        out['charge'] = charge
    if class_:
        out['class'] = int(class_[0])
    return out




def read_smiles(smiles):
    bond_to_order = {'-': 1, '=': 2, '#': 3, '$': 4, ':': 1.5, '.': 0}
    mol = nx.Graph()
    anchor = None
    idx = 0
    next_bond = 1
    branches = []
    ring_nums = {}
    for tokentype, token in tokenize(smiles):
        if tokentype == TokenType.ATOM:
            mol.add_node(idx, **parse_atom(token))
            if anchor is not None:
                mol.add_edge(anchor, idx, order=next_bond)
                next_bond = 1
            anchor = idx
            idx += 1
        elif tokentype == TokenType.BRANCH_START:
            branches.append(anchor)
        elif tokentype == TokenType.BRANCH_END:
            anchor = branches.pop()
        elif tokentype == TokenType.BOND_TYPE:
            # TODO: Bond type to ringnum
            next_bond = bond_to_order[token]
        elif tokentype == TokenType.RING_NUM:
            # idx is the index of the *next* atom we're adding. therefor: -1.
            if token in ring_nums:
                mol.add_edge(idx - 1, ring_nums[token], order=next_bond)
                next_bond = 1
                del ring_nums[token]
            else:
                ring_nums[token] = idx - 1
    return mol


if __name__ == '__main__':
    example = 'OCC(CCC)C(C(C)C)CCC'
    mol = read_smiles(example)
    mol = read_smiles('c1ccccc1')
    mol = read_smiles('C(C(C(C(C(C(C(C(C(C(C(C(C(C(C(C(C(C(C(C(C))))))))))))))))))))C')
    mol = read_smiles('N1CC2CCCC2CC1')
    mol = read_smiles('C%25CCCCC%25')
    mol = read_smiles('C1CCCCC1C1CCCCC1')
    mol = read_smiles('C1CC[13CH2]CC1C1CCCCC1')
    mol = read_smiles('[Rh-](Cl)(Cl)(Cl)(Cl)$[Rh-](Cl)(Cl)(Cl)Cl')
    mol2 = read_smiles('[15OH1-:4][HoH3]')

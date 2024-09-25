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
Contains helper functions for parsing and writing SMILES strings, as well as
some convenience functions for adding hydrogens, and detecting aromaticity.
"""

import logging
import re
import operator

import networkx as nx
from . import PTE

LOGGER = logging.getLogger(__name__)

ISOTOPE_PATTERN = r'(?P<isotope>[\d]+)?'
ELEMENT_PATTERN = r'(?P<element>b|c|n|o|s|p|as|se|\*|[A-Z][a-z]{0,2})'
STEREO_PATTERN = r'(?P<rs_isomer>@|@@|@TH[1-2]|@AL[1-2]|@SP[1-3]|@OH[\d]{1,2}|'\
                  r'@TB[\d]{1,2})?'
HCOUNT_PATTERN = r'(?P<hcount>H[\d]?)?'
CHARGE_PATTERN = r'(?P<charge>(-|\+)(\++|-+|[\d]{1,2})?)?'
CLASS_PATTERN = r'(?::(?P<class>[\d]+))?'
ATOM_PATTERN = re.compile(r'^\[' + ISOTOPE_PATTERN + ELEMENT_PATTERN +
                          STEREO_PATTERN + HCOUNT_PATTERN + CHARGE_PATTERN +
                          CLASS_PATTERN + r'\]$')

# Used to parse the electron config we get from the PTE. These configs look like
# [Rn]7s27p65f146d10(predicted) (for Og). We don't care about the core electrons
# ([Rn] in this case), nor whether it's predicted or measured. Note that the
# order of the orbitals is *not* consistent, but seems to be related to their
# energy. We capture the \d[spdf]\d+, and in particular the second number which
# is the number of electrons in that orbital. In addition, we capture the full
# electron config (from 7s...d10) in case we need to find the order of the
# orbitals later.
ELECTRON_CONFIG_PATTERN = (r'(?:\[[A-Z][a-z]?\])?((?:'
                           r'(?:[1-9]s(?P<s>[1-2]))|(?:[1-9]p(?P<p>[1-6]))|(?:[1-9]d(?P<d>10|[1-9]))|(?:[1-9]f(?P<f>1[0-4]|[1-9]))'
                           r'){1,4})(?:\(predicted\))?')
ELECTRON_CONFIG_PATTERN = re.compile(ELECTRON_CONFIG_PATTERN)

AROMATIC_ATOMS = "B C N O P S Se As *".split()


def parse_atom(atom):
    """
    Parses a SMILES atom token, and returns a dict with the information.

    Parameters
    ----------
    atom : str
        The atom string to interpret. Looks something like one of the
        following: "C", "c", "[13CH3-1:2]"

    Returns
    -------
    dict
        A dictionary containing at least 'element', 'aromatic', and 'charge'. If
        present, will also contain 'hcount', 'isotope', and 'class'.
    """
    defaults = {'charge': 0, 'hcount': 0, 'aromatic': False}
    if not atom.startswith('[') and not atom.endswith(']'):
        if atom != '*':
            # Don't specify hcount to signal we don't actually know anything
            # about it
            return {'element': atom.capitalize(), 'charge': 0,
                    'aromatic': atom.islower()}
        else:
            return defaults.copy()
    match = ATOM_PATTERN.match(atom)
    if match is None:
        raise ValueError('The atom {} is malformatted'.format(atom))
    out = defaults.copy()
    out.update({k: v for k, v in match.groupdict().items() if v is not None})

    if out.get('element', 'X').islower():
        out['aromatic'] = True

    parse_helpers = {
        'isotope': int,
        'element': str.capitalize,
        'rs_isomer': lambda x: x,
        'hcount': parse_hcount,
        'charge': parse_charge,
        'class': int,
        'aromatic': lambda x: x,
    }

    for attr, val_str in out.items():
        out[attr] = parse_helpers[attr](val_str)

    if out['element'] == '*':
        del out['element']

    if out.get('element') == 'H' and out.get('hcount', 0):
        raise ValueError("A hydrogen atom can't have hydrogens")

    if out.get('rs_isomer', False):
        out['rs_isomer'] = (out['rs_isomer'], [])

    return out


def format_atom(molecule, node_key, default_element='*'):
    """
    Formats a node following SMILES conventions. Uses the attributes `element`,
    `charge`, `hcount`, `rs_isomer`, `isotope` and `class`.

    Parameters
    ----------
    molecule : nx.Graph
        The molecule containing the atom.
    node_key : hashable
        The node key of the atom in `molecule`.
    default_element : str
        The element to use if the attribute is not present in the node.

    Returns
    -------
    str
        The atom as SMILES string.
    """
    node = molecule.nodes[node_key]
    name = node.get('element', default_element)
    charge = node.get('charge', 0)
    hcount = node.get('hcount', 0)
    stereo = node.get('rs_isomer', None)
    isotope = node.get('isotope', '')
    class_ = node.get('class', '')
    aromatic = node.get('aromatic', False)
    default_h = has_default_h_count(molecule, node_key)

    if stereo is not None or node.get('ez_isomer'):  # pragma: nocover
        LOGGER.warning("The SMILES writer does not write stereochemical information")

    if aromatic and name in AROMATIC_ATOMS:
        name = name.lower()

    if (stereo is None and isotope == '' and charge == 0 and default_h and class_ == '' and
            (name.lower() in 'b c n o p s *'.split() or name in 'F Cl Br I'.split())):
        return name

    if hcount:
        hcountstr = 'H'
        if hcount > 1:
            hcountstr += str(hcount)
    else:
        hcountstr = ''

    if charge > 0:
        chargestr = '+'
        if charge > 1:
            chargestr += str(charge)
    elif charge < 0:
        chargestr = '-'
        if charge < -1:
            chargestr += str(-charge)
    else:
        chargestr = ''

    if class_ != '':
        class_ = ':{}'.format(class_)
    fmt = '[{isotope}{name}{stereo}{hcount}{charge}{class_}]'
    return fmt.format(isotope=isotope, name=name, stereo='', hcount=hcountstr,
                      charge=chargestr, class_=class_)


def parse_hcount(hcount_str):
    """
    Parses a SMILES hydrogen count specifications.

    Parameters
    ----------
    hcount_str : str
        The hydrogen count specification to parse.

    Returns
    -------
    int
        The number of hydrogens specified.
    """
    if not hcount_str:
        return 0
    if hcount_str == 'H':
        return 1
    return int(hcount_str[1:])


def parse_charge(charge_str):
    """
    Parses a SMILES charge specification.

    Parameters
    ----------
    charge_str : str
        The charge specification to parse.

    Returns
    -------
    int
        The charge.
    """
    if not charge_str:
        return 0
    signs = {'-': -1, '+': 1}
    sign = signs[charge_str[0]]
    if len(charge_str) > 1 and charge_str[1].isdigit():
        charge = sign * int(charge_str[1:])
    else:
        charge = sign * charge_str.count(charge_str[0])
    return charge


def add_explicit_hydrogens(mol):
    """
    Adds explicit hydrogen nodes to `mol`, the amount is determined by the node
    attribute 'hcount'. Will remove the 'hcount' attribute.

    Parameters
    ----------
    mol : nx.Graph
        The molecule to which explicit hydrogens should be added. Is modified
        in-place.

    Returns
    -------
    None
        `mol` is modified in-place.
    """
    h_atom = parse_atom('[H]')
    if 'hcount' in h_atom:  # pragma: nocover; defensive
        del h_atom['hcount']
    for n_idx in list(mol.nodes):
        hcount = mol.nodes[n_idx].get('hcount', 0)
        idxs = range(max(mol) + 1, max(mol) + hcount + 1)
        # Get the defaults from parse_atom.
        mol.add_nodes_from(idxs, **h_atom.copy())
        mol.add_edges_from([(n_idx, jdx) for jdx in idxs], order=1)
        if 'hcount' in mol.nodes[n_idx]:
            del mol.nodes[n_idx]['hcount']
        if 'rs_isomer' in mol.nodes[n_idx]:
            # Replace the implicit hydrogen index  in the stereo definition for
            # the explicit one.
            mol.nodes[n_idx]['rs_isomer'] = tuple(n if n != n_idx else idxs[0] for n in mol.nodes[n_idx]['rs_isomer'])


def remove_explicit_hydrogens(mol):
    """
    Removes all explicit, simple hydrogens from `mol`. Simple means it is
    identical to the SMILES string "[H]", and has exactly one bond. Increments
    'hcount' where appropriate.

    Parameters
    ----------
    mol : nx.Graph
        The molecule whose explicit hydrogens should be removed. Is modified
        in-place.

    Returns
    -------
    None
        `mol` is modified in-place.
    """
    to_remove = set()
#    defaults = parse_atom('[H]')
    for n_idx in mol.nodes:
        node = mol.nodes[n_idx]
        neighbors = list(mol[n_idx])
        # TODO: get these defaults from parsing [H]. But do something smart
        #       with the hcount attribute.
        if (node.get('charge', 0) == 0 and node.get('element', '') == 'H' and
                'isotope' not in node and node.get('class', 0) == 0 and
                len(neighbors) == 1):
            neighbor = neighbors[0]
            if (mol.nodes[neighbor].get('element', '') == 'H' or
                    mol.edges[n_idx, neighbor].get('order', 1) != 1):
                # The molecule is H2, or the bond order is not 1.
                continue
            if mol.nodes[n_idx].get('ez_isomer'):
                anchor = mol.nodes[n_idx]['ez_isomer'][1]
                for ez_neighbor in mol[anchor]:
                    if (ez_neighbor != n_idx
                            and mol.edges[anchor, ez_neighbor].get('order', 1) == 1
                            and ez_neighbor not in to_remove):
                        # Found the other ligand that's connected to the anchor
                        # with a single bond
                        mol.nodes[ez_neighbor]['ez_isomer'] = mol.nodes[n_idx]['ez_isomer'].copy()
                        mol.nodes[ez_neighbor]['ez_isomer'][-1] = 'trans' if mol.nodes[ez_neighbor]['ez_isomer'][1] == 'cis' else 'cis'
                        mol.nodes[ez_neighbor]['ez_isomer'][0] = ez_neighbor

                        ez_partner = mol.nodes[n_idx]['ez_isomer'][3]
                        mol.nodes[ez_partner]['ez_isomer'][3] = ez_neighbor
                        mol.nodes[ez_partner]['ez_isomer'][-1] = mol.nodes[ez_neighbor]['ez_isomer'][-1]
                        break
                else:  # No break
                    # Don't remove the hydrogen, since we need it for the
                    # chirality
                    continue
            to_remove.add(n_idx)
            mol.nodes[neighbor]['hcount'] = mol.nodes[neighbor].get('hcount', 0) + 1
            if 'rs_isomer' in mol.nodes[neighbor]:
                # Replace the explicit hydrogen index for the implicit one
                mol.nodes[neighbor]['rs_isomer'] = tuple(n if n != n_idx else neighbor for n in mol.nodes[neighbor]['rs_isomer'])
    mol.remove_nodes_from(to_remove)
    for n_idx in mol.nodes:
        if 'hcount' not in mol.nodes[n_idx]:
            mol.nodes[n_idx]['hcount'] = 0


def fill_valence(mol, respect_hcount=True, respect_bond_order=True,
                 max_bond_order=3):
    """
    Sets the attribute 'hcount' on all nodes in `mol` that don't have it yet.
    The value to which it is set is based on the node's 'element', and the
    number of bonds it has. Default valences are determined based on the charge
    and the electron configuration from the PTE. If we can't determine the
    valence (for example for transition metals), then the 'hcount' will be set
    to 0.

    Parameters
    ----------
    mol : nx.Graph
        The molecule whose nodes should get a 'hcount'. Is modified in-place.
    respect_hcount : bool
        If True, don't change the hcount on nodes that already have it set.
    respect_bond_order : bool
        If False, first try to fill the valence by increasing bond orders, and
        add hydrogens after.
    max_bond_order : number
        Only meaningful if respect_bond_order is False. This is the highest
        bond order that will be set.

    Returns
    -------
    None
        `mol` is modified in-place.
    """
    if not respect_bond_order:
        increment_bond_orders(mol, max_bond_order=max_bond_order)
    for n_idx in mol:
        node = mol.nodes[n_idx]
        if ('hcount' in node and respect_hcount) or node.get('element') == 'H':
            continue
        missing = max(bonds_missing(mol, n_idx), 0)
        node['hcount'] = node.get('hcount', 0) + missing


def bonds_missing(mol, node_idx, use_order=True):
    """
    Returns how much the specified node is under valence. If use_order is
    False, treat all bonds as if they are order 1.

    Parameters
    ----------
    mol : nx.Graph
        The molecule.
    node_idx : hashable
        The node to look at. Should be in mol.
    use_order : bool
        If False, treat all bonds as single.

    Returns
    -------
    int
        The number of missing bonds.
    """
    bonds = _bonds(mol, node_idx, use_order)
    bonds += mol.nodes[node_idx].get('hcount', 0)

    val = valence(mol.nodes[node_idx])
    if not val:
        return 0
    val = [v for v in val if v >= bonds] or val[-1:]
    return int(val[0] - bonds)


def valence(atom):
    """
    Returns the valence of the atom. Since some elements can have
    multiple valences, the valence is returned as list.

    Parameters
    ----------
    atom: dict

    Returns
    -------
    list[int]
        The valences for the given atom. Returns an empty list if we don't know.
    """
    element = atom.get("element", '*')
    if element == '*':
        return []

    electrons = PTE[element.capitalize()]['AtomicNumber']
    electrons -= atom.get('charge', 0)
    if not electrons:
        return [0]
    try:
        pte_data = PTE[electrons]
    except KeyError as error:
        raise ValueError(f"Can't figure out an electron configuration for {atom}"
                         f" with {electrons} electrons") from error
    electron_config = pte_data['ElectronConfiguration']
    match = ELECTRON_CONFIG_PATTERN.match(electron_config)
    electron_config = {k: int(v) for k, v in match.groupdict(default=0).items()}
    # No idea how d or f electrons contribute to valency. The last case (p+s==0)
    # is there for Pd, which has config s0p0d10f0...
    if (electron_config['d'] not in (0, 10) or electron_config['f'] not in (0, 14)
            or (not electron_config['s']+electron_config['p'])):
        return []  # No clue what the valence should be or even means

    electrons = electron_config['s'] + electron_config['p']
    # Any sp-electrons we have we distribute over their orbitals. First 1
    # electron in each, then we start making pairs. The resulting valence will
    # be the number of unpaired electrons. Added bonus/complication: electrons
    # in pairs we can excite to higher shells, increasing the number of unpaired
    # electrons
    number_of_orbitals = 4  # 1*s + 3*p
    # Round 1: assign single electrons to orbitals
    single_electrons = min(number_of_orbitals, electrons)
    electrons -= single_electrons
    # Round 2: assign second electrons to orbitals to get electron pairs
    paired_electrons = min(number_of_orbitals, electrons)
    electrons -= paired_electrons
    # Of course all single electrons that got paired are no longer single
    single_electrons -= paired_electrons
    # Figure out how many orbitals are still empty. Maybe useful for later
    empty_orbitals = number_of_orbitals - single_electrons - paired_electrons
    assert electrons == 0, f"There are {electrons=}"

    # Excite any paired electrons to a higher orbital to deal with
    # multi-valency. Naturally, each electron pair consists of 2 electrons,
    # and results in 2 new bonding electrons
    # This is actually kind of wrong since spd hybridization just doesn't
    # happen. Produces the right answer though. Most of the time.
    # https://dx.doi.org/10.1021/ed084p783
    val = [single_electrons + 2*n for n in range(paired_electrons+1)]

    return val


def _bonds(mol, node_idx, use_order=True):
    """
    Returns how many explicit bonds the specified node has. If use_order is
    False, treat all bonds as if they are order 1.

    Parameters
    ----------
    mol : nx.Graph
        The molecule.
    node_idx : hashable
        The node to look at. Should be in mol.
    use_order : bool
        If False, treat all bonds as single.

    Returns
    -------
    int
        The number of bonds.
    """
    if use_order:
        bond_orders = map(operator.itemgetter(2),
                          mol.edges(nbunch=node_idx, data='order', default=1))
        bonds = sum(bond_orders)
    else:
        bonds = len(mol[node_idx])
    return bonds


def has_default_h_count(mol, node_idx, use_order=True):
    """
    Returns whether the hydrogen count for this atom is non-standard.

    Parameters
    ----------
    mol : nx.Graph
        The molecule.
    node_idx : hashable
        The node to look at. Should be in mol.
    use_order : bool
        If False, treat all bonds as single.

    Returns
    -------
    bool
    """
    bonds = _bonds(mol, node_idx, use_order)
    val = valence(mol.nodes[node_idx])
    val = [v for v in val if v >= bonds] or [0]
    hcount = mol.nodes[node_idx].get('hcount', 0)
    return max(val[0] - bonds, 0) == hcount


def _prune_nodes(nodes, mol):
    new_nodes = set()
    for node in nodes:
        # all wild card nodes are eligible
        if mol.nodes[node].get('element', '*') == '*':
            new_nodes.add(node)
            continue
        missing = bonds_missing(mol, node, use_order=True)
        if missing > 0:
            new_nodes.add(node)
    return new_nodes


def correct_aromatic_rings(mol, strict=True):
    """
    Aromaticity is defined here as regions that show delocalization induced
    molecular equivalence (DIME). In other words, regions where alternating
    single and double bonds can interchange without inducing any charges.
    This function will correctly identify and mark aromatic atoms within
    ``mol``, as well as kekulize regions that are marked as aromatic, but do not
    show DIME.

    Parameters
    ----------
    mol : nx.Graph
        The molecule.
    strict : bool
        Whether to raise a SyntaxError if the aromatic region of the molecule
        cannot be kekulized.

    Returns
    -------
    None
        `mol` is modified in-place.

    Raises
    ------
    SyntaxError
        If ``strict`` is True and the aromatic region of the molecule cannot
        be kekulized.
    """
    # set the aromatic attribute to False for all nodes, and swap all aromatic
    # bonds to single
    arom_atoms = {node for node, aromatic in mol.nodes(data='aromatic') if aromatic}
    nx.set_node_attributes(mol, False, 'aromatic')
    for edge in mol.edges:
        if mol.edges[edge].get('order') == 1.5:
            mol.edges[edge]['order'] = 1
    # prune all nodes from molecule that are eligible and have
    # full valency
    ds_graph = _prune_nodes(arom_atoms, mol)

    sub_ds_graph = mol.subgraph(ds_graph)
    max_match = nx.max_weight_matching(sub_ds_graph)
    # we check if a maximum matching exists and
    # if it is perfect. if it is not perfect,
    # this graph originates from a completely invalid
    # smiles and we raise an error
    if not nx.is_perfect_matching(sub_ds_graph, max_match):
        msg = "Your molecule is invalid and cannot be kekulized."
        if strict:
            raise SyntaxError(msg)
        else:
            LOGGER.warning(msg)

    # First we kekulize everything, then we dekekulize the aromatic parts.
    # Otherwise, c12ccccc1[nH]cc2 fails; you need to get rid of some nodes in
    # odd-numbered cycles (the N in this example), otherwise you can end up with
    # a suboptimal (but still maximal) matching
    for edge in max_match:
        mol.edges[edge]['order'] = 2

    dekekulize(mol)


def kekulize(mol):
    """
    Assigns alternating single and double bonds to all aromatic regions in
    ``mol``.

    Arguments
    ---------
    mol : nx.Graph
        The molecule.

    Returns
    -------
    None
        ``mol`` is modified in place.

    Raises
    ------
    ValueError
        If no alternating single and double bonds can be assigned to the
        aromatic region of ``mol``.
    """
    arom_nodes = {n for n in mol if mol.nodes[n].get('aromatic') and mol.nodes[n].get('element', '*') in AROMATIC_ATOMS}
    aromatic_mol = mol.subgraph(arom_nodes)
    matching = nx.max_weight_matching(aromatic_mol)
    if not nx.is_perfect_matching(aromatic_mol, matching):
        raise ValueError('Aromatic region cannot be kekulized.')
    for edge in aromatic_mol.edges:
        if edge in matching or edge[::-1] in matching:
            mol.edges[edge]['order'] = 2
        else:
            mol.edges[edge]['order'] = 1
        mol.nodes[edge[0]]['aromatic'] = False
        mol.nodes[edge[1]]['aromatic'] = False


def dekekulize(mol):
    """
    Finds all cycles in ``mol`` that consist of alternating single and double
    bonds, and marks them as aromatic.

    Arguments
    ---------
    mol : nx.Graph
        The molecule.

    Returns
    -------
    None
        ``mol`` is modified in place.
    """
    # The idea is the following:
    # 1) only nodes in cycles can be aromatic (in the sense that they show DIME)
    # 2) only nodes that have a double bond can be aromatic
    # Given those two requirements, the aromatic system is spanned by the
    # maximal matching
    cycles = nx.cycle_basis(mol)
    cycles_nodes = {n for cycle in cycles for n in cycle}
    double_bond_atoms = {n for edge in mol.edges for n in edge if mol.edges[edge].get('order', 1) == 2}
    maybe_aromatic = double_bond_atoms & cycles_nodes
    matching = nx.max_weight_matching(mol.subgraph(maybe_aromatic))
    aromatic_nodes = {n for e in matching for n in e}

    # The matching may extend into not fully aromatic cycles.
    for cycle in cycles:
        if set(cycle) <= aromatic_nodes:
            if not all(mol.nodes[n].get('element', '*') in AROMATIC_ATOMS for n in cycle):
                continue
            # Party!
            for node in cycle:
                mol.nodes[node]['aromatic'] = True
            for edge in zip(cycle, cycle[1:] + [cycle[0]]):
                mol.edges[edge]['order'] = 1.5


def increment_bond_orders(molecule, max_bond_order=3):
    """
    Increments bond orders up to what the atom's valence allows.

    Parameters
    ----------
    molecule : nx.Graph
        The molecule to process.
    max_bond_order : number
        The highest bond order allowed to make.

    Returns
    -------
    None
        molecule is modified in-place.
    """
    # Gather the number of open spots for all atoms beforehand, since some
    # might have multiple oxidation states (e.g. S). We don't want to change
    # oxidation state halfway through for some funny reason. It shouldn't be
    # necessary, but it can't hurt.
    missing_bonds = {}
    for idx in molecule:
        missing_bonds[idx] = max(bonds_missing(molecule, idx), 0)

    for idx, jdx in molecule.edges:
        missing_idx = missing_bonds[idx]
        missing_jdx = missing_bonds[jdx]
        edge_missing = min(missing_idx, missing_jdx)
        current_order = molecule.edges[idx, jdx].get("order", 1)
        if current_order == 1.5 or current_order >= max_bond_order:
            continue
        new_order = edge_missing + current_order
        new_order = min(new_order, max_bond_order)
        molecule.edges[idx, jdx]['order'] = new_order
        missing_bonds[idx] -= edge_missing
        missing_bonds[jdx] -= edge_missing


def _mark_chiral_atoms(molecule):
    """
    For all nodes tagged as chiral, figure out the three
    substituents and annotate the node with a tuple that
    has the order in which to rotate. This essentially
    corresponds to the definition of an improper dihedral
    angle centered on the chiral atom.
    """
    chiral_nodes = nx.get_node_attributes(molecule, 'rs_isomer')
    for node, (direction, rings) in chiral_nodes.items():
        # first the ring atoms in order
        # that they were connected then the
        # other neighboring atoms in order
        bonded_neighbours = sorted(molecule[node])
        neighbours = list(rings)
        for neighbour in bonded_neighbours:
            if neighbour not in neighbours:
                neighbours.append(neighbour)

        if molecule.nodes[node].get('hcount'):
            # We have an implicit hydrogen. This is the first atom to determine
            # the rotation. As a placeholder, we'll put the index of the parent
            # node.
            neighbours.insert(1, node)

        if len(neighbours) != 4:
            # FIXME Tetrahedral Allene-like Systems such as `NC(Br)=[C@]=C(O)C`
            msg = (f"Chiral node {node} has {len(neighbours)} neighbors, which "
                    "is different than the four expected for tetrahedral "
                    "chirality.")
            raise ValueError(msg)
        # the default is anti-clockwise sorting indicated by '@'
        # in this case the nodes are sorted with increasing
        # node index; however @@ means clockwise and the
        # order of nodes is reversed (i.e. with decreasing index)
        if direction == '@@':
            neighbours = [neighbours[0],  neighbours[1], neighbours[3], neighbours[2]]
        molecule.nodes[node]['rs_isomer'] = tuple(neighbours)


def _annotate_ez_isomers(molecule, ez_pairs):
    for first, second in ez_pairs:
        ligand_first, anchor_first, ez_first = first
        ligand_second, anchor_second, ez_second = second
        ez_isomer = None
        # here come all the ridiculous cases of EZ in smiles
        if ligand_first < anchor_first:
            # case 1 F/C=C/F is trans
            if ez_first == '/' and ez_second == '/':
                ez_isomer = 'trans'
            # case 2 F\C=C/F
            elif ez_first == '\\' and ez_second == '/':
                ez_isomer = 'cis'
            # case 3 F/C=C\F
            elif ez_first == '/' and ez_second == '\\':
                ez_isomer = 'cis'
            # case 4 F\C=C\F
            elif ez_first == '\\' and ez_second == '\\':  # pragma: no branch
                ez_isomer = 'trans'
        # as in C(/F)=C
        elif ligand_first > anchor_first:  # pragma: no branch
            # case 5 C(\F)=C/F is trans
            if ez_first == '\\' and ez_second == '/':
                ez_isomer = 'trans'
            # case 6 C(/F)=C/F
            elif ez_first == '/' and ez_second == '/':
                ez_isomer = 'cis'
            # case 7 C(/F)=C\F
            elif ez_first == '/' and ez_second == '\\':
                ez_isomer = 'trans'
            # case 8 C(\F)=C\F
            elif ez_first == '\\' and ez_second == '\\':  # pragma: no branch
                ez_isomer = 'cis'
        assert ez_isomer is not None
        # annotate ligands
        molecule.nodes[ligand_first]['ez_isomer'] = (ligand_first,
                                                     anchor_first,
                                                     anchor_second,
                                                     ligand_second,
                                                     ez_isomer)
        molecule.nodes[ligand_second]['ez_isomer'] = (ligand_second,
                                                      anchor_second,
                                                      anchor_first,
                                                      ligand_first,
                                                      ez_isomer)

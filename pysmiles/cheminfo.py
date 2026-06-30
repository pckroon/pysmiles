"""
Some simple cheminformatics tools.
"""
from collections import defaultdict
import numpy as np
import networkx as nx
from .smiles_helper import PTE, remove_explicit_hydrogens

def tanimoto_similarity(binarr1, binarr2):
    """
    Compute Tanimoto similarity between two binary arrays.

    Parameters
    ----------
    binarr1: np.array
    binarr2: np.array

    Returns
    -------
    float
        the tanimoto similarity
    """
    intersection = np.sum(binarr1 & binarr2)
    union = np.sum(binarr1) + np.sum(binarr2) - intersection
    return intersection / union if union != 0 else 0.0

def features_to_binary(feature_list, size=1024):
    """
    Converte a list of features to a binary array of
    given `size`.

    Parameters
    ----------
    feature_list: abc.iteratable[int]
        a iteratable with integer numbers
    size: int
        size of the binary array; default=1024

    Returns
    -------
    np.array
    """
    fp_array = np.zeros(size, dtype=int)
    nonzero = np.array(feature_list, dtype=int) % size
    fp_array[nonzero] = 1
    return fp_array

def _aggregate_features(graph, cycles=4, label='fp_invariant'):
    """
    Given atomic invariants in a molecule stored under the `label`
    keyword, iteratively aggregate those features for a given number
    of `cycles`. It returns a feature list, which is composed of the
    hashes of the aggregated features.

    Parameters
    ----------
    graph: networkx.Graph
    cycles: int
        number of cycles to run; default is 4
    label: abc.hashable
        the keyword under which the inital invariants are stored

    Returns
    -------
    list[int]
        list of the hashes of unique features
    """
    feature_list = list(nx.get_node_attributes(graph, label).values())
    substructures = set(nx.get_node_attributes(graph, 'subgraph').values())
    for idx in range(0, cycles):
        for node in graph.nodes:
            subgraph = set(graph.nodes[node]['subgraph'])
            feat_list_iter = [(idx, graph.nodes[node][label])]
            for neigh in graph.neighbors(node):
                order = graph.edges[(node, neigh)].get('order', 1)
                feat_list_iter.append((order, graph.nodes[neigh][label]))
                subgraph = subgraph | graph.nodes[neigh]['subgraph']
            graph.nodes[node]['subgraph'] = subgraph
            feat_list_iter = sorted(feat_list_iter)
            feat_tuple_iter = tuple(feat for sub_feat in feat_list_iter for feat in sub_feat)
            new_hash = hash(feat_tuple_iter)
            graph.nodes[node]['fp_invariant'] = new_hash
            if frozenset(subgraph) not in substructures:
                feature_list.append(new_hash)
                substructures.add(frozenset(subgraph))
    return feature_list
            
def _assign_invariants(graph, node_feats):
    """
    Assign inital atom invariants from the features
    provided in `node_feats`. Invarients get annotated
    under the `fp_invarient` keyword. It also stores all
    edges under the keyword `subraph`.

    Parameters
    ----------
    graph: networkx.graph
        molecule graph
    node_feats: list[abc.hashable]
        features used to assgin the invariants
    """
    for node in graph.nodes:
        feat_list = []
        for feat in node_feats:
            feat_list.append(graph.nodes[node][feat])
        _hash = hash(tuple(feat_list))
        graph.nodes[node]['fp_invariant'] = _hash
        graph.nodes[node]['subgraph'] = frozenset(graph.edges(node))

def _reduced_valency(graph):
    """
    Compute the reduced valency for fingerprint calculations.
    The reduced valency is given as the valency minus the
    hydrogen atoms.
    """
    red_valency = defaultdict(int)
    for e1, e2, order in graph.edges(data='order'):
        red_valency[e1] += order
        red_valency[e2] += order
    nx.set_node_attributes(graph, red_valency, 'red_valency')
    return graph

def annotate_rings(graph):
    """
    Find all nodes that are part of a ring and
    set the ring attribute to True otherwise 
    it is set to False.

    Parameters
    ----------
    graph: :class:`networkx.Graph`

    Returns
    -------
    :class:`networkx.Graph`
    """
    nx.set_node_attributes(graph, 0, 'ring')
    for cycle in nx.cycle_basis(graph):
        for node in cycle:
            graph.nodes[node]['ring'] = 1
    return graph

def _prepare_atomic_molecule(graph):
    """
    Set some attributes required for fingerprint calculations.
    """
    remove_explicit_hydrogens(graph)
    graph = _reduced_valency(graph)
    annotate_rings(graph)
    for node in graph.nodes:
        element = graph.nodes[node]['element']
        graph.nodes[node]['atomic_num'] = PTE[element.capitalize()]['AtomicNumber']
        graph.nodes[node]['mass'] = PTE[element.capitalize()]['AtomicMass']
        graph.nodes[node]['degree'] = graph.degree(node)
        
    return graph

def extended_connectivity_fp(graph,
                             node_feats=['red_valency',
                                         'hcount',
                                         'degree',
                                         'atomic_num',
                                         'mass',
                                         'charge']):
    """
    Computes an extended connectivity fingerprint from a networkx graph.
    Follows the workflow and explanation outliend in:
    https://chemicbook.com/2021/03/25/a-beginners-guide-for-understanding-extended-connectivity-fingerprints.html

    Parameters
    ----------
    graph: networkx.graph
        the molecule graph
    node_feats: list[abc.hashable]
        features to use for invariant computations

    """
    graph = _prepare_atomic_molecule(graph)
    _assign_invariants(graph, node_feats)
    feature_list = _aggregate_features(graph, cycles=4)
    fp_array = features_to_binary(feature_list)
    return fp_array

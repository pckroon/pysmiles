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
Contains a unittest.TestCase subclass that has an `assertEqualGraphs` method.
"""
from functools import partial
import operator
import unittest

import networkx as nx


def make_mol(node_data, edge_data):
    mol = nx.Graph()
    mol.add_nodes_from(node_data)
    mol.add_edges_from(edge_data)
    return mol


def _equal_dicts(dict1, dict2, excluded_keys=[]):
    dict1_keys = set(dict1.keys()) - set(excluded_keys)
    dict2_keys = set(dict2.keys()) - set(excluded_keys)
    return dict1_keys == dict2_keys and {k1: dict1[k1] for k1 in dict1_keys} == {k2: dict2[k2] for k2 in dict1_keys}


def assertEqualGraphs(graph1, graph2, excluded_node_attrs=[], excluded_edge_attrs=[]):  # pylint: disable=invalid-name
    """
    Asserts that `graph1` and `graph2` are equal up to isomorphism.
    """
    out = nx.is_isomorphic(graph1, graph2,
                           node_match=partial(_equal_dicts, excluded_keys=excluded_node_attrs),
                           edge_match=partial(_equal_dicts, excluded_keys=excluded_edge_attrs))
    if out:
        return

    matcher = nx.isomorphism.GraphMatcher(graph1, graph2)  # pylint: disable=no-member
    matches = list(matcher.isomorphisms_iter())
    if not matches:
        raise AssertionError("Graphs not isomorphic")

    # Now for the hard part. Find the isomorphism that matches best, and
    # then figure out where the attributes differ.
    scores = []
    for m_idx, match in enumerate(matches):
        score = 0
        for node_idx, node_jdx in match.items():
            node1 = graph1.nodes[node_idx]
            node2 = graph2.nodes[node_idx]
            for attr in (set(node1) | set(node2)) - set(excluded_node_attrs):
                if node1.get(attr) == node2.get(attr):
                    score += 1

        for idx1, jdx1 in graph1.edges:
            idx2, jdx2 = match[idx1], match[jdx1]
            edge1 = graph1.edges[idx1, jdx1]
            edge2 = graph2.edges[idx2, jdx2]
            for attr in (set(edge1) | set(edge2)) - set(excluded_edge_attrs):
                if edge1.get(attr) == edge2.get(attr):
                    score += 1
        scores.append((score, m_idx))

    # This one should be the best.
    _, m_idx = max(scores)
    match = matches[m_idx]
    assert dict(graph1.nodes) == {idx: graph2.nodes[match[idx]] for idx in graph1}
    bad_edges = []
    for idx1, jdx1 in graph1.edges:
        idx2, jdx2 = match[idx1], match[jdx1]
        edge1 = graph1.edges[idx1, jdx1]
        edge2 = graph2.edges[idx2, jdx2]
        if edge1 != edge2:
            bad_edges.append((idx1, jdx1, edge1, edge2))
    msg = ''
    fmt = 'Edge ({}, {}) between {} and {} is not equal: {} is not {};'
    for idx, jdx, edge1, edge2 in bad_edges:
        msg += fmt.format(idx, jdx, graph1.nodes[idx], graph2.nodes[jdx], edge1, edge2)
    raise AssertionError(msg)

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
import operator
import unittest

import networkx as nx


class GraphTest(unittest.TestCase):
    def assertEqualGraphs(self, graph1, graph2):
        out = nx.is_isomorphic(graph1, graph2,
                               node_match=operator.eq, edge_match=operator.eq)
        if out:
            return

        matcher = nx.isomorphism.GraphMatcher(graph1, graph2)
        matches = list(matcher.isomorphisms_iter())
        if not matches:
            raise AssertionError("Graphs not isomorphic")

        # Now for the hard part. Find the isomorphism that matches best, and
        # then figure out where the attributes differ.
        scores = []
        for m_idx, match in enumerate(matches):
            score = 0
            for node_idx, node_jdx in match.items():
                if graph1.nodes[node_idx] == graph2.nodes[node_jdx]:
                    score += 1
            for idx1, jdx1 in graph1.edges:
                idx2, jdx2 = match[idx1], match[jdx1]
                edge1 = graph1.edges[idx1, jdx1]
                edge2 = graph2.edges[idx2, jdx2]
                if edge1 == edge2:
                    score += 1
            scores.append((score, m_idx))

        # This one should be the best.
        _, m_idx = max(scores)
        match = matches[m_idx]
        for node_idx, node_jdx in match.items():
            self.assertEqual(graph1.nodes[node_idx], graph2.nodes[node_jdx])
        for idx1, jdx1 in graph1.edges:
            idx2, jdx2 = match[idx1], match[jdx1]
            edge1 = graph1.edges[idx1, jdx1]
            edge2 = graph2.edges[idx2, jdx2]
            if edge1 != edge2:
                fmt = 'Edge between {} and {} is not equal: {} is not {}'
                raise AssertionError(fmt.format(graph1.nodes[idx1],
                                                graph1.nodes[jdx1],
                                                edge1, edge2))

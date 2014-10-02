#include "transaction_minimizer.h"

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <limits>
#include <map>
#include <set>
#include <vector>

#include <boost/graph/copy.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/topological_sort.hpp>

using namespace boost;

using std::map;
using std::pair;
using std::set;
using std::vector;

namespace {
void removeCircularPayments(Graph *payments);
void computePaymentsSpan(const Graph &payments, Graph *paymentsSpan,
                         map<pair<Vertex, Vertex>, set<Vertex> > *routingTable);
void computePaymentAmounts(
        const Graph &payments, Graph *paymentsSpan,
        const map<pair<Vertex, Vertex>, set<Vertex> > &routingTable);

void verifyNodeGain(const Graph &left, const Graph &right);
} // namespace

void createMinimumTransactionsGraph(const Graph &payments,
                                    Graph *paymentsSpan) {
    // Early-abort for empty graphs.
    if (num_edges(payments) == 0) {
        copy_graph(payments, *paymentsSpan);
        return;
    }

    // Clones the payments graph into the payments DAG.
    Graph paymentsDAG;
    copy_graph(payments, paymentsDAG);

    // Removes circular payments by reducing the payment amounts on the circle
    // by the smallest amount and removing the edges that result in an amount of
    // zero.
    removeCircularPayments(&paymentsDAG);
    
    // Routing table that returns the next hops one can take at given source
    // node to reach given destination node.
    map<pair<Vertex, Vertex>, set<Vertex> > routingTable;

    // Computes the spanning tree of the payments graph s.t. the number of edges
    // is minimized while still keeping the same connected components in the
    // graph.
    computePaymentsSpan(paymentsDAG, paymentsSpan, &routingTable);

    // Computes the final payments along the edges in the payments span.
    computePaymentAmounts(paymentsDAG, paymentsSpan, routingTable);

    // Verifies the payment span.
    verifyNodeGain(payments, *paymentsSpan);
}

namespace {

class DFSVisitor : public default_dfs_visitor {
public:
    explicit DFSVisitor(bool *circleFound, vector<Vertex> *circle,
                        vector<bool> *visitedVertices)
        : m_circleFound(circleFound)
        , m_circle(circle)
        , m_visitedVertices(visitedVertices)
    {
        *m_circleFound = false;
        m_circle->clear();
    }

    void start_vertex(Vertex u, const Graph &) {
        if (*m_circleFound) return;

        assert(m_circle->empty());
        m_circle->push_back(u);
        m_verticesOnPath.insert(u);
        (*m_visitedVertices)[u] = true;
    }

    void examine_edge(Edge e, const Graph &g) {
        if (*m_circleFound) return;

        const Vertex to = target(e, g);
        if (m_verticesOnPath.find(to) != m_verticesOnPath.end()) {
            // We already have the target of this edge on the search path
            // -> we have found a circle.
            *m_circleFound = true;

            // Removes the vertices that are not part of the circle.
            while (m_circle->front() != to) {
                m_circle->erase(m_circle->begin());
            }
        }
    }

    void tree_edge(Edge e, const Graph &g) {
        if (*m_circleFound) return;

        // Adds the target vertex to the path.
        const Vertex to = target(e, g);
        m_circle->push_back(to);
        m_verticesOnPath.insert(to);
        (*m_visitedVertices)[to] = true;
    }

    void finish_vertex(Vertex u, const Graph &) {
        if (*m_circleFound) return;

        // Removes the vertex from the path.
        assert(m_circle->back() == u);
        m_circle->pop_back();
        m_verticesOnPath.erase(u);
    }

private:
    bool * const m_circleFound;
    vector<Vertex> * const m_circle;
    vector<bool> * const m_visitedVertices;

    set<Vertex> m_verticesOnPath;
};

void findCircle(const Graph &graph, bool *hasCircle, vector<Vertex> *circle) {
    vector<bool> visitedVertices(num_vertices(graph), false);

    // Searches for the sources of the graph.
    vector<Vertex> sources;
    Graph::vertex_iterator vit, vend;
    for (tie(vit, vend) = vertices(graph); vit != vend; ++vit) {
        if (in_degree(*vit, graph) == 0) {
            sources.push_back(*vit);
        }
    }
    if (sources.empty()) {
        sources.push_back(0); // just insert the first one.
    }

    // Runs over the sources of the graph.
    bool nonSourcesProcessed = false;
    while (!sources.empty()) {
        const Vertex source = sources.back();
        sources.pop_back();

        if (visitedVertices[source]) {
            continue;
        }

        DFSVisitor dfsVisitor(hasCircle, circle, &visitedVertices);
        
        vector<default_color_type> colorMap(num_vertices(graph));
        depth_first_visit(graph, source, dfsVisitor, &colorMap[0]);

        if (*hasCircle) {
            break;
        }

        if (sources.empty() && !nonSourcesProcessed) {
            // Checks if there are still some un-processed vertices.
            nonSourcesProcessed = true;

            for (tie(vit, vend) = vertices(graph); vit != vend; ++vit) {
                if (!visitedVertices[*vit]) {
                    sources.push_back(*vit);
                }
            }
        }
    }
}

void removeCircularPayments(Graph *payments) {
    do {
        bool hasCircle;
        vector<Vertex> circle;
        findCircle(*payments, &hasCircle, &circle);

        if (!hasCircle) {
            // No more circles -> we are done.
            break;
        }

        // Searches for the minimum payment along the circle.
        double minimumAmount = std::numeric_limits<double>::max();
        Vertex lastVertex = circle.back();
        for (vector<Vertex>::const_iterator it = circle.begin();
             it != circle.end(); ++it) {
            double amount = get(edge_capacity, *payments,
                                edge(lastVertex, *it, *payments).first);
            minimumAmount = std::min(minimumAmount, amount);

            lastVertex = *it;
        }

        // Reduces the amount on the edges and removes the edge with that
        // minimum amount.
        lastVertex = circle.back();
        for (vector<Vertex>::const_iterator it = circle.begin();
             it != circle.end(); ++it) {
            Edge e = edge(lastVertex, *it, *payments).first;

            double amount = get(edge_capacity, *payments, e);
            if (amount == minimumAmount) {
                // Removes the edge.
                remove_edge(e, *payments);
            } else {
                // Reduces the edge payment amount.
                put(edge_capacity, *payments, e, amount - minimumAmount);
            }

            lastVertex = *it;
        }
    } while (true);

#ifdef DEBUG
    // Prints the DAG.
    printf("Final DAG: \n");
    Graph::edge_iterator eit, eend;
    for (tie(eit, eend) = edges(*payments); eit != eend; ++eit) {
        const Vertex &u = source(*eit, *payments);
        const Vertex &v = target(*eit, *payments);
        const double &amount = get(edge_capacity, *payments, *eit);
        printf("%s --%.2f--> %s\n",
               get(boost::vertex_name, *payments, u).c_str(), amount,
               get(boost::vertex_name, *payments, v).c_str());
    }
    printf("\n");
#endif
}

/**
 * @brief Searches for the minimum set of edges that keep the DAG given in
 *        payments connected and writes this span into paymentsSpan. Moreover,
 *        a routing table is generated that tells for each (source, destination)
 *        touple, which next hops can be used.
 *
 * The search for the span is done by starting at the sinks which are marked as
 * reachable nodes. Now, in each round, the sinks are processed. Their reachable
 * nodes are set as the union of the reachable nodes of its children. Moreover,
 * the routing table gets filled by telling which nodes can be reached by using
 * which edge. Edges that only reach nodes that are also reachable by other
 * edges are removed.
 *
 * @param payments The DAG that contains the payments among all people
 * @param paymentsSpan (Out) The span with minimum edge size that keeps the
 *                           components in payments connected.
 * @param routingTable (Out) (source, dest) -> {next hops}. Multiple next hops
 *                           can appear if that edge is required anyways. One
 *                           can then e.g. distribute payments along these
 *                           paths to reduce the maximum transaction amount.
 */
void computePaymentsSpan(
        const Graph &payments, Graph *paymentsSpan,
        map<pair<Vertex, Vertex>, set<Vertex> > *routingTable) {
    const int n = num_vertices(payments);

    // Clones the payments DAG into the payments span.
    copy_graph(payments, *paymentsSpan);

    // The list of nodes that are reachable through a given node.
    map<Vertex, set<Vertex> > reachableNodes;

    // Searches for the sinks.
    queue<Vertex> sinks;
    for (int i = 0; i < n; ++i) {
        if (out_degree(i, payments) == 0) {
            reachableNodes[i].insert(i); // Adds itself as a reachable node.
            sinks.push(i);
        }
    }

    // Runs over the sinks and marks all incoming edges to it as useable. If by
    // that all out-edges of a node are marked as useable, it is processed and
    // added to the sinks of the next round.
    vector<Vertex> numUsableEdges(n, 0);
    do {
        queue<Vertex> nextSinks;
        while (!sinks.empty()) {
            const Vertex &t = sinks.front();
            sinks.pop();

            // Processes this sink node.
            // Adds itself as a reachable node.
            reachableNodes[t].insert(t);

            // Computes the reachable nodes of s by uniting the
            // reachable nodes of its children.
            vector<int> numPathsToNode(n, 0);
            Graph::out_edge_iterator eit2, eend2;
            for (tie(eit2, eend2) = out_edges(t, payments); eit2 != eend2;
                 ++eit2) {
                const Vertex &t2 = target(*eit2, payments);
                set<Vertex> &reachableNodesOfT2 = reachableNodes[t2];

                reachableNodes[t].insert(reachableNodesOfT2.begin(),
                                         reachableNodesOfT2.end());

                // Increases the path count for each reachable node of
                // t2.
                for (set<Vertex>::const_iterator it =
                     reachableNodesOfT2.begin(); it != reachableNodesOfT2.end();
                     ++it) {
                    ++numPathsToNode[*it];
                }
            }

            // Computes the routing table for node s. Every edge that
            // only reaches nodes that are also reachable by others
            // (numPathsToNode[i] > 1) are not used.
            for (tie(eit2, eend2) = out_edges(t, payments); eit2 != eend2;
                 ++eit2) {
                const Vertex &t2 = target(*eit2, payments);
                set<Vertex> &reachableNodesOfT2 = reachableNodes[t2];

                // Checks if the current edge is required.
                bool isRequiredEdge = false;
                for (set<Vertex>::const_iterator it =
                     reachableNodesOfT2.begin();
                     it != reachableNodesOfT2.end() && !isRequiredEdge; ++it) {
                    isRequiredEdge |= numPathsToNode[*it] == 1;
                }

                // Adds the edge to the routing table if it is required.
                if (isRequiredEdge) {
                    for (set<Vertex>::const_iterator it =
                         reachableNodesOfT2.begin();
                         it != reachableNodesOfT2.end(); ++it) {
                        // The next hop from t to *it is t2.
                        (*routingTable)[std::make_pair(t, *it)].insert(t2);
                    }
                } else {
                    // Removes the unused edge from the payments span.
                    remove_edge(t, t2, *paymentsSpan);
                }
            }

            // Marks the in-edges to t as useable.
            Graph::in_edge_iterator eit, eend;
            for (tie(eit, eend) = in_edges(t, payments); eit != eend; ++eit) {
                const Vertex &s = source(*eit, payments);
                ++numUsableEdges[s];

                if (out_degree(s, payments) == numUsableEdges[s]) {
                    // All out-edges of s are useable -> adds s to the
                    // next-round sinks.
                    nextSinks.push(s);
                }
            }
        }
        swap(sinks, nextSinks);
    } while (!sinks.empty());

#ifdef DEBUG
    printf("Payments span:\n");
    for (int i = 0; i < n; ++i) {
        printf("Vertex '%s'", get(vertex_name, *paymentsSpan, i).c_str());
        printf(" (rn: ");
        for (set<Vertex>::const_iterator it = reachableNodes[i].begin();
             it != reachableNodes[i].end(); ++it) {
            printf("%s%s",
                   it != reachableNodes[i].begin() ? ", " : "",
                   get(vertex_name, *paymentsSpan, *it).c_str());
        }
        printf(")\n");

        if (out_degree(i, *paymentsSpan) > 0) {
            printf ("  -> ");
            Graph::out_edge_iterator oeit, oeend;
            bool first = true;
            for (tie(oeit, oeend) = out_edges(i, *paymentsSpan); oeit != oeend;
                 ++oeit) {
                const Vertex &t = target(*oeit, *paymentsSpan);
                printf("%s%s",
                       get(vertex_name, *paymentsSpan, t).c_str(),
                       first ? "" : ", ");
                first = false;
            }
            printf("\n");
        }
    }
    printf("\n");
#endif
}

void computePaymentAmounts(
        const Graph &payments, Graph *paymentsSpan,
        const map<pair<Vertex, Vertex>, set<Vertex> > &routingTable) {
    // A per-vertex list that stores the amount that needs to be forwarded to
    // which destination node.
    map<Vertex, map<Vertex, double> > forwardPayments;

    // Resets the payment entries of all edges.
    Graph::edge_iterator eit, eend;
    for (tie(eit, eend) = edges(*paymentsSpan); eit != eend; ++eit) {
        put(edge_capacity, *paymentsSpan, *eit, 0);
    }

    // Sorts the vertices. If there is an edge (u, v) then v comes before u in
    // the ordering.
    vector<Vertex> sortedVertices;
    topological_sort(*paymentsSpan, std::back_inserter(sortedVertices));

    // Runs from the sources to the sinks over the payments span.
    for (vector<Vertex>::const_reverse_iterator vit = sortedVertices.rbegin();
         vit != sortedVertices.rend(); ++vit) {
        const Vertex from = *vit;

        map<Vertex, double> &paymentList = forwardPayments[from];

        // Adds all direct payments from the current node to its forward-list.
        Graph::out_edge_iterator oeit, oeend;
        for (tie(oeit, oeend) = out_edges(from, payments); oeit != oeend;
             ++oeit) {
            paymentList[target(*oeit, payments)] +=
                get(edge_capacity, payments, *oeit);
        }

        // Forwards the payments one hop by looking up the next hops for each
        // destination node and increasing their forwardPayments entries
        // accordingly.
        for (map<Vertex, double>::const_iterator it = paymentList.begin();
             it != paymentList.end(); ++it) {
            const Vertex &to = it->first;
            const double &amount = it->second;

            if (to == from) {
                // This is the payment entry to us. We keep that money.
                continue;
            }

            // Looks up the hops in the routing table.
            const set<Vertex> &hops = routingTable.at(std::make_pair(from, to));
            assert(!hops.empty());

            // Splits the amount evenly between the multiple paths.
            double amountPerPath = floor(amount / hops.size());
            for (set<Vertex>::const_iterator hopIt = hops.begin();
                 hopIt != hops.end(); ++hopIt) {
                const Vertex &hop = *hopIt;

                // We split the amount per path into integer parts. If we are
                // the first path, we cover for eventual rounding errors.
                double pathAmount = hopIt != hops.begin()
                        ? amountPerPath
                        : amount - (hops.size() - 1) * amountPerPath;

                // Increases the amount we send along the edge to the hop node.
                double &edgeAmount = get(edge_capacity, *paymentsSpan,
                                         edge(from, hop, *paymentsSpan).first);
                edgeAmount += pathAmount;

                // Increases the forward amount at the hop for the final target.
                forwardPayments[hop][to] += pathAmount;
            }
        }
    }
}

vector<double> computeNodeGain(const Graph &payments) {
    vector<double> nodeGain(num_vertices(payments), 0.0d);

    Graph::edge_iterator eit, eend;
    for (tie(eit, eend) = edges(payments); eit != eend; ++eit) {
        const Vertex from = source(*eit, payments);
        const Vertex to = target(*eit, payments);
        double amount = get(edge_capacity, payments, *eit);

        nodeGain[from] -= amount;
        nodeGain[to] += amount;
    }

    return nodeGain;
}

void verifyNodeGain(const Graph &left, const Graph &right) {
    const vector<double> &leftGain = computeNodeGain(left);
    const vector<double> &rightGain = computeNodeGain(right);

    bool error = false;
    for (size_t i = 0; i < leftGain.size(); ++i) {
        // Checks if the difference is bigger than one cent.
        if (std::fabs(leftGain[i] - rightGain[i]) > 0.01) {
            error = true;
            fprintf(stderr, "Gain-error at node '%s' (%.2f vs %.2f).\n",
                    get(boost::vertex_name, left, i).c_str(),
                    leftGain[i], rightGain[i]);
        }
    }

    if (error) {
        fprintf(stderr,
                "Please report this error with the payments you defined (if "
                "possible) at https://github.com/saemy/transaction_minimizer."
                "\n");
        std::exit(1);
    }
}

} // namespace

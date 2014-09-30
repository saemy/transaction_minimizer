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
#include <boost/numeric/ublas/matrix.hpp>

using namespace boost;

using boost::numeric::ublas::matrix;
using boost::numeric::ublas::zero_matrix;
using std::map;
using std::set;
using std::vector;

namespace {
void removeCircularPayments(Graph *payments);
void computePaymentsSpan(const Graph &payments, Graph *paymentsSpan,
                         map<Edge, set<Vertex> > *routingTable);
void computePaymentAmounts(const Graph &payments, Graph *paymentsSpan,
                           const map<Edge, set<Vertex> > &routingTable);
} // namespace

void reduceToMinimumTransactions(Graph *paymentsSpan) {
    // Early-abort for empty graphs.
    if (num_edges(*paymentsSpan) == 0) {
        return;
    }

    // Clones the input graph's vertices into the payments graph and removes the
    // edges from the input graph.
    Graph payments;
    copy_graph(*paymentsSpan, payments);
    while (num_edges(*paymentsSpan)) {
        remove_edge(*edges(*paymentsSpan).first, *paymentsSpan);
    }

    // Per-edge routing table that lists all the vertices that are reached when
    // using this edge in the payments span.
    map<Edge, set<Vertex> > routingTable;

    // Removes circular payments by reducing the payment amounts on the circle
    // by the smallest amount and removing the edges that result in an amount of
    // zero.
    removeCircularPayments(&payments);
    
    // Computes the spanning tree of the payments graph s.t. every node which
    // was connected to an other one still is, however, only over exactly one
    // path.
    computePaymentsSpan(payments, paymentsSpan, &routingTable);

    // Computes the final payments along the edges in the payments span.
    computePaymentAmounts(payments, paymentsSpan, routingTable);
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
}

// Assumes the incoming graph does not contain circles.
void computePaymentsSpan(const Graph &payments, Graph *paymentsSpan,
                         map<Edge, set<Vertex> > *routingTable) {
    // The number of vertices in the payments graph.
    const int n = num_vertices(payments);

    // Per vertex bitmap that flags nodes that are not reaching other nodes.
    vector<bool> longTimeNonReacher(n, true);

    // Constructs the adjacency matrix of the payments graph.
    matrix<bool> adjacencyMatrix = zero_matrix<bool>(n, n);

    Graph::edge_iterator eit, eend;
    for (tie(eit, eend) = edges(payments); eit != eend; ++eit) {
        const Vertex from = source(*eit, payments);
        const Vertex to = target(*eit, payments);

        adjacencyMatrix(from, to) = true;
        longTimeNonReacher[from] = false;
    }

    // Adds the sinks to the list of the unconnected non-reachers.
    set<Vertex> unconnectedNonReacher;
    for (int i = 0; i < n; ++i) {
        if (longTimeNonReacher[i]) {
            unconnectedNonReacher.insert(i);
        }
    }

#if DEBUG
        // prints the hop matrix.
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                std::cout << adjacencyMatrix(i, j) << ' ';
            }
            std::cout << '\n';
        }
        std::cout << '\n';
#endif

    // Computes the adjacency matrices for increasing hop count.
    matrix<bool> lastHopMatrix = adjacencyMatrix;
    do {

        // Calculates the next hop adjacency matrix.
        matrix<bool> hopMatrix = prod(lastHopMatrix, adjacencyMatrix);
        lastHopMatrix = hopMatrix;

#if DEBUG
        // prints the hop matrix.
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                std::cout << hopMatrix(i, j) << ' ';
            }
            std::cout << '\n';
        }
        std::cout << '\n';
#endif

        set<Vertex> nextUnconnectedNonReacher = unconnectedNonReacher;

        // Checks if there are still connections in the matrix and searches for
        // nodes that do not reach anything for the first time.
        bool hasConnection = false;
        for (int i = 0; i < n; ++i) {
            if (longTimeNonReacher[i]) {
                continue;
            }

            // Checks if the node j has a connection to another one.
            bool reacher = false;
            for (int j = 0; j < n && !reacher; ++j) {
                reacher |= hopMatrix(i, j);
            }

            // If this node does not reach anything (for the first time) -> adds
            // a link to the non-reacher of the last round (if existing in the
            // original graph).
            if (!reacher) {
                longTimeNonReacher[i] = true;

                for (set<Vertex>::const_iterator it =
                         unconnectedNonReacher.begin();
                     it != unconnectedNonReacher.end(); ++it) {
                    if (edge(i, *it, payments).second) {
                        assert((unsigned) i != *it);

                        // They are connected in the payments graph -> connects
                        // them in the span.
                        Edge e;
                        tie(e, tuples::ignore) =
                            add_edge(i, *it, *paymentsSpan);

                        nextUnconnectedNonReacher.erase(*it);

                        // Constructs the routing table entry for this edge.
                        set<Vertex> &rte = (*routingTable)[e];
                        rte.clear();
                        rte.insert(*it);
                        Graph::out_edge_iterator oeit, oeend;
                        for (tie(oeit, oeend) = out_edges(*it, *paymentsSpan);
                             oeit != oeend; ++oeit) {
                            const set<Vertex> &destRte = (*routingTable)[*oeit];

                            // Combines the rte with the outgoing rtes of the to
                            // vertex.
                            set<Vertex> tmpRte;
                            std::set_union(rte.begin(), rte.end(),
                                           destRte.begin(), destRte.end(),
                                           std::inserter(tmpRte, tmpRte.end()));
                            rte.swap(tmpRte);
                        }
                    }
                }
                nextUnconnectedNonReacher.insert(i);
            }
            hasConnection |= reacher;
        }

        unconnectedNonReacher = nextUnconnectedNonReacher;

        if (!hasConnection) {
            break;
        }
    } while (true);

#if DEBUG
    // Prints the payments span.
    for (tie(eit, eend) = edges(*paymentsSpan); eit != eend; ++eit) {
        printf("%lu -> %lu (rte:",
               source(*eit, *paymentsSpan), target(*eit, *paymentsSpan));
        for (set<Vertex>::const_iterator it = (*routingTable)[*eit].begin();
             it != (*routingTable)[*eit].end(); ++it) {
            std::cout << ' ' << *it;
        }
        std::cout << ")\n";
    }
#endif
}

void computePaymentAmounts(const Graph &payments, Graph *paymentsSpan,
                           const map<Edge, set<Vertex> > &routingTable) {
    // A per-vertex list that stores the amount that needs to be forwarded to
    // which destination node.
    map<Vertex, map<Vertex, double> > forwardPayments;

    vector<Vertex> sortedVertices;
    topological_sort(*paymentsSpan, std::back_inserter(sortedVertices));

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

        // Forwards the payments one hop by running over all outgoing edges and
        // adjusting the forwardPayments entries for these nodes.
        // It also computes and sets the total payment amount for these edges.
        for (tie(oeit, oeend) = out_edges(from, *paymentsSpan); oeit != oeend;
             ++oeit) {
            const Vertex to = target(*oeit, *paymentsSpan);

            double &edgeAmount = get(edge_capacity, *paymentsSpan, *oeit);
            edgeAmount = 0.0;

            // Runs over the nodes that are after the current edge and checks
            // if there are forward-payments that have to be written to the to's
            // forward payments list.
            for (set<Vertex>::const_iterator it =
                     routingTable.at(*oeit).begin();
                 it != routingTable.at(*oeit).end(); ++it) {
                const Vertex paymentTargetNode = *it;

                map<Vertex, double>::const_iterator paymentAmountIt =
                    paymentList.find(paymentTargetNode);
                if (paymentAmountIt != paymentList.end()) {
                    // Payments to paymentTargetNode are forwarded through this
                    // edge.
                    double forwardAmount = paymentAmountIt->second;

                    edgeAmount += forwardAmount;
                    forwardPayments[to][paymentTargetNode] += forwardAmount;
                }
            }
        }
    }

}

} // namespace

#include <cstdio>
#include <iostream>
#include <map>

#include "transaction_minimizer.h"

using std::map;

namespace {

Vertex mapNameToVertex(const std::string &name, Graph *graph,
                       map<std::string, Vertex> *nameToVertex) {
    Vertex ret;
    map<std::string, Vertex>::const_iterator it = nameToVertex->find(name);
    if (it == nameToVertex->end()) {
        ret = add_vertex(*graph);
        nameToVertex->insert(make_pair(name, ret));
        put(boost::vertex_name, *graph, ret, name);
    } else {
        ret = it->second;
    }

    return ret;
}

Graph readGraphFromStdIO() {
    Graph payments;

    map<std::string, Vertex> nameToVertex;
    while (true) {
        // Reads in the next edge.
        std::string fromName, toName;
        double amount;
        std::cin >> fromName >> toName >> amount;

        if (std::cin.eof()) {
            break;
        }

        // Maps the names to vertices.
        Vertex from = mapNameToVertex(fromName, &payments, &nameToVertex);
        Vertex to = mapNameToVertex(toName, &payments, &nameToVertex);

        // Adds the edge to the graph.
        Edge e;
        tie(e, boost::tuples::ignore) = add_edge(from, to, payments);
        put(boost::edge_capacity, payments, e, amount);
    }

    return payments;
}

void printPaymentsToStdIO(const Graph &payments) {
    // Searches for the longest name.
    size_t longestName = 0;
    int longestAmount = 0;
    Graph::edge_iterator eit, eend;
    for (tie(eit, eend) = edges(payments); eit != eend; ++eit) {
        longestName = std::max(longestName,
            get(boost::vertex_name, payments, source(*eit, payments)).length());
        longestName = std::max(longestName,
            get(boost::vertex_name, payments, target(*eit, payments)).length());
        longestAmount = std::max(longestAmount,
                                 (int) std::ceil(std::log10(
                                         get(boost::edge_capacity, payments,
                                             *eit))));
    }

    // Prints out the final payments graph.
    for (tie(eit, eend) = edges(payments); eit != eend; ++eit) {
        const std::string fromName =
            get(boost::vertex_name, payments, source(*eit, payments));
        const std::string toName =
            get(boost::vertex_name, payments, target(*eit, payments));
        double amount = get(boost::edge_capacity, payments, *eit);

        printf("%s%*s-> %s%*s: %*.2f\n",
               fromName.c_str(),
               (int) (fromName.length() - longestName - 1), " ",
               toName.c_str(),
               (int) (toName.length() - longestName - 1), " ",
               longestAmount+3, // +3 because of the point and precision
               amount);
    }
}

} // namespace

int main(int, char **) {
    // Reads in the payments graph.
    // Expects a format like the following:
    // A  B   100.00
    // B  A   20
    // C  A   55.20
    Graph payments = readGraphFromStdIO();

    // Computes the minimum transactions.
    reduceToMinimumTransactions(&payments);

    // Prints out the final payments.
    printPaymentsToStdIO(payments);
}


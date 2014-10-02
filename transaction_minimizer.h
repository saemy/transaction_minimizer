#ifndef TRANSACTION_MINIMIZER_H
#define TRANSACTION_MINIMIZER_H

#include <string>

#include <boost/graph/adjacency_list.hpp>

typedef boost::adjacency_list<boost::setS, boost::vecS, boost::bidirectionalS,
        boost::property<boost::vertex_name_t, std::string>,
        boost::property<boost::edge_capacity_t, double> > Graph;

typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
typedef boost::graph_traits<Graph>::edge_descriptor Edge;


void createMinimumTransactionsGraph(const Graph &payments, Graph *paymentsSpan);

#endif // TRANSACTION_MINIMIZER_H

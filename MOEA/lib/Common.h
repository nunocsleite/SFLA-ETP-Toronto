#ifndef COMMON_H
#define COMMON_H


#include <boost/graph/adjacency_list.hpp>


// Use a vector for vertices and a linked list for edges
typedef boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS> AdjacencyList;



#endif // COMMON_H

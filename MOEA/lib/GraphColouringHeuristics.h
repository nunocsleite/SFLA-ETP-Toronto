#ifndef GRAPH_COLOURING_HEURISTIC
#define GRAPH_COLOURING_HEURISTIC


#include "Common.h"
#include "Chromosome.h"

#include <boost/random.hpp>
#include <boost/graph/adjacency_list.hpp>
#include "Matrix.h"
#include "VertexPriorityQueue.h"


using namespace std;
using namespace boost;


class Chromosome;



class GCHeuristics {

public:
    static void SD(Matrix const& conflictMatrix, AdjacencyList const& graph, Chromosome& timetable);

    static void SD(Matrix const& conflictMatrix, AdjacencyList const& graph, Chromosome& timetable,
                   vector<int> const& fixedExams, unordered_map<int, int> const& examsPeriods, vector<int> & remainderExams);

    static void insertFixedExams(vector<unordered_set<int> >& periods, vector<bool>& removed_exams, int numPeriods,
                                 VertexPriorityQueue& pq, vector<unordered_set<int> >& examsAvailablePeriodsList,
                                 vector<int> const& fixedExams, unordered_map<int, int> const& examsPeriods,
                                 vector<int> &remainderExams, int N, AdjacencyList const& graph);
};

#endif

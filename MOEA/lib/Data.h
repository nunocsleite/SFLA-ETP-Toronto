#ifndef DATA
#define DATA

#include <vector>
#include <iostream>
#include "Matrix.h"
#include <boost/unordered_map.hpp>
#include "Common.h"

using namespace std;
using namespace boost;


class ProblemData {
    // Conflict matrix
    Matrix const* conflictMatrix;
    // Graph
    AdjacencyList const* graph;
    // Number of students
    int numStudents;
    // Number of exams
    int numExams;
    // Number of enrolments
    int numEnrolments;
    // Course student counts
    unordered_map<int, int> const* courseStudentCounts;

public:
	ProblemData() { }
    ProblemData(Matrix const* conflictMatrix, AdjacencyList const* graph,
                int numStudents, int numExams,
                int numEnrolments, unordered_map<int, int> const* courseStudentCounts)
        : conflictMatrix(conflictMatrix), graph(graph), numStudents(numStudents),
          numExams(numExams), numEnrolments(numEnrolments), courseStudentCounts(courseStudentCounts) { }
    Matrix const* getConflictMatrix() const { return conflictMatrix; }
    AdjacencyList const* getGraph() const { return graph; }

    int getNumStudents() const { return numStudents; }
    int getNumExams() const { return numExams; }
    int getNumEnrolments() const { return numEnrolments; }
    unordered_map<int, int> const* getCourseStudentCounts() const { return courseStudentCounts; }
};


class Data {
	//===
	// Algorithm data
	//===
	// Min number of periods.
    int minPeriods;
    // Max number of periods.
    int maxPeriods;
    // Period range
	vector<int> range;
    // Original number of Periods
    int originalNumPeriods;

	//===
	// Problem data
	//===
	ProblemData problData;
public:
    Data(int minPeriods, int maxPeriods, int originalNumPeriods)
        : minPeriods(minPeriods), maxPeriods(maxPeriods), originalNumPeriods(originalNumPeriods) {
		range.push_back(minPeriods);
		range.push_back(maxPeriods);
	}
	int getMinPeriods() const { return minPeriods; }
	int getMaxPeriods() const { return maxPeriods; }
	const vector<int>& getRange() const { return range; }
    int getOriginalNumPeriods() const { return originalNumPeriods; }
	void setProblemData(const ProblemData& problData) { this->problData = problData; }
	ProblemData const& getProblemData() const { return problData; }
};


#endif

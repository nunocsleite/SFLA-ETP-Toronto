#ifndef TESTSET_H
#define TESTSET_H

#include <iostream>
#include <fstream>
#include <vector>
#include <iterator>
#include <string>
#include <sstream>
#include <boost/tokenizer.hpp>
#include "GraphColouringHeuristics.h"
#include "Matrix.h"
#include <boost/unordered_map.hpp>


using namespace std;
using namespace boost;


std::string read_from_file(char const* infile);


///////////////////////////////////////////////////////////////////////////////
//  Class representing a test set
///////////////////////////////////////////////////////////////////////////////
class TestSet {

public:
    string getName() const { return name; }
    string getDesc() const { return desc; }
    int getNumStudents() const { return numStudents; }
    int getNumExams() const { return numExams; }
    int getNumEnrolments() const { return numEnrolments; }
    double getConflictDensity() const { return conflictDensity; }
    int getNumPeriods() const { return periods; }
    Data const& getData() const { return data; }
    // Get student counts
    unordered_map<int, int> const* getCourseStudentCounts() const { return &courseStudentCounts; }
    Matrix const* getConflictMatrix() const { return ptrConflictMatrix.get(); }
    AdjacencyList const* getGraph() const { return ptrGraph.get(); }
    friend ostream& operator<<(ostream& os, const TestSet& t);

protected:
    // Constructor
    TestSet(string _name, string _desc, Data& _data, string _rootDir)
        : name(_name), desc(_desc), periods(data.getOriginalNumPeriods()), rootDir(_rootDir), numStudents(0),
          numExams(0), numEnrolments(0), conflictDensity(0.0), data(_data)
    { }

    virtual void load(string _testSetName, string _rootDir) = 0;

//private:
protected:

    ////////////////////////
    // Fields
    ////////////////////////
    string name;
    string desc;
    int periods; // The minimum number of time slots of a feasible solution
    string rootDir;
    int numStudents;
    int numExams;
    int numEnrolments;
    double conflictDensity;
    Data& data;
    // Map to keep course student counts. Exams indexed from [1..numExams].
    unordered_map<int, int> courseStudentCounts;
    // Conflict matrix
    boost::shared_ptr<Matrix>  ptrConflictMatrix;
    // Graph representing exam relations
    boost::shared_ptr<AdjacencyList> ptrGraph;
};





#endif // TESTSET_H

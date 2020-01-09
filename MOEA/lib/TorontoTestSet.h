#ifndef TORONTOTESTSET_H
#define TORONTOTESTSET_H

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
#include "TestSet.h"

using namespace std;
using namespace boost;


//std::string read_from_file(char const* infile);



/////////////////////////////
// Toronto Dataset
//
/////////////////////////////
class TorontoTestSet : public TestSet {

public:
    // Constructor
    TorontoTestSet(string _name, string _desc, Data& _data, string _rootDir)
        : TestSet(_name, _desc, _data, _rootDir)  {

        cout << "TorontoTestSet ctor" << endl;

        // Load data
        load(_name, _rootDir);

    }


protected:
    void load(string _testSetName, string _rootDir);

private:
    void loadConflictMatrix();
    void loadCourseCounts();
    void computeStudentCount();
};



#endif // TORONTOTESTSET_H

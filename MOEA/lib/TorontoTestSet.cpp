
#include "TorontoTestSet.h"

///////////////////////////////////////////////////////////////////////////////
//  Helper function reading a file into a string
///////////////////////////////////////////////////////////////////////////////
std::string read_from_file(char const* infile)
{
    std::ifstream instream(infile);
    if (!instream.is_open()) {
        std::cerr << "Couldn't open file: " << infile << std::endl;
        exit(-1);
    }
    instream.unsetf(std::ios::skipws);      // No white space skipping!
    return std::string(std::istreambuf_iterator<char>(instream.rdbuf()),
                       std::istreambuf_iterator<char>());
}




////////////////////////////////////////////////////////////
// Toronto Dataset Methods
//
////////////////////////////////////////////////////////////

void TorontoTestSet::load(string _testSetName, string _rootDir)
{
    // Load course student counts from '.crs' file
    loadCourseCounts();
    // Compute student count from '.stu' file
    computeStudentCount();
    // Load conflict matrix and build adjacency list
    loadConflictMatrix();
    // Set Problem data
    boost::shared_ptr<ProblemData> ptrProblData( new ProblemData(this->getConflictMatrix(),
                                      this->getGraph(), this->getNumStudents(),
                                      this->getNumExams(),this->getNumEnrolments(), this->getCourseStudentCounts()) );
    data.setProblemData(*ptrProblData.get());
}


void TorontoTestSet::loadCourseCounts() {
    string filename = rootDir + "/" + name + ".crs";

//    cout << filename << endl;

    string wholeFile = read_from_file(filename.c_str());
    // Process whole file
    tokenizer<> tok(wholeFile);
    int idx, val;
    // Fill map
    for (tokenizer<>::iterator beg = tok.begin(); beg != tok.end();) {
        idx = atoi((*beg++).c_str());
        val = atoi((*beg++).c_str());
//        cout << idx << " " << val << endl;

        // Insert entry into map
        pair<unordered_map<int, int>::iterator, bool> p = courseStudentCounts.insert(pair<int, int>(idx, val));

        // Update number of enrolments
        numEnrolments += val;
        // Increment # exams
        numExams++;
    }
//    copy(courseStudentCounts.begin(), courseStudentCounts.end(), ostream_iterator<int>(cout, "\n"));
//    cin.get();
}


void TorontoTestSet::computeStudentCount() {
    string filename = rootDir + "/" + name + ".stu";

//    cout << filename << endl;

    string wholeFile = read_from_file(filename.c_str());
    // Process whole file escaping '\n'
    tokenizer<escaped_list_separator<char> > tok(wholeFile, escaped_list_separator<char>('\\', '\n'));
    int numLines = 0;
    for (tokenizer<escaped_list_separator<char> >::iterator beg = tok.begin(); beg != tok.end(); ++beg, ++numLines) { }
    // The number of students is equal to the number of records in '.stu' file
    numStudents = numLines;
//    cout << "numStudents = " << numStudents << endl;
}


// Load the conflict matrix for this test set
void TorontoTestSet::loadConflictMatrix() {
    string filename = rootDir + "/" + name + ".cm";
//    cout << filename << endl;
    string wholeFile = read_from_file(filename.c_str());
    // Get number of columns
    int ncols = 0;
    string firstLine;
    ifstream file(filename.c_str());
    getline(file, firstLine);
    file.close();
    tokenizer<> firstLineTok(firstLine);
    for (tokenizer<>::iterator beg = firstLineTok.begin(); beg != firstLineTok.end(); ++beg){
        ++ncols;
    }
    // Create shared ptr to manage Conflict Matrix
    boost::shared_ptr<Matrix> ptrMatrix(new Matrix(ncols, ncols));
    ptrConflictMatrix = ptrMatrix;
    // Obtain internal vector
    vector<int>& v = ptrConflictMatrix->getVec();
    // Instantiate graph with ncols vertices
    boost::shared_ptr<AdjacencyList> ptrGraphAux(new AdjacencyList(ncols+1));
    ptrGraph = ptrGraphAux;
    // Process whole file
    tokenizer<> tok(wholeFile);
    int v1 = 1, v2 = 1;
    int cost, nonZeroElements = 0;
    for (tokenizer<>::iterator beg = tok.begin(); beg != tok.end(); ++beg) {
        cost = atoi((*beg).c_str());
        // Also save arc's cost
        v.push_back(cost);
        if (cost != 0) {
            // Increment the number of non-zero elements
            ++nonZeroElements;
            if (v2 > v1) {
                add_edge(v1, v2, *ptrGraph.get());
            }
        }
        if (v2 == ncols) {
            ++v1;
            v2 = 1;
        }
        else
            ++v2;
    }
    // The ‘conflict’ density is the ratio of the number of non-zero elements
    // in the conflict matrix to the total number of conflict matrix elements.
    double matrixElements = ncols*ncols;
    conflictDensity = nonZeroElements / matrixElements;
}



ostream& operator<<(ostream& os, const TestSet& t) {
//    os << "Test set\t # students\t # exams\t # enrolments\t conflict density\t # time slots" << endl;

//    os << t.getName() << "\t " << t.getNumStudents() << "\t " << t.getNumExams() << "\t "
//       << t.getNumEnrolments() << "\t " << t.getConflictDensity() << "\t " << t.getNumPeriods() << endl;

    os << t.getName();

    return os;
}










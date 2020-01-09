#ifndef CHROMOSOME
#define CHROMOSOME

#include <boost/random.hpp>
#include <vector>
#include <boost/unordered_set.hpp>
#include <set>
#include <cmath>
#include "Data.h"
#include "Matrix.h"
#include <moeo>
#include <eo>
#include <es/eoRealInitBounded.h>
#include <es/eoRealOp.h>
#include <boost/unordered_map.hpp>
#include "Exam.h"
#include "Period.h"
#include "Common.h"


using namespace std;
using namespace boost;



// Chromosome implementation.
// A chromosome encodes a complete and feasible timetable.
// A period contains a list of exams to be held in the period.
// For each chromosome, a timetable with a random number of empty
// periods within the desired range is created.
////
class Chromosome {

public:
    Chromosome() {  }
    Chromosome(const Data& data) { init(data); }
    // Copy constructor
    Chromosome (Chromosome const& _chrom)
        : originalNumPeriods(_chrom.originalNumPeriods),
          range(_chrom.range),
          periods(_chrom.periods),
          proximityCost(_chrom.proximityCost),
          numClashes(_chrom.numClashes),
          conflictMatrix(_chrom.conflictMatrix),

          graph(_chrom.graph),

          numStudents(_chrom.numStudents),
          numExams(_chrom.numExams),
          numEnrolments(_chrom.numEnrolments),
          courseStudentCounts(_chrom.courseStudentCounts),
          courseNames(_chrom.courseNames)
    { }


    void init(Data const& data);

    vector<int> getFeasiblePeriods(int exam, int period) const;
    inline int getNumPeriods() const { return periods.size(); }
    inline int getOriginalNumPeriods() const { return originalNumPeriods; }
    inline vector<int> const& getRange() const { return range; }
    inline vector<unordered_set<int> >& getPeriods() { return periods; }
    inline vector<unordered_set<int> > const& getConstPeriods() const { return periods; }
    inline void setPeriods(vector<unordered_set<int> > const& _periods) {
        for (int i = 0; i < _periods.size(); ++i)
            periods[i] = _periods[i];
    }
    // Get Conflict matrix
    Matrix const* getConflictMatrix() { return conflictMatrix; }

    // Get graph
    AdjacencyList const* getGraph() { return graph; }

    inline double getProximityCost() const { return proximityCost; }


//    virtual void computeProximityCosts() = 0;


    int getNumClashesPeriod(int pi) const;
    int getNumStudentsPeriod(int pi) const;
    void pack();
    void expand();
    void validate() const;
    bool isFeasiblePeriodExam(int period, int exam) const;

    Chromosome const& getChromosome() const { return *this; }


    //	void setUCNames(char** names) { deetc_ucs = names; }
    //	char** getUCNames() { return deetc_ucs; }




    friend ostream& operator<<(ostream& os, const Chromosome& timetable);

    double getProximityCostExamPeriod(int exami, int pi) const;

//    void getInfoExamHighestProximityCost(int& examHighCost, int& period, double& highProxCost);
//    void getExamHighestProximityCost(int pi, int& exam, double& highProxCost);

//    void updateProximityCost(double proximityCost) {

////        this->proximityCost = proximityCost;

//        computeProximityCosts();

//        fitness(getProximityCost());

//    }

    typedef vector<unordered_set<int> >::iterator PeriodIterator;
    typedef vector<unordered_set<int> >::const_iterator ConstPeriodIterator;


//private:
protected:
    // Original number of Periods
    int originalNumPeriods;
    // Range
    vector<int> range;
    // Periods containing exams identifiers
    vector<unordered_set<int> > periods;
    // Number of conflicts
    double proximityCost;
    // Number of clashes
    int numClashes;
//    bool validNumClashes;
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
    boost::unordered_map<int, int> const* courseStudentCounts;
    // Course names
    vector<string> courseNames;
//	char** deetc_ucs;

    void baseComputeProximityCosts();

    int getMinPeriod();
    multiset<Period> getExamsAvailablePeriods(const unordered_set<int> &sourceExams, int sourcePeriod);
    void computeAvailablePeriods(unordered_set<int> const& sourceExams, int sourcePeriod, vector<Period>& availablePeriodsVec);
    int computeNumClashesPeriod(int exam, int p) const;
    int computeNumClashesExamPeriod(int exam, int p, int minPeriod) const;
    void removeEmptyPeriods();

};
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



/////////// EO CHROMOSOME /////////////////////////////////////////////////////////////////////////

// EO Chromosome implementation.
class eoChromosome : public Chromosome, public EO<double> {

public:


    /// TODO: OPTIMIZAR - USAR PONTEIRO PARA CROMOSOMA

    // Copy constructor
    eoChromosome (Chromosome const& _chrom)
        : Chromosome(_chrom)
    {
        computeProximityCosts();
    }

    eoChromosome () {  }
    eoChromosome (const Data& data) { init(data); }
    string className() const { return "EO Chromosome"; };
    void updateProximityCost(double proximityCost) {

//        this->proximityCost = proximityCost;

        computeProximityCosts();

//        fitness(getProximityCost());

    }
    /*virtual*/ void computeProximityCosts() {
        proximityCost = 0;

        baseComputeProximityCosts();
        // Update fitness
        fitness(getProximityCost());

//        cout << "fitness() = " << fitness() << endl;
    };

    friend ostream& operator<<(ostream& os, const eoChromosome& timetable);

};

///////////////////////////////////////////////////////////////////////////////////////////////////



///////// MOEO CHROMOSOME /////////////////////////////////////////////////////////////////////////

// The moeoObjectiveVectorTraits : minimizing 2 objectives
class ETTPObjectiveVectorTraits : public moeoObjectiveVectorTraits {
public:
    static bool minimizing (int i) { return true;  }
    static bool maximizing (int i) { return false; }
    static unsigned int nObjectives () { return 2; }
};

// Objective vector of real values
typedef moeoRealObjectiveVector<ETTPObjectiveVectorTraits> ETTPObjectiveVector;


// MOEO Chromosome implementation.
class moeoChromosome : public Chromosome, public MOEO<ETTPObjectiveVector> {

public:
    // Copy constructor
    moeoChromosome (Chromosome const& _chrom)
        : Chromosome(_chrom)
    {
//        computeProximityCosts(); // Added
    }

    moeoChromosome () {  }
    moeoChromosome (const Data& data) { init(data); }
    string className() const { return "MOEO Chromosome"; };

    void updateProximityCost(double proximityCost) {

//        this->proximityCost = proximityCost;

        computeProximityCosts();

//        fitness(getProximityCost());

    }
//    virtual void computeProximityCosts() {
    void computeProximityCosts() {
        proximityCost = 0;
        baseComputeProximityCosts();
        // Update fitness
        fitness(getProximityCost());
    };

    friend ostream& operator<<(ostream& os, const moeoChromosome& timetable);

};

///////////////////////////////////////////////////////////////////////////////////////////////////




//// Chromosome implementation.
//// A chromosome encodes a complete and feasible timetable.
//// A period contains a list of exams to be held in the period.
//// For each chromosome, a timetable with a random number of empty
//// periods within the desired range is created.
//////
//class Chromosome : public MOEO<ETTPObjectiveVector> {

//public:
//    Chromosome() {  }
//    Chromosome(const Data& data) { init(data); }
//    void init(Data const& data);
//    string className() const { return "ETTP"; };

//    vector<int> getFeasiblePeriods(int exam, int period) const;
//    inline int getNumPeriods() const { return periods.size(); }
//    inline int getOriginalNumPeriods() const { return originalNumPeriods; }
//    inline vector<int> const& getRange() const { return range; }
////	inline vector<vector<int> >& getPeriods() { return periods; }
//    inline vector<unordered_set<int> >& getPeriods() { return periods; }
////	inline vector<vector<int> > const& getConstPeriods() const { return periods; }
//    inline vector<unordered_set<int> > const& getConstPeriods() const { return periods; }

//    inline double getProximityCost() const { return proximityCost; }
//    // Get Conflict matrix
//    Matrix const* getConflictMatrix() { return conflictMatrix; }
//    void computeProximityCosts();
//    int getNumClashesPeriod(int pi) const;
//    int getNumStudentsPeriod(int pi) const;
//    void pack();
//    void expand();
//    void validate() const;
//    //	void setUCNames(char** names) { deetc_ucs = names; }
//    //	char** getUCNames() { return deetc_ucs; }

//    friend ostream& operator<<(ostream& os, const Chromosome& timetable);

//    bool invalidNumClashes() const { return !validNumClashes; }
////    int getNumClashes() const {

//////        if (invalidNumClashes())
//////        {
//////            throw std::runtime_error("invalid number of clashes in chromosome");
//////        }
////        return numClashes;
////    }
//    void invalidateNumClashes() { validNumClashes = false; }
////    void computeNumClashes();

////    int getNumClashesExamPeriod(int exami, int pi) const;
//    double getProximityCostExamPeriod(int exami, int pi) const;

////    void updateNumClashes(int numClashes) {
////        this->numClashes = numClashes;
////        validNumClashes = true;
////        // Set fitness
////        fitness(numClashes);

//////        cout << numClashes << endl;
//////        cin.get();
////    }

//    void updateProximityCost(double proximityCost) {
//        this->proximityCost = proximityCost;
////        validNumClashes = true;
//        // Set fitness
//        fitness(proximityCost);

////        cout << proximityCost << endl;
////        cin.get();
//    }


////    typedef vector<vector<int> >::iterator PeriodIterator;
////    typedef vector<vector<int> >::const_iterator ConstPeriodIterator;
//    typedef vector<unordered_set<int> >::iterator PeriodIterator;
//    typedef vector<unordered_set<int> >::const_iterator ConstPeriodIterator;

//private:
//    // Original number of Periods
//    int originalNumPeriods;
//    // Range
//    vector<int> range;
//    // Periods containing exams identifiers
////    vector<vector<int> > periods;

//    vector<unordered_set<int> > periods;
//    // Number of conflicts
//    double proximityCost;
//    // Number of clashes
//    int numClashes;
//    bool validNumClashes;
//    // Conflict matrix
//    Matrix const* conflictMatrix;
//    // Number of students
//    int numStudents;
//    // Number of exams
//    int numExams;
//    // Number of enrolments
//    int numEnrolments;
//    // Course student counts
//    unordered_map<int, int> const* courseStudentCounts;
//    // Course names
//    vector<string> courseNames;
////	char** deetc_ucs;

//    bool isFeasiblePeriodExam(int period, int exam) const;
//    int getMinPeriod();
//    multiset<Period> getExamsAvailablePeriods(const unordered_set<int> &sourceExams, int sourcePeriod);
//    void computeAvailablePeriods(unordered_set<int> const& sourceExams, int sourcePeriod, vector<Period>& availablePeriodsVec);
//    int computeNumClashesPeriod(int exam, int p) const;
//    int computeNumClashesExamPeriod(int exam, int p, int minPeriod) const;
//    void removeEmptyPeriods();
//};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


#endif

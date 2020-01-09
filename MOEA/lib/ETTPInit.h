#ifndef ETTPINIT_H
#define ETTPINIT_H

#include "GraphColouringHeuristics.h"
#include "Matrix.h"
#include "Data.h"
#include "TestSet.h"
#include "eoETTPEval.h"

using namespace std;
using namespace boost;


template <typename EOT>
class ETTPInit : public eoInit<EOT> {

//class ETTPInit : public eoInit<Chromosome> {

public:
    const Data& data;
    const Matrix& conflictMatrix;
    const AdjacencyList& graph;

    ETTPInit(Data const& data, Matrix const& conflictMatrix, AdjacencyList const& graph)
        : fixedRange(false), numPeriods(0), data(data), conflictMatrix(conflictMatrix), graph(graph) {}

//    virtual void operator()(Chromosome& _chrom) {
        virtual void operator()(EOT& _chrom) {

    // For each chromosome, a timetable with a random number of empty
        // periods within the desired range is created.
//        _chrom.init(data);

//        cout << "After chrom.init" << endl;

        do {
            // For each chromosome, a timetable with a random number of empty
            // periods within the desired range is created.
            _chrom.init(data);

            // Exams are then inserted into randomly selected periods in the order
            // determined by the graph colouring heuristic, depending on the version
            // of the MOEA. Like the mutation operator, when it is not possible to
            // schedule an exam without violating any of the hard constraints,
            // a new period will be created at the end of the timetable to
            // accommodate the exam.
            GCHeuristics::SD(conflictMatrix, graph, _chrom);

        } while (!isFeasible(_chrom));


        /// ???

        _chrom.computeProximityCosts();
        // Set fitness
//        _chrom.fitness(_chrom.getProximityCost());

//        cout << "Got feasible solution -> fitness = " << _chrom.fitness() << endl << endl;

//        cout << "Got feasible solution -> fitness = " << _chrom.getProximityCost() << endl << endl;


//        cout << "Fitness = " << _chrom.fitness() << endl;
//        chrom.computeNumClashes();

//        moeoETTPEval eval;
//        eval(chrom);
    }

    void setFixedRange(int _numPeriods) {
        fixedRange = true;
        numPeriods = _numPeriods;
    }

private:

    bool fixedRange;
    int numPeriods;

//    bool isFeasible(Chromosome const& _chrom) const {
    bool isFeasible(EOT const& _chrom) const {
        if (!fixedRange) {
            if (_chrom.getNumPeriods() > _chrom.getRange()[1]) {
//                cout << "Infeasible chromosome: Number of periods = " << _chrom.getNumPeriods()
//                     << " above upper limit (" << _chrom.getRange()[1] << ")" << endl;
                return false;
            }

//            if (_chrom.getNumPeriods() < _chrom.getRange()[0]) {
//                cout << "Infeasible chromosome: Number of periods = " << _chrom.getNumPeriods()
//                     << " below lower limit (" << _chrom.getRange()[0] << ")" << endl;
//                return false;
//            }
            return true;
        }
        else { // Fixed unitary range
            if (_chrom.getNumPeriods() > numPeriods) {
                cout << "Infeasible chromosome: Number of periods = " << _chrom.getNumPeriods()
                     << " above upper limit (" << numPeriods << ")" << endl;
                return false;
            }

            if (_chrom.getNumPeriods() < numPeriods) {
                cout << "Infeasible chromosome: Number of periods = " << _chrom.getNumPeriods()
                     << " below lower limit (" << numPeriods << ")" << endl;
                return false;
            }
            return true;
        }
    }
};

#endif // ETTPINIT_H















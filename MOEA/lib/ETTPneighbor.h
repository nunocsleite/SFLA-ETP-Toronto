#ifndef ETTPNEIGHBOR_H
#define ETTPNEIGHBOR_H

#include <neighborhood/moNeighbor.h>
#include "Chromosome.h"
#include "ETTPKempeChainHeuristic.h"


/// TODO: ??
/// use different EO instead of Chromosome because of different fitness function used in SHC
/// (only considers min # clashes)


/**
 * Neighbor computation for ETTP problem (used in local search)
 */
//template<class EOT, class Fitness = typename EOT::Fitness>
//class ETTPneighbor: virtual public moNeighbor<EOT, Fitness> {
class ETTPneighbor: virtual public moNeighbor<eoChromosome, eoChromosome::Fitness> {

public:

//    using moNeighbor<Chromosome, Fitness>::fitness;
    using moNeighbor<eoChromosome, eoChromosome::Fitness>::fitness;

    /**
     * Return the class Name
     * @return the class name as a std::string
     */
    virtual std::string className() const {
        return "ETTPneighbor";
    }


    void setSolution(eoChromosome& _timetable) {
        timetable = _timetable; /// Operator= hum

        // Kempe Chain neighbourhood
        ETTPKempeChainHeuristic kempeChain;
        kempeChain(timetable);
        // Set fitness
        fitness(timetable.fitness());
    }


    /**
     * Apply move on a solution
     * @param _sol the solution to move
     */
    virtual void move(eoChromosome& _timetable) {

        ///
        /// ALTERAR
        ///

        _timetable = timetable;


        /*

//        cout << "ETTPneighbor::move" << endl;
//        cout << "_timetable fitness = " << _timetable.fitness() << endl;
//        cout << "_timetable.numClashes = " << _timetable.getNumClashes() << endl;

        // Move selected neighbor this to _timetable
//        vector<vector<int> >& periods = _timetable.getPeriods();
//        vector<int>& sourcePeriodExams = periods[clashPeriod];
//        vector<int>& targetPeriodExams = periods[periodToMove];
        vector<unordered_set<int> >& periods = _timetable.getPeriods();
        unordered_set<int>& sourcePeriodExams = periods[clashPeriod];
        unordered_set<int>& targetPeriodExams = periods[periodToMove];
        // Remove clashed exam from clashed period
//        vector<int>::iterator removeItr = find_if(sourcePeriodExams.begin(), sourcePeriodExams.end(),
//                                                                bind2nd(equal_to<int>(), clashedExam));

        /// OPTIMIZAR - Done

//        unordered_set<int>::iterator removeItr = find_if(sourcePeriodExams.begin(), sourcePeriodExams.end(),
//                                                                bind2nd(equal_to<int>(), clashedExam));
//        if (removeItr == sourcePeriodExams.end())
//            throw std::runtime_error("Could not remove clashed exam");

//        sourcePeriodExams.erase(removeItr);

        if (sourcePeriodExams.erase(clashedExam) == 0)
            throw std::runtime_error("Could not remove clashed exam");


        // Insert clashed exam into new period
//        targetPeriodExams.push_back(clashedExam);
        targetPeriodExams.insert(clashedExam);
        // Update number of clashes
//        _timetable.updateNumClashes(numClashes);
        _timetable.updateProximityCost(proximityCost);


//        cout << "Proximity Cost = " << proximityCost << endl;
//        cout << "Fitness = " << _timetable.fitness() << endl;


        _timetable.computeProximityCosts();

//        cout << "After computing: Proximity Cost = " << _timetable.fitness() << endl;


//        cout << "After update -> _timetable fitness = " << _timetable.fitness() << endl;
//        cout << "_timetable.numClashes = " << _timetable.getNumClashes() << endl;

//        cin.get();

*/




    }

    /**
     * Compute move cost
     * @param _timetable the solution
     * @param _clashPeriod the source period
     * @param _clashedExam the source exam
     * @param _periodToMove the period where the exam is to move
     */
    void computeMoveCost(eoChromosome const& _timetable, int _clashPeriod, int _clashedExam, int _periodToMove) {

        clashPeriod = _clashPeriod;
        clashedExam = _clashedExam;
        periodToMove = _periodToMove;

        eoChromosome aux = _timetable;

        vector<unordered_set<int> >& periods = aux.getPeriods();
        unordered_set<int>& sourcePeriodExams = periods[clashPeriod];
        unordered_set<int>& targetPeriodExams = periods[periodToMove];

        // Remove clashed exam from clashed period
        if (sourcePeriodExams.erase(clashedExam) == 0)
            throw std::runtime_error("[computeMoveCost method] Could not remove clashed exam");


        // Insert clashed exam into new period
        targetPeriodExams.insert(clashedExam);
        // Update number of clashes
        aux.computeProximityCosts();

        // Number of clashes after move
        proximityCost = aux.fitness();
        // Set fitness to number of clashes
        fitness(proximityCost);

//        cout << "[computeMoveCost method] Fitness = " << fitness() << endl;


        /// SEE ERROR WHEN COMPUTING DELTA EVALUATION
        ////////////////////////////////////////////////////


        /*       clashPeriod = _clashPeriod;
        clashedExam = _clashedExam;
        periodToMove = _periodToMove;
        // Compute neighbor number of clashes as a result of
        // moving '_clashedExam' from '_clashPeriod' to '_periodToMove'
        // It's assumed that solution timetable has already the number of clashes computed
//        numClashes = _timetable.getNumClashes();
//        proximityCost = _timetable.getProximityCost();

        proximityCost = _timetable.fitness();

        // Update number of clashes obtained if we move solution to neighbor:
        //   a. Get number of clashes of clashed exam at the clash period
//        int n1 = _timetable.getNumClashesExamPeriod(clashedExam, clashPeriod);
        double n1 = _timetable.getProximityCostExamPeriod(clashedExam, clashPeriod);
        //   b. Get number of clashes of clashed exam at the new period
//        int n2 = _timetable.getNumClashesExamPeriod(clashedExam, periodToMove);
        double n2 = _timetable.getProximityCostExamPeriod(clashedExam, periodToMove);

        //   Number of clashes after move
//        numClashes = numClashes - n1 + n2;
        proximityCost = proximityCost - n1 + n2;
        // Set fitness to number of clashes
//        fitness(numClashes);
        fitness(proximityCost);

//        cout << "[computeMoveCost]" << endl;
//        cout << "clashPeriod = " << clashPeriod << ", clashedExam = " << clashedExam << " periodToMove = " << periodToMove << endl;
//        cout << "n1 = " << n1 << ", n2 = " << n2 << ", proximityCost = " << proximityCost << endl;
*/






    }

//    int getNumClashes() const { return numClashes; }
    double getProximityCost() const { return proximityCost; }

//    function moveInfo = computeMoveClashes(Data, timetable, clashPeriod, ...
            //                                        clashedExam, period)
            //        moveInfo.timetable = [];
            //        moveInfo.numClashes = -1;
            //        % Verify if clashedExam can move to period (if it's feasible)
            //        % If it can move, then compute the resulting timetable and number of
            //        % clashes.
            //        if (isFeasibleExamPeriod(Data, timetable, period, clashedExam))
            //            % Add clashed exam to period
            //            timetable.Periods{period} = [timetable.Periods{period} clashedExam];
            //            % Remove clashed exam from original period
            //            timetable.Periods{clashPeriod} = remove(clashedExam, timetable.Periods{clashPeriod});
            //            % Compute number of clashes of new timetable
            //            moveInfo.timetable = computeNumClashes(Data, timetable);
            //            moveInfo.numClashes = moveInfo.timetable.NumClashes;
            //        end
            //    end

private:
    int clashPeriod;
    int clashedExam;
    int periodToMove;
    // Neighbor number of clashes
//    int numClashes;
    double proximityCost;


    eoChromosome timetable;
};

#endif // ETTPNEIGHBOR_H












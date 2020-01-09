#ifndef EOEVOLUTIONOPERATOR_H
#define EOEVOLUTIONOPERATOR_H

#include <eoEvalFunc.h>
#include "Chromosome.h"
#include "eoSFLA.h"
#include "ETTPKempeChainHeuristic.h"


//////////////////////////////////////////////////////////////
//
// eoEvolutionOperator 1 - With variation operators
//
//////////////////////////////////////////////////////////////
class eoSFLAEvolOperator_1 : public eoBinOp<eoChromosome> {

    typedef eoChromosome POT;

public:
    bool operator()(POT& _Pw, POT const& _Pb) {
        double originalPbFitness = _Pb.fitness();
        double originalPwFitness = _Pw.fitness();

        auto newPw = _Pb;

        if (rng.uniform() < 0.1) {

            POT& frog1 = newPw;
            // Insert exams
            int randNumPeriods = 1;

            for (int i = 1; i <= randNumPeriods; ++i) {
                // Get Frog1 periods
                vector<unordered_set<int> >& frog1Periods = frog1.getPeriods();
                // Generate a random index in frogs
                int randPeriod1 = rng.random(frog1.getNumPeriods());
                int randPeriod2 = rng.random(frog1.getNumPeriods());
                // Swap period exams
                auto aux = frog1Periods[randPeriod1];
                frog1Periods[randPeriod1] = frog1Periods[randPeriod2];
                frog1Periods[randPeriod2] = aux;
            }

    //            double T1 = 10; // 39.18, pr = 0.01, with period swap, 3 periods
    //            // Stochastic Hill Climber algorithm instance
    //            //
    //            // Neighborhood
    //            ETTPneighborhood neighborhood;
    //            // Full evaluation function
    //            ProximityCostEval<POT> fullEval;
    //            // Neighbor evaluation function
    //            ETTPneighborEval neighEval;
    //            moSimpleCoolingSchedule<eoChromosome> cool(T1, 0.1, 5, T1); // SHC
    //            moSA<ETTPneighbor> sa(neighborhood, fullEval, neighEval, cool);
    //            sa(frog1);

        }

        vector<unordered_set<int> >& pwPeriods = _Pw.getPeriods();
        // Get newPw periods
        vector<unordered_set<int> >& newPwPeriods = newPw.getPeriods();

    //        //////////////////
    //        // Get best from frog Pw
    //        int bestDayIdx = bestDay(_Pw);
    //        // Generate a random index in Pb
    //        int randPeriod = rng.random(_Pb.getNumPeriods());
    //        // Pw exams
    //        unordered_set<int>& pwExams = pwPeriods[bestDayIdx];
    //        // NewPw exams
    //        unordered_set<int>& newPwExams = newPwPeriods[randPeriod];
    //        //////////////////



        // Generate random number of periods to insert from Pw
        int randNumPeriods = 3;
//        int randNumPeriods = 1;

        for (int i = 1; i <= randNumPeriods; ++i) {
            // Get a random period from Pw
            int randPwPeriod = rng.random(_Pw.getNumPeriods());
            // Insert Pw random period exams (those that are feasible) into the new frog and remove duplicated exams

            // Generate a random index in the new frog
            int randNewPwPeriod = rng.random(newPw.getNumPeriods());
            // Pw exams from best day
            unordered_set<int>& pwExams = pwPeriods[randPwPeriod];
            // newPw exams from rand period
            unordered_set<int>& newPwExams = newPwPeriods[randNewPwPeriod];

            // Insert Pw period into newPw into that random position
            for (auto it = pwExams.begin(); it != pwExams.end(); ++it) {
                // If the exam does not exist and is feasible then insert it
                if (newPwExams.find(*it) == newPwExams.end() && newPw.isFeasiblePeriodExam(randNewPwPeriod, *it)) {
    //                if (newPwExams.find(*it) == newPwExams.end() && newPw.isFeasiblePeriodExam(bestDayIdx, *it)) {

                newPwExams.insert(*it);
    //                ++numInsertedExams;
    //                cout << "Insert: " << *it << endl;

                    /// TODO / SEE POSSIBLE OPTIMIZATIONS

                    // Remove duplicated exams
                    for (int i = 0; i < newPwPeriods.size(); ++i) {
                        if (i != randNewPwPeriod) {
    //                        if (i != bestDayIdx) {
                            unordered_set<int>& exams = newPwPeriods[i];
                            if (exams.find(*it) != exams.end()) {
                                exams.erase(*it);
    //                            ++numRemovedExams;
    //                            cout << "Remove: " << *it << endl;
                                break;
                            }
                        }
                    }

                }
    //            else {
    //                notInserted.insert(*it);
    //                ++numNotInsertedExams;
    //            }

            }
        }

        if (rng.uniform() < 0.05) {
            POT& frog1 = newPw;
            // Insert exams
            // Get Frog1 periods
            vector<unordered_set<int> >& frog1Periods = frog1.getPeriods();
            // Generate a random index in frogs
            int randPeriod1 = rng.random(frog1.getNumPeriods());
            int randPeriod2 = rng.random(frog1.getNumPeriods());
            // Swap period exams
            auto aux = frog1Periods[randPeriod1];
            frog1Periods[randPeriod1] = frog1Periods[randPeriod2];
            frog1Periods[randPeriod2] = aux;
        }

        // Copy newPw to _Pw
        _Pw = newPw;
        _Pw.computeProximityCosts();

    //        _Pw.validate();

    //        cout << "originalPbFitness = " << originalPbFitness << endl;
//        cout << "originalPwFitness = " << originalPwFitness << endl;
//        cout << "[After] Pw Fitness " << _Pw.fitness() << endl;
    }
};

/*

//////////////////////////////////////////////////////////////
//
// eoEvolutionOperator 3 - With variation operators
//
//////////////////////////////////////////////////////////////
class eoSFLAEvolOperator_3 : public eoBinOp<eoChromosome> {

    typedef eoChromosome POT;

public:
    bool operator()(POT& _Pw, POT const& _Pb) {
        double originalPbFitness = _Pb.fitness();
        double originalPwFitness = _Pw.fitness();

        auto newPw = _Pb;

//        if (rng.uniform() < 0.5) {
        if (rng.uniform() < 0.1) {
//        if (rng.uniform() < 0) {
            POT& frog1 = newPw;
            // Insert exams
            int randNumPeriods = 1;

            for (int i = 1; i <= randNumPeriods; ++i) {
                // Get Frog1 periods
                vector<unordered_set<int> >& frog1Periods = frog1.getPeriods();
                // Generate a random index in frogs
                int randPeriod1 = rng.random(frog1.getNumPeriods());
                int randPeriod2 = rng.random(frog1.getNumPeriods());
                // Swap period exams
                auto aux = frog1Periods[randPeriod1];
                frog1Periods[randPeriod1] = frog1Periods[randPeriod2];
                frog1Periods[randPeriod2] = aux;
            }

    //            double T1 = 10; // 39.18, pr = 0.01, with period swap, 3 periods
    //            // Stochastic Hill Climber algorithm instance
    //            //
    //            // Neighborhood
    //            ETTPneighborhood neighborhood;
    //            // Full evaluation function
    //            ProximityCostEval<POT> fullEval;
    //            // Neighbor evaluation function
    //            ETTPneighborEval neighEval;
    //            moSimpleCoolingSchedule<eoChromosome> cool(T1, 0.1, 5, T1); // SHC
    //            moSA<ETTPneighbor> sa(neighborhood, fullEval, neighEval, cool);
    //            sa(frog1);

        }

        vector<unordered_set<int> >& pwPeriods = _Pw.getPeriods();
        // Get newPw periods
        vector<unordered_set<int> >& newPwPeriods = newPw.getPeriods();

    //        //////////////////
    //        // Get best from frog Pw
    //        int bestDayIdx = bestDay(_Pw);
    //        // Generate a random index in Pb
    //        int randPeriod = rng.random(_Pb.getNumPeriods());
    //        // Pw exams
    //        unordered_set<int>& pwExams = pwPeriods[bestDayIdx];
    //        // NewPw exams
    //        unordered_set<int>& newPwExams = newPwPeriods[randPeriod];
    //        //////////////////



        // Generate random number of periods to insert from Pw
//        int randNumPeriods = 3;
        int randNumPeriods = 3;

        for (int i = 1; i <= randNumPeriods; ++i) {
            // Get a random period from Pw
            int randPwPeriod = rng.random(_Pw.getNumPeriods());
            // Insert Pw random period exams (those that are feasible) into the new frog and remove duplicated exams

            // Generate a random index in the new frog
            int randNewPwPeriod = rng.random(newPw.getNumPeriods());
            // Pw exams from best day
            unordered_set<int>& pwExams = pwPeriods[randPwPeriod];
            // newPw exams from rand period
            unordered_set<int>& newPwExams = newPwPeriods[randNewPwPeriod];

            // Insert Pw period into newPw into that random position
            for (auto it = pwExams.begin(); it != pwExams.end(); ++it) {
                // If the exam does not exist and is feasible then insert it
                if (newPwExams.find(*it) == newPwExams.end() && newPw.isFeasiblePeriodExam(randNewPwPeriod, *it)) {
    //                if (newPwExams.find(*it) == newPwExams.end() && newPw.isFeasiblePeriodExam(bestDayIdx, *it)) {

                newPwExams.insert(*it);
    //                ++numInsertedExams;
    //                cout << "Insert: " << *it << endl;

                    /// TODO / SEE POSSIBLE OPTIMIZATIONS

                    // Remove duplicated exams
                    for (int i = 0; i < newPwPeriods.size(); ++i) {
                        if (i != randNewPwPeriod) {
    //                        if (i != bestDayIdx) {
                            unordered_set<int>& exams = newPwPeriods[i];
                            if (exams.find(*it) != exams.end()) {
                                exams.erase(*it);
    //                            ++numRemovedExams;
    //                            cout << "Remove: " << *it << endl;
                                break;
                            }
                        }
                    }

                }
    //            else {
    //                notInserted.insert(*it);
    //                ++numNotInsertedExams;
    //            }

            }
        }

//        if (rng.uniform() < 0) {
//        if (rng.uniform() < 0.5) {
        if (rng.uniform() < 0.05) {
            POT& frog1 = newPw;
            // Insert exams
            // Get Frog1 periods
            vector<unordered_set<int> >& frog1Periods = frog1.getPeriods();
            // Generate a random index in frogs
            int randPeriod1 = rng.random(frog1.getNumPeriods());
            int randPeriod2 = rng.random(frog1.getNumPeriods());
            // Swap period exams
            auto aux = frog1Periods[randPeriod1];
            frog1Periods[randPeriod1] = frog1Periods[randPeriod2];
            frog1Periods[randPeriod2] = aux;
        }

        // Copy newPw to _Pw
        _Pw = newPw;
        _Pw.computeProximityCosts();

//        _Pw = _Pb;
//        POT& frog1 = _Pw;

//    if (rng.uniform() < 0) {

        if (rng.uniform() < 0.3) {

//  if (rng.uniform() < 0.05) {
        POT& frog1 = _Pw;

        double T1 = 5;
        // Stochastic Hill Climber algorithm instance
        //
        // Neighborhood
        ETTPneighborhood neighborhood;
        // Full evaluation function
        ProximityCostEval<POT> fullEval;
        // Neighbor evaluation function
        ETTPneighborEval neighEval;
        moSimpleCoolingSchedule<eoChromosome> cool(T1, 0.0001, 5, 0.1); // 39.53 p=0.1
//        moSimpleCoolingSchedule<eoChromosome> cool(T1, 0.0001, 5, 1); // SA

        moSA<ETTPneighbor> sa(neighborhood, fullEval, neighEval, cool);
        sa(frog1);

        _Pw.computeProximityCosts();
}
    //        _Pw.validate();

    //        cout << "originalPbFitness = " << originalPbFitness << endl;
//        cout << "originalPwFitness = " << originalPwFitness << endl;
//        cout << "[After] Pw Fitness " << _Pw.fitness() << endl;
    }
};
*/


///////////////////////////////////////////////////////////////////////////////
typedef eoVector<double, int> Individual;


class PeriodPermutationEval : public eoEvalFunc<Individual>
{
public:
    // Ctor
    PeriodPermutationEval(eoChromosome& _chrom) : chrom(_chrom) {}

    /**
     * Count the number of 1 in a bitString
     * @param _sol the solution to evaluate
     */
    void operator() (Individual& _solution) {
        // Evaluate the period ordering according to the current solution
        auto& periods = chrom.getPeriods();
        eoChromosome newChrom = chrom;
        auto& newPeriods = newChrom.getPeriods();
        // Fitness
        double fit = 0;
        for (unsigned int i = 0; i < _solution.size(); i++)
             newPeriods[i] = periods[_solution[i]];

        newChrom.computeProximityCosts();
        fit = newChrom.fitness();
        _solution.fitness(fit);

//        cout << "Eval ok, fitness = " << fit << endl;
    }

private:
    eoChromosome& chrom;
};



class eoSFLAPeriodEvolOperator : public eoBinOp<eoVector<double, int> > {

    typedef eoVector<double, int> EOT;

public:
    // Ctor
//    eoSFLAPeriodEvolOperator(eoChromosome& _chrom) : chrom(_chrom) {}

    bool operator()(EOT& _Pw, EOT const& _Pb) {
        // Within each submemeplex, the worst performance frog
        // is updated according to the following rule:
        //  the qth frog in a submemeplex U(q) is obtained
        //  by randomly selecting a subsequence in Pb to replace
        //  the corresponding position in Pw, while keeping the
        //  other positions in Pw unchanged or if violating the
        //  feasibility constraint, just randomly relocate the
        //  remain positions to form a new feasible solution.
        ////

//        cout << "eoSFLAPeriodEvolOperator" << endl;

//        cout << "_Pb: ";
//        copy(_Pb.begin(), _Pb.end(), ostream_iterator<int>(cout, " "));
//        cout << endl;

//        cout << "Initial _Pw: ";
//        copy(_Pw.begin(), _Pw.end(), ostream_iterator<int>(cout, " "));
//        cout << endl;


        // Get random sequence in Pb
        int low = rng.random(_Pb.size());
        int high = rng.random(_Pb.size());
        if (low > high)
            swap(low, high);


//        cout << "_Pb.size() = " << _Pb.size() << endl;
//        cout << "low = " << low << ", high = " << high << endl;

        // Replace the corresponding position in Pw, while keeping the
        // other positions in Pw unchanged or if violating the
        // feasibility constraint, just randomly relocate the
        // remain positions to form a new feasible solution.
        int sizeSequence = high-low+1;

//        cout << "sizeSequence = " << sizeSequence << endl;

        // Copy elements to be replaced in Pw
        vector<int> replacedIdxs(_Pw.begin()+low, _Pw.begin()+low+sizeSequence);

//        cout << "replacedIdxs: ";
//        copy(replacedIdxs.begin(), replacedIdxs.end(), ostream_iterator<int>(cout, " "));
//        cout << endl;

        // Verify feasibility and relocate if necessary the remainder elements
        // 1. Remove duplicated elements
        for (int i = low; i < low + sizeSequence; ++i) {
            // If exist in replacedIdxs, remove it
            replacedIdxs.erase(remove_if(replacedIdxs.begin(), replacedIdxs.end(), bind2nd(equal_to<int>(), _Pb[i])), replacedIdxs.end());
            // Mark as -1 in _Pw
            auto it = find_if(_Pw.begin(), _Pw.end(), bind2nd(equal_to<int>(), _Pb[i]));
            if (it != _Pw.end())
                *it = -1;
            else
                throw runtime_error("Error: Cannot find _Pw element in class eoSFLAPeriodEvolOperator");
        }

//        cout << "After removing from replacedIdxs: ";
//        copy(replacedIdxs.begin(), replacedIdxs.end(), ostream_iterator<int>(cout, " "));
//        cout << endl;
//        cout << "replacedIdxs.size() = " << replacedIdxs.size() << endl;


        // 2. Replace elements in Pw
        for (int i = low; i < low + sizeSequence; ++i)
            _Pw[i] = _Pb[i];

        // 3. Randomly select elements to insert
        random_shuffle(replacedIdxs.begin(), replacedIdxs.end());
        for (int i1 = 0, i2 = 0; i1 < _Pw.size(); ++i1) {
            if (_Pw[i1] == -1) {
                _Pw[i1] = replacedIdxs[i2++];
            }
        }

//        cout << "Final _Pw: ";
//        copy(_Pw.begin(), _Pw.end(), ostream_iterator<int>(cout, " "));
//        cout << endl;

        // Verify chromosome integrity
        int numPeriods = _Pb.size();
        while (--numPeriods >= 0) {
            auto it = _Pw.begin();
            bool found = false;
            do {
                it = find_if(it, _Pw.end(), bind2nd(equal_to<int>(), numPeriods));
                if (found == true && it != _Pw.end() ) {
                    std::ostringstream error_msg;
                    error_msg << "Error: Duplicated _Pw element " << *it << " - in class eoSFLAPeriodEvolOperator";
                    throw runtime_error(error_msg.str());
                }
                else if (found == false && it != _Pw.end())
                    found = true;
                else if (found == false && it == _Pw.end()) {
                    std::ostringstream error_msg;
                    error_msg << "Error: Missing period " << numPeriods << " in _Pw - in class eoSFLAPeriodEvolOperator";
                    throw runtime_error(error_msg.str());
                }
            }
            while (found == true && it == _Pw.end());
        }

//        cout << "eoSFLAPeriodEvolOperator Done" << endl;
    }
};


/*
 *
 *
////////////////////////////////////////////
// Code Paper ECTA 2013 LNCS
//
///////////////////////////////////////////


//////////////////////////////////////////////////////////////
//
// eoEvolutionOperator 2 - With variation operators + Micro SFLA
//
//////////////////////////////////////////////////////////////
class eoSFLAEvolOperator_2 : public eoBinOp<eoChromosome> {

    typedef eoChromosome POT;

public:
    bool operator()(POT& _Pw, POT const& _Pb) {
        double originalPbFitness = _Pb.fitness();
        double originalPwFitness = _Pw.fitness();

        auto newPw = _Pb;
//        POT& frog1 = newPw;


//        auto newPw = _Pw;

        cout << endl << endl << "Initial fitness = " << newPw.fitness() << endl;
//        cout << endl << endl << "Best  frog fitness = " << _Pb.fitness() << endl;
//        cout << "Worst frog fitness = " << _Pw.fitness() << endl;


        if (rng.uniform() < 0) {
//        if (rng.uniform() < 0.1) { // LAST
//        if (rng.uniform() < 0.3) {
//        if (rng.uniform() < 1) {

//            ETTPneighbor<

            // Kempe Chain neighbourhood
            ETTPKempeChainHeuristic kempeChain;
            kempeChain(newPw);

        }



         if (rng.uniform() < 0) {
//        if (rng.uniform() < 0.1) { // LAST
//        if (rng.uniform() < 0.3) {
//        if (rng.uniform() < 1) {

            // Insert exams
            int randNumPeriods = 1;

            for (int i = 1; i <= randNumPeriods; ++i) {
                // Get Frog1 periods
                vector<unordered_set<int> >& frog1Periods = newPw.getPeriods();
                // Generate a random index in frogs
                int randPeriod1 = rng.random(newPw.getNumPeriods());
                int randPeriod2 = rng.random(newPw.getNumPeriods());
                // Swap period exams
                auto aux = frog1Periods[randPeriod1];
                frog1Periods[randPeriod1] = frog1Periods[randPeriod2];
                frog1Periods[randPeriod2] = aux;
            }
        }





        if (rng.uniform() < 0) {
//         if (rng.uniform() < 0.05) {
//        if (rng.uniform() < 0.1) { // LAST
//            if (rng.uniform() < 0.3) { // best
//          if (rng.uniform() < 1) {

                // Neighborhood
                ETTPneighborhood neighborhood;
                // Full evaluation function
                ProximityCostEval<POT> fullEval;
                // Neighbor evaluation function
                ETTPneighborEval neighEval;

                // Yor83
//                moSimpleCoolingSchedule<eoChromosome> cool(0.1, 0.001, 2, 0.000001); // Yor83
//                moSimpleCoolingSchedule<eoChromosome> cool(0.01, 0.0001, 2, 0.000001); // Yor83
                moSimpleCoolingSchedule<eoChromosome> cool(0.002, 0.000001, 2, 0.001); // Yor83



//                moSimpleCoolingSchedule<eoChromosome> cool(0.01, 0.000001, 2, 0.001); 10.10
//                moSimpleCoolingSchedule<eoChromosome> cool(0.005, 0.000001, 2, 0.001); 10.07


//                moSimpleCoolingSchedule<eoChromosome> cool(0.002, 0.000001, 2, 0.001); // SA - 10.03 - hec
                // kfu - 14.15 - good diversity
                // lse - 11.17
                // rye - 9.44
                // sta - 157.28
                // tre - 8.24
                // uta - 3.92
                // ute - 25.28
                // yor - 34.96 (começou em 35.69)

                //-------------



//                moSimpleCoolingSchedule<eoChromosome> cool(1, 0.0001, 0, 0.00001); // 4.97

//                moSimpleCoolingSchedule<eoChromosome> cool(1, 0.001, 5, 0.0001); // Good
//                moSimpleCoolingSchedule<eoChromosome> cool(0.5, 0.01, 3, 0.001); // Good

                // SHC
//                moSimpleCoolingSchedule<eoChromosome> cool(0.1, 0.01, 0, 0.1);

//                moSimpleCoolingSchedule<eoChromosome> cool(0.1, 0.01, 0, 0.1); // 40.46
//                moSimpleCoolingSchedule<eoChromosome> cool(0.2, 0.01, 2, 0.02); // 40.98
//                moSimpleCoolingSchedule<eoChromosome> cool(50, 0.01, 0, 50); // 40.83

                moSA<ETTPneighbor> sa(neighborhood, fullEval, neighEval, cool);
//                cout << endl << "fitness before SA " << frog1.fitness() << endl;
                sa(newPw);
//                cout << "fitness after SA " << frog1.fitness() << endl;

        }
//        else {
//        if (true) { // LAST
        if (false) {
//        if (rng.uniform() < 0.1) {
            vector<unordered_set<int> >& pwPeriods = _Pw.getPeriods();
            // Get newPw periods
            vector<unordered_set<int> >& newPwPeriods = newPw.getPeriods();

            // Generate random number of periods to insert from Pw
            int randNumPeriods = 3;
//              int randNumPeriods = 10; // LAST
//            int randNumPeriods = 30;
//            int randNumPeriods = 1;
            for (int i = 1; i <= randNumPeriods; ++i) {
                // Get a random period from Pw
                int randPwPeriod = rng.random(_Pw.getNumPeriods());
                // Insert Pw random period exams (those that are feasible) into the new frog and remove duplicated exams

                // Generate a random index in the new frog
                int randNewPwPeriod = rng.random(newPw.getNumPeriods());
                // Pw exams from best day
                unordered_set<int>& pwExams = pwPeriods[randPwPeriod];

                // newPw exams from rand period
                unordered_set<int>& newPwExams = newPwPeriods[randNewPwPeriod];

                // Insert Pw period into newPw into that random position
                for (auto it = pwExams.begin(); it != pwExams.end(); ++it) {
                    // If the exam does not exist and is feasible then insert it
                    if (newPwExams.find(*it) == newPwExams.end() && newPw.isFeasiblePeriodExam(randNewPwPeriod, *it)) {
//                        if (newPwExams.find(*it) == newPwExams.end() && newPw.isFeasiblePeriodExam(bestDayIdx, *it)) {


                    newPwExams.insert(*it);
        //                ++numInsertedExams;
//                        cout << "Insert: " << *it << endl;
//                        cin.get();


                        /// TODO / SEE POSSIBLE OPTIMIZATIONS

                        // Remove duplicated exams
                        for (int i = 0; i < newPwPeriods.size(); ++i) {
                            if (i != randNewPwPeriod) {
        //                        if (i != bestDayIdx) {
                                unordered_set<int>& exams = newPwPeriods[i];
                                if (exams.find(*it) != exams.end()) {
                                    exams.erase(*it);
        //                            ++numRemovedExams;
        //                            cout << "Remove: " << *it << endl;
                                    break;
                                }
                            }
                        }

                    }
        //            else {
        //                notInserted.insert(*it);
        //                ++numNotInsertedExams;
        //            }

                }
            }
            newPw.computeProximityCosts();

            cout << "After Pw merging = " << newPw.fitness() << endl;
        }



//        if (true) {
//            newPw = _Pw;

//            vector<unordered_set<int> >const& pbPeriods = _Pb.getConstPeriods();
//            // Get newPw periods
//            vector<unordered_set<int> >& newPwPeriods = newPw.getPeriods();

//            // Generate random number of periods to insert from Pb
////            int randNumPeriods = 3;
////              int randNumPeriods = 10;
////            int randNumPeriods = 30;
//            int randNumPeriods = 1;
//            for (int i = 1; i <= randNumPeriods; ++i) {
//                // Get a random period from Pb
//                int randPbPeriod = rng.random(_Pb.getNumPeriods());
//                // Insert Pb random period exams (those that are feasible) into the new frog and remove duplicated exams
//                // Generate a random index in the frog
//                int randNewPwPeriod = rng.random(newPw.getNumPeriods());
//                unordered_set<int>const& pbExams = pbPeriods[randPbPeriod];
//                // newPw exams from rand period
//                unordered_set<int>& newPwExams = newPwPeriods[randNewPwPeriod];

//                // Insert Pb period into newPw into that random position
//                for (auto it = pbExams.begin(); it != pbExams.end(); ++it) {
//                    // If the exam does not exist and is feasible then insert it
//                    if (newPwExams.find(*it) == newPwExams.end() && newPw.isFeasiblePeriodExam(randNewPwPeriod, *it)) {

//                    newPwExams.insert(*it);
//        //                ++numInsertedExams;
////                        cout << "Insert: " << *it << endl;
////                        cin.get();


//                        /// TODO / SEE POSSIBLE OPTIMIZATIONS

//                        // Remove duplicated exams
//                        for (int i = 0; i < newPwPeriods.size(); ++i) {
//                            if (i != randNewPwPeriod) {
//        //                        if (i != bestDayIdx) {
//                                unordered_set<int>& exams = newPwPeriods[i];
//                                if (exams.find(*it) != exams.end()) {
//                                    exams.erase(*it);
//        //                            ++numRemovedExams;
//        //                            cout << "Remove: " << *it << endl;
//                                    break;
//                                }
//                            }
//                        }

//                    }
//        //            else {
//        //                notInserted.insert(*it);
//        //                ++numNotInsertedExams;
//        //            }

//                }
//            }
//            newPw.computeProximityCosts();

//            cout << "After Pw merging = " << newPw.fitness() << endl;
//        }


        if (rng.uniform() < 0) {
//        if (rng.uniform() < 0.1) { // LAST
//        if (rng.uniform() < 0.3) {
//        if (rng.uniform() < 1) {

           // Insert exams
           int randNumPeriods = 1;
//           int randNumPeriods = 3; // LAST

           for (int i = 1; i <= randNumPeriods; ++i) {
               // Get Frog1 periods
               vector<unordered_set<int> >& frog1Periods = newPw.getPeriods();
               // Generate a random index in frogs
               int randPeriod1 = rng.random(newPw.getNumPeriods());
               int randPeriod2 = rng.random(newPw.getNumPeriods());
               // Swap period exams
               auto aux = frog1Periods[randPeriod1];
               frog1Periods[randPeriod1] = frog1Periods[randPeriod2];
               frog1Periods[randPeriod2] = aux;
           }

           newPw.computeProximityCosts();

           cout << "Swap periods - fitness = " << newPw.fitness() << endl;

       }

        if (rng.uniform() < 0) {
//        if (rng.uniform() < 1) {  // With Kempe chain (twist) and SA, Yor83 - 35.59
            // Kempe Chain neighbourhood
            ETTPKempeChainHeuristic kempeChain;
            kempeChain(newPw);
        }


//        if (rng.uniform() < 0) {
//        if (rng.uniform() < 0.01) {
//        if (rng.uniform() < 0.1) {
//        if (rng.uniform() < 0.2) {
//        if (rng.uniform() < 0.5) {
        if (rng.uniform() < 1) {
            //
            // Neighborhood
            ETTPneighborhood neighborhood;
            // Full evaluation function
            ProximityCostEval<POT> fullEval;
            // Neighbor evaluation function
            ETTPneighborEval neighEval;
// LAST
//moSimpleCoolingSchedule<eoChromosome> cool2(0.1, 0.00001, 5, 0.0000001); // 35.09 from 51.916 + pw merging 51.78
//                                                                         // Yor83 - 34.79, with pr=0.1 Swap, m=10, N=3

//            moSimpleCoolingSchedule<eoChromosome> cool2(0.1, 0.01, 5, 0);

//            moSimpleCoolingSchedule<eoChromosome> cool2(0.1, 0.9, 5, 0.0000001); // Yor - 40 so Kempe modificado

//            moSimpleCoolingSchedule<eoChromosome> cool2(0.1, 0.001, 5, 0.0000001); // Hec - 10.3

            moSimpleCoolingSchedule<eoChromosome> cool2(0.1, 0.9, 5, 0); // Yor83 - 37.19, Hec - 10.16
                                                                         // Yor83 - 35.59 with Kempe twist
                                                                         // Yor83 - 35.77 so Kempe modificado

            moSA<ETTPneighbor> sa2(neighborhood, fullEval, neighEval, cool2);

            cout << "Before SA2 - fitness = " << newPw.fitness() << endl;
            cout << "SA2" << endl;
            sa2(newPw);
            cout << "After SA2 - fitness = " << newPw.fitness() << endl << endl << endl;

        }



        _Pw = newPw;
        _Pw.computeProximityCosts();


//        _Pw.validate();

        return true;
    }
*/


//////////////////////////////////////////////////////////////
//
// eoEvolutionOperator 2 - With variation operators + Micro SFLA
//
//////////////////////////////////////////////////////////////
class eoSFLAEvolOperator_2 : public eoBinOp<eoChromosome> {

    typedef eoChromosome POT;


    int temp;


public:

    void setTemperature(int _temp) {
        temp = _temp;
    }

    bool operator()(POT& _Pw, POT const& _Pb) {
        double originalPbFitness = _Pb.fitness();
        double originalPwFitness = _Pw.fitness();

        auto newPw = _Pb;

//        cout << endl << endl << "Initial fitness = " << newPw.fitness() << endl;
//        cout << endl << endl << "Best  frog fitness = " << _Pb.fitness() << endl;
//        cout << "Worst frog fitness = " << _Pw.fitness() << endl;



//         if (rng.uniform() < 0) {
////        if (rng.uniform() < 0.1) { // LAST
////        if (rng.uniform() < 0.3) {
////        if (rng.uniform() < 1) {

//            // Insert exams
//            int randNumPeriods = 1;

//            for (int i = 1; i <= randNumPeriods; ++i) {
//                // Get Frog1 periods
//                vector<unordered_set<int> >& frog1Periods = newPw.getPeriods();
//                // Generate a random index in frogs
//                int randPeriod1 = rng.random(newPw.getNumPeriods());
//                int randPeriod2 = rng.random(newPw.getNumPeriods());
//                // Swap period exams
//                auto aux = frog1Periods[randPeriod1];
//                frog1Periods[randPeriod1] = frog1Periods[randPeriod2];
//                frog1Periods[randPeriod2] = aux;
//            }
//        }




////        if (true) { // LAST
//        if (false) {
////        if (rng.uniform() < 0.1) {
//            vector<unordered_set<int> >& pwPeriods = _Pw.getPeriods();
//            // Get newPw periods
//            vector<unordered_set<int> >& newPwPeriods = newPw.getPeriods();

//            // Generate random number of periods to insert from Pw
//            int randNumPeriods = 3;
////              int randNumPeriods = 10; // LAST
////            int randNumPeriods = 30;
////            int randNumPeriods = 1;
//            for (int i = 1; i <= randNumPeriods; ++i) {
//                // Get a random period from Pw
//                int randPwPeriod = rng.random(_Pw.getNumPeriods());
//                // Insert Pw random period exams (those that are feasible) into the new frog and remove duplicated exams

//                // Generate a random index in the new frog
//                int randNewPwPeriod = rng.random(newPw.getNumPeriods());
//                // Pw exams from best day
//                unordered_set<int>& pwExams = pwPeriods[randPwPeriod];

//                // newPw exams from rand period
//                unordered_set<int>& newPwExams = newPwPeriods[randNewPwPeriod];

//                // Insert Pw period into newPw into that random position
//                for (auto it = pwExams.begin(); it != pwExams.end(); ++it) {
//                    // If the exam does not exist and is feasible then insert it
//                    if (newPwExams.find(*it) == newPwExams.end() && newPw.isFeasiblePeriodExam(randNewPwPeriod, *it)) {
////                        if (newPwExams.find(*it) == newPwExams.end() && newPw.isFeasiblePeriodExam(bestDayIdx, *it)) {


//                    newPwExams.insert(*it);
//        //                ++numInsertedExams;
////                        cout << "Insert: " << *it << endl;
////                        cin.get();


//                        /// TODO / SEE POSSIBLE OPTIMIZATIONS

//                        // Remove duplicated exams
//                        for (int i = 0; i < newPwPeriods.size(); ++i) {
//                            if (i != randNewPwPeriod) {
//        //                        if (i != bestDayIdx) {
//                                unordered_set<int>& exams = newPwPeriods[i];
//                                if (exams.find(*it) != exams.end()) {
//                                    exams.erase(*it);
//        //                            ++numRemovedExams;
//        //                            cout << "Remove: " << *it << endl;
//                                    break;
//                                }
//                            }
//                        }

//                    }
//        //            else {
//        //                notInserted.insert(*it);
//        //                ++numNotInsertedExams;
//        //            }

//                }
//            }
//            newPw.computeProximityCosts();

//            cout << "After Pw merging = " << newPw.fitness() << endl;
//        }



//        if (rng.uniform() < 0) {
////        if (rng.uniform() < 0.1) { // LAST
////        if (rng.uniform() < 0.3) {
////        if (rng.uniform() < 1) {

//           // Insert exams
//           int randNumPeriods = 1;
////           int randNumPeriods = 3; // LAST

//           for (int i = 1; i <= randNumPeriods; ++i) {
//               // Get Frog1 periods
//               vector<unordered_set<int> >& frog1Periods = newPw.getPeriods();
//               // Generate a random index in frogs
//               int randPeriod1 = rng.random(newPw.getNumPeriods());
//               int randPeriod2 = rng.random(newPw.getNumPeriods());
//               // Swap period exams
//               auto aux = frog1Periods[randPeriod1];
//               frog1Periods[randPeriod1] = frog1Periods[randPeriod2];
//               frog1Periods[randPeriod2] = aux;
//           }

//           newPw.computeProximityCosts();

//           cout << "Swap periods - fitness = " << newPw.fitness() << endl;

//       }


        // Neighborhood
        ETTPneighborhood neighborhood;
        // Full evaluation function
        ProximityCostEval<POT> fullEval;
        // Neighbor evaluation function
        ETTPneighborEval neighEval;
        // Cooling Schedule
        moSimpleCoolingSchedule<eoChromosome> cool(temp, 0.1, 1, temp-0.5);

        moSA<ETTPneighbor> sa(neighborhood, fullEval, neighEval, cool);

//            cout << "Before SA2 - fitness = " << newPw.fitness() << endl;
//            cout << "SA2" << endl;
//        sa(newPw);
//            cout << "After SA2 - fitness = " << newPw.fitness() << endl << endl << endl;


        // Kempe Chain neighbourhood
        ETTPKempeChainHeuristic kempeChain;
        kempeChain(newPw);


        if (rng.uniform() < 0) {
//                if (rng.uniform() < 0.1) { // LAST
        //        if (rng.uniform() < 0.3) {
        //        if (rng.uniform() < 1) {

                   // Insert exams
                   int randNumPeriods = 1;
        //           int randNumPeriods = 3; // LAST

                   for (int i = 1; i <= randNumPeriods; ++i) {
                       // Get Frog1 periods
                       vector<unordered_set<int> >& frog1Periods = newPw.getPeriods();
                       // Generate a random index in frogs
                       int randPeriod1 = rng.random(newPw.getNumPeriods());
                       int randPeriod2 = rng.random(newPw.getNumPeriods());
                       // Swap period exams
                       auto aux = frog1Periods[randPeriod1];
                       frog1Periods[randPeriod1] = frog1Periods[randPeriod2];
                       frog1Periods[randPeriod2] = aux;
                   }

                   newPw.computeProximityCosts();

//                   cout << "Swap periods - fitness = " << newPw.fitness() << endl;

               }



        _Pw = newPw;
        _Pw.computeProximityCosts();


//        _Pw.validate();

        return true;
    }




////////////////////////////////////////////////////////////////
////
//// eoEvolutionOperator 2 -
////
////////////////////////////////////////////////////////////////
//class eoSFLAEvolOperator_2 : public eoBinOp<eoChromosome> {

//    typedef eoChromosome POT;
//public:
//    bool operator()(POT& _Pw, POT const& _Pb) {
//        double originalPbFitness = _Pb.fitness();
//        double originalPwFitness = _Pw.fitness();

//        auto newPw = _Pb;
////        POT& frog1 = newPw;


////        auto newPw = _Pw;

//        cout << endl << endl << "Initial fitness = " << newPw.fitness() << endl;
////        cout << endl << endl << "Best  frog fitness = " << _Pb.fitness() << endl;
////        cout << "Worst frog fitness = " << _Pw.fitness() << endl;


///*
//        if (rng.uniform() < 0.1) {
//            vector<unordered_set<int> >& pwPeriods = _Pw.getPeriods();
//            // Get newPw periods
//            vector<unordered_set<int> >& newPwPeriods = newPw.getPeriods();

//            // Generate random number of periods to insert from Pw
//            int randNumPeriods = 3;
////              int randNumPeriods = 10; // LAST
////            int randNumPeriods = 30;
////            int randNumPeriods = 1;
//            for (int i = 1; i <= randNumPeriods; ++i) {
//                // Get a random period from Pw
//                int randPwPeriod = rng.random(_Pw.getNumPeriods());
//                // Insert Pw random period exams (those that are feasible) into the new frog and remove duplicated exams

//                // Generate a random index in the new frog
//                int randNewPwPeriod = rng.random(newPw.getNumPeriods());
//                // Pw exams from best day
//                unordered_set<int>& pwExams = pwPeriods[randPwPeriod];

//                // newPw exams from rand period
//                unordered_set<int>& newPwExams = newPwPeriods[randNewPwPeriod];

//                // Insert Pw period into newPw into that random position
//                for (auto it = pwExams.begin(); it != pwExams.end(); ++it) {
//                    // If the exam does not exist and is feasible then insert it
//                    if (newPwExams.find(*it) == newPwExams.end() && newPw.isFeasiblePeriodExam(randNewPwPeriod, *it)) {
////                        if (newPwExams.find(*it) == newPwExams.end() && newPw.isFeasiblePeriodExam(bestDayIdx, *it)) {


//                    newPwExams.insert(*it);
//        //                ++numInsertedExams;
////                        cout << "Insert: " << *it << endl;
////                        cin.get();


//                        /// TODO / SEE POSSIBLE OPTIMIZATIONS

//                        // Remove duplicated exams
//                        for (int i = 0; i < newPwPeriods.size(); ++i) {
//                            if (i != randNewPwPeriod) {
//        //                        if (i != bestDayIdx) {
//                                unordered_set<int>& exams = newPwPeriods[i];
//                                if (exams.find(*it) != exams.end()) {
//                                    exams.erase(*it);
//        //                            ++numRemovedExams;
//        //                            cout << "Remove: " << *it << endl;
//                                    break;
//                                }
//                            }
//                        }

//                    }
//        //            else {
//        //                notInserted.insert(*it);
//        //                ++numNotInsertedExams;
//        //            }

//                }
//            }
//            newPw.computeProximityCosts();

//            cout << "After Pw merging = " << newPw.fitness() << endl;
//        }
//*/


//        if (true) {
//            newPw = _Pw;

//            vector<unordered_set<int> >const& pbPeriods = _Pb.getConstPeriods();
//            // Get newPw periods
//            vector<unordered_set<int> >& newPwPeriods = newPw.getPeriods();

//            // Generate random number of periods to insert from Pb
////            int randNumPeriods = 3;
////              int randNumPeriods = 10;
////            int randNumPeriods = 30;
//            int randNumPeriods = 1;
//            for (int i = 1; i <= randNumPeriods; ++i) {
//                // Get a random period from Pb
//                int randPbPeriod = rng.random(_Pb.getNumPeriods());
//                // Insert Pb random period exams (those that are feasible) into the new frog and remove duplicated exams
//                // Generate a random index in the frog
//                int randNewPwPeriod = rng.random(newPw.getNumPeriods());
//                unordered_set<int>const& pbExams = pbPeriods[randPbPeriod];
//                // newPw exams from rand period
//                unordered_set<int>& newPwExams = newPwPeriods[randNewPwPeriod];

//                // Insert Pb period into newPw into that random position
//                for (auto it = pbExams.begin(); it != pbExams.end(); ++it) {
//                    // If the exam does not exist and is feasible then insert it
//                    if (newPwExams.find(*it) == newPwExams.end() && newPw.isFeasiblePeriodExam(randNewPwPeriod, *it)) {

//                    newPwExams.insert(*it);
//        //                ++numInsertedExams;
////                        cout << "Insert: " << *it << endl;
////                        cin.get();


//                        /// TODO / SEE POSSIBLE OPTIMIZATIONS

//                        // Remove duplicated exams
//                        for (int i = 0; i < newPwPeriods.size(); ++i) {
//                            if (i != randNewPwPeriod) {
//        //                        if (i != bestDayIdx) {
//                                unordered_set<int>& exams = newPwPeriods[i];
//                                if (exams.find(*it) != exams.end()) {
//                                    exams.erase(*it);
//        //                            ++numRemovedExams;
//        //                            cout << "Remove: " << *it << endl;
//                                    break;
//                                }
//                            }
//                        }

//                    }
//        //            else {
//        //                notInserted.insert(*it);
//        //                ++numNotInsertedExams;
//        //            }

//                }
//            }
//            newPw.computeProximityCosts();

//            cout << "After Pw merging = " << newPw.fitness() << endl;
//        }



////        if (rng.uniform() < 0) {
//        if (rng.uniform() < 0.1) { // LAST
////        if (rng.uniform() < 0.3) {
////        if (rng.uniform() < 1) {

//           // Insert exams
//           int randNumPeriods = 1;
////           int randNumPeriods = 3; // LAST

//           for (int i = 1; i <= randNumPeriods; ++i) {
//               // Get Frog1 periods
//               vector<unordered_set<int> >& frog1Periods = newPw.getPeriods();
//               // Generate a random index in frogs
//               int randPeriod1 = rng.random(newPw.getNumPeriods());
//               int randPeriod2 = rng.random(newPw.getNumPeriods());
//               // Swap period exams
//               auto aux = frog1Periods[randPeriod1];
//               frog1Periods[randPeriod1] = frog1Periods[randPeriod2];
//               frog1Periods[randPeriod2] = aux;
//           }

//           newPw.computeProximityCosts();

//           cout << "Swap periods - fitness = " << newPw.fitness() << endl;

//       }


////        if (rng.uniform() < 0) {
////        if (rng.uniform() < 0.01) {
//        if (rng.uniform() < 0.1) {
////        if (rng.uniform() < 0.2) {
////        if (rng.uniform() < 0.5) {
////        if (rng.uniform() < 1) {
//            //
//            // Neighborhood
//            ETTPneighborhood neighborhood;
//            // Full evaluation function
//            ProximityCostEval<POT> fullEval;
//            // Neighbor evaluation function
//            ETTPneighborEval neighEval;

////            moSimpleCoolingSchedule<eoChromosome> cool2(0.1, 0.00001, 5, 0.0000001); // 35.09 from 51.916 + pw merging 51.78
//                                                                            // Yor83 - 34.79, with pr=0.1 Swap, m=10, N=3

//            moSimpleCoolingSchedule<eoChromosome> cool2(0.1, 0.001, 5, 0.0000001);

//            moSA<ETTPneighbor> sa2(neighborhood, fullEval, neighEval, cool2);
//            cout << "Before SA2 - fitness = " << newPw.fitness() << endl;
//            cout << "SA2" << endl;
//            sa2(newPw);
//            cout << "After SA2 - fitness = " << newPw.fitness() << endl << endl << endl;

//        }

//        if (rng.uniform() < 0) {
////        if (rng.uniform() < 0.1) { // LAST
////        if (rng.uniform() < 0.3) {
////        if (rng.uniform() < 1) {

//            // Kempe Chain neighbourhood
//            ETTPKempeChainHeuristic kempeChain;
//            kempeChain(newPw);
//        }


//        _Pw = newPw;
//        _Pw.computeProximityCosts();

//        _Pw.validate();

//        return true;
//    }


private:

    void createChromosomeFromPeriodPermutation(eoChromosome& _chrom, Individual& _solution) {
        auto& periods = _chrom.getPeriods();
        eoChromosome newChrom = _chrom;
        auto& newPeriods = newChrom.getPeriods();
        for (unsigned int i = 0; i < _solution.size(); i++)
             newPeriods[i] = periods[_solution[i]];

        _chrom.setPeriods(newPeriods);
        _chrom.computeProximityCosts();
    }

};


#endif // EOEVOLUTIONOPERATOR_H

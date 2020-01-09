#ifndef EOSFLA_H
#define EOSFLA_H


#include "MOEA.h"
#include "moeoLocalSearchUpdater.h"
#include "ETTPneighborhood.h"
#include "ETTPneighborEval.h"
//#include "algo/moTS.h"
#include "algo/moSA.h"



/**
 * Shuffled Frog-Leaping Algorithm
 *
 * Global exploration
 *
 *  Step 0 Initialize. Select m and n, where m is the number of memeplexes and n is the number of frogs in each
 *  memeplex. Therefore, the total sample size, F, in the swamp is given by F = mn.
 *
 *  Step 1 Generate a virtual population. Sample F virtual frogs, U(1),...,U(F).
 *  Compute the performance value f(i) for each frog U.
 *
 *  Step 2 Rank frogs. Sort the F frogs in order of decreasing performance value. Store them in an array X={U(i),
 *  f(i), i=1,...,F} so that i=1 represents the frog with the best performance value. Record the best frog’s
 *  position, PX, in the entire population (F frogs; where PX =U(1)).
 *
 *  Step 3 Partition frogs into memeplexes. Partition array X into m memeplexes Y1, Y2,...,Ym, each containing
 *  n frogs. E.g., for m=3, rank 1 goes to memeplex 1, rank 2 goes to memeplex 2, rank 3 goes to memeplex 3,
 *  rank 4 goes to memeplex 1, and so on (Fig. 2).
 *
 *  Step 4 Memetic evolutions within each memeplex. Evolve each memeplex Y(l), l=1,..., m.
 *
 *  After partitioning frogs to m memeplexes, evolve each memeplex and each of them should iterate N times. After
 *  the memeplexes have been evolved, the algorithm returns to the global exploration for shuffling.
 *
 *  Local exploration: frog-leaping algorithm
 *  Step 4-0 Set im=0, where im counts the number of memeplexes and will be compared with
 *  the total number m of memeplexes. Set iN=0, where iN counts the number of evolutionary
 *  steps and will be compared with the maximum number N of steps to be completed within each
 *  memeplex. Within each memeplex (Fig. 3b), the frogs with the best and the worst fitness are
 *  identified as Pb and Pw, respectively. Also, the frog with the global best fitness is identified as
 *  Pg. Then, an evolution process is applied to improve only the frog with the worst fitness (i.e.,
 *  not all frogs) in each cycle.
 *  Step 4-1 Set im=im+1.
 *  Step 4-2 Set iN=iN+1.
 *  Step 4-3 Improve the worst frog’s position.
 *  The position of the frog with the worst fitness is adjusted as follows:
 *     Change in frog position (Di) = rand() x (Pb - Pw)     (5)
 *
 *     New position Pw = current position Pw + Di;
 *                          (Dmax geq Di geq Dmax)           (6)
 *
 *  where rand() is a random number between 0 and 1 and Dmax is the maximum allowed change in a frog’s position.
 *  Step 4-4 If this process produces a better frog (solution), it replaces the worst frog. Otherwise, the
 *  calculations in Eqs. 5 and 6 are repeated with respect to the global best frog (i.e., Pg replaces Pb).
 *  Step 4-5 If no improvement becomes possible in this latter case, then a new solution is
 *  randomly generated to replace the worst frog with another frog having any arbitrary fitness (as shown in Fig. 3b).
 *  Step 4-6 If iN < N, go to step 4-2.
 *  Step 4-7 If im<m, go to step 4-1. Otherwise, return to the global search to shuffle memeplexes.
 *
 *  Step 5 Shuffle memeplexes. After a defined number of memetic evolutionary steps within each memeplex,
 *  replace Y1,...,Ym into X such that X={Yk, k = 1..., m}. Sort X in order of decreasing performance value.
 *  Update the population the best frog’s position PX.
 *
 *  Step 6 Check convergences. If the convergence criteria are satisfied, stop. Otherwise, return to step 3.
 *  Typically, the decision on when to stop is made by a prespecified number of consecutive time loops when
 *  at least one frog carries the “best memetic pattern” without change. Alternatively, a maximum total number of
 *  function evaluations can be defined (Fig. 4).
 */

#include <eoAlgo.h>
#include "eoSFLAInitializer.h"
#include <eoContinue.h>
#include <eoEvalFunc.h>
#include <vector>
#include <boost/shared_ptr.hpp>
#include "Chromosome.h"


using namespace std;
using namespace boost;



/**
    This is a generic class for Shuffled Frog-Leaping algorithms. There
    is only one operator defined, which takes a population and does stuff to
    it. It needn't be a complete algorithm, can be also a step of an
    algorithm. Only used for mono-objective cases.

    @ingroup Algorithms
*/

////
//// TODO: GENERIC CLASS eoSFLA
////

//template <class POT>
//class eoSFLA : public eoAlgo<POT>
class eoSFLA : public eoAlgo<Chromosome>
{
public:

    typedef Chromosome POT;

    /** Constructor
     * @param n - number of memeplexes
     * @param N - number of frogs in each memeplex
     * @param F - total sample size in the swamp and given by F = mN
     * @param _init - An eoInit that initializes each frog (solution)
     * @param _continuator - An eoContinue that manages the stopping criterion and the checkpointing system
     * @param _eval - An eoEvalFunc: the evaluation performer
     */
//    eoSFLA(int _m, int _N, int _F, eoInit<POT>& _init, eoContinue<POT>& _continuator, eoEvalFunc<POT>& _eval) :
        eoSFLA(int _m, int _N, int _F, ETTPInit& _init, eoContinue<POT>& _continuator, eoEvalFunc<POT>& _eval) :
               m(_m), N(_N), F(_F), init(_init), continuator(_continuator), eval(_eval)
    {}

    // Apply SFLA to the population
    virtual void operator() (eoPop<POT> &_pop)
    {
        try
        {
            // Step 0 and 1 - Initialization and generation of virtual population of frogs
            //   The population is already initialized by the 'init' object
            //
            // Step 2 - Frog ranking
            rankFrogs(_pop);
            // Create memeplexe structure
            vector<vector<POT> > memplexes = createMemeplexes();

            do
            {
                // Step 3 Partition frogs into memeplex
                partitionFrogsIntoMemplexes(_pop, memplexes);
                // Step 4 Frog-leaping algorithm for local search
                frogLeapingLocalSearch(memplexes);
                // Step 5 Shuffle memeplexes. After a defined number of memetic evolutionary steps within each memeplex,
                // replace Y1,...,Ym into X such that X={Yk, k = 1..., m}. Sort X in order of decreasing performance value.
                // Update the population the best frog’s position PX.
                shuffleMemplexes(memplexes, _pop);
                // Step 6 Check convergences. If the convergence criteria are satisfied, stop. Otherwise, return to step 3.
                // Typically, the decision on when to stop is made by a prespecified number of consecutive time loops when
                // at least one frog carries the “best memetic pattern” without change. Alternatively, a maximum total number of
                // function evaluations can be defined (Fig. 4).
            }
            while (continuator(_pop));
        }
        catch (std::exception& _e)
        {
            std::string s = _e.what();
            s.append (" in eoSFLA");
            throw std::runtime_error(s);
        }

        cout << "Best frog overall fitness = " << Pg.fitness() << endl;

    }


    POT getBestFrog() {
        return Pg;
    }

protected:
    // m is the number of memeplexes
    int m;
    // N is the number of frogs in each memeplex
    int N;
    // The total sample size, F, in the swamp is given by F = mN
    int F;
    // SFLA Initializer
    //eoSFLAInitializer<POT>& init;
    //eoInit<POT>& init;
    ETTPInit& init;

    // Continuator
    eoContinue<POT>& continuator;
    // Evaluation function
    eoEvalFunc<POT>& eval;
    // The frog with the global best fitness is identified as Pg.
    POT Pg;


    /////////////////////////////////////////////////////////////////////////////
    // Auxiliary methods
    /////////////////////////////////////////////////////////////////////////////

    void updateGlobalBestFrog(POT _chrom) {
        // Pre-condition: vector _pop is sorted by fitness
        // Global best frog is located at index 0 in the population vector
        Pg = _chrom;
    }
    POT getMemeplexBestFrog(vector<POT>& _memeplex) {
        // Pre-condition: vector _memeplex is sorted by fitness
        // Memeplex best frog is located at index 0 in the memeplex vector
        return _memeplex[0];
    }
    POT getMemeplexWorstFrog(vector<POT>& _memeplex) {
        // Pre-condition: vector _memeplex is sorted by fitness
        // Memeplex worst frog is located at the last index in the memeplex vector
        return _memeplex[_memeplex.size()-1];
    }

    struct FitnessCmp {
        bool operator()(Chromosome const& _chrom1, Chromosome const& _chrom2) {
            return _chrom1.fitness() < _chrom2.fitness(); // Ascending order because we're minimizing
        }
    };

    void rankFrogs(eoPop<POT>& _pop) {

        cout << "_pop.size() = " << _pop.size() << endl;
        // Step 2 Rank frogs. Sort the F frogs in order of decreasing performance value. Store them in an array X={U(i),
        // f(i), i=1,...,F} so that i=1 represents the frog with the best performance value. Record the best frog’s
        // position, PX, in the entire population (F frogs; where PX =U(1)).
        sort(_pop.begin(), _pop.end(), FitnessCmp());
        // Record the best frog’s position, PX, in the entire population (F frogs; where PX =U(1)).
        //   The best frog is just the 0 index frog
        updateGlobalBestFrog(_pop[0]);

        cout << "After ranking" << endl;

        cout << "Best frog fitness: " << _pop[0] << endl;

//        cin.get();
    }

    /// TODO: Optimizar copia do vector no retorno

    vector<vector<POT> > createMemeplexes() {

        cout << "createMemeplexes" << endl;
        cout << "m = " << m << endl;

        vector<vector<POT> > memeplexes;
        for (int i = 0; i < m; ++i)
            memeplexes.push_back(vector<POT>(N));

//        cout << memplexes.size() << endl;
//        for (int i = 0; i < m; ++i)
//            cout << "memplex size = " << memplexes[i].size() << endl;

//        cin.get();

        return memeplexes;
    }

    void partitionFrogsIntoMemplexes(eoPop<POT>& _pop, vector<vector<POT> >& _memeplexes) {
        // Step 3 Partition frogs into memeplex. Partition array X into m memeplexes Y1, Y2,...,Ym, each containing
        // n frogs. E.g., for m=3, rank 1 goes to memeplex 1, rank 2 goes to memeplex 2, rank 3 goes to memeplex 3,
        // rank 4 goes to memeplex 1, and so on (Fig. 2).

//        cout << "partitionFrogsIntoMemeplexes" << endl;

        vector<int> memeplexIndexes(m);
        for (int i = 0; i < _pop.size(); ++i) {
            // Each entry in each memeplex references a solution contained in population vector
            _memeplexes[i%m][memeplexIndexes[i%m]] = _pop[i];
            // Increment memeplex
            ++memeplexIndexes[i%m];
        }

//        cout << "List Memeplexes" << endl;

//        // Iterate over the memeplexes
//        for (int im = 0; im < m; ++im) {
//            cout << "Memeplexe " << im << endl;

//            // Perform N evolutionary steps
//            for (int iN = 0; iN < N; ++iN) {
//                cout << "\t" << _memeplexes[im][iN].fitness() << endl;
//            }
//        }

//        cin.get();

    }

    void shuffleMemplexes(vector<vector<POT> >& _memeplexes, eoPop<POT>& _pop) {
        // Step 5 Shuffle memeplexes. After a defined number of memetic evolutionary steps within each memeplex,
        // replace Y1,...,Ym into X such that X={Yk, k = 1..., m}. Sort X in order of decreasing performance value.
        // Update the population the best frog’s position PX.

//        cout << "shuffleMemeplexes" << endl;

        // Perhaps Missing...

        int k = 0;
        // Iterate over the memeplexes
        for (int im = 0; im < m; ++im) {
            // Perform N evolutionary steps
            for (int iN = 0; iN < N; ++iN) {
                _pop[k++] = _memeplexes[im][iN];
            }
        }


        // As the memeplexes contain references to population we don't have to combine mememplxes into population vector.
        // We just need to sort the population vector.
        sort(_pop.begin(), _pop.end(), FitnessCmp());
        // Record the best frog’s position, PX, in the entire population (F frogs; where PX =U(1)).
        //   The best frog is just the 0 index frog
        updateGlobalBestFrog(_pop[0]);
    }

    void frogLeapingLocalSearch(vector<vector<POT> >& _memeplexes) {
        // Step 4 Memetic evolutions within each memeplex. Evolve each memeplex Y(l), l=1,..., m.
        // After partitioning frogs to m memeplexes, evolve each memeplex and each of them should iterate N times. After
        // the memeplexes have been evolved, the algorithm returns to the global exploration for shuffling.
        //
        // Local exploration: frog-leaping algorithm
        //
        // Step 4-0 Set im = 0, where im counts the number of memeplexes and will be compared with the total number m of memeplexes.
        // Set iN = 0, where iN counts the number of evolutionary steps and will be compared with the maximum number N of steps
        // to be completed within each memeplex. Within each memeplex (Fig. 3b), the frogs with the best and the worst fitness are
        // identified as Pb and Pw, respectively. Also, the frog with the global best fitness is identified as Pg.
        // Then, an evolution process is applied to improve only the frog with the worst fitness (i.e., not all frogs) in each cycle.
        // Step 4-1 Set im = im+1.
        // Step 4-2 Set iN = iN+1.
        // Step 4-3 Improve the worst frog’s position.
        // The position of the frog with the worst fitness is adjusted as follows:
        //    Change in frog position (Di) = rand() x (Pb - Pw)     (5)
        //
        //    New position Pw = current position Pw + Di;
        //                         (Dmax geq Di geq Dmax)           (6)
        //
        // where rand() is a random number between 0 and 1 and Dmax is the maximum allowed change in a frog’s position.
        // Step 4-4 If this process produces a better frog (solution), it replaces the worst frog. Otherwise, the
        // calculations in Eqs. 5 and 6 are repeated with respect to the global best frog (i.e., Pg replaces Pb).
        // Step 4-5 If no improvement becomes possible in this latter case, then a new solution is
        // randomly generated to replace the worst frog with another frog having any arbitrary fitness (as shown in Fig. 3b).
        // Step 4-6 If iN < N, go to step 4-2.
        // Step 4-7 If im<m, go to step 4-1. Otherwise, return to the global search to shuffle memeplexes.
        //////

//        cout << "frogLeapingLocalSearch" << endl;

        // Iterate over the memeplexes
        for (int im = 0; im < m; ++im) {
            // Perform N evolutionary steps
            for (int iN = 0; iN < N; ++iN) {
                // Identify as Pb and Pw, respectively, the frogs with the best and the worst fitness.
                // Also, the frog with the global best fitness is identified as Pg.
                // Then, an evolution process is applied to improve only the frog with
                // the worst fitness (i.e., not all frogs) in each cycle.
                ////


//                cin.get();

                POT Pb = getMemeplexBestFrog(_memeplexes[im]);
                POT Pw = getMemeplexWorstFrog(_memeplexes[im]);
                // Step 4-3 Improve the worst frog’s position.
                POT newFrog = improveWorstFrogPosition(Pb, Pw);
                // Step 4-4 If this process produces a better frog (solution), it replaces the worst frog.
                // Otherwise, the calculations in Eqs. 5 and 6 are repeated with respect to the global best frog (i.e., Pg replaces Pb).
                // Step 4-5 If no improvement becomes possible in this latter case, then a new solution is
                // randomly generated to replace the worst frog with another frog having any arbitrary fitness (as shown in Fig. 3b).
                if (newFrog.fitness() < Pw.fitness()) {
                    replaceWorstFrog(newFrog, _memeplexes[im]);
                }
                else {
                    // Repeat computations with respect to the global best frog (i.e., Pg replaces Pb).
                    POT newFrog = improveWorstFrogPosition(Pg, Pw);
                    if (newFrog.fitness() < Pw.fitness()) {
                        replaceWorstFrog(newFrog, _memeplexes[im]);
                    }
                    else {
                        // Randomly generate a solution to replace the worst frog with another frog having any arbitrary fitness
                        POT newFrogRandom = generateRandomFrogPosition();
//                        int l1 = rng.random(Pg.getNumPeriods());
//                        int l2 = rng.random(Pg.getNumPeriods());
//                        int randFrog = rng.random(_memeplexes[im].size());
//                        POT newFrogRandom = _memeplexes[im][randFrog];

//                        auto aux = newFrogRandom.getPeriods()[l1];
//                        newFrogRandom.getPeriods()[l1] = newFrogRandom.getPeriods()[l2];
//                        newFrogRandom.getPeriods()[l2] = aux;

                        // 15 x 60 x 100 - 5.316
                        // 20 x 100 x 300 - hum 4.48 e saturou


                        replaceWorstFrog(newFrogRandom, _memeplexes[im]);
                    }
                }
                // Update Global frog
                POT bestMemeplexFrog = getMemeplexBestFrog(_memeplexes[im]);
                if (bestMemeplexFrog.fitness() < Pg.fitness())
                    updateGlobalBestFrog(bestMemeplexFrog);

            }
        }
    }

    POT generateRandomFrogPosition() {
//        boost::shared_ptr<POT> ptrNewFrogRandom(new POT());
        POT newFrogRandom;

//        cout << "improveWorstFrogPosition" << endl;

        init(newFrogRandom);

//        cout << "[improveWorstFrogPosition] after init random solution" << endl;

//        cin.get();


        return newFrogRandom;
    }

    void replaceWorstFrog(POT& newFrog, vector<POT>& _memeplex) {

//        cout << "replaceWorstFrog" << endl;

//        cout << "Before replace" << endl;

//        listMemeplex(_memeplex);

        // Replace worst frog by new one and then orderly relocte new frog in the memeplex vector
        _memeplex[_memeplex.size()-1] = newFrog;
        int i;
        for (i = _memeplex.size()-2; i >= 0; --i) {
            if (_memeplex[i].fitness() > _memeplex[i+1].fitness()) {
                // Swap entries
                POT aux = _memeplex[i];
                _memeplex[i] = _memeplex[i+1];
                _memeplex[i+1] = aux;
            }
            else // Stop. The memeplex is ordered.
                break;
        }
        if (_memeplex[0].fitness() < Pg.fitness())
            updateGlobalBestFrog(_memeplex[0]);


//        cout << _memeplex[0].fitness() << endl;
//        cout << "Best frog overall fitness = " << Pg.fitness() << endl;

//        cout << "After replace" << endl;

        listMemeplex(_memeplex);

//        cin.get();
    }

    void listMemeplex(vector<POT> const& _memeplex) {

//        cout << "List Memeplex" << endl;

        // Perform N iterations
        for (int iN = 0; iN < N; ++iN) {
            cout << "\t" << _memeplex[iN].fitness() << endl;
        }
    }


    POT improveWorstFrogPosition(POT _Pb, POT _Pw) {

//        cout << "improveWorstFrogPosition" << endl;

//        POT newFrog;

//        // Just for testing

//        init(newFrog);
//        cin.get();
//        return newFrog;



        examsExchangeOperator(_Pb, _Pw);




//        // Stochastic Hill Climber algorithm instance
//        //
//        // Neighborhood
//        ETTPneighborhood neighborhood;
//        // Full evaluation function
//    //    NumClashesEval fullEval;
//        ProximityCostEval fullEval;

//        // Neighbor evaluation function
//        ETTPneighborEval neighEval;
//        // tmax
//        int tmax = 10;
//    //    int tmax = 3;
//        // Temperature T
//        double tempT = 0.0001;
//        //moSHC(Neighborhood& _neighborhood, eoEvalFunc<EOT>& _fullEval, moEval<Neighbor>& _eval,
//        //double _tmax, double _tempT, ProblSense _sense)
//        moSHC<ETTPneighbor> shc(neighborhood, fullEval, neighEval, tmax, tempT, minimize);

//        cout << "Fitness before " << _Pw.fitness() << endl;
//        shc(_Pw);
//        cout << "Fitness after " << _Pw.fitness() << endl;

//        // Stochastic Hill Climber algorithm instance
//        //
//        // Neighborhood
//        ETTPneighborhood neighborhood;
//        // Full evaluation function
//        ProximityCostEval fullEval;

//        // Neighbor evaluation function
//        ETTPneighborEval neighEval;
//        // tmax
////                int tmax = 10;
//        int tmax = 3;
//        // Temperature T
//        double tempT = 0.0001;
//        //moSHC(Neighborhood& _neighborhood, eoEvalFunc<EOT>& _fullEval, moEval<Neighbor>& _eval,
//        //double _tmax, double _tempT, ProblSense _sense)
//        moSHC<ETTPneighbor> shc(neighborhood, fullEval, neighEval, tmax, tempT, minimize);

////        cout << "Fitness before " << _memeplexes[im][0].fitness() << endl;

//        shc(_Pw);
//        _Pw.computeProximityCosts();
//        _Pw.fitness(_Pw.getProximityCost());

//        cout << "Fitness after " << _memeplexes[im][0].fitness() << endl;


        return _Pw;

    }


//    void examsExchangeOperator(POT const& _Pb, POT& _Pw) {
//        cout << "Pb Fitness " << _Pb.fitness() << endl;
//        cout << "Pw Fitness " << _Pw.fitness() << endl;

//        double originalPbFitness = _Pb.fitness();
//        double originalPwFitness = _Pw.fitness();
//        int originalPwNumPeriods = _Pw.getNumPeriods();

//        // Generate random period interval in frog _Pb
//        int low, high;
//        int i = 0;
//        do {
//            low = rng.random(_Pb.getNumPeriods());
////            high = rng.random(_Pb.getNumPeriods());
//            high = low+1;

//            ++i;
////        } while (low > high);
//        } while (!(high < _Pb.getNumPeriods()));


//        cout << "i = " << i << endl;
//        cout << "low = " << low << ", high = " << high << endl;
//        // Get Pb periods
//        vector<unordered_set<int> > const& pbPeriods = _Pb.getConstPeriods();
//        // Get Pw periods
//        vector<unordered_set<int> >& pwPeriods = _Pw.getPeriods();
//        // Backup Pw original exams and replace them by Pb periods in the specified interval
//        vector<unordered_set<int> > pwReplacedPeriods;
////        for (int i = low; i <= high; ++i) {
////            pwReplacedPeriods.push_back(pwPeriods[i]);
////            pwPeriods[i] = pbPeriods[i]; // Copy periods from Pb
////        }

//        for (int i = 0; i < pwPeriods.size(); ++i) {
//            pwReplacedPeriods.push_back(pwPeriods[i]);
//            pwPeriods[i] = unordered_set<int>(0);
//        }

//        for (int i = low; i <= high; ++i) {
////            pwReplacedPeriods.push_back(pwPeriods[i]);
//            pwPeriods[i] = pbPeriods[i]; // Copy periods from Pb
//        }

//        int size = 0;
//        cout << "pbPeriods in interval [" << low << "," << high << "]" << endl;
//        for (int i = low; i <= high; ++i) {
//            for (unordered_set<int>::const_iterator it = pbPeriods[i].begin(); it != pbPeriods[i].end(); ++it) {
//                cout << *it << " ";
//                ++size;
//            }
//            cout << endl;
//        }
//        cout << "pbPeriods #exams in interval " << size << endl;


//        size = 0;
//        cout << "pwReplacedPeriods" << endl;
//        for (i = 0; i < pwReplacedPeriods.size(); ++i) {
//            for (unordered_set<int>::iterator it = pwReplacedPeriods[i].begin(); it != pwReplacedPeriods[i].end(); ++it) {
//                cout << *it << " ";
//                ++size;
//            }
//            cout << endl;
//        }
//        cout << "pwReplacedPeriods #exams " << size << endl;

//        size = 0;
//        cout << "pwPeriods in interval [" << low << "," << high << "]" << endl;
//        for (int i = low; i <= high; ++i) {
//            for (unordered_set<int>::iterator it = pwPeriods[i].begin(); it != pwPeriods[i].end(); ++it) {
//                cout << *it << " ";
//                ++size;
//            }
//            cout << endl;
//        }
//        cout << "pwPeriods #exams in interval " << size << endl;


//        int removed = 0;


//        // To ensure the feasibility of chromosomes after the
//        // crossover, duplicated exams are deleted. These exams are
//        // removed from the original periods, while the newly inserted
//        // periods are left intact.
//        for (int j = low; j <= high; ++j) {
//            unordered_set<int> const& examsToRemove = pbPeriods[j];
//            // Remove each exam in examsToRemove vector
//            for (unordered_set<int>::const_iterator itRem = examsToRemove.begin(); itRem != examsToRemove.end(); ++itRem) {
//                int i;
//                for (i = 0; i < pwPeriods.size(); ++i) {
//                    if (i < low || i > high) {
//                        unordered_set<int>& exams = pwPeriods[i];
//                        if (exams.find(*itRem) != exams.end()) {
//                            exams.erase(*itRem);
//                            ++removed;
//                            break;
//                        }

//                    }
//                }
//                // There will be some non removed exams

////                if (i == pwPeriods.size()) {
////                    cout << endl << "numPeriods = " << pwPeriods.size() << endl;
////                    cout << endl << "Couldn't find exam " << *itRem << endl;
//////                    cin.get();
////                }
//            }
//        }
//        cout << "Removed Exams from pwPeriods = " << removed << endl;

//        removed = 0;
//        // Now, remove exams from pwReplacedPeriods
//        for (int j = low; j <= high; ++j) {
//            unordered_set<int> const& examsToRemove = pbPeriods[j];
//            // Remove each exam in examsToRemove vector
//            for (unordered_set<int>::const_iterator itRem = examsToRemove.begin(); itRem != examsToRemove.end(); ++itRem) {
//                int i;
//                for (i = 0; i < pwReplacedPeriods.size(); ++i) {
//                    unordered_set<int>& exams = pwReplacedPeriods[i];
//                    if (exams.find(*itRem) != exams.end()) {
//                        exams.erase(*itRem);
//                        ++removed;
//                        break;
//                    }
//                }
//                // There will be some non removed exams. Those must be reinserted again in Pw

////                if (i == pwPeriods.size()) {
////                    cout << endl << "numPeriods = " << pwPeriods.size() << endl;
////                    cout << endl << "Couldn't find exam " << *itRem << endl;
//////                    cin.get();
////                }
//            }
//        }
//        cout << "Removed Exams from pwReplacedPeriods = " << removed << endl;


//        cout << "Exams remaining to insert" << endl;


//        // Reinsert remaining exams in Pw
//        size = 0;
//        for (int i = 0; i < pwReplacedPeriods.size(); ++i) {
//            unordered_set<int>& exams = pwReplacedPeriods[i];
//            size += exams.size();

////            cout << "exams[" << i << "].size() = " << exams.size() << endl;

////            for (unordered_set<int>::iterator it = exams.begin(); it != exams.end(); ++it) {
////                cout << *it << " ";
////            }
//        }
//        cout << "#Exams from pwReplacedPeriods = " << size << endl;


//        vector<int> fixedExams;
//        unordered_map<int, int> examsPeriods;
//        for (int i = 0; i < pwPeriods.size(); ++i) {
//            unordered_set<int> const& sourceFixedExams = pwPeriods[i];
//            // Add exam and respective period to fixedExams and examsPeriods vectors
//            for (unordered_set<int>::const_iterator it = sourceFixedExams.begin(); it != sourceFixedExams.end(); ++it) {
//                fixedExams.push_back(*it);
////                examsPeriods.push_back(i+1);
//                examsPeriods.insert(pair<int, int>(*it, i+1));
////                pair<unordered_map<int, int>::iterator, bool> p =
////                        exams.insert(pair<int, int>(*examIt, (*courseStudentCounts).at(*examIt)));
//            }
//        }
//        cout << "fixedExams.size() = " << fixedExams.size() << endl;

//        cout << "examsPeriods:" << endl;

//        for (unordered_map<int, int>::const_iterator it = examsPeriods.begin(); it != examsPeriods.end(); ++it) {
//            cout << it->first << " " << it->second << " " << endl;
//        }
//        cout << "examsPeriods.size() = " << examsPeriods.size() << endl;

//        vector<int> remainderExams;
//        for (int i = 0; i < pwReplacedPeriods.size(); ++i) {
//            unordered_set<int> const& sourceRemainderExams = pwReplacedPeriods[i];
//            // Add remainder exam to remainderExams vector
//            for (unordered_set<int>::const_iterator it = sourceRemainderExams.begin(); it != sourceRemainderExams.end(); ++it) {
//                remainderExams.push_back(*it);
//            }
//        }

//        cout << "remainderExams.size() = " << remainderExams.size() << endl;


////        _chrom.init(data);

//        POT copyPw = _Pw;

//        do {
//            _Pw = copyPw;

//            GCHeuristics::SD(init.conflictMatrix, init.graph, _Pw, fixedExams, examsPeriods, remainderExams);

//            _Pw.validate();

//            // Actualize proximity costs
//            _Pw.computeProximityCosts();
//            // Set fitness
//            _Pw.fitness(_Pw.getProximityCost());

//            cout << "[After] Pw Fitness " << _Pw.fitness() << endl;
//            cout << "[After] Pw # Periods " << _Pw.getNumPeriods() << endl;

//            _Pw.pack();

//            // Actualize proximity costs
//            _Pw.computeProximityCosts();
//            // Set fitness
//            _Pw.fitness(_Pw.getProximityCost());

//            cout << "[After Packing] Pw Fitness " << _Pw.fitness() << endl;
//            cout << "[After Packing] Pw # Periods " << _Pw.getNumPeriods() << endl;


//            cout << "Original Pb Fitness " << originalPbFitness << endl;
//            cout << "Original Pw Fitness " << originalPwFitness << endl;

////            if (_Pw.getNumPeriods() != originalPwNumPeriods)
////                cin.get();

//        }
//        while (_Pw.getNumPeriods() != originalPwNumPeriods);
//    }


//    void examsExchangeOperator(POT const& _Pb, POT& _Pw) {

    void examsExchangeOperator(POT const& Pb, POT& _Pw) {

        POT _Pb = _Pw;
        _Pw = Pb;

        double originalPbFitness = _Pb.fitness();
        double originalPwFitness = _Pw.fitness();

        // Get best from frog Pb
        int bestDayIdx = bestDay(_Pb);

        // Get Pb periods
        vector<unordered_set<int> > const& pbPeriods = _Pb.getConstPeriods();
        // Get Pw periods
        vector<unordered_set<int> >& pwPeriods = _Pw.getPeriods();

        // Insert best day in the worst frog and remove duplicated and infeasible exams

        // Generate a random index in Pw
        int randPeriod = rng.random(_Pw.getNumPeriods());
        // Pw exams
        unordered_set<int>& pwExams = pwPeriods[randPeriod];
        unordered_set<int> notInserted;
        int numInsertedExams = 0;
        int numNotInsertedExams = 0;


        // Insert Pb best period into Pw into that random position
        for (auto it = pbPeriods[bestDayIdx].begin(); it != pbPeriods[bestDayIdx].end(); ++it) {
            if (pwExams.find(*it) == pwExams.end() && _Pw.isFeasiblePeriodExam(randPeriod, *it)) {
                pwExams.insert(*it);
                ++numInsertedExams;
            }
            else {
                notInserted.insert(*it);
                ++numNotInsertedExams;
            }

        }
        // Remove duplicated exams

        // To ensure the feasibility of chromosomes after the
        // crossover, duplicated exams are deleted. These exams are
        // removed from the original periods, while the newly inserted
        // periods are left intact.
        unordered_set<int> const& examsToRemove = pbPeriods[bestDayIdx];
        int numRemovedExams = 0;

        for (auto itRem = examsToRemove.begin(); itRem != examsToRemove.end(); ++itRem) {
            int i;
            for (i = 0; i < pwPeriods.size(); ++i) {
                if (i != randPeriod) {
                    unordered_set<int>& exams = pwPeriods[i];
                    if (exams.find(*itRem) != exams.end() && notInserted.find(*itRem) == notInserted.end()) {
                        exams.erase(*itRem);
                        ++numRemovedExams;
                        break;
                    }
                }
            }
        }



        // Sem SHC: 5.62


        // Stochastic Hill Climber algorithm instance
        //
        // Neighborhood
        ETTPneighborhood neighborhood;
        // Full evaluation function
        ProximityCostEval fullEval;

        // Neighbor evaluation function
        ETTPneighborEval neighEval;

//        moSimpleCoolingSchedule<Chromosome> cool(10, 0.9, 100, 0.01);
//        moSimpleCoolingSchedule<Chromosome> cool(0.1, 0.01, 5, 0.001); // 5.62  5 x 20 x 30 (sem SA faz 6.73)

//        moSA<ETTPneighbor> sa(neighborhood, fullEval, neighEval, cool); // 5.54 10 x 50 x 30
                                                                        // 5.73  20 x 25 x 30
                                                                        // 5.80  5 x 30 x 100

        moSimpleCoolingSchedule<Chromosome> cool(0.1, 0.01, 3, 0.001);
        moSA<ETTPneighbor> sa(neighborhood, fullEval, neighEval, cool); // 5.47 15 x 60 x 100
                                                                        // 5.465 15 x 60 x 100 (com random neighbor escolhido no memeplex)

                                                                        // 5.103 15 x 60 x 300

//        sa(_Pw);




//        // tmax
//        int tmax = 10;
////        int tmax = 3;
//        // Temperature T
// //       double tempT = 0.0001; //
////        double tempT = 0.001; // 5.61

////        double tempT = 0.01; // tmax = 3,  5.51
////        double tempT = 0.01; // tmax = 10, m = 10, N = 20, L = 100  5.47
//        double tempT = 0.01; // tmax = 10, m = 15, N = 60, L = 300  5.259


//        //moSHC(Neighborhood& _neighborhood, eoEvalFunc<EOT>& _fullEval, moEval<Neighbor>& _eval,
//        //double _tmax, double _tempT, ProblSense _sense)
//        moSHC<ETTPneighbor> shc(neighborhood, fullEval, neighEval, tmax, tempT, minimize);

////        cout << "Fitness before " << _memeplexes[im][0].fitness() << endl;

//        shc(_Pw);


//        _Pw.validate();

        _Pw.computeProximityCosts();
        _Pw.fitness(_Pw.getProximityCost());

        cout << "[After] Pw Fitness " << _Pw.fitness() << endl;



    }


    // Compute the "Best day" to exchange between chromosomes.
    // The best day consist of three periods (we exclude saturdays because
    // saturdays are always clash-free) and is the day with the lowest number
    // of clashes per student.
    int bestDay(POT const& chromosome) {
        int numPeriods = chromosome.getNumPeriods();
        int indexBestDay = 0, numClashesDay = 0, numStudentsDay = 0;
        double numClashesPerStudent = 0.0, bestNumClashes = 0.0;

        for (int i = 0; i < numPeriods; ++i) {
            // Get number of clashes of the current "day" (or period)
            numClashesDay = chromosome.getNumClashesPeriod(i);
            // Get number of students of the current "day" (or period)
            numStudentsDay = chromosome.getNumStudentsPeriod(i);
            // Compute number of clashes per student
            numClashesPerStudent = (double)numClashesDay / numStudentsDay;

    //        cout << "numClashesDay = " << numClashesDay << endl;
    //        cout << "numStudentsDay = " << numStudentsDay << endl;
    //        cout << "numClashesPerStudent = " << numClashesPerStudent << endl;
    //        cin.get();

            // Actualize best day
            if (i == 0 || numClashesPerStudent < bestNumClashes) {
                bestNumClashes = numClashesPerStudent;
                indexBestDay = i;
            }
        }
        return indexBestDay;
    }

////////////////////////////////////////////////////////////////////////////////////////////
//    bool Crossover::operator()(Chromosome& _chromosome1, Chromosome& _chromosome2)
//    {
//    //    cout << "Crossover" << endl;

//        bool oneAtLeastIsModified = true;
//        //
//        // Day-exchange crossover
//        //
//        // In day-exchange crossover, only the best days (excluding
//        // Saturdays, since exams scheduled on Saturdays are always
//        // clash-free) of chromosomes, selected based on the crossover
//        // rate, are eligible for exchange. The best day consists of three
//        // periods and is the day with the lowest number of clashes per student.

//        // Computation of the offspring
//        generateOffspring(_chromosome1, _chromosome2);
//        generateOffspring(_chromosome2, _chromosome1);


//    //    // does at least one genotype has been modified ?
//    //    if ((_chromosome1 != offspring1) || (_chromosome2 != offspring2))
//    //    {
//    //        // update
//    //        _chromosome1.value(offspring1);
//    //        _chromosome2.value(offspring2);
//    //        // at least one genotype has been modified
//    //        oneAtLeastIsModified = true;
//    //    }
//    //    else
//    //    {
//    //        // no genotype has been modified
//    //        oneAtLeastIsModified = false;
//    //    }
//        // return 'true' if at least one genotype has been modified
//        return oneAtLeastIsModified;
//    }


//    void Crossover::generateOffspring(Chromosome& _parent1, Chromosome& _parent2)
//    {
//        // We change directly _parent1 to be the new offspring because it's already a copy
//        // of the original parent.
//        Chromosome& offspring = _parent1;
//        // Compute the "Best day" of the second parent.
//        int day2 = bestDay(_parent2);
//        // Insert best day in the first parent and remove duplicated exams
//        insertDay(offspring, _parent2, day2);
//        // Actualize proximity costs
//        offspring.computeProximityCosts();

//    //  std::vector<bool> taken_values(result.size(), false);
//    //  if (_point1 > _point2)
//    //    std::swap(_point1, _point2);
//    //  /* first parent */
//    //  for (unsigned int i=0 ; i<=_point1 ; i++)
//    //    {
//    //      // result[i] == _parent1[i]
//    //      taken_values[_parent1[i]] = true;
//    //    }
//    //  for (unsigned int i=_point2 ; i<result.size() ; i++)
//    //    {
//    //      // result[i] == _parent1[i]
//    //      taken_values[_parent1[i]] = true;
//    //    }
//    //  /* second parent */
//    //  unsigned int i = _point1+1;
//    //  unsigned int j = 0;
//    //  while (i<_point2 && j<_parent2.size())
//    //    {
//    //      if (! taken_values[_parent2[j]])
//    //        {
//    //          result[i] = _parent2[j];
//    //          i++;
//    //        }
//    //      j++;
//    //    }
//    //  return result;
//    }

//    // Compute the "Best day" to exchange between chromosomes.
//    // The best day consist of three periods (we exclude saturdays because
//    // saturdays are always clash-free) and is the day with the lowest number
//    // of clashes per student.
//    int Crossover::bestDay(Chromosome const& chromosome) {
//        int numPeriods = chromosome.getNumPeriods();
//        int indexBestDay = 0, numClashesDay = 0, numStudentsDay = 0;
//        double numClashesPerStudent = 0.0, bestNumClashes = 0.0;

//        for (int i = 0; i < numPeriods; ++i) {
//            // Get number of clashes of the current "day" (or period)
//            numClashesDay = chromosome.getNumClashesPeriod(i);
//            // Get number of students of the current "day" (or period)
//            numStudentsDay = chromosome.getNumStudentsPeriod(i);
//            // Compute number of clashes per student
//    //        numClashesPerStudent = (double)numClashesDay / numStudentsDay;
//            numClashesPerStudent = numStudentsDay;

//    //        cout << "numClashesDay = " << numClashesDay << endl;
//    //        cout << "numStudentsDay = " << numStudentsDay << endl;
//    //        cout << "numClashesPerStudent = " << numClashesPerStudent << endl;
//    //        cin.get();

//            // Actualize best day
//            if (i == 0 || numClashesPerStudent < bestNumClashes) {
//                bestNumClashes = numClashesPerStudent;
//                indexBestDay = i;
//            }
//        }
//        return indexBestDay;
//    }

//    class IsInExamsToRemove {
//        vector<int> const& examsToRemove;
//    public:
//        IsInExamsToRemove(vector<int> const& examsToRemove) : examsToRemove(examsToRemove) { }
//        bool operator()(int exam) {
//            return find_if(examsToRemove.begin(), examsToRemove.end(), bind2nd(equal_to<int>(), exam)) != examsToRemove.end();
//        }
//    };


//    void Crossover::insertDay(Chromosome& _parent1, Chromosome& _parent2, int day2) {

//    //    cout << endl << endl << "Before inserting day" << endl;
//        _parent1.validate();

//        // The newly inserted day takes the place of a randomly chosen day
//        // which is pushed to the end of the timetable.
//    //    int numPeriods = _parent1.getNumPeriods();
//        // Get periods
//    //    vector<vector<int> >& periods1 = _parent1.getPeriods();
//    //    vector<vector<int> >& periods2 = _parent2.getPeriods();
//        vector<unordered_set<int> >& periods1 = _parent1.getPeriods();
//        vector<unordered_set<int> >& periods2 = _parent2.getPeriods();

//        // Generate random period
//    //    int randPeriod = rng.random(numPeriods);
//        int randPeriod = rng.random(periods1.size());


//    //    cout << "randPeriod = " << randPeriod << endl;
//    //    cout << "periods1 exams before inserting day from periods2" << endl;
//    //    int i = 1;
//    //    for (vector<vector<int> >::iterator it = periods1.begin(); it != periods1.end(); ++it) {
//    //        cout << endl << "Period " << i << endl;
//    //        for (vector<int>::const_iterator examit = (*it).begin(); examit != (*it).end(); ++examit) {
//    //            cout << *examit << " ";
//    //        }
//    //        cout << endl;
//    //        ++i;
//    //    }

//        // Push existing day to the end of the timetable
//        periods1.push_back(periods1[randPeriod]);
//        // Insert new day from second parent
//        periods1[randPeriod] = periods2[day2];
//        // Increment number of periods
//    //    ++numPeriods;
//    //    _parent1.setNumPeriods(numPeriods);

//    //    cout << "periods1 exams after inserting day from periods2" << endl;
//    //    i = 1;
//    //    for (vector<vector<int> >::iterator it = periods1.begin(); it != periods1.end(); ++it) {
//    //        cout << endl << "Period " << i << endl;
//    //        for (vector<int>::const_iterator examit = (*it).begin(); examit != (*it).end(); ++examit) {
//    //            cout << *examit << " ";
//    //        }
//    //        cout << endl;
//    //        ++i;
//    //    }



//        // To ensure the feasibility of chromosomes after the
//        // crossover, duplicated exams are deleted. These exams are
//        // removed from the original periods, while the newly inserted
//        // periods are left intact.
//    //    vector<int> const& examsToRemove = periods2[day2];
//        unordered_set<int> const& examsToRemove = periods2[day2];

//    //    cout << "periods2 exams to remove:" << endl;
//    //    for (vector<int>::const_iterator examit = examsToRemove.begin(); examit != examsToRemove.end(); ++examit) {
//    //        cout << *examit << " ";
//    //    }

//    //    IsInExamsToRemove pred(examsToRemove);

//        int numRemovedExams = 0;

//        /// TODO: USAR unordored_map<INT> PARA GUARDAR EXAMES

//        // Remove each exam in examsToRemove vector
//    //    for (vector<int>::const_iterator itRem = examsToRemove.begin(); itRem != examsToRemove.end(); ++itRem) {
//        for (unordered_set<int>::const_iterator itRem = examsToRemove.begin(); itRem != examsToRemove.end(); ++itRem) {

//            int i;
//    //        for (i = 0; i < numPeriods; ++i) {
//            for (i = 0; i < periods1.size(); ++i) {
//                if (i != randPeriod) {
//    //                vector<int>& exams = periods1[i];
//                    unordered_set<int>& exams = periods1[i];
//        //            vector<int>::iterator newEnd = remove_if(exams.begin(), exams.end(), pred);
//        //            numRemovedExams += (exams.end()-newEnd);
//        //            exams.erase(newEnd, exams.end());


//                    if (exams.find(*itRem) != exams.end()) {
//                        exams.erase(*itRem);
//                        break;
//                    }



//    //                vector<int>::iterator newEnd = remove_if(exams.begin(), exams.end(), bind2nd(equal_to<int>(), *itRem));

//                    // ESTA LINHA NAO FUNCIONA
//    //                unordered_set<int>::iterator newEnd = remove_if(exams.begin(), exams.end(),
//    //                                                                     bind2nd(equal_to<int>(), *itRem));

//    //                if (newEnd != exams.end()) {
//    //                    // Remove exam
//    //                    exams.erase(newEnd, exams.end());
//    //                    ++numRemovedExams;
//    //                    break;
//    //                }
//    //                int numExamsToRemove = (exams.end()-newEnd);
//    //                cout << endl << "numExamsToRemove = " << numExamsToRemove << endl;
//    //                cout << endl << "exams before remove" << endl;
//    //                copy(exams.begin(), exams.end(), ostream_iterator<int>(cout, " "));
//    //                cout << endl;
//    //                numRemovedExams += numExamsToRemove;
//    //                exams.erase(newEnd, exams.end());
//    //                cout << endl << "exams after remove" << endl;
//    //                copy(exams.begin(), exams.end(), ostream_iterator<int>(cout, " "));
//    //                cout << endl;
//                }
//            }
//            if (i == periods1.size()) {
//                cout << endl << "numPeriods = " << periods1.size() << endl;
//                cout << endl << "Couldn't find exam " << *itRem << endl;
//                cin.get();
//            }
//        }

//    //    cout << endl << endl << "numRemovedExams = " << numRemovedExams << endl;
//    //    cout << endl << endl << "examsToRemove.size() = " << examsToRemove.size() << endl;

//    //    cout << endl << endl << "After inserting day and remove duplicated exams" << endl;

//    //    cout << "periods1 exams" << endl;
//    //    i = 1;
//    //    for (vector<vector<int> >::iterator it = periods1.begin(); it != periods1.end(); ++it) {
//    //        cout << endl << "Period " << i << endl;
//    //        for (vector<int>::const_iterator examit = (*it).begin(); examit != (*it).end(); ++examit) {
//    //            cout << *examit << " ";
//    //        }
//    //        cout << endl;
//    //        ++i;
//    //    }

//        _parent1.validate();
//    }

//    /////////////////////////////////////////////////////////////////////////////////////////////////




};

#endif // EOSFLA_H












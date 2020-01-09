#ifndef EOSFLA_H
#define EOSFLA_H



#include "Chromosome.h"
#include "MOEA.h"
#include "ETTPneighborhood.h"
#include "ETTPneighborEval.h"
#include "algo/moSA.h"
#include "eoEvolutionOperator.h"



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

//#include <eoAlgo.h>
//#include "eoSFLAInitializer.h"
//#include <eoContinue.h>
//#include <eoEvalFunc.h>
//#include <vector>
//#include <boost/shared_ptr.hpp>
//#include "Chromosome.h"


using namespace std;
using namespace boost;



/**
    This is a generic class for Shuffled Frog-Leaping algorithms. There
    is only one operator defined, which takes a population and does stuff to
    it. It needn't be a complete algorithm, can be also a step of an
    algorithm. Only used for single-objective cases.

    @ingroup Algorithms
*/


template <class POT>
class eoSFLA : public eoAlgo<POT>
{
public:

    /** Constructor
     * @param n - number of memeplexes
     * @param N - number of frogs in each memeplex
     * @param F - total sample size in the swamp and given by F = mN
     * @param _init - An eoInit that initializes each frog (solution)
     * @param _continuator - An eoContinue that manages the stopping criterion and the checkpointing system
     * @param _eval - An eoEvalFunc: the evaluation performer
     */
    eoSFLA<POT>(int _m, int _N, int _q, int _F, long _L, eoInit<POT>& _init, eoContinue<POT>& _continuator,
                    eoEvalFunc<POT>& _eval, eoBinOp<POT>& _chromEvolOperator) :
            m(_m), N(_N), q(_q), F(_F), L(_L), init(_init), continuator(_continuator),
            eval(_eval), evolution(L*m+1), counter(0), chromEvolOperator(_chromEvolOperator)
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

            // Update evolution registry
            evolution[counter++] = Pg.fitness();



            ////////////////////////////////
            // Added: 7.Feb.2014
            //
            t = 0;
//            initT = 1, r = 0.000000000001; // Yor 37.76 and stagnate

//            initT = 1, r = 0.00001; // 38.81 with swap

            initT = 0.1, r = 0.00000001; // 38.1


//            initT = 10, r = 0.0000001;// Yor 38.56

//            initT = 1, r = 0.0000001; // Yor 37.37 stagnates
//            initT = 0.1, r = 0.000001;// Yor 37.76 and stagnate

            _temp = initT;
            ////////////////////////////////



            do
            {

                // Step 3 Partition frogs into memeplex
                partitionFrogsIntoMemplexes(_pop, memplexes);

                printMemplexes(memplexes);


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



                ////////////////////////////////
                // Added: 7.Feb.2014
                //
                t = t+1; // Increment rate variable
                _temp = Temp(t, initT, r);
                ////////////////////////////////

                cout << "Temp = " << _temp << endl;

            }
            while (continuator(_pop));
        }
        catch (std::exception& _e)
        {
            std::string s = _e.what();
            s.append (" in eoSFLA operator()");
            throw std::runtime_error(s);
        }

        cout << "Best frog overall fitness = " << Pg.fitness() << endl;

        ofstream outFile("Evolution.txt");
        outFile << "[";
        for (auto it = evolution.begin(); it != evolution.end(); ++it)
            outFile << *it << ", ";
        outFile << "]" << endl;
        outFile.close();
    }


    POT getBestFrog() {
        return Pg;
    }

protected:
    // m is the number of memeplexes
    int m;
    // N is the number of frogs in each memeplex
    int N;
    // q is the number of frogs in each submemeplex (q < N)
    int q;
    // The total sample size, F, in the swamp is given by F = mN
    int F;
    // Number of time loops
    long L;
    // Create a vector of dimension L*m
    vector<double> evolution;
    // Iteration counter
    int counter;
    // SFLA Initializer
    //eoSFLAInitializer<POT>& init;
    eoInit<POT>& init;
    // Continuator
    eoContinue<POT>& continuator;
    // Evaluation function
    eoEvalFunc<POT>& eval;
    // Chromosome evolution operator
    eoBinOp<POT>& chromEvolOperator;


    // The frog with the global best fitness is identified as Pg.
    POT Pg;


    /////////////////////////////////////////////////////////////////////////////
    // Auxiliary methods
    /////////////////////////////////////////////////////////////////////////////

    void printMemplexes(vector<vector<POT> >& _memeplexes) {

        // Iterate over the memeplexes
        for (int im = 0; im < m; ++im) {
            cout << "Memplex #" << im << ":" << endl;
            // Iterate memeplex
            for (int iN = 0; iN < N; ++iN) {
                cout << '\t' << _memeplexes[im][iN].fitness() << endl;
            }
        }
        cout << endl;
    }


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
        bool operator()(POT const& _chrom1, POT const& _chrom2) {
            return _chrom1.fitness() < _chrom2.fitness(); // Ascending order because we're minimizing
        }
    };

    void rankFrogs(eoPop<POT>& _pop) {
//        cout << "_pop.size() = " << _pop.size() << endl;
        // Step 2 Rank frogs. Sort the F frogs in order of decreasing performance value. Store them in an array X={U(i),
        // f(i), i=1,...,F} so that i=1 represents the frog with the best performance value. Record the best frog’s
        // position, PX, in the entire population (F frogs; where PX =U(1)).
        sort(_pop.begin(), _pop.end(), FitnessCmp());
        // Record the best frog’s position, PX, in the entire population (F frogs; where PX =U(1)).
        //   The best frog is just the 0 index frog
        updateGlobalBestFrog(_pop[0]);

//        cout << "After ranking" << endl;
//        cout << "Best frog fitness: " << _pop[0].getChromosome() << endl;
    }

    /// TODO: Optimizar copia do vector no retorno

    vector<vector<POT> > createMemeplexes() {
//        cout << "createMemeplexes" << endl;
//        cout << "m = " << m << endl;

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

        // As the memeplexes contain references to population we don't have to combine memeplexes into population vector.
        // We just need to sort the population vector.
        sort(_pop.begin(), _pop.end(), FitnessCmp());
        // Record the best frog’s position, PX, in the entire population (F frogs; where PX =U(1)).
        //   The best frog is just the 0 index frog
        updateGlobalBestFrog(_pop[0]);

        ////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////


        cout << "Best frog fitness = " << getBestFrog().fitness() << endl;
    }


    //
    // Added: 7.Feb.2014
    //
    // Temperature actualization
    double Temp(double t, double Tmax, double R) {
        double newtemp = Tmax*exp(-R*t);
        return newtemp;
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
//            vector<bool> marked(N);
//            auto resultPair = createSubMemeplex(_memeplexes[im], q, marked);
//            auto subMemeplex = resultPair.first;
//            auto randomIdxs = resultPair.second;

//            sort(subMemeplex.begin(), subMemeplex.end(), FitnessCmp());

            // CHANGED 31-Out-2013
auto& subMemeplex  = _memeplexes[im];

//            cout << "Sorted submemeplex fitnesses" << endl;

//            for (int i = 0; i < subMemeplex.size(); ++i) {
//                cout << subMemeplex[i].fitness() << endl;
//            }

            // Perform N evolutionary steps
//            for (int iN = 0; iN < N; ++iN) {
            for (int iN = 0; iN < q; ++iN) { // hum
//            for (int iN = 0; iN < 10; ++iN) { // hum

                // Identify as Pb and Pw, respectively, the frogs with the best and the worst fitness.
                // Also, the frog with the global best fitness is identified as Pg.
                // Then, an evolution process is applied to improve only the frog with
                // the worst fitness (i.e., not all frogs) in each cycle.
                ////

                // Optimize submemeplex
                POT Pb = getMemeplexBestFrog(subMemeplex);
                POT Pw = getMemeplexWorstFrog(subMemeplex);
                // Step 4-3 Improve the worst frog’s position.

// CHANGED 31-Out-2013

                POT Prandom = subMemeplex[rng.random(subMemeplex.size()/2)];
                POT newFrog = improveWorstFrogPosition(Prandom, Pw);
//                POT newFrog = improveWorstFrogPosition(Pb, Pw);

                // Step 4-4 If this process produces a better frog (solution), it replaces the worst frog.
                // Otherwise, the calculations in Eqs. 5 and 6 are repeated with respect to the global best frog (i.e., Pg replaces Pb).
                // Step 4-5 If no improvement becomes possible in this latter case, then a new solution is
                // randomly generated to replace the worst frog with another frog having any arbitrary fitness (as shown in Fig. 3b).

                double fit1;
                double fit2;
                fit1 = (double)newFrog.fitness();
                fit2 = (double)Pw.fitness();

                double alpha = exp((fit2 - fit1) / (fit2*_temp) );
                double r = rng.uniform();
                bool isAccept = (r < alpha);


                if (isAccept)
                    replaceWorstFrog(newFrog, subMemeplex);

//                    replaceWorstFrog(newFrog, subMemeplex);
                /*if (newFrog.fitness() < Pw.fitness()) {
                    replaceWorstFrog(newFrog, subMemeplex);
                }
                else {
                    // Repeat computations with respect to the global best frog (i.e., Pg replaces Pb).
                    POT newFrog = improveWorstFrogPosition(Pg, Pw);
                    if (newFrog.fitness() < Pw.fitness()) {
                        replaceWorstFrog(newFrog, subMemeplex);
                    }
                    else {
                        // Randomly generate a solution to replace the worst frog with another frog having any arbitrary fitness
                        POT newFrogRandom = generateRandomFrogPosition();
                        replaceWorstFrog(newFrogRandom, subMemeplex);
                    }
                }
                */
                // Update Global frog
                POT bestMemeplexFrog = getMemeplexBestFrog(subMemeplex);
                if (bestMemeplexFrog.fitness() < Pg.fitness())
                    updateGlobalBestFrog(bestMemeplexFrog);

//                listMemeplex(subMemeplex);
            }


            if (rng.uniform() < 0) {
//            if (rng.uniform() < 0.1) {
//              if (rng.uniform() < 0.3) {
//            if (rng.uniform() < 1) {
/*
                int randIdx = rng.random(subMemeplex.size());
                POT& frog1 = subMemeplex[randIdx];
//                int randIdx = rng.random(_memeplexes[im].size());
//                POT& frog1 = _memeplexes[im][randIdx];
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

*/

////                double T1 = 10; // 39.18, pr = 0.01, with period swap, 3 periods
////                double T1 = 2; // Very good on car92
////                double T1 = 0.1;
////                double T1 = 0.01;
//                double T1 = 0.001;

//                // Stochastic Hill Climber algorithm instance
//                //
//                // Neighborhood
//                ETTPneighborhood neighborhood;
//                // Full evaluation function
//                ProximityCostEval<POT> fullEval;
////                ProximityCostEval<eoChromosome> fullEval;
//                // Neighbor evaluation function
//                ETTPneighborEval neighEval;
////                moSimpleCoolingSchedule<eoChromosome> cool(T1, 0.1, 10, T1); // SHC
//                moSimpleCoolingSchedule<eoChromosome> cool(T1, 1, 3, T1); // SHC


////                moSimpleCoolingSchedule<eoChromosome> cool(1, 0.0001, 0, 0.00001); // 4.97

////                moSimpleCoolingSchedule<eoChromosome> cool(1, 0.0001, 0, 0.01); // Good
////                moSimpleCoolingSchedule<eoChromosome> cool(1, 0.0001, 0, 0.1); // Good


//                moSA<ETTPneighbor> sa(neighborhood, fullEval, neighEval, cool);
//                sa(frog1);

            }

//            cout << "[After improvement] Sorted submemeplex fitnesses" << endl;

//            for (int i = 0; i < subMemeplex.size(); ++i) {
//                cout << subMemeplex[i].fitness() << endl;
//            }

//            cin.get();

//            cout << "After iterating submemeplex" << endl;
//            copy(marked.begin(), marked.end(), ostream_iterator<bool>(cout, " "));
//            cout << endl;
//            copy(randomIdxs.begin(), randomIdxs.end(), ostream_iterator<int>(cout, " "));
//            cout << endl;


            /*
             // CHANGED 31-Out-2013


            // Upgrade the memeplex
            int n = N, k = 0;
//            for (int i = 1; i <= n; ++i) {
            for (int i = 1; i <= q; ++i) {
//                if (marked[i-1]) {
                    int frogIdx = randomIdxs[i-1];
//                    cout << "frogIdx = " << frogIdx << endl;
                    _memeplexes[im][frogIdx] = subMemeplex[k++];
//                    _memeplexes[im][frogIdx-1] = subMemeplex[k++];
//                }
            }
*/



            sort(_memeplexes[im].begin(), _memeplexes[im].end(), FitnessCmp()); // don't know if necessary....

//            cout << "k = " << k << endl;
//            cin.get();

            // Update evolution registry
            evolution[counter++] = Pg.fitness();

        }
    }

    double triangular(double a,double b,double c) {
//       double U = rand() / (double) RAND_MAX;
        double U = rng.uniform();
        double F = (c - a) / (b - a);

        // cout << "U = " << U << "F = " << F << endl;

        if (U <= F)
          return a + sqrt(U * (b - a) * (c - a));
        else
          return b - sqrt((1 - U) * (b - a) * (b - c));
    }


    pair<vector<POT>, vector<int> > createSubMemeplex(vector<POT>& _memeplex, int _q, vector<bool>& marked) {
        int n = _memeplex.size();
        // Create a vector of n random indexes
        vector<int> randomIdxs;
//        for (int i = 1; i <= n; ++i)
//            randomIdxs[i-1] = i;
//        random_shuffle(randomIdxs.begin(), randomIdxs.end());

////        copy(randomIdxs.begin(), randomIdxs.end(), ostream_iterator<int>(cout, " "));


//        // Create a boolean vector in order to mark chosen positions
////        vector<bool> marked(n);
//        // Create submemeplex
        vector<POT> submemeplex(_q);
//        int k = 0;
//        double sum = 0;
//        // Select q random frogs to form the submemeplex
//        for (int i = 1; i <= n; ++i) {
////            int idx = randomIdxs[i-1];
////            // Index of frog in original memeplex
////            int j = idx;
//////            cout << "j = " << j << ", n = " << n << endl;


////            // Individual probability of being chosen
////            double pj = (double)(2*(n+1-j)) / (n*(n+1));
////            sum += pj;
////            cout << "pj = " << pj << endl;
////            double r = rng.uniform();
////            cout << "r = " << r << endl;


//            if (r < pj) {
//                // Select individual
//                submemeplex[k++] = _memeplex[j-1];
//                marked[i-1] = true;
//            }
//            if (k == _q)
//                break; // We have selected q frogs
//        }

////        cout << "sum = " << sum << endl;
////        cout << "triangular = " << triangular(1, n, (n-10)) << endl;



////cin.get();

//        // If there weren't chosen q frogs select the remainder manually
////        if (k < _q) {
////            for (int i = 1; i <= n; ++i) {
////                if (marked[i-1]) {
////                    int idx = randomIdxs[i-1];
////                    int j = idx;
////                    // Select individual
////                    submemeplex[k++] = _memeplex[j-1];
////                    marked[i-1] = true;
////                }
////                if (k == _q)
////                    break; // We have selected q frogs
////            }
////        }

////        cout << "k = " << k << ", q = " << _q << endl;

////        cout << "k == q: " << (k == _q) << endl;


        int numSelected = 0;
        int k = 0;
        while (numSelected < _q) {
            int idx = (int)triangular(1, n, 1);

            // cout << "idx = " << idx << endl;

            if (!marked[idx]) {
                marked[idx] = true;
                ++numSelected;
                randomIdxs.push_back(idx);
//                submemeplex[k++] = _memeplex[idx-1]; // humm
                submemeplex[k++] = _memeplex[idx];
            }
        }
//        cout << "done " << endl;
//        copy(marked.begin(), marked.end(), ostream_iterator<bool>(cout, " "));
//        cout << endl;

        sort(submemeplex.begin(), submemeplex.end(), FitnessCmp());

        return pair<vector<POT>, vector<int> >(submemeplex, randomIdxs);
    }


    POT generateRandomFrogPosition() {
//        boost::shared_ptr<POT> ptrNewFrogRandom(new POT());
        POT newFrogRandom = POT();

        init(newFrogRandom);
        eval(newFrogRandom);

        ////
        // Neighborhood
//        ETTPneighborhood neighborhood;
//        // Full evaluation function
//        ProximityCostEval<POT> fullEval;
//        // Neighbor evaluation function
//        ETTPneighborEval neighEval;
//        moSimpleCoolingSchedule<POT> cool(0.01, 0.001, 10, 0.0000001);
//        moSA<ETTPneighbor> sa(neighborhood, fullEval, neighEval, cool);

//        cout << "Random Frog improved with SA" << endl;
//        sa(newFrogRandom);
//        cout << newFrogRandom.fitness() << endl;
        ////


        /*
        ////
        POT newFrogRandom = getBestFrog();

        // Neighborhood
        ETTPneighborhood neighborhood;
        // Full evaluation function
        ProximityCostEval<POT> fullEval;
//        ProximityCostEval<eoChromosome> fullEval;
        // Neighbor evaluation function
        ETTPneighborEval neighEval;
//        moSimpleCoolingSchedule<eoChromosome> cool1(0.1, 0.000001, 3, 0.0999);
        moSimpleCoolingSchedule<POT> cool1(0.1, 0.000001, 3, 0.0999);
        moSA<ETTPneighbor> sa1(neighborhood, fullEval, neighEval, cool1);

        cout << "Random Frog with SA" << endl;
        sa1(newFrogRandom);
//        eval(newFrogRandom);
        cout << newFrogRandom.fitness() << endl;
        ////

//        eval(newFrogRandom);

*/

        return newFrogRandom;
    }


    void replaceWorstFrog(POT& newFrog, vector<POT>& _memeplex) {
//        cout << "replaceWorstFrog" << endl;
//        cout << "Before replace" << endl;
//        listMemeplex(_memeplex);

        // Replace worst frog by new one and then orderly relocate new frog in the memeplex vector
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
    }


    void listMemeplex(vector<POT> const& _memeplex) {

        cout << "List Memeplex" << endl;

        // Perform N iterations
        for (int iN = 0; iN < _memeplex.size(); ++iN) {
            cout << "\t" << _memeplex[iN].fitness() << endl;
        }
    }


    POT improveWorstFrogPosition(POT const& _Pb, POT const& _Pw) {
        POT improvedFrog = _Pw;
        // Added: 7.Feb.2014
//        chromEvolOperator.setTemperature(_temp);
        eoSFLAEvolOperator_2* evol = dynamic_cast<eoSFLAEvolOperator_2*>(&chromEvolOperator);
        if (evol != nullptr)
            evol->setTemperature(_temp);
        else
            std::cerr << "This object is not of type eoSFLAEvolOperator_2\n";

        chromEvolOperator(improvedFrog, _Pb);
        // Evaluate new frog
        eval(improvedFrog);
        return improvedFrog;
    }



/*
    // Compute the "Best day" to exchange between chromosomes.
    // The best day consist of three periods (we exclude saturdays because
    // saturdays are always clash-free) and is the day with the lowest number
    // of clashes per student.
    int bestDay(POT const& chromosome) {
        int numPeriods = chromosome.getNumPeriods();
        int indexBestDay = 0, numClashesDay = 0, numStudentsDay = 0;
        double numClashesPerStudent = 0.0, bestNumClashes = 0.0;


        /// ELE ESCOLHE SEMPRE O ULTIMO PERIODO

        for (int i = 0; i < numPeriods-1; ++i) { /// ALTERADO
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
    */


private:
    double _temp;
    ////////////////////////////////
    // Added: 7.Feb.2014
    //
    int t;
    double initT, r;
    ////////////////////////////////


};



#endif // EOSFLA_H












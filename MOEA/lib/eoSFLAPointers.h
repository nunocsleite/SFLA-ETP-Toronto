#ifndef EOSFLA_H
#define EOSFLA_H

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
    eoSFLA(int _m, int _N, int _F, eoInit<POT>& _init, eoContinue<POT>& _continuator, eoEvalFunc<POT>& _eval) :
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
            vector<vector<POT*> > memplexes = createMemeplexes();

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
    eoInit<POT>& init;
    // Continuator
    eoContinue<POT>& continuator;
    // Evaluation function
    eoEvalFunc<POT>& eval;
    // The frog with the global best fitness is identified as Pg.
    POT* Pg;


    /////////////////////////////////////////////////////////////////////////////
    // Auxiliary methods
    /////////////////////////////////////////////////////////////////////////////

    void updateGlobalBestFrog(POT* _chrom) {
        // Pre-condition: vector _pop is sorted by fitness
        // Global best frog is located at index 0 in the population vector
        Pg = _chrom;
    }
    POT* getMemeplexBestFrog(vector<POT*>& _memeplex) {
        // Pre-condition: vector _memeplex is sorted by fitness
        // Memeplex best frog is located at index 0 in the memeplex vector
        return _memeplex[0];
    }
    POT* getMemeplexWorstFrog(vector<POT*>& _memeplex) {
        // Pre-condition: vector _memeplex is sorted by fitness
        // Memeplex worst frog is located at the last index in the memeplex vector
        return _memeplex[_memeplex.size()-1];
    }

    struct FitnessCmp {
        bool operator()(Chromosome const& _chrom1, Chromosome const& _chrom2) {

            /// MELHORAR

//            bool b = true;

//            if (_chrom1.getNumPeriods() > _chrom1.getRange()[1]) {
//                cout << "Infeasible chromosome 1: Number of periods = " << _chrom1.getNumPeriods()
//                     << " above upper limit (" << _chrom1.getRange()[1] << ")" << endl;
//                // throw excp.
//                b = false;
//            }
//            else if (_chrom1.getNumPeriods() < _chrom1.getRange()[0]) {
//                cout << "Infeasible chromosome 1: Number of periods = " << _chrom1.getNumPeriods()
//                     << " below lower limit (" << _chrom1.getRange()[0] << ")" << endl;
//                // throw excp.
//                b = false;
//            }

//            if (_chrom2.getNumPeriods() > _chrom2.getRange()[1]) {
//                cout << "Infeasible chromosome 2: Number of periods = " << _chrom2.getNumPeriods()
//                     << " above upper limit (" << _chrom2.getRange()[1] << ")" << endl;
//                // throw excp.
//                b = false;
//            }
//            else if (_chrom2.getNumPeriods() < _chrom2.getRange()[0]) {
//                cout << "Infeasible chromosome 2: Number of periods = " << _chrom2.getNumPeriods()
//                     << " below lower limit (" << _chrom2.getRange()[0] << ")" << endl;
//                // throw excp.
//                b = false;

//            }

//            if (b)
                return _chrom1.fitness() < _chrom2.fitness(); // Ascending order because we're minimizing
//            else
//                return false;
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
        updateGlobalBestFrog(&_pop[0]);

        cout << "After ranking" << endl;

        cout << "Best frog fitness: " << _pop[0] << endl;

//        cin.get();
    }

    /// TODO: Optimizar copia do vector no retorno

    vector<vector<POT*> > createMemeplexes() {

        cout << "createMemeplexes" << endl;
        cout << "m = " << m << endl;

        vector<vector<POT*> > memeplexes;
        for (int i = 0; i < m; ++i)
            memeplexes.push_back(vector<POT*>(N));

//        cout << memplexes.size() << endl;
//        for (int i = 0; i < m; ++i)
//            cout << "memplex size = " << memplexes[i].size() << endl;

//        cin.get();

        return memeplexes;
    }

    void partitionFrogsIntoMemplexes(eoPop<POT>& _pop, vector<vector<POT*> >& _memeplexes) {
        // Step 3 Partition frogs into memeplex. Partition array X into m memeplexes Y1, Y2,...,Ym, each containing
        // n frogs. E.g., for m=3, rank 1 goes to memeplex 1, rank 2 goes to memeplex 2, rank 3 goes to memeplex 3,
        // rank 4 goes to memeplex 1, and so on (Fig. 2).

        cout << "partitionFrogsIntoMemeplexes" << endl;

        vector<int> memeplexIndexes(m);
        for (int i = 0; i < _pop.size(); ++i) {
            // Each entry in each memeplex references a solution contained in population vector
            _memeplexes[i%m][memeplexIndexes[i%m]] = &_pop[i];
            // Increment memeplex
            ++memeplexIndexes[i%m];
        }

        cout << "List Memeplexes" << endl;

        // Iterate over the memeplexes
        for (int im = 0; im < m; ++im) {
            cout << "Memeplexe " << im << endl;

            // Perform N evolutionary steps
            for (int iN = 0; iN < N; ++iN) {
                cout << "\t" << (*_memeplexes[im][iN]).fitness() << endl;
            }
        }

//        cin.get();

    }

    void shuffleMemplexes(vector<vector<POT*> >& _memeplexes, eoPop<POT>& _pop) {
        // Step 5 Shuffle memeplexes. After a defined number of memetic evolutionary steps within each memeplex,
        // replace Y1,...,Ym into X such that X={Yk, k = 1..., m}. Sort X in order of decreasing performance value.
        // Update the population the best frog’s position PX.

        cout << "shuffleMemeplexes" << endl;


        // As the memeplexes contain references to population we don't have to combine mememplxes into population vector.
        // We just need to sort the population vector.
        sort(_pop.begin(), _pop.end(), FitnessCmp());
        // Record the best frog’s position, PX, in the entire population (F frogs; where PX =U(1)).
        //   The best frog is just the 0 index frog
        updateGlobalBestFrog(&_pop[0]);
    }

    void frogLeapingLocalSearch(vector<vector<POT*> >& _memeplexes) {
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

        cout << "frogLeapingLocalSearch" << endl;

        // Iterate over the memeplexes
        for (int im = 0; im < m; ++im) {
            // Perform N evolutionary steps
            for (int iN = 0; iN < N; ++iN) {
                // Identify as Pb and Pw, respectively, the frogs with the best and the worst fitness.
                // Also, the frog with the global best fitness is identified as Pg.
                // Then, an evolution process is applied to improve only the frog with
                // the worst fitness (i.e., not all frogs) in each cycle.
                ////
                POT* Pb = getMemeplexBestFrog(_memeplexes[im]);
                POT* Pw = getMemeplexWorstFrog(_memeplexes[im]);
                // Step 4-3 Improve the worst frog’s position.
                POT newFrog = improveWorstFrogPosition(Pb, Pw);
                // Step 4-4 If this process produces a better frog (solution), it replaces the worst frog.
                // Otherwise, the calculations in Eqs. 5 and 6 are repeated with respect to the global best frog (i.e., Pg replaces Pb).
                // Step 4-5 If no improvement becomes possible in this latter case, then a new solution is
                // randomly generated to replace the worst frog with another frog having any arbitrary fitness (as shown in Fig. 3b).
                if (newFrog.fitness() < Pw->fitness()) {
                    replaceWorstFrog(newFrog, _memeplexes[im]);
                }
                else {
                    // Repeat computations with respect to the global best frog (i.e., Pg replaces Pb).
                    ptrNewFrog = improveWorstFrogPosition(Pg, Pw);
                    if (ptrNewFrog->fitness() < Pw->fitness()) {
                        replaceWorstFrog(ptrNewFrog, _memeplexes[im]);
                    }
                    else {
                        // Randomly generate a solution to replace the worst frog with another frog having any arbitrary fitness
                        boost::shared_ptr<POT> ptrNewFrogRandom = generateRandomFrogPosition();
                        replaceWorstFrog(ptrNewFrogRandom, _memeplexes[im]);
                    }
                }
//                // Update Global frog
//                POT* bestMemeplexFrog = getMemeplexBestFrog(_memeplexes[im]);
//                if (bestMemeplexFrog->fitness() < Pg->fitness())
//                    updateGlobalBestFrog(bestMemeplexFrog);

            }
        }
    }

    POT generateRandomFrogPosition() {
//        boost::shared_ptr<POT> ptrNewFrogRandom(new POT());
        POT newFrogRandom;

        cout << "improveWorstFrogPosition" << endl;

        init(newFrogRandom);

        cout << "[improveWorstFrogPosition] after init random solution" << endl;

//        cin.get();

        return newFrogRandom;
    }

    void replaceWorstFrog(POT& newFrog, vector<POT*>& _memeplex) {

        cout << "replaceWorstFrog" << endl;

        cout << "Before replace" << endl;

        listMemeplex(_memeplex);

        // Replace worst frog by new one and then orderly relocte new frog in the memeplex vector
        _memeplex[_memeplex.size()-1] = &newFrog;
        for (int i = _memeplex.size()-2; i >= 0; --i) {
            if (_memeplex[i]->fitness() > _memeplex[i+1]->fitness()) {
                // Swap entries
                POT* aux = _memeplex[i];
                _memeplex[i] = _memeplex[i+1];
                _memeplex[i+1] = aux;
            }
            else // Stop. The memeplex is ordered.
                break;
        }

        cout << "After replace" << endl;

        listMemeplex(_memeplex);

//        cin.get();
    }

    void listMemeplex(vector<POT*> const& _memeplex) {

        cout << "List Memeplex" << endl;

        // Perform N iterations
        for (int iN = 0; iN < N; ++iN) {
            cout << "\t" << (*_memeplex[iN]).fitness() << endl;
        }
    }


    boost::shared_ptr<POT> improveWorstFrogPosition(POT* _Pb, POT* _Pw) {

        cout << "improveWorstFrogPosition" << endl;

        boost::shared_ptr<POT> ptrNewFrog(new POT());

        // Just for testing

        init(*ptrNewFrog);



//        cin.get();
        return ptrNewFrog;
    }



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












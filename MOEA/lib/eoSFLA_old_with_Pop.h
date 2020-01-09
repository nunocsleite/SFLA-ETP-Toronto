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
 *  Step 3 Partition frogs into memeplex. Partition array X into m memeplexes Y1, Y2,...,Ym, each containing
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


//#include <eoContinue.h>
//#include <eoPSO.h>
//#include <eoVelocity.h>
//#include <eoFlight.h>

///** An easy-to-use particle swarm algorithm.
//* Use any particle, any flight, any topology...
//*
//*   The main steps are :
//*         (The population  is expected to be already evaluated)
//*        - for each generation and each particle pi
//*        - evaluate the velocities
//*               -- perform the fligth of pi
//*       -- evaluate pi
//*       -- update the neighborhoods
//*
//* @ingroup Algorithms
//*/
//template < class POT > class eoEasyPSO:public eoPSO < POT >
//{
//public:

//    /** Full constructor
//    * @param _init - An eoInitializerBase that initializes the topology, velocity, best particle(s)
//    * @param _continuator - An eoContinue that manages the stopping criterion and the checkpointing system
//    * @param _eval - An eoEvalFunc: the evaluation performer
//    * @param _velocity - An eoVelocity that defines how to compute the velocities
//    * @param _flight - An eoFlight that defines how to make the particle flying: that means how
//    * to modify the positions according to the velocities
//    */
//    eoEasyPSO (
//        eoInitializerBase <POT> &_init,
//        eoContinue < POT > &_continuator,
//        eoEvalFunc < POT > &_eval,
//        eoVelocity < POT > &_velocity,
//        eoFlight < POT > &_flight):
//            init(_init),
//            continuator (_continuator),
//            eval (_eval),
//            velocity (_velocity),
//            flight (_flight)
//    {}


//    /** Constructor without eoFlight. For special cases when the flight is performed withing the velocity.
//       * @param _init - An eoInitializerBase that initializes the topology, velocity, best particle(s)
//       * @param _continuator - An eoContinue that manages the stopping criterion and the checkpointing system
//       * @param _eval - An eoEvalFunc: the evaluation performer
//       * @param _velocity - An eoVelocity that defines how to compute the velocities
//    */
//    eoEasyPSO (
//        eoInitializerBase <POT> &_init,
//        eoContinue < POT > &_continuator,
//        eoEvalFunc < POT > &_eval,
//        eoVelocity < POT > &_velocity):
//            init(_init),
//            continuator (_continuator),
//            eval (_eval),
//            velocity (_velocity),
//            flight (dummyFlight)
//    {}


//        /* Constructor without eoInitializerBase. Assume the initialization is done before running the algorithm
//    * @param _continuator - An eoContinue that manages the stopping criterion and the checkpointing system
//    * @param _eval - An eoEvalFunc: the evaluation performer
//    * @param _velocity - An eoVelocity that defines how to compute the velocities
//    * @param _flight - An eoFlight that defines how to make the particle flying: that means how
//    * to modify the positions according to the velocities
//    */
//    eoEasyPSO (
//        eoContinue < POT > &_continuator,
//        eoEvalFunc < POT > &_eval,
//        eoVelocity < POT > &_velocity,
//        eoFlight < POT > &_flight):
//            init(dummyInit),
//            continuator (_continuator),
//            eval (_eval),
//            velocity (_velocity),
//            flight (_flight)
//    {}


//    /** Constructor without eoFlight nor eoInitializerBase. For special cases when the flight is performed withing the velocity.
//       * @param _continuator - An eoContinue that manages the stopping criterion and the checkpointing system
//       * @param _eval - An eoEvalFunc: the evaluation performer
//       * @param _velocity - An eoVelocity that defines how to compute the velocities
//    */
//    eoEasyPSO (
//        eoContinue < POT > &_continuator,
//        eoEvalFunc < POT > &_eval,
//        eoVelocity < POT > &_velocity):
//            init(dummyInit),
//            continuator (_continuator),
//            eval (_eval),
//            velocity (_velocity),
//            flight (dummyFlight)
//    {}

//    /// Apply a few iteration of flight to the population (=swarm).
//    virtual void operator  () (eoPop < POT > &_pop)
//    {
//        try
//        {
//            // initializes the topology, velocity, best particle(s)
//            init();
//            do
//            {
//                // loop over all the particles for the current iteration
//                for (unsigned idx = 0; idx < _pop.size (); idx++)
//                {
//                    // perform velocity evaluation
//                    velocity (_pop[idx],idx);

//                    // apply the flight
//                    flight (_pop[idx]);

//                    // evaluate the position
//                    eval (_pop[idx]);

//                    // update the topology (particle and local/global best(s))
//                    velocity.updateNeighborhood(_pop[idx],idx);
//                }

//            }
//            while (continuator (_pop));

//        }
//        catch (std::exception & e)
//        {
//            std::string s = e.what ();
//            s.append (" in eoEasyPSO");
//            throw std::runtime_error (s);
//        }

//    }

//protected:
//    eoInitializerBase <POT> &init;
//    eoContinue < POT > &continuator;
//    eoEvalFunc < POT > &eval;
//    eoVelocity < POT > &velocity;
//    eoFlight < POT > &flight;

//         // if the flight does not need to be used, use the dummy flight instance
//        class eoDummyFlight:public eoFlight < POT >
//        {
//         public:
//        eoDummyFlight () {}
//        void operator  () (POT & _po) {}
//        }dummyFlight;

//        // if the initializer does not need to be used, use the dummy one instead
//        class eoDummyInitializer:public eoInitializerBase < POT >
//        {
//         public:
//        eoDummyInitializer () {}
//        void operator  () (POT & _po) {}
//        }dummyInit;

//};
///**
// * @example t-eoEasyPSO.cpp
// * Example of a test program building a PSO algorithm.
// */





/**
    This is a generic class for Shuffled Frog-Leaping algorithms. There
    is only one operator defined, which takes a population and does stuff to
    it. It needn't be a complete algorithm, can be also a step of an
    algorithm. Only used for mono-objective cases.

    @ingroup Algorithms
*/
template <class POT>
class eoSFLA : public eoAlgo<POT>
{
public:

    /** Constructor
     * @param _init - An eoInitializerBase that initializes the SFLA data structures
     * @param _continuator - An eoContinue that manages the stopping criterion and the checkpointing system
     * @param _eval - An eoEvalFunc: the evaluation performer
     */
    eoSFLA (
        eoSFLAInitializer<POT>& _init,
        eoContinue<POT>& _continuator,
        eoEvalFunc<POT>& _eval) :
            init(_init),
            continuator(_continuator),
            eval(_eval)
    {}

    // Apply SFLA to the population.
    virtual void operator() (eoPop<POT> &_pop)
    {
        try
        {
            // Step 0, 1, and 2 - Initialization, generation of virtual population of frogs, and ranking
            init();
            //
            // Step 3 Partition frogs into memeplex
            shuffleMemplexes(pop);

            do
            {
                //
                // Step 4 Frog-leaping algorithm for local search
                frogLeapingLocalSearch(pop);
                //
                // Step 5 Shuffle memeplexes. After a defined number of memetic evolutionary steps within each memeplex,
                // replace Y1,...,Ym into X such that X={Yk, k = 1..., m}. Sort X in order of decreasing performance value.
                // Update the population the best frog’s position PX.
                shuffleMemplexes(pop);
                //
                // Step 6 Check convergences. If the convergence criteria are satisfied, stop. Otherwise, return to step 3.
                // Typically, the decision on when to stop is made by a prespecified number of consecutive time loops when
                // at least one frog carries the “best memetic pattern” without change. Alternatively, a maximum total number of
                // function evaluations can be defined (Fig. 4).
            }
            while (continuator(_pop));
        }
        catch (std::exception& e)
        {
            std::string s = e.what();
            s.append (" in eoSFLA");
            throw std::runtime_error(s);
        }
    }


protected:
    // m is the number of memeplexes
    int m;
    // N is the number of frogs in each memeplex.
    int N;
    // The total sample size, F, in the swamp is given by F = mN.
    int F;
    // SFLA Initializer
    eoSFLAInitializer<POT>& init;
    // Continuator
    eoContinue<POT>& continuator;
    // Evaluation function
    eoEvalFunc<POT>& eval;

    /////////////////////////////////////////////////////////////////////////////
    // Auxiliary methods
    /////////////////////////////////////////////////////////////////////////////

    void shuffleMemplexes() {
        // Step 3 Partition frogs into memeplex. Partition array X into m memeplexes Y1, Y2,...,Ym, each containing
        // n frogs. E.g., for m=3, rank 1 goes to memeplex 1, rank 2 goes to memeplex 2, rank 3 goes to memeplex 3,
        // rank 4 goes to memeplex 1, and so on (Fig. 2).


    }

    void frogLeapingLocalSearch() {
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

        // Iterate over the memeplexes
        for (int im = 0; im < m; ++im) {
            // Perform N evolutionary steps
            for (int iN = 0; iN < N; ++iN) {
                // Identify as Pb and Pw, respectively, the frogs with the best and the worst fitness.
                // Also, the frog with the global best fitness is identified as Pg.
                // Then, an evolution process is applied to improve only the frog with
                // the worst fitness (i.e., not all frogs) in each cycle.
                ////

                // Step 4-3 Improve the worst frog’s position.

                // Step 4-4 If this process produces a better frog (solution), it replaces the worst frog.
                // Otherwise, the calculations in Eqs. 5 and 6 are repeated with respect to the global best frog (i.e., Pg replaces Pb).

                // Step 4-5 If no improvement becomes possible in this latter case, then a new solution is
                // randomly generated to replace the worst frog with another frog having any arbitrary fitness (as shown in Fig. 3b).


            }
        }
    }

};

#endif // EOSFLA_H












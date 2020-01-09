#ifndef MOEOSFLA_H
#define MOEOSFLA_H

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
template <class MOEOT>
class moeoSFLA: public moeoEA<MOEOT>
{
public:



    /**
     * Apply a the algorithm to the population _pop until the stopping criteria is satified.
     * @param _pop the population
     */
    virtual void operator () (eoPop<MOEOT> &_pop)
    {
//        eoPop < MOEOT > offspring, empty_pop;
//        popEval (empty_pop, _pop);	// a first eval of _pop
//        // evaluate fitness and diversity
//        fitnessAssignment(_pop);
//        diversityAssignment(_pop);
//        do
//        {
//            // generate offspring, worths are recalculated if necessary
//            breed (_pop, offspring);
//            // eval of offspring
//            popEval (_pop, offspring);
//            // after replace, the new pop is in _pop. Worths are recalculated if necessary
//            replace (_pop, offspring);
//        }
//        while (continuator (_pop));


        // Global exploration
        //
        // Step 0 Initialize. Select m and n, where m is the number of memeplexes and n is the number of frogs in each
        // memeplex. Therefore, the total sample size, F, in the swamp is given by F = mn.
        //
        // Step 1 Generate a virtual population. Sample F virtual frogs, U(1),...,U(F).
        // Compute the performance value f(i) for each frog U.
        eoPop<MOEOT> empty_pop;
        popEval(empty_pop, _pop);	// a first eval of _pop
        // evaluate fitness and diversity
        fitnessAssignment(_pop);
        //diversityAssignment(_pop);

        // Step 2 Rank frogs. Sort the F frogs in order of decreasing performance value. Store them in an array X={U(i),
        // f(i), i=1,...,F} so that i=1 represents the frog with the best performance value. Record the best frog’s
        // position, PX, in the entire population (F frogs; where PX =U(1)).
        //
        // Step 3 Partition frogs into memeplex. Partition array X into m memeplexes Y1, Y2,...,Ym, each containing
        // n frogs. E.g., for m=3, rank 1 goes to memeplex 1, rank 2 goes to memeplex 2, rank 3 goes to memeplex 3,
        // rank 4 goes to memeplex 1, and so on (Fig. 2).
        //
        // Step 4 Memetic evolutions within each memeplex. Evolve each memeplex Y(l), l=1,..., m.
        //
        // After partitioning frogs to m memeplexes, evolve each memeplex and each of them should iterate N times. After
        // the memeplexes have been evolved, the algorithm returns to the global exploration for shuffling.
        //
        // Local exploration: frog-leaping algorithm
        // Step 4-0 Set im=0, where im counts the number of memeplexes and will be compared with
        // the total number m of memeplexes. Set iN=0, where iN counts the number of evolutionary
        // steps and will be compared with the maximum number N of steps to be completed within each
        // memeplex. Within each memeplex (Fig. 3b), the frogs with the best and the worst fitness are
        // identified as Pb and Pw, respectively. Also, the frog with the global best fitness is identified as
        // Pg. Then, an evolution process is applied to improve only the frog with the worst fitness (i.e.,
        // not all frogs) in each cycle.
        // Step 4-1 Set im=im+1.
        // Step 4-2 Set iN=iN+1.
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
        //
        // Step 5 Shuffle memeplexes. After a defined number of memetic evolutionary steps within each memeplex,
        // replace Y1,...,Ym into X such that X={Yk, k = 1..., m}. Sort X in order of decreasing performance value.
        // Update the population the best frog’s position PX.
        //
        // Step 6 Check convergences. If the convergence criteria are satisfied, stop. Otherwise, return to step 3.
        // Typically, the decision on when to stop is made by a prespecified number of consecutive time loops when
        // at least one frog carries the “best memetic pattern” without change. Alternatively, a maximum total number of
        // function evaluations can be defined (Fig. 4).





    }


protected:
    // m is the number of memeplexes
    int m;
    // n is the number of frogs in each memeplex.
    int n;
    // The total sample size, F, in the swamp is given by F = mn.
    int F;


};



#endif // MOEOSFLA_H










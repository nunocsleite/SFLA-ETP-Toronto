#ifndef MOEOLOCALSEARCHUPDATER_H
#define MOEOLOCALSEARCHUPDATER_H

#include <eoPop.h>
#include <utils/eoUpdater.h>
#include "moSHC.h"
#include "Chromosome.h"
#include "moeoETTPEval.h"
#include "ETTPInit.h"
#include "eoSFLA.h"


/**
* This class allows to improve some solutions contained in the
* main population at each generation.
*/
//template <class MOEOT>
template <class MOEOT, class Neighbor>
class moeoLocalSearchUpdater : public eoUpdater
{
public:

    /**
     * Ctor
     * @param _pop the main population
     */
//    moeoLocalSearchUpdater(eoPop<MOEOT>& _pop, moLocalSearch<Neighbor> & _localSearchAlg)
//        : pop(_pop), localSearchAlg(_localSearchAlg), iterationIdx(1) { }

    moeoLocalSearchUpdater(ETTPInit<eoChromosome>& init, eoPop<MOEOT>& _pop, moLocalSearch<Neighbor> & _localSearchAlg)
        : eoinit(init), pop(_pop), localSearchAlg(_localSearchAlg), iterationIdx(1) { }


    /**
     * Improve some solutions contained in the main population
     */
    void operator()()
    {
//        cout << "Local search updater"  << endl;
//        cout << pop.size() << endl;
//        pop.printOn(cout);
//        cin.get();


//        function pop = improvePopulation(S, Data, pop, HCAlgorithm)
//            % The local search operator is applied to chromosomes
//            % selected based on a tournament selection scheme,
//            % where all the chromosomes in the population are randomly
//            % grouped into fours and from each group, the chromosome
//            % with the smallest rank is selected. Only a quarter of the
//            % population will undergo local exploitation.
//            popSize = length(pop);
//            p = randperm(popSize);

//            for i = 1 : 4 : popSize
//                ranks = [S.ri(i), S.ri(i+1), S.ri(i+2), S.ri(i+3)];
//                [M, I] = min(ranks);
//                improvIdx = i+I-1;
//                % Improve with Local search
//                pop{improvIdx} = HCAlgorithm(pop{improvIdx}, Data);
//            end
//        end

        // Extract first front of the final population using an moeoArchive (this is the output of nsgaII)
        moeoUnboundedArchive<moeoChromosome> arch;
        arch(pop);
        // Printing of the final archive
        cout << endl << "Archive before improving:" << endl;
        arch.sortedPrintOn(cout);
        cout << endl;

        // The local search operator is applied to chromosomes
        // selected based on a tournament selection scheme,
        // where all the chromosomes in the population are randomly
        // grouped into fours and from each group, the chromosome
        // with the smallest rank is selected. Only a quarter of the
        // population will undergo local exploitation.
//        int popSize = pop.size();
        int popSize = arch.size();

//        for i = 1 : 4 : popSize
//            ranks = [S.ri(i), S.ri(i+1), S.ri(i+2), S.ri(i+3)];
//            [M, I] = min(ranks);
//            improvIdx = i+I-1;
//            % Improve with Local search
//            pop{improvIdx} = HCAlgorithm(pop{improvIdx}, Data);
//        end

//        cout << popSize << endl;
        cout << "Iteration: " << iterationIdx++ << endl;

        moeoETTPEval eval;

//        for (int i = 0; i < popSize; ++i) {

//            cout << "Improve pop[" << i << "]" << endl;

////            cout << "pop[" << i << "].getNumPeriods() = " << pop[i].getNumPeriods() << endl;

//            /// TODO: SEE HOW TO MANIPULATE THE Rt POPULATION

//////            pop[i].invalidateObjectiveVector();

            //eval(pop[i]);

//           int i = rng.random(pop.size());

        for (int i = 0; i < arch.size(); ++i) {


//            eoChromosome sol(pop[i].getChromosome());
            eoChromosome sol(arch[i].getChromosome());


//            cout << "Before Local Search" << endl;
//            cout << "pop [" << i << "] num periods = " << pop[i].getNumPeriods() << endl;
//            cout << "pop [" << i << "] fitness = " << pop[i].fitness() << endl;

//            cout << "sol num periods = " << sol.getNumPeriods() << endl;
//            cout << "sol fitness = " << sol.fitness() << endl;

//            localSearchAlg(sol);

            //////////////////////////////////////////////////////////////////////

            // Simulated Annealing
            // Neighborhood
            ETTPneighborhood neighborhood;
            // Full evaluation function
            ProximityCostEval<eoChromosome> fullEval;
            // Neighbor evaluation function
            ETTPneighborEval neighEval;
            // Cooling schedule
//            moSimpleCoolingSchedule<eoChromosome> cool(5, 0.9, 10, 0.00001); //
//            moSimpleCoolingSchedule<eoChromosome> cool(5, 0.9, 10, 0); // 10.65 - Hec, 37.877 Yor
            moSimpleCoolingSchedule<eoChromosome> cool(0.01, 0.001, 10, 0.0000001); //  Yor 36.71 - 21 periods, 39.48 - 20 periods

//            moSimpleCoolingSchedule<eoChromosome> cool(0.002, 0.0001, 10, 0.001);
            //  Yor83
            // cost = 44.8002 - periods = 19
            // cost = 40.678 - periods = 20
            // cost = 37.0213 - periods = 21

            //            moSimpleCoolingSchedule<eoChromosome> cool(0.1, 0.001, 10, 0.001); //  Yor

//            moSimpleCoolingSchedule<eoChromosome> cool(0.1, 0.001, 3, 0.0001); //  Yor

//            moSimpleCoolingSchedule<eoChromosome> cool(0.0002, 0.00001, 3, 0.0001); //  Yor lento...
//            moSimpleCoolingSchedule<eoChromosome> cool(0.002, 0.0001, 3, 0.001); //  Yor


//                        moSimpleCoolingSchedule<eoChromosome> cool(0.01, 0.0001, 10, 0.0000001); // 35.67 Yor
//                        moSimpleCoolingSchedule<eoChromosome> cool(0.000001, 0.0001, 10, 0.000001); //  Yor

            // SA object
            moSA<ETTPneighbor> sa(neighborhood, fullEval, neighEval, cool);

        //    cout << "SA" << endl;
        //    sa(frog1);
            sa(sol);

            //////////////////////////////////////////////////////////////////////


//            cout << "After Local Search" << endl;
//            cout << "pop [" << i << "] num periods = " << pop[i].getNumPeriods() << endl;
//            cout << "pop [" << i << "] fitness = " << pop[i].fitness() << endl;

//            cout << "sol num periods = " << sol.getNumPeriods() << endl;
//            cout << "sol fitness = " << sol.fitness() << endl;

//            cin.get();

//            pop[i] = MOEOT(sol.getChromosome());
            arch[i] = MOEOT(sol.getChromosome());


            // Evaluate solution in order to validate the objective vector
//            eval(pop[i]);

            eval(arch[i]);

//            int j = rng.random(pop.size());
//            pop[j] = arch[i];

            pop[i] = arch[i];




//            cout << endl << endl << pop << endl;

//            pop[i].validate();



        }




        /*
        ////////////////// SFLA ///////////////////////////////////////////////////
        unordered_map<int, eoPop<eoChromosome> > solutions;
        // Organize solutions by number of periods
        for (int i = 0; i < popSize; ++i) {
            // Create instance of eoChromosome (by copy)
            eoChromosome sol(pop[i].getChromosome());
            sol.computeProximityCosts();
            // (by copy)
            solutions[pop[i].getNumPeriods()].push_back(pop[i]);
        }

//        // Print map
//        cout << "Print map" << endl;

//        for (auto it = solutions.begin(); it != solutions.end(); ++it) {
//            cout << "# periods: " << it->first << endl;
//            cout << "solutions: " << endl;
//            for (auto itSol = it->second.begin(); itSol != it->second.end(); ++itSol) {
//                cout << itSol->fitness() << endl;
//            }
//        }

        cout << "SFLA" << endl;

        for (auto it = solutions.begin(); it != solutions.end(); ++it) {
//            // Dont't optimize infeasible solutions
//            if (it->first < 28 || it->first > 36)
//                continue;
            // Set range of chromosome initializer to fixed number of periods
            ETTPInit<eoChromosome> chromInit = eoinit; // (by copy)
            chromInit.setFixedRange(it->first);
            runSFLA(chromInit, it->second);
//            eoPop<eoChromosome> copyPop = it->second;
//            runSFLA(chromInit, copyPop);
//            //
//            it->second[0] = copyPop[0];
        }

        // Copy solutions to pop
        int i = 0;
        moeoETTPEval eval;
        for (auto it = solutions.begin(); it != solutions.end(); ++it) {
            for (auto itSol = it->second.begin(); itSol != it->second.end(); ++itSol) {
                pop[i] = MOEOT(itSol->getChromosome());
                eval(pop[i]);
                ++i;
            }
        }


        ///////////////////////////////////////////////////////////////////////////

//        cin.get();

*/
        /** fitness assignment used in NSGA */
        moeoDominanceDepthFitnessAssignment<MOEOT> fitnessAssignment;
        /** diversity assignment used in NSGA-II */
        moeoFrontByFrontCrowdingDiversityAssignment<MOEOT> diversityAssignment;

        // Evaluate fitness and diversity
        fitnessAssignment(pop);
        diversityAssignment(pop);
    }


private:
    /** Solution initializer */
    ETTPInit<eoChromosome>& eoinit;
    /** the main population */
    eoPop<MOEOT> & pop;
    /** the local search algorithm */
    moLocalSearch<Neighbor> & localSearchAlg;

    int iterationIdx;

    void runSFLA(ETTPInit<eoChromosome>& init, eoPop<eoChromosome>& _pop) {
/* TO DO
        ////////////////// SFLA ///////////////////////////////////////////////////

        if (_pop.size() == 1)
            return;


        const int F = _pop.size(); // The total sample size, F, in the swamp is given by F = mN.
        const int m = 1; // 3;  // m is the number of memeplexes
        const int N = F/m;  // N is the number of frogs in each memeplex.
        // Objective function evaluation
        eoETTPEval eval;
        // Number of consecutive time loops
//        int L = 100;
        int L = 10;


        eoGenContinue<eoChromosome> terminator(L); // Terminate after concluding L time loops
        eoCheckPoint<eoChromosome> checkpoint(terminator);
        // Build SFLA
        eoSFLA<eoChromosome> sfla(m, N, F, L, init, checkpoint, eval);
        // Run the algorithm
        sfla(_pop);

        // Print best frog
        //...

        cout << "Fim" << endl;

    //    outFile << sfla.getBestFrog().fitness() << " - " << sfla.getBestFrog().getNumPeriods() << endl;

        cout << "SFLA Best Frog: " << sfla.getBestFrog().getChromosome() << endl;

//        _pop[0] = sfla.getBestFrog();

        //////////////////////////////////////////////////////////////////////////

        */
    }


};


#endif // MOEOLOCALSEARCHUPDATER_H

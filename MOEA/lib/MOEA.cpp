
#include "MOEA.h"
#include "Repairing.h"
#include "moeoLocalSearchUpdater.h"
#include "ETTPneighborhood.h"
#include "ETTPneighborEval.h"


NSGAII::NSGAII(unsigned int popSize, unsigned int maxGen, double crossProb, double mutProb, double reInsertionRate,
               ETTPInit<moeoChromosome>& ettpInit, ETTPInit<eoChromosome>& eoinit)
    : popSize(popSize), maxGen(maxGen), crossProb(crossProb), mutProb(mutProb),
      reInsertionRate(reInsertionRate), ettpInit(ettpInit), eoinit(eoinit)
{

}


void NSGAII::run(ofstream& outFile) {

    cout << "NSGAII::Run" << endl;

    // Generate initial population
    eoPop<moeoChromosome> pop(popSize, ettpInit);

    cout << pop << endl;

//    outFile << "Initial Pop" << endl;
//    pop.sortedPrintOn(outFile);
//    outFile << endl;

//    cin.get();

    // Build NSGA-II
//    moeoNSGAII<Chromosome> nsgaII(maxGen, eval, crossover, crossProb, mutation, mutProb);

    // Timetable packing and Continuation
    eoRepair<moeoChromosome> repair(pop);
    eoGenContinue<moeoChromosome> terminator(maxGen);
    eoCheckPoint<moeoChromosome> checkpoint(terminator);
    checkpoint.add(repair);

    // Use a hybridization at the checkpointing step.
    // We have to implement an "eoUpdater" which applies a local search.
    // This updater should be added in a "eoCheckpoint".
//    eoCheckPoint<FlowShop> checkpoint (term);
//    moeoArchiveUpdater<FlowShop> updater(arch, pop);
//    checkpoint.add(updater);
    /// TODO: Apply LS updater to both pop and offspring by implementing a custom popeval

    // Stochastic Hill Climber algorithm instance
    //
    // Neighborhood
//    ETTPneighborhood neighborhood;
    // Full evaluation function
//    NumClashesEval fullEval;
//    ProximityCostEval<eoChromosome> fullEval;

    // Neighbor evaluation function
//    ETTPneighborEval neighEval;
    // tmax
//    int tmax = 5;
//    int tmax = 3;
    // Temperature T
//    double tempT = 0.0001;
    //moSHC(Neighborhood& _neighborhood, eoEvalFunc<EOT>& _fullEval, moEval<Neighbor>& _eval,
    //double _tmax, double _tempT, ProblSense _sense)
//    moSHC<ETTPneighbor> shc(neighborhood, fullEval, neighEval, tmax, tempT, minimize);

    // Simulated Annealing
    // Neighborhood
    ETTPneighborhood neighborhood;
    // Full evaluation function
    ProximityCostEval<eoChromosome> fullEval;
    // Neighbor evaluation function
    ETTPneighborEval neighEval;
    // Cooling schedule
    moSimpleCoolingSchedule<eoChromosome> cool(5, 0.9, 10, 0.00001); //
//    moSimpleCoolingSchedule<eoChromosome> cool(5, 0.9, 10, 0); // 10.65 - Hec
//    moSimpleCoolingSchedule<eoChromosome> cool(0.01, 0.0001, 10, 0.0000001); // 35.67 Yor
    // SA object
    moSA<ETTPneighbor> sa(neighborhood, fullEval, neighEval, cool);


//    moeoLocalSearchUpdater<moeoChromosome, ETTPneighbor> updater(eoinit, pop, shc);
    moeoLocalSearchUpdater<moeoChromosome, ETTPneighbor> updater(eoinit, pop, sa);

    checkpoint.add(updater);

    // Breeding
    eoSGATransform<moeoChromosome> transform(crossover, crossProb, mutation, mutProb);
    // Build NSGA-II
    moeoNSGAII<moeoChromosome> nsgaII(checkpoint, eval, transform);

    // indicator
//    moeoAdditiveEpsilonBinaryMetric < ETTPObjectiveVector > indicator;
//    moeoIBEA<moeoChromosome> nsgaII(checkpoint, eval, transform, indicator);


    // Run the algorithm
    nsgaII(pop);

    // Extract first front of the final population using an moeoArchive (this is the output of nsgaII)
    moeoUnboundedArchive<moeoChromosome> arch;
    arch(pop);

    // Printing of the final archive
    cout << endl << endl << "Final Archive" << endl;
    arch.sortedPrintOn(cout);
    cout << endl;

//    cout << pop << endl;

//    outFile << "Final Archive" << endl;
//    arch.sortedPrintOn(outFile);
//    outFile << endl;


    // IBEA results

    //cool(0.002, 0.0001, 3, 0.001);
//Iteration: 25
//Final Archive

//Timetable

// cost = 44.9617 - periods = 19


//Timetable

// cost = 40.6206 - periods = 20


//Timetable

// cost = 37.1382 - periods = 21




}



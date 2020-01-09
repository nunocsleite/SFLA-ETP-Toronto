#ifndef _MOEA_H_
#define _MOEA_H_


#include <iostream>
#include <moeo>
#include <es/eoRealInitBounded.h>
#include <es/eoRealOp.h>
#include "Crossover.h"
#include "Mutation.h"
#include "Chromosome.h"
#include "ETTPInit.h"
#include "moeoETTPEval.h"

using namespace std;



// Single-objective (minimize number of clashes) evaluation
// Used in Local Search exploration
//class ProximityCostEval : public eoEvalFunc<Chromosome> {
template <typename EOT>
class ProximityCostEval : public eoEvalFunc<EOT> {
public:
//    void operator () (EOT& chrom) {

    void operator () (eoChromosome& chrom) {

        chrom.computeProximityCosts(); // shall include this, hum
        // Evaluate proximity cost
        chrom.fitness(chrom.getProximityCost());
     }
};




class NSGAII {
    unsigned int popSize;
    unsigned int maxGen;
    double crossProb;
    double mutProb;
    double reInsertionRate;
    // Chromosome initializer
//    ETTPInit& ettpInit;
    ETTPInit<moeoChromosome>& ettpInit;
    ETTPInit<eoChromosome>& eoinit;

    // Objective functions evaluation
    moeoETTPEval eval;
    // Crossover and mutation
    Mutation mutation;
    Crossover crossover;

public:
    NSGAII(unsigned int popSize, unsigned int maxGen, double crossProb, double mutProb,
//           double reInsertionRate, ETTPInit& init);
           double reInsertionRate, ETTPInit<moeoChromosome>& init, ETTPInit<eoChromosome>& eoinit);

    void run(ofstream &outFile);
};


#endif

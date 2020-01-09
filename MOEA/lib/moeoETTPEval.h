#ifndef MOEOETTPEVAL_H
#define MOEOETTPEVAL_H

#include "Chromosome.h"
#include <moeo>



// Evaluation of objective functions
class moeoETTPEval : public moeoEvalFunc<moeoChromosome> {
public:
    void operator () (moeoChromosome& chrom) {
        ETTPObjectiveVector objVector;
        objVector[0] = proximityCost(chrom);
        objVector[1] = numPeriods(chrom);
        chrom.objectiveVector(objVector);
//        chrom.
    }
    inline double proximityCost(moeoChromosome& chrom) {
        chrom.computeProximityCosts();
        chrom.fitness(chrom.getProximityCost()); ///////////////// VER


//        if (chrom.getNumPeriods() > chrom.getRange()[1])
//            return chrom.getNumClashes() + (chrom.getNumPeriods() - chrom.getRange()[1])*1000; // Penalize solution
//        else
//            return chrom.getNumClashes();

        /// TODO: HOW TO PENALIZE SOLUTION?

        if (chrom.getNumPeriods() > chrom.getRange()[1])
            return chrom.getProximityCost() + (chrom.getNumPeriods() - chrom.getRange()[1])*1000; // Penalize solution
        else
            return chrom.getProximityCost();
    }
    inline double numPeriods(moeoChromosome& chrom) {
        if (chrom.getNumPeriods() > chrom.getRange()[1])
            return chrom.getNumPeriods() + (chrom.getNumPeriods() - chrom.getRange()[1])*10; // Penalize solution
        else
            return chrom.getNumPeriods();
    }
};


#endif // MOEOETTPEVAL_H

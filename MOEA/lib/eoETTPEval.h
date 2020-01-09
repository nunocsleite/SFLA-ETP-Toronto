#ifndef EOETTPEVAL_H
#define EOETTPEVAL_H

#include "Chromosome.h"


// Evaluation of objective function
class eoETTPEval : public eoEvalFunc<eoChromosome> {
public:
    void operator () (eoChromosome& _chrom) {

//        _chrom.fitness(proximityCost(_chrom));

        proximityCost(_chrom);
    }

    inline double proximityCost(eoChromosome& _chrom) {
        _chrom.computeProximityCosts();

//        if (chrom.getNumPeriods() > chrom.getRange()[1])
//            return chrom.getNumClashes() + (chrom.getNumPeriods() - chrom.getRange()[1])*1000; // Penalize solution
//        else
//            return chrom.getNumClashes();

        /// TODO: HOW TO PENALIZE SOLUTION?

//        if (_chrom.getNumPeriods() > _chrom.getRange()[1])
//            return _chrom.getProximityCost() + (_chrom.getNumPeriods() - _chrom.getRange()[1])*0.2; // Penalize solution
//        else
            return _chrom.getProximityCost();
    }
};


#endif // EOETTPEVAL_H

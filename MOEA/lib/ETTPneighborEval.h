#ifndef ETTPNEIGHBOREVAL_H
#define ETTPNEIGHBOREVAL_H

#include <eval/moEval.h>
#include "ETTPneighbor.h"

class ETTPneighborEval : public moEval<ETTPneighbor> {

    virtual void operator()(eoChromosome & timetable, ETTPneighbor & neighbor) {
//        cout << endl << "ETTPneighborEval" << endl;

//        cout << "timetable.fitness() = " << timetable.fitness()
//             << ", neighbor.fitness() = " << neighbor.fitness() << endl;

        // Nothing to do because neighbor has already a fitness value

    }
};

#endif // ETTPNEIGHBOREVAL_H

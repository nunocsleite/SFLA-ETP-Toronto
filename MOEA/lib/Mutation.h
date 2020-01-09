#ifndef MUTATION_H
#define MUTATION_H

#include "Chromosome.h"


using namespace std;


class Mutation : public eoMonOp<moeoChromosome>
{

public:

    /**
    * the class name (used to display statistics)
    */
    string className() const;

    /**
    * eoMon mutation - _chromosome is the chromosome to mutate.
    * @param _chromosome the chromosome
    */
    bool operator()(moeoChromosome& _chromosome);

private:

    /**
    * generation of an offspring
    * @param _parent1 the first parent
    * @param _parent2 the second parent
    */
//    void generateOffspring(Chromosome &_parent1, Chromosome &_parent2);

};




#endif // MUTATION_H

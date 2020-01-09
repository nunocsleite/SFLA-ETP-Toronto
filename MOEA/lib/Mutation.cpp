#include "Mutation.h"



string Mutation::className() const
{
    return "ETTP Mutation";
}


bool Mutation::operator()(moeoChromosome& chromosome)
{
//    cout << "Mutation" << endl;

    bool chromosomeIsModified = true;




    // return 'true' if at least one genotype has been modified
    return chromosomeIsModified;
}



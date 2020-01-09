
#include <cmath>
#include <iostream>

using namespace std;


extern void testMOEA();

//extern void testSFLA();

// Temperature actualization
double Temp(double t, double Tmax, double R) {
    double newtemp = Tmax*exp(-R*t);
    return newtemp;
}


int getSANumberEvaluations(double tmax, double r, double k, double tmin) {
    double t = 0, temp = tmax;
    long numberEvaluations = 0;
    do {
        for (int i = 1; i <= k; ++i)
            ++numberEvaluations;
        // Actualize temperature
        ++t;
        temp = Temp(t, tmax, r);
    } while (temp >= tmin);

    return numberEvaluations;
}



int main(int argc, char* argv[])
{
/*
//    cout << getSANumberEvaluations(0.0001, 0.001, 5, 0.0001) << endl;
//    cout << getSANumberEvaluations(1, 0.001, 1, 0.0001) << endl;
//    cout << getSANumberEvaluations(1, 0.001, 5, 0.0001) << endl;
//    cout << getSANumberEvaluations(3, 0.001, 1, 0.0001) << endl;
//    cout << getSANumberEvaluations(3, 0.001, 5, 0.0001) << endl;

    // moSimpleCoolingSchedule<eoChromosome> cool2(0.1, 0.00001, 5, 0.0000001);
    double tmax = 0.1, r = 0.00001, k = 5, tmin = 0.0000001;

    cout << getSANumberEvaluations(tmax, r, k, tmin) << endl;
// number of evaluations: 6 907 760
*/

    testMOEA();

    // Shuffled Leap Frog Algorithm
////    testSFLA();

	return 0;
}












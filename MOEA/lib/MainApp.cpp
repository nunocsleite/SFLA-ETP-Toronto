#include <iostream>
#include <vector>
#include <fstream>
#include "TestSetDescription.h"
#include "Data.h"
#include "ETTPInit.h"
#include "MOEA.h"
#include "eoSFLA.h"
#include "eoETTPEval.h"
#include "Repairing.h"
#include "eoEvolutionOperator.h"
#include "TorontoTestSet.h"
#include "ITC07TestSet.h"

#include <string>
#include <stdio.h>
#include <time.h>

//#include "boost/graph/sequential_vertex_coloring.hpp"

using namespace std;


// Get current date/time, format is YYYY-MM-DD.HH:mm:ss
const std::string currentDateTime() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
    // for more information about date/time format
    strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);

    return buf;
}



// These functions are defined below
void runSFLA(ofstream& outFile, TestSet const& testSet);
void runNSGAII(ofstream& outFile, TestSet const& testSet);
void runTorontoDatasets();
void runITC2007Datasets();



void testMOEA() {
    runTorontoDatasets();
//    runITC2007Datasets();

//    typedef adjacency_list<listS, vecS, undirectedS> Graph;
//    typedef graph_traits<Graph>::vertex_descriptor vertex_descriptor;
//    typedef graph_traits<Graph>::vertices_size_type vertices_size_type;
//    typedef property_map<Graph, vertex_index_t>::const_type vertex_index_map;

//    typedef std::pair<int, int> Edge;
//    enum nodes {A, B, C, D, E, n};
//    Edge edge_array[] = { Edge(A, C), Edge(B, B), Edge(B, D), Edge(B, E),
//                        Edge(C, B), Edge(C, D), Edge(D, E), Edge(E, A),
//                        Edge(E, B) };
//    int m = sizeof(edge_array) / sizeof(Edge);
//    Graph g(edge_array, edge_array + m, n);

//    // Test with the normal order
//    std::vector<vertices_size_type> color_vec(num_vertices(g));
//    iterator_property_map<vertices_size_type*, vertex_index_map>
//    color(&color_vec.front(), get(vertex_index, g));
//    vertices_size_type num_colors = sequential_vertex_coloring(g, color);
//    cout << "Number of colors = " << num_colors << endl;

    /*
    /////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////

    // Time measurement
    clock_t begin = clock();
    // Seed of random generator
    srand (time(NULL));

    vector<TestSetDescription> torontoTestSet;

    torontoTestSet.push_back(TestSetDescription("car-f-92", "Carleton University", 32));
    torontoTestSet.push_back(TestSetDescription("car-s-91", "Carleton University", 35));
    torontoTestSet.push_back(TestSetDescription("ear-f-83", "Earl Haig Collegiate", 24));
    torontoTestSet.push_back(TestSetDescription("hec-s-92", "Ecole des Hautes Etudes Commerciales", 18));
    torontoTestSet.push_back(TestSetDescription("kfu-s-93", "King Fahd University", 20));
    torontoTestSet.push_back(TestSetDescription("lse-f-91", "London School of Economics", 18));
    torontoTestSet.push_back(TestSetDescription("pur-s-93", "Purdue University", 42));
    torontoTestSet.push_back(TestSetDescription("rye-s-93", "Ryerson University", 23));
    torontoTestSet.push_back(TestSetDescription("sta-f-83", "St. Andrews High school", 13));
    torontoTestSet.push_back(TestSetDescription("tre-s-92", "Trent University", 23));
    torontoTestSet.push_back(TestSetDescription("uta-s-92", "University of Toronto, Arts & Science", 35));
    torontoTestSet.push_back(TestSetDescription("ute-s-92", "University of Toronto, Engineering", 10));
    torontoTestSet.push_back(TestSetDescription("yor-f-83", "York Mills Collegiate", 21));
    /////////////////////////////////////////////////////////////////////////////////////////////////
    //copy(torontoTestSet.begin(), torontoTestSet.end(), ostream_iterator<TestSet>(cout, "\n"));

    ofstream outFile("Results.txt");

    cout << "Test set\t # students\t # exams\t # enrolments\t conflict density\t # time slots" << endl;

//    outFile << "Test set\t # students\t # exams\t # enrolments\t conflict density\t # time slots" << endl;

//    string rootDir = "./../../Toronto Testprob/all_file_with_II_plus_Corrected IIc";

    // Version I of Toronto benchmarks
    string rootDir = "./../../Toronto Testprob/Toronto";

//    for (vector<TestSetDescription>::iterator it = torontoTestSet.begin(); it != torontoTestSet.end(); ++it) {

    // Car-f-92
//    for (vector<TestSetDescription>::iterator it = torontoTestSet.begin(); it != torontoTestSet.begin()+1; ++it) {

    // Car91
//    for (vector<TestSetDescription>::iterator it = torontoTestSet.begin()+1; it != torontoTestSet.begin()+2; ++it) {

    // Ear-f-83
//    for (vector<TestSetDescription>::iterator it = torontoTestSet.begin()+2; it != torontoTestSet.begin()+3; ++it) {

    // Hec92
//    for (vector<TestSetDescription>::iterator it = torontoTestSet.begin()+3; it != torontoTestSet.begin()+4; ++it) {


    // Kfu93
//    for (vector<TestSetDescription>::iterator it = torontoTestSet.begin()+4; it != torontoTestSet.begin()+5; ++it) {


    // Lse91
//    for (vector<TestSetDescription>::iterator it = torontoTestSet.begin()+5; it != torontoTestSet.begin()+6; ++it) {

    //Pur93
//    for (vector<TestSetDescription>::iterator it = torontoTestSet.begin()+6; it != torontoTestSet.begin()+7; ++it) {


    // Rye92
//    for (vector<TestSetDescription>::iterator it = torontoTestSet.begin()+7; it != torontoTestSet.begin()+8; ++it) {


    // Sta83
//    for (vector<TestSetDescription>::iterator it = torontoTestSet.begin()+8; it != torontoTestSet.begin()+9; ++it) {


    // Tre92
//    for (vector<TestSetDescription>::iterator it = torontoTestSet.begin()+9; it != torontoTestSet.begin()+10; ++it) {


    // Uta92
//    for (vector<TestSetDescription>::iterator it = torontoTestSet.begin()+10; it != torontoTestSet.begin()+11; ++it) {

    // Ute92
//    for (vector<TestSetDescription>::iterator it = torontoTestSet.begin()+11; it != torontoTestSet.begin()+12; ++it) {


    typedef adjacency_list<listS, vecS, undirectedS> Graph;
    typedef graph_traits<Graph>::vertex_descriptor vertex_descriptor;
    typedef graph_traits<Graph>::vertices_size_type vertices_size_type;
    typedef property_map<Graph, vertex_index_t>::const_type vertex_index_map;

    // Yor83
//    for (vector<TestSetDescription>::iterator it = torontoTestSet.begin()+12; it != torontoTestSet.begin()+13; ++it) {
    for (vector<TestSetDescription>::iterator it = torontoTestSet.begin(); it != torontoTestSet.end(); ++it) {

        // Periods range definition.
        int periods = (*it).getPeriods();

        // Specify MinPeriods and MaxPeriods

        Data data(periods, periods, periods); // Fixed length timetables

        // Create TestSet instance
        TestSet testSet((*it).getName(), (*it).getDescription(), data, rootDir);

        cout << testSet << endl;
        // Run test set
//        runNSGAII(outFile, testSet);
//        runSFLA(outFile, testSet);

        // Get graph corresponding
        Graph g = testSet.getGraph();



        // Test with the normal order
        std::vector<vertices_size_type> color_vec(num_vertices(g));
        iterator_property_map<vertices_size_type*, vertex_index_map>
        color(&color_vec.front(), get(vertex_index, g));
        vertices_size_type num_colors = sequential_vertex_coloring(g, color);
        cout << "Number of colors = " << num_colors << endl;


    }

    // Print time elapsed
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    cout << "time elapsed: " << elapsed_secs << " seconds" << endl;
*/

}



void runTorontoDatasets() {

//    // Time measurement
//    clock_t begin = clock();
//    // Seed of random generator
//    srand (time(NULL));


    vector<TestSetDescription> torontoTestSet;

    /// TODO    COLOCAR NUM FICHEIRO DEETC.CRS


    //torontoTestSet.push_back(TestSet("isel", "ISEL", 30));
    //torontoTestSet.push_back(TestSet("DEETC", "DEETC ISEL", 18, 1000));
    //		char* deetc_ucs[] = { "SET", "AA", "AC2", "ACir", "ACp", "AED", "AIEC", "ALGA", "AM1", "AM2", "Ant", "ASI", "AVE", "BD", "CAC", "CCD", "CEE", "CG", "CGAV",
    //						"CIA", "CMov", "Com", "Cpl", "CSDin", "CSDist", "CSI", "CSM", "CSO", "E1", "E2", "EA", "EGP", "Elctr", "ES", "F1", "F2", "FAE", "FIA",
    //						"FT", "GSI", "IRS", "IS", "ITI", "LC", "LIC", "LSD", "M1", "M2", "MAT", "MNO", "MSr", "OE", "OGE", "PC", "PCI", "PCM", "PDSr", "PE",
    //						"PF", "Pg", "PI", "PICC/CPg", "PIV", "POO", "PR", "PSC", "PSTR", "RCom", "RCp", "RDC", "RI", "RM", "RSCM", "SCDig", "SCDist", "SCO",
    //						"SD", "SE1", "SEAD1", "SEAD2", "SEADI", "SEAS", "SG", "SI", "SI1", "SI2", "SOi", "SOt", "SS", "ST", "STBL", "STDS" };

    //timetable.setUCNames(deetc_ucs);

    torontoTestSet.push_back(TestSetDescription("car-f-92", "Carleton University", 32));
    torontoTestSet.push_back(TestSetDescription("car-s-91", "Carleton University", 35));
    torontoTestSet.push_back(TestSetDescription("ear-f-83", "Earl Haig Collegiate", 24));
    torontoTestSet.push_back(TestSetDescription("hec-s-92", "Ecole des Hautes Etudes Commerciales", 18));
    torontoTestSet.push_back(TestSetDescription("kfu-s-93", "King Fahd University", 20));
    torontoTestSet.push_back(TestSetDescription("lse-f-91", "London School of Economics", 18));
    torontoTestSet.push_back(TestSetDescription("pur-s-93", "Purdue University", 42));
    torontoTestSet.push_back(TestSetDescription("rye-s-93", "Ryerson University", 23));
    torontoTestSet.push_back(TestSetDescription("sta-f-83", "St. Andrews High school", 13));
    torontoTestSet.push_back(TestSetDescription("tre-s-92", "Trent University", 23));
    torontoTestSet.push_back(TestSetDescription("uta-s-92", "University of Toronto, Arts & Science", 35));
    torontoTestSet.push_back(TestSetDescription("ute-s-92", "University of Toronto, Engineering", 10));
    torontoTestSet.push_back(TestSetDescription("yor-f-83", "York Mills Collegiate", 21));
    /////////////////////////////////////////////////////////////////////////////////////////////////
    //copy(torontoTestSet.begin(), torontoTestSet.end(), ostream_iterator<TestSetDescription>(cout, "\n"));

//    ofstream outFile( + ".txt");

    cout << "Test set\t # students\t # exams\t # enrolments\t conflict density\t # time slots" << endl;

//    outFile << "Test set\t # students\t # exams\t # enrolments\t conflict density\t # time slots" << endl;

//    string rootDir = "./../../Toronto Testprob/all_file_with_II_plus_Corrected IIc";

    // Version I of Toronto benchmarks
    string rootDir = "./../../Toronto Testprob/Toronto";

//    for (vector<TestSetDescription>::iterator it = torontoTestSet.begin(); it != torontoTestSet.end(); ++it) {

    // Car-f-92
//    for (vector<TestSetDescription>::iterator it = torontoTestSet.begin(); it != torontoTestSet.begin()+1; ++it) {

    // Car91
//    for (vector<TestSetDescription>::iterator it = torontoTestSet.begin()+1; it != torontoTestSet.begin()+2; ++it) {

    // Ear-f-83
//    for (vector<TestSetDescription>::iterator it = torontoTestSet.begin()+2; it != torontoTestSet.begin()+3; ++it) {

    // Hec92
//    for (vector<TestSetDescription>::iterator it = torontoTestSet.begin()+3; it != torontoTestSet.begin()+4; ++it) {

    // Kfu93
//    for (vector<TestSetDescription>::iterator it = torontoTestSet.begin()+4; it != torontoTestSet.begin()+5; ++it) {


    // Lse91
//    for (vector<TestSetDescription>::iterator it = torontoTestSet.begin()+5; it != torontoTestSet.begin()+6; ++it) {

    //Pur93
//    for (vector<TestSetDescription>::iterator it = torontoTestSet.begin()+6; it != torontoTestSet.begin()+7; ++it) {


    // Rye92
//    for (vector<TestSetDescription>::iterator it = torontoTestSet.begin()+7; it != torontoTestSet.begin()+8; ++it) {


    // Sta83
    for (vector<TestSetDescription>::iterator it = torontoTestSet.begin()+8; it != torontoTestSet.begin()+9; ++it) {


    // Tre92
//    for (vector<TestSetDescription>::iterator it = torontoTestSet.begin()+9; it != torontoTestSet.begin()+10; ++it) {


    // Uta92
//    for (vector<TestSetDescription>::iterator it = torontoTestSet.begin()+10; it != torontoTestSet.begin()+11; ++it) {

    // Ute92
//    for (vector<TestSetDescription>::iterator it = torontoTestSet.begin()+11; it != torontoTestSet.begin()+12; ++it) {


    // Yor83
//    for (vector<TestSetDescription>::iterator it = torontoTestSet.begin()+12; it != torontoTestSet.begin()+13; ++it) {

        // Periods range definition.
        int periods = (*it).getPeriods();

        // Specify MinPeriods and MaxPeriods

        Data data(periods, periods, periods); // Fixed length timetables
//        Data data(periods-1, periods+1, periods); // Variable length timetables
//        Data data(periods-4, periods+4, periods);

        // LAST USED IN NSGA-II
        //Data data(periods, periods+4, periods);


        // Create TestSet instance
        TorontoTestSet testSet((*it).getName(), (*it).getDescription(), data, rootDir);

        ofstream outFile((*it).getName() + ".txt");

        // Run test set

//        runNSGAII(outFile, testSet);
        runSFLA(outFile, testSet);


    }




//    // Print time elapsed
//    clock_t end = clock();
//    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

//    cout << "time elapsed: " << elapsed_secs << " seconds" << endl;
}



void runITC2007Datasets()
{
    // Time measurement
    clock_t begin = clock();
    // Seed of random generator
    srand (time(NULL));

    vector<TestSetDescription> testSet;

    testSet.push_back(TestSetDescription("exam_comp_set1.exam", "Dataset 1", 54));
    testSet.push_back(TestSetDescription("exam_comp_set2.exam", "Dataset 2", 40));
    testSet.push_back(TestSetDescription("exam_comp_set3.exam", "Dataset 3", 36));
    testSet.push_back(TestSetDescription("exam_comp_set4.exam", "Dataset 4", 21));
    testSet.push_back(TestSetDescription("exam_comp_set5.exam", "Dataset 5", 42));
    testSet.push_back(TestSetDescription("exam_comp_set6.exam", "Dataset 6", 16));
    testSet.push_back(TestSetDescription("exam_comp_set7.exam", "Dataset 7", 80));
    testSet.push_back(TestSetDescription("exam_comp_set8.exam", "Dataset 8", 80));
    /////////////////////////////////////////////////////////////////////////////////////////////////
    copy(testSet.begin(), testSet.end(), ostream_iterator<TestSetDescription>(cout, "\n"));

    ofstream outFile("Results.txt");

    cout << "Test set\t # students\t # exams\t # enrolments\t conflict density\t # time slots" << endl;

//    outFile << "Test set\t # students\t # exams\t # enrolments\t conflict density\t # time slots" << endl;

    // ITC 2007 benchmarks
    string rootDir = "./../../ITC2007";

    // Dataset 1
    for (vector<TestSetDescription>::iterator it = testSet.begin(); it != testSet.begin()+1; ++it) {

    // Dataset 2
//    for (vector<TestSetDescription>::iterator it = itc2007TestSet.begin()+1; it != itc2007TestSet.begin()+2; ++it) {

    // Dataset 3
//    for (vector<TestSetDescription>::iterator it = itc2007TestSet.begin()+2; it != itc2007TestSet.begin()+3; ++it) {

    // Dataset 4
//    for (vector<TestSetDescription>::iterator it = itc2007TestSet.begin()+3; it != itc2007TestSet.begin()+4; ++it) {


    // Dataset 5
//    for (vector<TestSetDescription>::iterator it = itc2007TestSet.begin()+4; it != itc2007TestSet.begin()+5; ++it) {


    // Dataset 6
//    for (vector<TestSetDescription>::iterator it = itc2007TestSet.begin()+5; it != itc2007TestSet.begin()+6; ++it) {

    // Dataset 7
//    for (vector<TestSetDescription>::iterator it = itc2007TestSet.begin()+6; it != itc2007TestSet.begin()+7; ++it) {

    // Dataset 8
//    for (vector<TestSetDescription>::iterator it = testSet.begin()+7; it != testSet.begin()+8; ++it) {

        // Periods range definition.
        int periods = (*it).getPeriods();

        // Specify MinPeriods and MaxPeriods
        Data data(periods, periods, periods); // Fixed length timetables

        cout << "periods: " << periods << endl;

        // Create ITC07TestSet instance
        ITC07TestSet testSet((*it).getName(), (*it).getDescription(), data, rootDir);

        // Run test set
//        runSFLA(outFile, testSet);
    }

    // Print time elapsed
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

    cout << "time elapsed: " << elapsed_secs << " seconds" << endl;
}



void runSFLA(ofstream& outFile, TestSet const& testSet) {

    cout << "Start Date/Time = " << currentDateTime() << endl;
    outFile << "Start Date/Time = " << currentDateTime() << endl;


    cout << "Run SFLA" << endl;

    cout << testSet << endl;
    outFile << testSet << endl;

    // SFLA parameters

//    const int m = 3;   // m is the number of memeplexes
//    const int N = 10;  // N is the number of frogs in each memeplex
//    const int q = 2;  // q is the number of frogs in each submemeplex

    // LAST
//    const int m = 3;   // m is the number of memeplexes
//    const int N = 3;  // N is the number of frogs in each memeplex
//    const int q = 3;  // q is the number of frogs in each submemeplex

//    const int m = 10;   // m is the number of memeplexes
//    const int N = 3;  // N is the number of frogs in each memeplex
//    const int q = 3;  // q is the number of frogs in each submemeplex

//    const int m = 5;   // m is the number of memeplexes
//    const int N = 15;  // N is the number of frogs in each memeplex
//    const int q = 15;  // q is the number of frogs in each submemeplex


// PAPER

//    const int m = 10;   // m is the number of memeplexes
//    const int N = 5;  // N is the number of frogs in each memeplex
//    const int q = 5;  // q is the number of frogs in each submemeplex

//    const int m = 10;   // m is the number of memeplexes
//    const int N = 20;  // N is the number of frogs in each memeplex
//    const int q = 20;  // q is the number of frogs in each submemeplex


    const int m = 10;   // m is the number of memeplexes
    const int N = 5;  // N is the number of frogs in each memeplex
    const int q = 5;  // q is the number of frogs in each submemeplex





//    const int m = 30;   // m is the number of memeplexes
//    const int N = 300;  // N is the number of frogs in each memeplex.
//    const int q = 300;


//    const int m = 30;   // m is the number of memeplexes
//    const int N = 50;  // N is the number of frogs in each memeplex.
//    const int q = 50;

//    const int m = 10;   // m is the number of memeplexes
//    const int N = 25;  // N is the number of frogs in each memeplex.
//    const int q = 25;



    const int F = m*N; // The total sample size, F, in the swamp is given by F = mN.

    // Solution initializer
    ETTPInit<eoChromosome> init(testSet.getData(), *(testSet.getConflictMatrix()), *(testSet.getGraph()));
    // Generate initial population
    eoPop<eoChromosome> pop(F, init);
    // Objective function evaluation
    eoETTPEval eval;
    // Timetable packing and Continuation
//    eoRepair<eoChromosome> repair(pop);

    // Number of consecutive time loops
//    int L = 1;

//    int L = 3; // PAPER

//    int L = 15;
//    int L = 50;
      int L = 1000000;

    eoGenContinue<eoChromosome> terminator(L); // Terminate after concluding L time loops
    eoCheckPoint<eoChromosome> checkpoint(terminator);
//    checkpoint.add(repair);
    // Chromosome evolution operator
//    eoSFLAEvolOperator_1<eoChromosome> evolutionOperator; // generic hum...
//    eoSFLAEvolOperator_1 evolutionOperator; // generic hum...

    // The best
    eoSFLAEvolOperator_2 evolutionOperator; // generic hum...


    //    eoSFLAEvolOperator_3 evolutionOperator; // generic hum...

    // Build SFLA
    eoSFLA<eoChromosome> sfla(m, N, q, F, L, init, checkpoint, eval, evolutionOperator);

    // Run the algorithm
    sfla(pop);

    cout << "Fim" << endl;

//    outFile << sfla.getBestFrog().fitness() << " - " << sfla.getBestFrog().getNumPeriods() << endl;

    // Print best frog
    outFile << sfla.getBestFrog().getChromosome() << endl;

    cout << "End Date/Time = " << currentDateTime() << endl;
    outFile << "End Date/Time = " << currentDateTime() << endl;


//    cin.get();
}



void runNSGAII(ofstream& outFile, TestSet const& testSet) {

    cout << testSet << endl;
    outFile << testSet << endl;

    // NSGA-II parameters
//    const unsigned int POP_SIZE = 50; // Population size
//    const unsigned int POP_SIZE = 10; // Population size
//    const unsigned int POP_SIZE = 30; // Population size
        const unsigned int POP_SIZE = 100; // Population size
//    const unsigned int POP_SIZE = 200; // Population size


    //    const unsigned int MAX_GEN = 100; // Maximum number of generations
//    const unsigned int MAX_GEN = 50; // Maximum number of generations
    const unsigned int MAX_GEN = 200; // Maximum number of generations
//    const unsigned int MAX_GEN = 5; // Maximum number of generations


        const double P_CROSS = 0.8; // Crossover probability
//    const double P_CROSS = 1; // Crossover probability

//    const double P_CROSS = 0; // Crossover probability

    /// Com 0.8 nao consegue apanhar a frente completa..

    const double P_MUT = 0.15; // Mutation probability
//    const double P_MUT = 0; // Mutation probability


    const double REINSERTION_RATE = 0.02; // Rate at which exams are reinserted for each chromosome selected on mutation rate.

    ETTPInit<moeoChromosome> init(testSet.getData(), *(testSet.getConflictMatrix()), *(testSet.getGraph()));

    // eoChromosome solution initializer
    ETTPInit<eoChromosome> eoinit(testSet.getData(), *(testSet.getConflictMatrix()), *(testSet.getGraph()));

    NSGAII nsgaii(POP_SIZE, MAX_GEN, P_CROSS, P_MUT, REINSERTION_RATE, init, eoinit);

    // Run NSGA-II
    nsgaii.run(outFile);
}




#ifndef TESTSETDESCRIPTION_H
#define TESTSETDESCRIPTION_H

#include <string>

using namespace std;


class TestSetDescription {
    string name;
    string description;
    int numPeriods;
public:
    TestSetDescription(string name, string description, int numPeriods) :
        name(name), description(description), numPeriods(numPeriods) { }
    string getName() const { return name; }
    string getDescription() const { return description; }
    int getPeriods() const { return numPeriods; }
    friend ostream& operator<<(ostream& os, const TestSetDescription& t);
};


#endif // TESTSETDESCRIPTION_H

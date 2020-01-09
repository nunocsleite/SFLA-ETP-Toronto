#ifndef ITC07TESTSET_H
#define ITC07TESTSET_H

#include "TestSet.h"
#include <string>


using namespace std;


///////////////////////////////////////
// Unary and Binary hard constraint
// definition
///////////////////////////////////////

class UnaryHardConstraint : public eoUF<int, bool>
{};

class BinaryHardConstraint : public eoBF<int, int, bool>
{};


// E.g.
// [PeriodHardConstraints]
// 11, AFTER, 10
class AfterConstraint : public BinaryHardConstraint {
    bool operator()(int e1, int e2) {
        // TODO
    }
};

// [PeriodHardConstraints]
// 425, EXAM_COINCIDENCE, 426
class ExamCoincidenceConstraint : public BinaryHardConstraint {
    bool operator()(int e1, int e2) {
        // TODO
    }
};

// [PeriodHardConstraints]
// 100, EXCLUSION, 120
class ExclusionConstraint : public BinaryHardConstraint {
    bool operator()(int e1, int e2) {
        // TODO
    }
};

// [RoomHardConstraints]
// 2, ROOM_EXCLUSIVE
class RoomExclusiveConstraint : public UnaryHardConstraint {
    bool operator()(int e) {
        // TODO
    }
};

//////////////////////////////////////////




class ITC07TestSet : public TestSet {

public:
    ITC07TestSet(string _testSetName, string _description, Data& _data, string _rootDir)
        : TestSet(_testSetName, _description, _data, _rootDir) {

        cout << "ITC07TestSet ctor" << endl;

        // Load data
        load(_testSetName, _rootDir);
    }


protected:
    void load(string _testSetName, string _rootDir);

private:
    // Read exams and students
    void readExams(tokenizer<escaped_list_separator<char> >::iterator& it,
                   tokenizer<escaped_list_separator<char> > const& tok);
    // Read periods
    void readPeriods(tokenizer<escaped_list_separator<char> >::iterator& it,
                     tokenizer<escaped_list_separator<char> > const& tok);
    // Read rooms
    void readRooms(tokenizer<escaped_list_separator<char> >::iterator& it,
                   tokenizer<escaped_list_separator<char> > const& tok);
    // Read constraints and weightings
    void readConstraints(tokenizer<escaped_list_separator<char> >::iterator& it,
                         tokenizer<escaped_list_separator<char> > const& tok);


    ///////////////////////////
    // Inner classes
    ///////////////////////////
    struct SeqLineInfo {
        // Period info
        int day;
        int month;
        int year;

        int hour;
        int minute;
        int second;

        int duration;
        int penalty;

        // Room info
        int roomCapacity;
        int roomPenalty;
    };


    ///////////////////////////
    // Auxiliary methods
    ///////////////////////////
    boost::shared_ptr<AdjacencyList> buildExamGraph(boost::shared_ptr<Matrix> const& ptrConflictMatrix);
    int matchPeriods(tokenizer<escaped_list_separator<char> >::iterator& it);
    int matchPeriodSequenceLine(tokenizer<escaped_list_separator<char> >::iterator& it, SeqLineInfo& seqLineInfo);

    int matchRooms(tokenizer<escaped_list_separator<char> >::iterator& it);
    int matchRoomSequenceLine(tokenizer<escaped_list_separator<char> >::iterator& it, SeqLineInfo& seqLineInfo);

    void readPeriodHardConstraints(tokenizer<escaped_list_separator<char> >::iterator& it,
                                   vector<BinaryHardConstraint> &periodConstraints);
    bool readPeriodConstraint(tokenizer<escaped_list_separator<char> >::iterator& it,
                              vector<BinaryHardConstraint>& periodConstraints);


    ///////////////////////////
    // Fields
    ///////////////////////////


};


#endif // ITC07TESTSET_H

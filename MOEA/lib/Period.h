#ifndef PERIOD_H
#define PERIOD_H


#include <boost/shared_ptr.hpp>
#include <vector>
#include "Exam.h"


using namespace std;
using namespace boost;


class Period {
    vector<boost::shared_ptr<Exam> > examsToSchedule;
    int id;

public:
    Period() : examsToSchedule(0), id(0) { }
    Period(int id) : examsToSchedule(0), id(id) { }
    int getId() const { return id; }
    void setId(int id) { this->id = id; }
    vector<boost::shared_ptr<Exam> > const& getExams() const { return examsToSchedule; }
    int getCapacity() const { return examsToSchedule.size(); }
    void addExam(boost::shared_ptr<Exam> const& exam) {
        examsToSchedule.push_back(exam);
    }
    bool operator<(Period const& period) const {
        return getCapacity() < period.getCapacity();
    }

    friend ostream& operator<<(ostream& out, Period const& period) {
        out << "Id = " << period.getId() << endl;
        out << "capacity = " << period.getCapacity() << endl;
        out << "Exams = ";
        out << "Exams size = " << period.getExams().size() << endl;
        for (vector<boost::shared_ptr<Exam> >::const_iterator it = period.getExams().begin(); it != period.getExams().end(); ++it) {
            out << (*it).get()->getId() << endl;
        }
        return out;
    }
};



//class Period {
////    vector<shared_ptr<Exam> > examsToSchedule;
//    vector<Exam*> examsToSchedule;
//    int capacity;
//    int id;

//public:
//    Period() : examsToSchedule(0), capacity(0), id(0) { }
//    Period(int id) : examsToSchedule(0), capacity(0), id(id) { }
//    ~Period() {
//        for (vector<Exam*>::iterator it = examsToSchedule.begin(); it != examsToSchedule.end(); ++it)
//            delete *it;
//    }
//    Period(Period const& period) {
//        for (vector<Exam*>::const_iterator it = period.examsToSchedule.begin(); it != period.examsToSchedule.end(); ++it)
//            examsToSchedule.push_back(new Exam((*it)->getId()));
//        capacity = period.capacity;
//        id = period.id;
//    }
//    Period& operator=(Period const& period) {
//        examsToSchedule.clear();
//        for (vector<Exam*>::const_iterator it = period.examsToSchedule.begin(); it != period.examsToSchedule.end(); ++it)
//            examsToSchedule.push_back(new Exam((*it)->getId()));
//        capacity = period.capacity;
//        id = period.id;
//    }
//    int getId() const { return id; }
//    void setId(int id) { this->id = id; }
////    vector<shared_ptr<Exam> > const& getExams() const { return examsToSchedule; }
//    vector<Exam*> const& getExams() const { return examsToSchedule; }
//    int getCapacity() const { return capacity; }
////    void addExam(shared_ptr<Exam> const& exam) {
//    void addExam(Exam* exam) {
//        examsToSchedule.push_back(exam);
//        ++capacity;

////        cout << "addExam -> " << *this << endl;
//    }
//    bool operator<(Period const& period) const {
//        return capacity < period.capacity;
//    }

//    friend ostream& operator<<(ostream& out, Period const& period) {
//        out << "Id = " << period.getId() << endl;
//        out << "capacity = " << period.getCapacity() << endl;
//        out << "Exams = ";
//        out << "Exams size = " << period.getExams().size() << endl;
////        for (vector<shared_ptr<Exam> >::const_iterator it = period.getExams().begin(); it != period.getExams().end(); ++it) {
////            out << (*it).get()->getId() << endl;
////        }
//        for (vector<Exam*>::const_iterator it = period.getExams().begin(); it != period.getExams().end(); ++it) {
//            out << (*it)->getId() << endl;
//        }
//        return out;
//    }
//};

#endif // PERIOD_H

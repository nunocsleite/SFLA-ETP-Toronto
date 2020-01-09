
#include "TestSetDescription.h"


ostream& operator<<(ostream& os, const TestSetDescription& t) {
//    os << "Test set\t # students\t # exams\t # enrolments\t conflict density\t # time slots" << endl;

//    os << t.getName() << "\t " << t.getNumStudents() << "\t " << t.getNumExams() << "\t "
//       << t.getNumEnrolments() << "\t " << t.getConflictDensity() << "\t " << t.getNumPeriods() << endl;

    os << t.getName();

    return os;
}


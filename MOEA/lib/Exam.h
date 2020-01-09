#ifndef EXAM_H
#define EXAM_H


class Exam {
    int id;
    bool scheduled;
    int period;
public:
    Exam(int id, int p) : id(id), scheduled(false), period(p) { }
    Exam(int id) : id(id), scheduled(false) { }
    int getId() const { return id; }
    bool isAlreadyScheduled() { return scheduled; }
    void schedule() { scheduled = true; }
    int getPeriod() { return period; }
//    void setPeriod(int p) { period = p; }
};


#endif // EXAM_H

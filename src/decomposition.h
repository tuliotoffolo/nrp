/**
 * This class creates and keeps the sub-problems derived from the
 * decomposition approaches, ensuring each one is solved at most once
 * unless an improvement may be obtained.
 *
 * @author Tulio Toffolo
 */

#ifndef DECOMPOSITION_H_INCLUDED
#define DECOMPOSITION_H_INCLUDED

#include "instance.h"
#include "solution.h"

// decompositions
#define DEC_DAYS   1
#define DEC_NURSES 2
#define DEC_SHIFTS 3


struct Subproblem
{
    int index;
    int sortedIndex;
    char name[20];
    std::vector<int> varsToFix;
};


class Decomposition
{
protected:
    int nSubproblems;
    Instance &instance;

    Decomposition(int index, Instance &instance)
            : index(index), instance(instance)
    {}

    int getVarIndex(int nurse, int shift, int day);

public:
    const int index;

    int getNSubproblems()
    { return nSubproblems; }

    virtual Subproblem *getSubproblem()
    { return nullptr; }

    virtual void reset(bool improved)
    {}

    virtual void reset(bool improved, Subproblem *ignore)
    {}
};


class DaysDecomposition : public Decomposition
{
private:
    int nDays, daysStep, maxDays;
    std::vector<Subproblem *> subproblems, sortedSubproblems;

    void makeSubproblems();

public:
    DaysDecomposition(Instance &instance, int nDays, int daysStep, int maxDays);
    Subproblem *getSubproblem();
    void reset(bool improved);
    void reset(bool improved, Subproblem *ignore);
};


class NursesDecomposition : public Decomposition
{
private:
    //int nNurses, maxNurses;
    //vector<int *> combinations;
    //Subproblem *refSubproblem;

    int nNurses, nursesStep, maxNurses;
    std::vector<Subproblem *> subproblems, sortedSubproblems;

    void makeSubproblems();


public:
    NursesDecomposition(Instance &instance, int nNurses, int nursesStep, int maxNurses);
    Subproblem *getSubproblem();
    void reset(bool improved);
    void reset(bool improved, Subproblem *ignore);
};


class ShiftsDecomposition : public Decomposition
{
private:
    std::vector<Subproblem *> subproblems, sortedSubproblems;

public:
    ShiftsDecomposition(Instance &instance);
    Subproblem *getSubproblem();
    void reset(bool improved);
    void reset(bool improved, Subproblem *ignore);
};


class DecompositionsManager
{
private:
    std::vector<Decomposition *> decompositions;
    std::vector<Decomposition *> sortedDecompositions;

    bool improved;
    int lastDecomposition;
    Subproblem *lastSubproblem;

public:
    DecompositionsManager(Instance &instance, int nDays, int daysStep, int maxDays, int nNurses, int nursesStep, int maxNurses);
    Subproblem *getSubproblem(ostream &output);
    void reset();
};

#endif //DECOMPOSITION_H_INCLUDED

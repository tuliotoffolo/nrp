/**
 * This class represents a solution for the NRP.
 *
 * @author Tulio Toffolo
 */

#ifndef SOLUTION_H_INCLUDED
#define SOLUTION_H_INCLUDED

#include <vector>
#include <set>
#include <map>
#include <string>
#include <cstring>
#include <cstdlib>
#include <cfloat>

using namespace std;

#include "instance.h"

#ifdef CPLEX

#include "ilcplex/cplex.h"
#include "ilcplex/ilocplex.h"

#endif

class Solution
{
public:
    Instance &instance;
    vector<double> x;
    vector<int> idx;
    double objval;
    double totaltime;

    Solution(Instance &_instance) : instance(_instance)
    {
        objval = DBL_MAX;
        totaltime = 0;
        int size = instance.nNurses * instance.nDays * instance.nShifts;
        for (int i = 0; i < size; i++) {
            x.push_back(0);
            idx.push_back(-1);
        }
    }

    Solution(Instance &_instance, vector<int> &_idx) : instance(_instance)
    {
        objval = DBL_MAX;
        totaltime = 0;
        int size = instance.nNurses * instance.nDays * instance.nShifts;
        for (int i = 0; i < size; i++) {
            x.push_back(0);
            idx.push_back(_idx[i]);
        }
    }

    Solution(Solution *sol) : instance(sol->instance)
    {
        this->objval = sol->objval;
        this->totaltime = sol->totaltime;

        this->x.reserve(sol->x.size());
        this->idx.reserve(sol->idx.size());
        for (int i = 0; i < sol->x.size(); i++) {
            this->x.push_back(sol->x[i]);
            this->idx.push_back(sol->idx[i]);
        }
    }

    static bool compare(Solution *s1, Solution *s2);

    void read(ifstream &input);
    void read(double *vars, double _objval);

#ifdef CPLEX
    void read(IloNumArray &cols, double _objval);
#endif

    void save(const char *fileName, double time);
    void save(const char *fileName);
    void saveXML(const char *fileName, const time_t startDate);

    double getGap(double LB);

    int getNurse(int index);
    int getNurseBegin(int nurse);
    int getNurseEnd(int nurse);

    int getX(int index);
    int getX(int nurse, int shift, int day);
    int getX(int nurse, string shift, int day);

    int getIdx(int nurse, int shift, int day);
    int getIdx(int nurse, string shift, int day);

    void setX(int nurse, int shift, int day, int value);
    void setX(int nurse, string shift, int day, int value);

    void setIdx(int nurse, int shift, int day, int value);
    void setIdx(int nurse, string shift, int day, int value);

    void printSolution(ostream &output);

    double getVariacao(Solution *sol);

    double getTime()
    { return totaltime; }

private:
    int getIndex(int nurse, int shift, int day);
    int getIndex(int nurse, string shift, int day);
};

class NurseWeight
{
public:
    int nurse;
    double weight;

    NurseWeight()
    {
        this->nurse = -1;
        this->weight = 0;
    }

    NurseWeight(int nurse, double weight)
    {
        this->nurse = nurse;
        this->weight = weight;
    }

    static bool compare(NurseWeight s1, NurseWeight s2);
};

#endif // SOLUTION_H_INCLUDED

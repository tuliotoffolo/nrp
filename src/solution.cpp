/**
 * This class represents a solution for the NRP.
 *
 * @author Tulio Toffolo
 */

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <cmath>
#include "solution.h"

bool Solution::compare(Solution *s1, Solution *s2)
{
    return s1->objval < s2->objval;
}

void Solution::read(ifstream &input)
{
    char obj[10], time[10], eq[10];
    input >> obj >> eq >> objval;
    for (int i = 0; i < x.size(); i++)
        x[i] = 0;

    int shift;
    for (int nurse = 0; nurse < instance.nNurses; nurse++) {
        for (int day = 0; day < instance.nDays; day++) {
            input >> shift;
            if (shift >= 0)
                setX(nurse, shift, day, 1);
        }
    }
    if (strcmp(obj, "Obj.") == 0) {
        input >> time >> eq;
        if (strcmp(time, "Time") == 0)
            input >> totaltime;
    }
}

void Solution::read(double *cols, double _objval)
{
    this->objval = _objval;
    for (int i = 0; i < x.size(); i++)
        x[i] = cols[idx[i]];
}

#ifdef CPLEX

void Solution::read(IloNumArray &cols, double _objval)
{
    this->objval = _objval;
    for (int i = 0; i < x.size(); i++)
        x[i] = cols[idx[i]];
}

#endif

void Solution::printSolution(ostream &output)
{
    for (int nurse = 0; nurse < instance.nNurses; nurse++)
        for (int day = 0; day < instance.nDays; day++)
            for (int shift = 0; shift < instance.nShifts; shift++)
                if (round(getX(nurse, shift, day)) == 1)
                    output << "X(N" << nurse << "," << instance.shiftNames[shift] << "," << day << ")"
                           << " = " << getX(nurse, shift, day) << endl;
}

void Solution::save(const char *fileName, double time)
{
    ofstream output(fileName, ios::out);
    output << "Obj. = " << objval << endl;
    int s;
    for (int nurse = 0; nurse < instance.nNurses; nurse++) {
        for (int day = 0; day < instance.nDays; day++) {
            s = -1;
            for (int shift = 0; shift < instance.nShifts; shift++)
                if (getX(nurse, shift, day) == 1)
                    s = shift;
            output << s << " ";
        }
        output << endl;
    }
    double ctime = time > 600 ? 600 : time;
    output << "Time = " << ctime << endl;
    output.close();
}

void Solution::save(const char *fileName)
{
    save(fileName, totaltime);
}

void Solution::saveXML(const char *fileName, const time_t startDate)
{
    int s;

    FILE *f = fopen(fileName, "w");
    if (!f) {
        fprintf(stderr, "Could not open file %s", fileName);
        exit(EXIT_FAILURE);
    }

    fprintf(f, "<Solution>\n");
    fprintf(f, "    <SchedulingPeriodID>%s</SchedulingPeriodID>\n", instance.problemName);
    fprintf(f, "    <Competitor>Tulio A. M. Toffolo</Competitor>\n");
    for (int nurse = 0; nurse < instance.nNurses; nurse++) {
        for (int day = 0; day < instance.nDays; day++) {
            s = -1;
            for (int shift = 0; shift < instance.nShifts; shift++) {
                if (getX(nurse, shift, day) == 1) {
                    fprintf(f, "    <Assignment>\n");
                    time_t date = Instance::dateForDay(&startDate, day);
                    struct tm btime;
                    localtime_r(&date, &btime);
                    fprintf(f, "        <Date>%d-%d-%d</Date>\n", btime.tm_year + 1900, btime.tm_mon + 1, btime.tm_mday);
                    fprintf(f, "        <Employee>%d</Employee>\n", nurse);
                    fprintf(f, "        <ShiftType>%s</ShiftType>\n", instance.shiftNames[shift].c_str());
                    fprintf(f, "    </Assignment>\n");
                }
            }
        }
    }
    fprintf(f, "</Solution>\n");

    fclose(f);
}

double Solution::getGap(double LB)
{
    if (LB == 0) return objval;
    return (objval - LB) / objval;
}

int Solution::getNurse(int index)
{
    return index / instance.nNurses;
}

int Solution::getNurseBegin(int nurse)
{
    return nurse * instance.nShifts * instance.nDays;
}

int Solution::getNurseEnd(int nurse)
{
    return nurse * instance.nShifts * instance.nDays + instance.nShifts * instance.nDays;
}

int Solution::getX(int index)
{
    return x[index];
}

int Solution::getX(int nurse, int shift, int day)
{
    return x[getIndex(nurse, shift, day)];
}

int Solution::getX(int nurse, string shift, int day)
{
    return x[getIndex(nurse, shift, day)];
}

int Solution::getIdx(int nurse, int shift, int day)
{
    return idx[getIndex(nurse, shift, day)];
}

int Solution::getIdx(int nurse, string shift, int day)
{
    return idx[getIndex(nurse, shift, day)];
}

void Solution::setX(int nurse, int shift, int day, int value)
{
    x[getIndex(nurse, shift, day)] = value;
}

void Solution::setX(int nurse, string shift, int day, int value)
{
    x[getIndex(nurse, shift, day)] = value;
}

void Solution::setIdx(int nurse, int shift, int day, int value)
{
    idx[getIndex(nurse, shift, day)] = value;
}

void Solution::setIdx(int nurse, string shift, int day, int value)
{
    idx[getIndex(nurse, shift, day)] = value;
}

int Solution::getIndex(int nurse, int shift, int day)
{
    int index = nurse * instance.nShifts * instance.nDays
                + shift * instance.nDays
                + day;
#ifdef DEBUG
    if (index > x.size()) {
        cerr << "Index error concerning vector X within class Solution:" << endl;
        cerr << "  size = " << x.size() << " and requested index = " << index << endl;
        exit(0);
    }
#endif
    return index;
}

int Solution::getIndex(int nurse, string shift, int day)
{
    return getIndex(nurse, instance.getShiftIndex(shift), day);
}

double Solution::getVariacao(Solution *sol)
{
    int contDifs = 0, contTotal = 0;
    for (int i = 0; i < x.size(); i++) {
        if (this->x[i] == 1) {
            contTotal++;
            if (this->x[i] != sol->x[i])
                contDifs++;
        }
    }

    return (double) contDifs / (double) contTotal;
}

bool NurseWeight::compare(NurseWeight s1, NurseWeight s2)
{
    return s1.weight < s2.weight;
}

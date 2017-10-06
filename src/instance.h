/**
 * This class represents an instance of the NRP.
 *
 * @author Haroldo Santos
 * @author Tulio Toffolo
 */

#ifndef INSTANCE_H_INCLUDED
#define INSTANCE_H_INCLUDED

#include <vector>
#include <set>
#include <map>
#include <string>
#include <time.h>

using namespace std;

class Instance
{
public:
    Instance(const char *fileName);

    /* reads the start of the time horizon from a txt instance file */
    static time_t readDateFromTXT(const char *fileName);

    int nDays;       /* days numbered from 0..nDays-1 */
    int nWeekends;   /* number of weekends in the month */
    int nNurses;     /* nurses numbered from 0 */
    int nShifts;

    char problemName[256];

    /* original shift names */
    vector<string> shiftNames;
    vector<char> shiftC;   /* shiftC[i] returns the character which
                                 uniquely identifies shift C */

    int nContracts;

    /* number of weekends considered for each contract */
    vector<int> nWE;

    vector<vector<vector<int> > > WE;

    vector<int> nurseContract;  /* stores the index of the contract of a nurse */

    vector<vector<bool> > NurseShift; /* NurseShift[i,j] returns true if nurse i can work at shift j */

    /* at a given day and shift the required demand of nurses */
    vector<vector<int> > demand;

    /* days and shift where a nurse does not wants to work:
       nx[i,j,k]  returns 0 if nuse i does WANTS to work at shift j at day k
       or greater than zero (a penalty) for undesired days/shifts */

    vector<vector<vector<int> > > nx;

    /* maximum number of assignments per contract */
    vector<int> na_Max;

    /* maximum number of assignments per contract */
    vector<int> wna_Max;

    /* minimum number of assignments per contract */
    vector<int> na_Min;

    /* minimum number of assignments per contract */
    vector<int> wna_Min;

    /* maximum number of consecutive working days */
    vector<int> nw_Max;
    vector<int> wnw_Max;

    /* minimum number of consecutive working days */
    vector<int> nw_Min;
    vector<int> wnw_Min;

    /* maximum number of consecutive free days */
    vector<int> nf_Max;
    vector<int> wnf_Max;

    /* minimum number of consecutive free days */
    vector<int> nf_Min;
    vector<int> wnf_Min;

    /* maximum number of consecutive working weekends */
    vector<int> nwec_Max;
    vector<int> wnwec_Max;

    /* minimum number of consecutive working weekends */
    vector<int> nwec_Min;
    vector<int> wnwec_Min;

    /* Complete weekends */
    vector<int> fw;
    vector<int> wfw;

    /* Weight of identical shift types during weekend */
    vector<int> ssw;
    vector<int> wssw;

    static time_t dateForDay(const time_t *date, const int nDay);

    int getShiftIndex(const string str) const;
private:
    /* will read at most 1024 cols of size 1024 per line in the dat file */
    char **cols;
    int readCols(FILE *f);

    /* remove ; . */
    void removeGarbage(char *str);

    /* removes non numeric characters */
    void removeAlpha(char *str);

    /* check contents of current line */
    bool isParam(const char *param);
    bool isSet(const char *setName);
    /* checks if it is a set whose name starts with setname */
    bool isSetStart(const char *setName);

    /* returns true if there is already one shift identified by this char */
    bool shiftCUsed(char c);

    /* for each char stores pre-computed the shift index, -1 for characters which
       do not correspont to shifts  */
    vector<int> idxShift;

    map<string, int> shiftIndex;

};

#endif // INSTANCE_H_INCLUDED

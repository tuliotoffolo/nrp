/**
 * This class represents an instance of the NRP.
 *
 * @author Haroldo Santos
 * @author Tulio Toffolo
 */

#include <cstdio>
#include <string>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <string.h>
#include <ctype.h>
#include "instance.h"

/**
 * splits a string usign a delimiter char
 * returns how many columns were found
 * multDel = 1 indicates that multiple delimiters
 * grouped together represent only one
 **/
int splitString(char **columns, const char *str, const char delimiter,
                const int maxColumns, const int columnSize, const char multDel);

Instance::Instance(const char *fileName)
{
    FILE *f = fopen(fileName, "r");
    if (!f) {
        fprintf(stderr, "Could not open file %s", fileName);
        exit(EXIT_FAILURE);
    }

    // initialization of cols
    cols = new char *[1024];
    cols[0] = new char[1024 * 1024];
    unsigned int it1D;
    for (it1D = 0; (it1D < 1024); it1D++)
        cols[it1D] = &(cols[0][(it1D * 1024)]);

    strcpy(problemName, fileName);

    // intial value, actual value comes after reading param r
    nDays = 512;

    // removing path
    while (strstr(problemName, "/")) {
        char str[256];
        strcpy(str, strstr(problemName, "/") + 1);
        strcpy(problemName, str);
    }
    while (strstr(problemName, "\\")) {
        char str[256];
        strcpy(str, strstr(problemName, "\\") + 1);
        strcpy(problemName, str);
    };
    // and extension
    char *s;
    if ((s = strstr(problemName, "."))) {
        if (s[1] == 'd')
            *s = '\0';
    }

    if ((s = strstr(problemName, "_"))) {
        if (s[1] == 'j')
            *s = '\0';
        else if ((s = strstr(strstr(problemName, "_") + 1, "_")))
            if (s[1] == 'j')
                *s = '\0';
    }


    nContracts = 0;

    int ncols;

    // every possible ascii character
    idxShift.resize(256, -1);

    while ((ncols = readCols(f)) != -1) {
        if (ncols == 0)
            continue;

        if (isParam("nW")) {
            removeGarbage(cols[3]);
            this->nWeekends = atoi(cols[3]);
        }
        else if (isSet("S")) {
            ncols = readCols(f);
            assert(ncols >= 1);
            for (int i = 0; (i < ncols); ++i) {
                removeGarbage(cols[i]);
                this->shiftNames.push_back(string(cols[i]));

                char c = cols[i][0];

                while (shiftCUsed(c))
                    ++c;

                shiftC.push_back(c);

                idxShift[c] = shiftC.size() - 1;

                shiftIndex[cols[i]] = shiftC.size() - 1;

                //printf( "%s %c\n", shiftNames[i].c_str(), shiftC[i] );
            }

            nShifts = shiftC.size();
        }
        else if (isSet("C")) {
            this->nContracts = 0;
            do {
                ncols = readCols(f);
                assert(ncols >= 1);
                if ((ncols == 1) && (cols[0][0] == ';'))
                    break;

                this->nContracts += ncols;
            }
            while (1);

            this->nWE.resize(this->nContracts, 0);
            this->WE.resize(this->nContracts, vector<vector<int> >(10, vector<int>(4, -1)));

        }
        else if (isParam("nWE")) {
            do {
                ncols = readCols(f);
                assert(ncols >= 1);
                if ((ncols == 1) && (cols[0][0] == ';'))
                    break;
                if (ncols == 2) {
                    removeGarbage(cols[0]);
                    removeGarbage(cols[1]);
                    removeAlpha(cols[0]);
                    int nContract = atoi(cols[0]);
                    int nwe = atoi(cols[1]);
                    this->nWE[nContract] = nwe;
                    this->WE[nContract].resize(nwe, vector<int>(4, -1));
                }
            }
            while (1);
        }
        else if (isSet("N")) {
            this->nNurses = 0;
            do {
                ncols = readCols(f);
                assert(ncols >= 1);
                if ((ncols == 1) && (cols[0][0] == ';'))
                    break;

                this->nNurses += ncols;
            }
            while (1);
        }
        if (isParam("r")) {
            ncols = readCols(f);
            assert(ncols == (int) shiftC.size() + 1);

            vector<int> shiftI(shiftC.size(), -1);

            for (int i = 0; i < (int) shiftC.size(); ++i) {
                removeGarbage(cols[i]);
                int idx = getShiftIndex(cols[i]);
                shiftI[i] = idx;
                assert((idx >= 0) && (idx < (int) shiftC.size()));
            }
            do {
                ncols = readCols(f);
                assert(ncols >= 1);
                if ((ncols == 1) && (cols[0][0] == ';'))
                    break;
                assert(ncols > (int) shiftC.size());

                demand.push_back(vector<int>(shiftC.size(), 0));
                vector<int> &v = demand[demand.size() - 1];
                for (int i = 1; i <= (int) shiftC.size(); ++i) {
                    removeGarbage(cols[i]);

                    int shiftIdx = shiftI[i - 1];

                    v[shiftIdx] = atoi(cols[i]);
                }
            }
            while (1);

            nDays = demand.size();

            if (nNurses > (int) nx.size())
                nx.resize(nNurses, vector<vector<int> >(nShifts, vector<int>(nDays, 0)));

        }
        else if (isSetStart("NS[")) {
            //printf("%s %s\n", cols[0], cols[1] );
            char *p = strstr(cols[1], "[");
            removeGarbage(p);
            removeAlpha(p);

            //printf("N : %d\n", atoi(p) );

            int nIdx = atoi(p);

            if (nIdx + 1 > (int) NurseShift.size()) {
                NurseShift.resize(nIdx + 1, vector<bool>((int) shiftC.size(), false));
            }

            vector<bool> &vNurse = NurseShift[nIdx];
            do {
                ncols = readCols(f);
                assert(ncols >= 1);
                if ((ncols == 1) && (cols[0][0] == ';'))
                    break;

                for (int i = 0; (i < ncols); ++i) {
                    removeGarbage(cols[i]);
                    int sIdx = getShiftIndex(string(cols[i]));
                    assert((sIdx >= 0) && (sIdx < (int) shiftC.size()));
                    vNurse[sIdx] = true;
                }
            }
            while (1);

            //            printf("nurse shift compatibility:\n");
            //            for ( int i=0 ; (i<NurseShift.size()) ; ++i  )
            //            {
            //                for ( int j=0 ; (j<shiftC.size() ) ; ++j )
            //                {
            //                    if (NurseShift[i][j])
            //                        printf("1 ");
            //                    else
            //                        printf("0 ");
            //
            //                }
            //                printf("\n");
            //            }

            nNurses = NurseShift.size();
        }
        else if (isParam("nx")) {
            do {
                ncols = readCols(f);
                if ((ncols == 1) && (cols[0][0] == ';'))
                    break;

                assert(ncols >= 1);

                int nurse, shiftIdx, day;
                removeGarbage(cols[0]);
                removeGarbage(cols[1]);

                //printf("%s %s\n", cols[0], cols[1]);
                char *sptr, *tk;
                tk = strtok_r(cols[0], " ", &sptr);
                assert(tk);
                nurse = atoi(tk + 1);
                tk = strtok_r(NULL, " ", &sptr);
                assert(tk);
                shiftIdx = getShiftIndex(tk);
                tk = strtok_r(NULL, " ", &sptr);
                assert(tk);
                day = atoi(tk);

                //printf("%d %d %d   %d\n", nurse, shiftIdx, day, atoi(cols[1]) );

                /* since days is specified latter, giving a slack */
                if (nurse + 1 > (int) nx.size())
                    nx.resize(max(nurse + 1, nNurses), vector<vector<int> >(nShifts, vector<int>(max(nDays, 512), 0)));

                assert((nurse < (int) nx.size()) && (nurse >= 0));
                assert((shiftIdx < (int) nx[nurse].size()) && (shiftIdx >= 0));
                assert((day - 1 < (int) nx[nurse][shiftIdx].size()) && (day - 1 >= 0));

                nx[nurse][shiftIdx][day - 1] = atoi(cols[1]); //TODO verificar -1
            }
            while (1);
        }
        else if (isParam("cn")) {
            do {
                ncols = readCols(f);
                if ((ncols == 1) && (cols[0][0] == ';'))
                    break;

                removeGarbage(cols[0]);
                removeGarbage(cols[1]);

                int n = atoi(cols[0] + 1);
                int c = atoi(cols[1] + 1);

                //printf( "%d %d\n", n, c );

                if (n + 1 > (int) nurseContract.size())
                    nurseContract.resize(n + 1, 0);
                nurseContract[n] = c;

                nContracts = max(c + 1, nContracts);
            }
            while (1);
        }
        else if (isParam("na_Max")) {
            assert(nContracts > 0);
            assert(nDays > 0);
            na_Max.resize(nContracts, nDays);
            do {
                ncols = readCols(f);


                if ((ncols == 1) && (cols[0][0] == ';'))
                    break;

                removeGarbage(cols[0]);
                removeGarbage(cols[1]);

                int c = atoi(cols[0] + 1);
                int m = atoi(cols[1]);

                na_Max[c] = m;
            }
            while (1);
        }
        else if (isParam("na_Min")) {
            assert(nContracts > 0);
            assert(nDays > 0);
            na_Min.resize(nContracts, 0);
            do {
                ncols = readCols(f);

                if ((ncols == 1) && (cols[0][0] == ';'))
                    break;

                removeGarbage(cols[0]);
                removeGarbage(cols[1]);

                int c = atoi(cols[0] + 1);
                int m = atoi(cols[1]);

                na_Min[c] = m;
            }
            while (1);
        }
        else if (isSetStart("WE[")) {
            //         printf("%s --- ", cols[1]);
            char *p = strstr(cols[1], "[") + 1;
            removeGarbage(p);
            removeAlpha(p);
            //       printf("%d %s\n ", atoi(p), cols[2] );
            removeGarbage(cols[2]);
            removeAlpha(cols[2]);

            int c = atoi(p);
            int nWeekend = atoi(cols[2]) - 1;

            do {
                ncols = readCols(f);

                if ((ncols == 1) && (cols[0][0] == ';'))
                    break;
                WE[c][nWeekend].clear();
                for (int i = 0; (i < ncols); ++i) {
                    removeGarbage(cols[i]);
                    removeAlpha(cols[i]);
                    WE[c][nWeekend].push_back(atoi(cols[i])); //TODO -1 adicionado
                }
            }
            while (1);
        }
        else if (isParam("nw_Max")) {
            assert(nContracts > 0);
            assert(nDays > 0);
            nw_Max.resize(nContracts, nDays);
            do {
                ncols = readCols(f);


                if ((ncols == 1) && (cols[0][0] == ';'))
                    break;

                removeGarbage(cols[0]);
                removeGarbage(cols[1]);

                int c = atoi(cols[0] + 1);
                int m = atoi(cols[1]);

                nw_Max[c] = m;
            }
            while (1);
        }
        else if (isParam("nw_Min")) {
            assert(nContracts > 0);
            assert(nDays > 0);
            nw_Min.resize(nContracts, 0);
            do {
                ncols = readCols(f);


                if ((ncols == 1) && (cols[0][0] == ';'))
                    break;

                removeGarbage(cols[0]);
                removeGarbage(cols[1]);

                int c = atoi(cols[0] + 1);
                int m = atoi(cols[1]);

                nw_Min[c] = m;
            }
            while (1);
        }
        else if (isParam("wna_Max")) {
            assert(nContracts > 0);
            assert(nDays > 0);
            wna_Max.resize(nContracts, 1);
            do {
                ncols = readCols(f);


                if ((ncols == 1) && (cols[0][0] == ';'))
                    break;

                removeGarbage(cols[0]);
                removeGarbage(cols[1]);

                int c = atoi(cols[0] + 1);
                int m = atoi(cols[1]);

                wna_Max[c] = m;
            }
            while (1);
        }
        else if (isParam("wna_Min")) {
            assert(nContracts > 0);
            assert(nDays > 0);
            wna_Min.resize(nContracts, 1);
            do {
                ncols = readCols(f);

                if ((ncols == 1) && (cols[0][0] == ';'))
                    break;

                removeGarbage(cols[0]);
                removeGarbage(cols[1]);

                int c = atoi(cols[0] + 1);
                int m = atoi(cols[1]);

                wna_Min[c] = m;
            }
            while (1);
        }
        else if (isParam("wnw_Max")) {
            assert(nContracts > 0);
            assert(nDays > 0);
            wnw_Max.resize(nContracts, 1);
            do {
                ncols = readCols(f);


                if ((ncols == 1) && (cols[0][0] == ';'))
                    break;

                removeGarbage(cols[0]);
                removeGarbage(cols[1]);

                int c = atoi(cols[0] + 1);
                int m = atoi(cols[1]);

                wnw_Max[c] = m;
            }
            while (1);
        }
        else if (isParam("wnw_Min")) {
            assert(nContracts > 0);
            assert(nDays > 0);
            wnw_Min.resize(nContracts, 1);
            do {
                ncols = readCols(f);


                if ((ncols == 1) && (cols[0][0] == ';'))
                    break;

                removeGarbage(cols[0]);
                removeGarbage(cols[1]);

                int c = atoi(cols[0] + 1);
                int m = atoi(cols[1]);

                wnw_Min[c] = m;
            }
            while (1);
        }
        else if (isParam("nf_Min")) {
            assert(nContracts > 0);
            assert(nDays > 0);
            nf_Min.resize(nContracts, 1);
            do {
                ncols = readCols(f);


                if ((ncols == 1) && (cols[0][0] == ';'))
                    break;

                removeGarbage(cols[0]);
                removeGarbage(cols[1]);

                int c = atoi(cols[0] + 1);
                int m = atoi(cols[1]);

                nf_Min[c] = m;
            }
            while (1);
        }
        else if (isParam("nf_Max")) {
            assert(nContracts > 0);
            assert(nDays > 0);
            nf_Max.resize(nContracts, 1);
            do {
                ncols = readCols(f);


                if ((ncols == 1) && (cols[0][0] == ';'))
                    break;

                removeGarbage(cols[0]);
                removeGarbage(cols[1]);

                int c = atoi(cols[0] + 1);
                int m = atoi(cols[1]);

                nf_Max[c] = m;
            }
            while (1);
        }
        else if (isParam("wnf_Min")) {
            assert(nContracts > 0);
            assert(nDays > 0);
            wnf_Min.resize(nContracts, 1);
            do {
                ncols = readCols(f);


                if ((ncols == 1) && (cols[0][0] == ';'))
                    break;

                removeGarbage(cols[0]);
                removeGarbage(cols[1]);

                int c = atoi(cols[0] + 1);
                int m = atoi(cols[1]);

                wnf_Min[c] = m;
            }
            while (1);
        }
        else if (isParam("wnf_Max")) {
            assert(nContracts > 0);
            assert(nDays > 0);
            wnf_Max.resize(nContracts, 1);
            do {
                ncols = readCols(f);


                if ((ncols == 1) && (cols[0][0] == ';'))
                    break;

                removeGarbage(cols[0]);
                removeGarbage(cols[1]);

                int c = atoi(cols[0] + 1);
                int m = atoi(cols[1]);

                wnf_Max[c] = m;
            }
            while (1);
        }
        else if (isParam("nwec_Min")) {
            assert(nContracts > 0);
            assert(nDays > 0);
            nwec_Min.resize(nContracts, 0);
            do {
                ncols = readCols(f);


                if ((ncols == 1) && (cols[0][0] == ';'))
                    break;

                removeGarbage(cols[0]);
                removeGarbage(cols[1]);

                int c = atoi(cols[0] + 1);
                int m = atoi(cols[1]);

                nwec_Min[c] = m;
            }
            while (1);
        }
        else if (isParam("nwec_Max")) {
            assert(nContracts > 0);
            assert(nDays > 0);
            nwec_Max.resize(nContracts, 0);
            do {
                ncols = readCols(f);


                if ((ncols == 1) && (cols[0][0] == ';'))
                    break;

                removeGarbage(cols[0]);
                removeGarbage(cols[1]);

                int c = atoi(cols[0] + 1);
                int m = atoi(cols[1]);

                nwec_Max[c] = m;
            }
            while (1);
        }
        else if (isParam("wnwec_Min")) {
            assert(nContracts > 0);
            assert(nDays > 0);
            wnwec_Min.resize(nContracts, 0);
            do {
                ncols = readCols(f);


                if ((ncols == 1) && (cols[0][0] == ';'))
                    break;

                removeGarbage(cols[0]);
                removeGarbage(cols[1]);

                int c = atoi(cols[0] + 1);
                int m = atoi(cols[1]);

                wnwec_Min[c] = m;
            }
            while (1);
        }
        else if (isParam("wnwec_Max")) {
            assert(nContracts > 0);
            assert(nDays > 0);
            wnwec_Max.resize(nContracts, 0);
            do {
                ncols = readCols(f);


                if ((ncols == 1) && (cols[0][0] == ';'))
                    break;

                removeGarbage(cols[0]);
                removeGarbage(cols[1]);

                int c = atoi(cols[0] + 1);
                int m = atoi(cols[1]);

                wnwec_Max[c] = m;
            }
            while (1);
        }
        else if (isParam("fw")) {
            assert(nContracts > 0);
            assert(nDays > 0);
            fw.resize(nContracts, 1);
            do {
                ncols = readCols(f);


                if ((ncols == 1) && (cols[0][0] == ';'))
                    break;

                removeGarbage(cols[0]);
                removeGarbage(cols[1]);

                int c = atoi(cols[0] + 1);
                int m = atoi(cols[1]);

                fw[c] = m;
            }
            while (1);
        }
        else if (isParam("wfw")) {
            assert(nContracts > 0);
            assert(nDays > 0);
            wfw.resize(nContracts, 1);
            do {
                ncols = readCols(f);


                if ((ncols == 1) && (cols[0][0] == ';'))
                    break;

                removeGarbage(cols[0]);
                removeGarbage(cols[1]);

                int c = atoi(cols[0] + 1);
                int m = atoi(cols[1]);

                wfw[c] = m;
            }
            while (1);
        }
        else if (isParam("ssw")) {
            assert(nContracts > 0);
            assert(nDays > 0);
            ssw.resize(nContracts, 1);
            do {
                ncols = readCols(f);


                if ((ncols == 1) && (cols[0][0] == ';'))
                    break;

                removeGarbage(cols[0]);
                removeGarbage(cols[1]);

                int c = atoi(cols[0] + 1);
                int m = atoi(cols[1]);

                ssw[c] = m;
            }
            while (1);
        }
        else if (isParam("wssw")) {
            assert(nContracts > 0);
            assert(nDays > 0);
            wssw.resize(nContracts, 1);
            do {
                ncols = readCols(f);


                if ((ncols == 1) && (cols[0][0] == ';'))
                    break;

                removeGarbage(cols[0]);
                removeGarbage(cols[1]);

                int c = atoi(cols[0] + 1);
                int m = atoi(cols[1]);

                wssw[c] = m;
            }
            while (1);
        }
    }

    fclose(f);
}

/* returns -1 if file ended.
   if not returns number of columns (strings or numbers which are separated by spaces) */
int Instance::readCols(FILE *f)
{
    char line[1024];

    if (!fgets(line, 1024, f))
        return -1;

    return splitString(cols, line, ' ', 1024, 1024, 1);
}

void Instance::removeGarbage(char *str)
{
    char *p;

    // removing spaces at the end */
    p = str + strlen(str) - 1;
    while (p >= str) {
        if ((*p == ' ') || (*p == '\r') || (*p == '\n'))
            *p = '\0';
        else
            break;

        --p;
    }


    p = str;

    while (*p != '\0') {
        if (*p == '\n') {
            *p = '\0';
            return;
        }
        if (!isalnum(*p))
            *p = ' ';
        ++p;
    }

}

bool Instance::isParam(const char *param)
{
    if (strcasecmp(cols[0], "param"))
        return false;

    if (strcasecmp(cols[1], param))
        return false;

    return true;
}

bool Instance::isSet(const char *setName)
{
    if (strcasecmp(cols[0], "set"))
        return false;

    if (strcasecmp(cols[1], setName))
        return false;

    return true;
}

bool Instance::isSetStart(const char *setName)
{
    if (strcasecmp(cols[0], "set"))
        return false;

    if (!strcasestr(cols[1], setName))
        return false;

    return true;
}

bool Instance::shiftCUsed(char c)
{
    for (int i = 0; (i < (int) shiftC.size()); ++i)
        if (c == shiftC[i])
            return true;

    return false;
}

void Instance::removeAlpha(char *str)
{
    char *p = str;
    while (*p != '\0') {
        if (*p == '\n') {
            *p = '\0';
            return;
        }
        if (isalpha(*p))
            *p = ' ';
        ++p;
    }
}

int Instance::getShiftIndex(const string str) const
{
    map<string, int>::const_iterator mIt = shiftIndex.find(str);
    if (mIt == shiftIndex.end()) {
        printf("Error on shift index\n");
        exit(EXIT_FAILURE);
    }

    return mIt->second;
}

time_t Instance::readDateFromTXT(const char *fileName)
{
    FILE *f = fopen(fileName, "r");
    if (!f) {
        fprintf(stderr, "Could not open file %s", fileName);
        exit(EXIT_FAILURE);
    }

    char line[256];
    while (fgets(line, 256, f)) {
        if (strstr(line, "SchedulingPeriod")) {
            fgets(line, 256, f);
            assert(strlen(line));

            char *s = strstr(line, ">");
            assert(s);
            s += 1;
            char strDate[256];
            strcpy(strDate, s);
            *(strstr(strDate, "<")) = '\0';
            printf("date: %s\n", strDate);

            char strYear[256];
            strcpy(strYear, strDate);
            *(strstr(strYear, "-")) = '\0';

            char strMonth[256];
            strcpy(strMonth, strstr(strDate, "-") + 1);
            *(strstr(strMonth, "-")) = '\0';

            char strDay[256];
            strcpy(strDay, strstr(strstr(strDate, "-") + 1, "-") + 1);

            //printf("day %s month %s year %s\n", strDay, strMonth, strYear);

            int day = atoi(strDay);
            int month = atoi(strMonth);
            int year = atoi(strYear);

            struct tm btime;

            btime.tm_year = year - 1900;
            btime.tm_mon = month - 1;
            btime.tm_mday = day;

            fclose(f);

            return mktime(&btime);
        }
    }

    fprintf(stderr, "Could not read date from instance file.\n");
    exit(EXIT_FAILURE);

    // so that the compiler will not complain
    return mktime(NULL);


    fclose(f);
}

time_t Instance::dateForDay(const time_t *date, const int nDay)
{
    struct tm btime;
    btime = *localtime(date);
    btime.tm_mday += nDay;
    return mktime(&btime);
}

int splitString(char **columns, const char *str, const char delimiter,
                const int maxColumns, const int columnSize, const char multDel)
{
    int sizeColumn;
    int ncolumn = 0;
    const char *send = str + strlen(str);
    const char *s = str;
    if (str[0] == '\0')
        return 0;
    const char *ns = s;
    PROCESS_COLUMN:
    if (ncolumn + 1 == maxColumns)
        return ncolumn;

    // finds the next delimiter
    FIND_DELIMITER:
    if (ns == send)
        goto FOUND_COLUMN;
    if (*ns != delimiter) {
        ns++;
        goto FIND_DELIMITER;
    }
    FOUND_COLUMN:
    sizeColumn = ns - s;
    if ((!multDel) || (sizeColumn > 0)) {
        if (sizeColumn)
            memcpy(columns[ncolumn], s, sizeColumn);
        columns[ncolumn][sizeColumn] = '\0';
        ncolumn++;
    }
    if (ns == send)
        return ncolumn;
    ++ns;
    s = ns;
    if (ns != send)
        goto PROCESS_COLUMN;

    return ncolumn;
}

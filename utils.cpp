//
// Created by 12859 on 2021/8/25.
//

#include "utils.h"
#include <cstring>
#include <fstream>

using namespace std;

void helpmsg(const char* pgm)
{
    cout << "Suggested arguments:" << endl;
    cout << "> " << pgm << " ";
    cout << "-k 5 -d 4 -f U400K4.dat -m CSA -h 5" << endl;
    cout << "explanations:" << endl;

    cout << "-d: data dimensionality, e.g. 2" << endl;
    cout << "-f: a tab or space delimited file, format:" << endl;
    cout << "    id xmin ymin [zmin] xmax ymax [zmax]" << endl;
    cout << "-k: top-k" << endl;
    cout << "-m: method, e.g., CSA, CSA+, MDA, HMDA" << endl;
    cout << "-h: the height of quad-tree, 1 to ?" << endl;
}

const char* read(int a_argc, const char** a_argv, const char* a_param, const char* a_def)
{
    for (int i = 0; i < a_argc; i++)
    {
        if (strcmp(a_argv[i], a_param) == 0)
        {
            if (i + 1 == a_argc)
                return "";
            else
                return a_argv[i + 1];
        }
    }
    return "";
}

void read_options(const char* datafile){
    fstream fpdata;
    fpdata.open(datafile, ios::in);

    while(!fpdata.eof()){

        if(fpdata.eof()){
            break;
        }
    }
    fpdata.close();
}
#include <iostream>
#include "utils.h"
#include "cell.h"
#include "CSA.h"
#include <cstring>
using namespace std;

int main(const int argc, const char** argv) {
    cout.precision(6);
    cout << "hmp" << endl;
    clock_t at, ad;

    cout << "Parse Parameters" << endl;
    if (argc == 1)
    {
        helpmsg(argv[0]);
        return -1;
    }

    const int k = atoi(read(argc, argv, "-k", ""));
    int dim = atoi(read(argc, argv, "-d", ""));
    const char* methodName = read(argc, argv, "-m", "");
    const char* datafile = read(argc, argv, "-f", "");  // option file name
    int h = atoi(read(argc, argv, "-h", ""));

    vector<vector<double>> P0=read_options(datafile, dim);

    vector<int> kskyband;
    kskyband_nortree(kskyband, P0, k);

    vector<vector<double>> P;
    P.reserve(kskyband.size());
    for (int i:kskyband) {
        P.push_back(P0[i]);
    }

    vector<double> b;
    b.reserve(2*dim);
    for (int j = 0; j < dim; ++j) {
        b.push_back(0.0);
        b.push_back(1.0);
    }
    string m=string(methodName);
    cell *root_ptr=new cell(b, 0, h, P.size(), k, m);
    if (strcmp(methodName, "CSA") == 0) {
        cout<<"CSA begin"<<endl;
        CSA(*root_ptr, P);
        cout<<"CSA end"<<endl;
    }


    return 0;
}

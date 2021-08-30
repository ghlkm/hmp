#include <iostream>
#include "utils.h"
#include "cell.h"
#include "hmp.h"
#include <cstring>
#include <chrono>

using namespace std;

vector<int> cell_debug;


int m2m(const char* s){
    int ret=1;
    if(strcmp(s, "CSA") == 0){
        return mCSA;
    }else if(strcmp(s, "CSA+") == 0){
        return mCSAp;
    }else if(strcmp(s, "MDA") == 0){
        return mMDA;
    }else if(strcmp(s, "MDA+") == 0){
        return mMDAp;
    }else{
        return ret;
    }
}

int main(const int argc, const char** argv) {
    cout.precision(6);
    cout << "hmp" << endl;

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
    vector<vector<int>> w;
    string s=string(datafile);
    s+=".kskyband";
//    kskyband_write(P0, k, s, w);
//    exit(0);

    vector<vector<int>> r;
    kskyband_read(s, r);
    for (int ik=0;ik<k;++ik) {
        for (int &i:r[ik]) {
            kskyband.push_back(i);
        }
    }
//    kskyband=vector<int>(kskyband.rbegin(), kskyband.rend());
//    assert(w.size()==r.size());
//    for (int l = 0; l < w.size(); ++l) {
//        assert(w[l].size()==r[l].size());
//        for (int i = 0; i < w[l].size(); ++i) {
//            assert(w[l][i]==r[l][i]);
//        }
//    }
//    kskyband_nortree(kskyband, P0, k);

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
    int method=m2m(methodName);
    cell_debug=vector<int>(h+1, 0);
    cout<<"k-ksyband size: "<<P.size()<<endl;
    cout<<sizeof(cell)<<endl;
    cell *root_ptr=new cell(b, 0, h, P.size(), k, method);
//    exit(0);
    auto ab = chrono::steady_clock::now();
    if (method==mCSA) {
        cout<<"CSA begin"<<endl;
        CSA(*root_ptr, P);
        cout<<"CSA end"<<endl;
    }else if(method==mCSAp){
        cout<<"CSA+ begin"<<endl;
        CSAp(*root_ptr, P);
        cout<<"CSA+ end"<<endl;
    }else if(method==mMDA){
        cout<<"MDA begin"<<endl;
        MDA(*root_ptr, P);
        cout<<"MDA end"<<endl;
    }else if(method==mMDAp){
        cout<<"MDA+ begin"<<endl;
        MDAp(*root_ptr, P);
        cout<<"MDA+ end"<<endl;
    }
    cout<<cell_debug<<endl;

    auto ae = chrono::steady_clock::now();
    chrono::duration<double> elapsed_seconds= ae-ab;
    cout << "Total time cost: " << elapsed_seconds.count() << endl;

    return 0;
}

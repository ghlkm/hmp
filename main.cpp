#include <iostream>

#include "utils.h"
#include "cell.h"
#include "hmp.h"
#include <cstring>
#include <chrono>
#include <fstream>

using namespace std;

vector<int> cell_debug;

vector<int> vt_debug;

extern long rsky_c; // CSA, CSA+, leaf
extern long dmc_c;  // CSA, CSA+, all
extern long rdo_g_c; // CSA+, dominate number
extern long s_rsky_c; // MDA, all
extern long s_rsky_p_c; // MDA+, max k
extern vector<long> dmc_p_c;
extern vector<long> dg_p_c;
extern long rtest_c;
extern long score_size; // baseline


int m2m(const char* s){
    int ret=1;
    if(strcmp(s, "BS") == 0){
        return mBASELINE;
    }else if(strcmp(s, "UTK") == 0){
        return mUTK;
    }else if(strcmp(s, "CSA") == 0){
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
    int method=m2m(methodName);

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
    if(method!=mUTK){
        P0.clear();
        vector<vector<double>> EMPTY;
        P0.swap(EMPTY);
    }

    vector<double> b;
    b.reserve(2*dim);
    for (int j = 0; j < dim; ++j) {
        b.push_back(0.0);
        b.push_back(1.0);
    }
    cell_debug=vector<int>(h+1, 0);
    vt_debug=vector<int>(h+1, 0);
    cout<<"k-ksyband size: "<<P.size()<<endl;
    cout<<sizeof(cell)<<endl;
    cell *root_ptr=nullptr;
    if(method!=mUTK){
        root_ptr=new cell(b, 0, h, P0.size(), k, method);
    }else{
        root_ptr=new cell(b, 0, h, P.size(), k, method);
    }
//    exit(0);
    auto ab = chrono::steady_clock::now();
    if(method==mBASELINE){
        cout<<"Baseline begin"<<endl;
        Baseline(*root_ptr, P);
        cout<<"Baseline end"<<endl;
    }else if(method==mUTK){
        cout<<"Baseline utk begin"<<endl;
        Baseline2(*root_ptr, P0);
        cout<<"Baseline utk end"<<endl;
    }else if (method==mCSA) {
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
    cout<<vt_debug<<endl;

    auto ae = chrono::steady_clock::now();
    chrono::duration<double> elapsed_seconds= ae-ab;
    cout << "Total time cost: " << elapsed_seconds.count() << endl;

    cal_mem(*root_ptr);

    cout << "rsky_c: " << rsky_c <<endl;
    cout << "dmc_c: " << dmc_c <<endl;
    cout << "rdo_g_c: " << rdo_g_c <<endl;
    cout << "s_rsky_c: " << s_rsky_c <<endl;
    cout << "s_rsky_p_c: " << s_rsky_p_c<<endl;
    cout << "score_size: " << score_size <<endl;
    cout << dmc_p_c << endl;
    cout << dg_p_c << endl;
    cout << rtest_c << endl;

    fstream log;
    string path=string(datafile);
    string df=path.substr(path.rfind('/'), path.size());
    df=df.substr(0, df.rfind('.'));
    string filename=string("~/data/log/")+df+
            string("_k")+to_string(k)+
            string ("_h")+to_string(h)+
            string (methodName)+string (".log");
    s+=".log";
    log.open(filename, ios::out);
    log << "Total time cost: " << elapsed_seconds.count() << endl;
    log << "rsky_c: " << rsky_c <<endl;
    log << "dmc_c: " << dmc_c <<endl;
    log << "rdo_g_c: " << rdo_g_c <<endl;
    log << "s_rsky_c: " << s_rsky_c <<endl;
    log << "s_rsky_p_c: " << s_rsky_p_c<<endl;
    log << "score_size: " << score_size <<endl;
    log << dmc_p_c << endl;
    log << dg_p_c << endl;
    log << rtest_c << endl;
    log.close();
    return 0;
}

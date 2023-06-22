#include <iostream>

#include "utils.h"
#include "cell.h"
#include "hmp.h"
#include <cstring>
#include <chrono>
#include <fstream>

using namespace std;
#define RTREE_IO
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
unordered_map<long int, RtreeNode*> ramTree; // load Rtree to main-memory
Rtree* rtree;
std::size_t page_access=0;

int m2m(const char* s){
    int ret=1;
    if(strcmp(s, "BS") == 0){
        return mBASELINE;
    }else if(strcmp(s, "BL") == 0){
        return mUTK;
    }else if(strcmp(s, "CSA") == 0) {
        return mCSA;
    }else if(strcmp(s, "MDA") == 0){
        return mMDA;
    }else if(strcmp(s, "MDA+") == 0){
        return mMDAp;
    }else if(strcmp(s, "TOPK")==0){
        return mTOPK1;
    }else if(strcmp(s, "TOPK+")==0){
        return mTOPK2;
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
    vector<RtreeNodeEntry*> p;
    vector<int> kskyband;
#ifndef RTREE_IO
    vector<vector<double>> P0=read_options(datafile, dim);
#endif
#ifdef RTREE_IO
    vector<vector<double>> P0;
//    string indexf=string("./")+string(datafile)+string("index.txt");
//    const char* indexfile = indexf.c_str();
    const char* indexfile="index4.txt";
    fstream fpdata;
    fpdata.open(datafile, ios::in);
    int id;
    vector<vector<float>> PointSet;
    PointSet.emplace_back(dim*2);
    P0.emplace_back(dim);
    int objCnt=0;
    while (!fpdata.eof()){
        fpdata >> id;
        if (fpdata.eof())
            break;
        P0.emplace_back(dim);
        PointSet.emplace_back(2*dim);
        for (int d = 0; d < dim; d++){
            fpdata >> PointSet.back()[d];
        }
        for (int d = 0; d < dim; d++){
            fpdata >> PointSet.back()[d+dim];
        }
        for (int d = 0; d < dim; d++){
            P0.back()[d]=(PointSet.back()[d]+PointSet.back()[d+dim])/2.0;
        }

        Hypercube hc(dim, &PointSet[objCnt + 1][0], &PointSet[objCnt + 1][dim]);

        p.push_back(new RtreeNodeEntry(id, hc));
        objCnt++;
        //log information
        if (objCnt % 1000 == 0)
            cout << ".";
        if (objCnt % 10000 == 0)
            cout << objCnt << " objects loaded" << endl;
    }
    fpdata.close();
    cout << "Total number of options: " << objCnt << endl;
    // build rtree
    cout << "Bulkloading R-tree..." << endl;
    const int maxChild = (PAGESIZE - RtreeNode::size()) / RtreeNodeEntry::size(dim);
    FileMemory mem(PAGESIZE, indexfile, RtreeNodeEntry::fromMem, true);
    cout << "[Rtree allocate mem done]" << endl;
    rtree = TGS::bulkload(mem, dim, maxChild, maxChild, (int)maxChild*0.3, (int)maxChild*0.3, p.data(), objCnt, false);
    cout << "[Rtree build done]" << endl;
    // in-memory rtree
    cout << "cache R-tree into memory" << endl;
    rtreeRAM(*rtree, ramTree);
    // aggregate rtree
    aggregateRecords(*rtree);
    cout << "[Aggregate Rtree done]" << endl;
//    void kskyband_rtree(const int dimen, Rtree& a_rtree, vector<long int>& kskyband, float* PG[], const int k)
//    auto kskyb = chrono::steady_clock::now();
//    kskyband_rtree(dim, *rtree, kskyband, PointSet, k); //
//    auto kskye = chrono::steady_clock::now();
//    chrono::duration<double> kskyt= kskye-kskyb;
//    cout << "kskyband time cost: " << kskyt.count() << endl;
//    cout<<"kskyband"<<kskyband.size()<<endl;
//    cout<<"page_access: "<<page_access<<endl;
//    exit(1);
#endif
//    vector<vector<int>> w;
    string s=string(datafile);
    s+=".kskyband";


//    kskyband_write(P0, k, s, w);
//    exit(0);

    vector<vector<int>> r;
    kskyband_read(s, r);
    for (int ik=0;ik<k;++ik) {
        if(ik>=r.size()){
            break;
        }
        for (int &i:r[ik]) {
            kskyband.push_back(i+1);
//            cout<<i<<": "<<P0[i]<<endl;
//            exit(0);
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
        root_ptr=new cell(b, 0, h, P.size(), k, method);
    }else{
        root_ptr=new cell(b, 0, h, P0.size(), k, method);
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
    }else if(method==mMDA){
        cout<<"MDA begin"<<endl;
        MDA(*root_ptr, P);
        cout<<"MDA end"<<endl;
    }else if(method==mMDAp){
        cout<<"MDA+ begin"<<endl;
        MDAp(*root_ptr, P);
        cout<<"MDA+ end"<<endl;
    }else if(method==mTOPK1){
        // usual rtree top-k
//        unordered_map<long int, RtreeNode *> empty_now; // not empty later
        topk_multi(*root_ptr, k,  P, 100, nullptr, false);
    }else if(method==mTOPK2){
        // usual rtree top-k and initialize bound as MaxMinK
        // 1. first find weight vector in which cell
        // 2. then find the cell's MaxMin_k
        // 3. find topk using MaxMin_k as a bound to prune
//        unordered_map<long int, RtreeNode *> empty_now; // not empty later
        topk_multi(*root_ptr, k,  P, 100, nullptr, true);
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

#ifdef DEBUG
    fstream hlog;
    string path=string(datafile);
    string df=path.substr(path.rfind('/'), path.size());
    string folder=path.substr(0, path.rfind('/'));
    df=df.substr(0, df.rfind('.'));
    string filename=folder+"/../"+string("/log/")+df+
            string("_k")+to_string(k)+
            string ("_h")+to_string(h)+
            string (methodName)+string (".heat");
    s+=".log";
    hlog.open(filename, ios::out);
    if (hlog.is_open()){
        cout<<filename<<endl;
        std::cout << "Output operation successfully performed\n";
    }else{
        cout<<filename<<endl;
        std::cout << "Error opening file";
        exit(-1);
    }
    for (auto &i: wheat) {
        for (double j:i) {
            hlog<<j<<" ";
        }
        hlog<<"\n";
    }
    hlog.close();
    cout<<"wheat: "<<wheat.size()<<endl;
#endif
//    fstream log;
//    string path=string(datafile);
//    string df=path.substr(path.rfind('/'), path.size());
//    string folder=path.substr(0, path.rfind('/'));
//    df=df.substr(0, df.rfind('.'));
//    string filename=folder+"/../"+string("/log/")+df+
//            string("_k")+to_string(k)+
//            string ("_h")+to_string(h)+
//            string (methodName)+string (".log");
//    s+=".log";
//    log.open(filename, ios::out);
//    if (log.is_open())
//    {
//        cout<<filename<<endl;
//        std::cout << "Output operation successfully performed\n";
//    }
//    else
//    {
//        cout<<filename<<endl;
//        std::cout << "Error opening file";
//        exit(-1);
//    }
//    log << "Total time cost: " << elapsed_seconds.count() << endl; // TODO add memory usage
//    log << "rsky_c: " << rsky_c <<endl;
//    log << "dmc_c: " << dmc_c <<endl;
//    log << "rdo_g_c: " << rdo_g_c <<endl;
//    log << "s_rsky_c: " << s_rsky_c <<endl;
//    log << "s_rsky_p_c: " << s_rsky_p_c<<endl;
//    log << "score_size: " << score_size <<endl;
//    log << dmc_p_c << endl;
//    log << dg_p_c << endl;
//    log << rtest_c << endl;
//    log.close();
    return 0;
}

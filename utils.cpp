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

inline bool v1_dominate_v2(vector<double>& v1, vector<double>& v2){
    assert(v1.size()==v2.size());
    for (int i = 0; i < v1.size(); ++i) {
        if(v1[i] < v2[i]){
            return false;
        }
    }
    return true;
}

vector<vector<double>> read_options(const char* datafile, int dim){
    fstream fpdata;
    fpdata.open(datafile, ios::in);
    int id;
    double tmp;
    vector<vector<double>> P0;
    while(!fpdata.eof()){
        fpdata >> id;
        if(fpdata.eof()){
            break;
        }
        vector<double> option(dim);
        for (int i = 0; i < dim; ++i) {
            fpdata>>option[i];
        }
        for (int j = 0; j <dim; ++j) {
            fpdata>>tmp;
            option[j]+=tmp;
            option[j]/=2.0;
        }
        P0.push_back(option);
        if(P0.size()%1000==0){
            cout << P0.size() << " objects loaded" << endl;
        }
    }
    fpdata.close();
    cout << "Total number of options: " << P0.size() << endl;
    return P0;
}

void kskyband_nortree(vector<int> &ret, vector<vector<double>> &data,int k) {
    /*
     * the k-skyband contains thoes records that are dominated by fewer than k others
     */
    vector<int> do_cnt(data.size(), 0);
    vector<int> k_cnt(k, 0);
    for (auto i = 0; i < data.size(); ++i) {
        for (auto j = i + 1; j < data.size(); ++j) {
            if (do_cnt[i] >= k) {
                break;
            }
            if (v1_dominate_v2(data[i], data[j])) {
                ++do_cnt[j];
            } else if (v1_dominate_v2(data[j], data[i])) {
                ++do_cnt[i];
            }
        }
        if (do_cnt[i] < k) {
            ret.push_back(i);
            k_cnt[do_cnt[i]]+=1;
        }
    }
    cout<<"for k-skyband, # of options:"<<endl;
    for (int l = 0; l <k ; ++l) {
        cout<<"@k="<<l+1<<":num="<<k_cnt[l]<<endl;
    }
}

int r_dominate(std::vector<std::vector<double>> &vs, std::vector<double> &p1, std::vector<double> &p2){
    int ret=NONDO;
    double p1s, p2s;
    for (auto &v: vs) {
        p1s=p1*v;
        p2s=p2*v;
        if(p1s>p2s){
            if(ret==DOMINATED){
                return NONDO;
            }else{
                ret=DOMINATE;
            }
        }else if(p1s<p2s){
            if(ret==DOMINATE){
                return NONDO;
            }else{
                ret=DOMINATED;
            }
        }
    }
    return ret;
}
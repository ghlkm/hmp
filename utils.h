

#ifndef HEATMAP_UTILS_H
#define HEATMAP_UTILS_H

#include <vector>
#include <rtree.h>
#include <hypercube.h>
#include <filemem.h>
#include <tgs.h>
#include "vector_operator.h"
#include <rentry.h>
#include <rnode.h>

extern long rtest_c;

#define DOMINATE 1
#define DOMINATED -1
#define NONDO 0

void helpmsg(const char* pgm);

const char* read(int a_argc, const char** a_argv, const char* a_param, const char* a_def);

vector<vector<double>> read_options(const char* datafile, int dim);

inline bool v1_dominate_v2(vector<double>& v1, vector<double>& v2);

void kskyband_nortree(vector<int> &ret, vector<vector<double>> &data,int k);

bool r_dominate(std::vector<double> &s1, std::vector<double> &s2);

void kskyband_write(vector<vector<double>> &data, int k, string &filename, vector<vector<int>> &ret);

void kskyband_read (const string &filename, vector<vector<int>> &ret);

int closeTo(double d);

double sum(vector<double> &v);

void rtreeRAM(const Rtree& rt, unordered_map<long int, RtreeNode*>& ramTree);

int countRecords(Rtree& a_rtree, int pageID);

void aggregateRecords(Rtree& a_rtree);

template<typename VVF>
void build_rtree(Rtree* &rtree_rt, unordered_map<long int, RtreeNode*>& ramTree,VVF &data){
    if(data.empty()){
        return;
    }
    int dim=data[0].size();
    RtreeNodeEntry** p = new RtreeNodeEntry*[data.size()];
    for (int id=0;id<data.size();++id)
    {
        float* cl = new float[dim];
        float* cu = new float[dim];
        for (int i = 0; i <dim ; ++i) {
            cl[i]=data[id][i]-SIDELEN;
            cu[i]=data[id][i]+SIDELEN;
        }
        Hypercube hc(dim, cl, cu);
        p[id] = new RtreeNodeEntry(id, hc);
    }
    // build rtree
    const int maxChild = (PAGESIZE - RtreeNode::size()) / RtreeNodeEntry::size(dim);
    //FileMemory mem(PAGESIZE, "./result/index.txt", RtreeNodeEntry::fromMem, true);
    string indexf_name="index"+to_string(time(nullptr)%1000)+"txt"; // TODO specific by user
    FileMemory *mem=new FileMemory(PAGESIZE, indexf_name.c_str(), RtreeNodeEntry::fromMem, true);
    rtree_rt = TGS::bulkload(*mem, dim, maxChild, maxChild, (int)maxChild*0.3, (int)maxChild*0.3, p, data.size(), false);
    rtreeRAM(*rtree_rt, ramTree);
    aggregateRecords(*rtree_rt);
}


template<typename FLOAT>
FLOAT utk_orderScore(vector<FLOAT>& pivot, vector<FLOAT>& entry){
    FLOAT ret=0;
    for (int i = 0; i < pivot.size(); i++){
        ret += pivot[i] * entry[i];
    }
    return ret;
}

bool r_dominate(vector<vector<float>>& vs, vector<float> &v1, vector<float>& v2);


int r_dominate(std::vector<std::vector<double>> &vs, std::vector<double> &p1, std::vector<double> &p2);

template<typename FLOAT>
bool utk_countRegionDominator(int dimen, vector<FLOAT> &pt, vector<int>& rskyband, vector<vector<FLOAT>> &PG, vector<vector<FLOAT>>& vs, const int k)
{
    // pt, dim
    // PG, dim
    // vs, region vertexes, dim-1
    if(rskyband.size()<k){
        return true;
    }
    vector<FLOAT> record(dimen,0);
    int count = 0;
    for (int i : rskyband){
        rtest_c+=1;
        if(r_dominate(vs, PG[i], pt)) {
            count++;
            if (count >= k) {
                return false;
            }
        }
    }
    return count<k;
}

template<typename FLOAT>
long utk_rskyband(vector<vector<FLOAT>>& region_v, const int dimen, Rtree& a_rtree, vector<int>& rskyband,
                  vector<vector<FLOAT>> &PG, unordered_map<long int, RtreeNode *> &ramTree, int k=1){
    // region_v: dim
    // dimen:  dim
    rskyband.clear();
    vector<FLOAT> pivot(dimen, 0);
    for (auto &v: region_v) {
        for (int i = 0; i <v.size() ; ++i) {
            pivot[i]+=v[i];
        }
    }
    for (FLOAT & l : pivot) {
        l/=region_v.size();
    }
    RtreeNode* node;
    priority_queue<pair<FLOAT, int>> heap;
    int NegPageid;
    FLOAT maxscore;
    int pageID;
    FLOAT tmpScore;
    unordered_set<long int> dominators;
    heap.emplace(INFINITY, a_rtree.m_memory.m_rootPageID);
    vector<FLOAT> pt(dimen, 0);
    long ret=0;
    while (!heap.empty()){
        ret=heap.size()>ret?heap.size():ret;
        tmpScore = heap.top().first;
        pageID = heap.top().second;
        heap.pop();
        if (pageID >= MAXPAGEID){
            if (utk_countRegionDominator(dimen, PG[pageID - MAXPAGEID], rskyband, PG, region_v, k)){
                rskyband.push_back(pageID - MAXPAGEID);
            }
        }else{
            node = ramTree[pageID];
            if (node->isLeaf()){
                for (int i = 0; i < node->m_usedspace; i++){
                    if (utk_countRegionDominator(dimen, PG[node->m_entry[i]->m_id], rskyband, PG, region_v, k)){
                        maxscore = utk_orderScore(pivot, PG[node->m_entry[i]->m_id]);
                        heap.push(make_pair(maxscore, node->m_entry[i]->m_id + MAXPAGEID));
                    }
                }
            }
            else{
                for (int i = 0; i < node->m_usedspace; i++){
                    for (int j = 0; j < dimen; j++){
                        pt[j] = node->m_entry[i]->m_hc.getUpper()[j];
                    }
                    if (utk_countRegionDominator(dimen, pt, rskyband, PG, region_v, k)){
                        maxscore = utk_orderScore(pivot, pt);
                        heap.push(make_pair(maxscore, node->m_entry[i]->m_id));
                    }
                }
            }
        }
    }
    return ret;
}

#endif //HEATMAP_UTILS_H

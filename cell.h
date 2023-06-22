
#ifndef HEATMAP_CELL_H
#define HEATMAP_CELL_H

#include <vector>
#include <string>
#include "utils.h"
#include <unordered_set>
#include <cmath>
#include <algorithm>
#include <set>

#define mUTK -1
#define mBASELINE 0
#define mCSA  1
#define mCSAp 2
#define mMDA  3
#define mMDAp 4
#define mTOPK1 5
#define mTOPK2 6
#define DEBUG
#ifdef DEBUG
extern vector<vector<double>> wheat;
#endif
extern vector<int> cell_debug;

extern vector<int> vt_debug;

extern long s_rsky_p_c;
extern long rsky_c;
extern long rtest_c;
extern vector<long> dmc_p_c;
extern vector<long> dg_p_c;
extern long score_size; // baseline
extern int debug_vcnt;

class cell{
public:
    std::vector<cell*> children;
    std::vector<int> rkskyband;
    int dim;
    int k;
    std::vector<double> bounds;
    std::vector<std::vector<double>> vertexes;
    std::vector<int> dmc; // only used for CSA, "dominated" count array, dmc[i]=j means option i dominated by totally j options
    double uHeat; // MaxMin_k, the k_th largest minimum score
    double MaxKMaxUHeat;
    double AvgMaxUHeat;
    double rHeat; // No. of RSK
    double dHeat; // No. of r-dominate relationships
    int cur_l; // current level
    int tar_l; // target level
    int method; // method, which cell

//    vector<unordered_set<int>> rdo_graph; // only used for CSA+, a "dominate" graph
    vector<set<int>> rdo_graph; // only used for CSA+, a "dominate" graph

    vector<double> heap;  // only used for MDA, size k, a min heap

    // only used for MDA and MDA+, a superset of rkskyband, not empty if this cell is a leaf node
    // each element: <id, <lb, ub>>
    vector<pair<int, pair<double, double>>> s_rskyband;

    vector<int> mdap_s_rksyband;

    vector<double> rskyband_lb_MDA; // each element: <id, lb>

    cell()= default;

    cell(std::vector<double> &b, int cur_level, int tar_level, int card, int K, int m=mMDAp){
        /*
         * \para b, bounds
         * \para cur_level, this cell would be in which level of the quad-tree
         * \para tar_level, the quad-tree will build how many level use
         * \para card, cardinality of the option dataset, used for init dmc size and so on
         */
        dim=b.size()/2;
        bounds=b;
        uHeat=0.0;
        MaxKMaxUHeat=0.0;
        AvgMaxUHeat=0.0;
        rHeat=0.0;
        dHeat=0.0;
        cur_l=cur_level;
        tar_l=tar_level;
        method=m;
        this->k=K;
        // generate this cell's vertexes
        this->get_vertexes();

        if(m==mCSAp){
//            this->rdo_graph=vector<unordered_set<int>>(card);
            this->rdo_graph=vector<set<int>>(card);
        }
        if(m==mMDA || m==mMDAp){
            this->heap=vector<double>(k, 0.0);
        }else if(m==mCSA || m==mCSAp){
            dmc=std::vector<int>(card, 0);
        }
        if(m==mCSA || m==mCSAp || m==mMDA || cur_level==0){ // root
            if(cur_level<tar_level){
                this->get_children(card, m);
            }
        }
    }


    void get_children(int card, int m){
        double u=1.0;
        double dis=(bounds[1]-bounds[0])/2;  // child should be divided by 2
        double l=1.0-dis*dim;
        unsigned int child_p_num=(1<<dim);
        for(int iter=0; iter<child_p_num;++iter){
            vector<double> lv(dim);
            unsigned iter_cp=iter;
            for (int j = 0; j < dim; ++j) {
                if(iter_cp%2==0){
                    lv[j]=bounds[j*2];
                }else{
                    lv[j]=(bounds[j*2]+bounds[j*2+1])/2.0;
                }
                iter_cp=iter_cp>>1;
            }
            double slv=sum(lv);

            if(u>slv && slv>l){
                vector<double> child_b(dim*2);
                for (int i = 0; i <dim ; ++i) {
                    child_b[i*2]=lv[i];
                    child_b[i*2+1]=lv[i]+dis;
                }
                this->children.push_back(new cell(child_b, cur_l+1, tar_l, card, k, m));
            }
        }
        cell_debug[cur_l+1]+=this->children.size();
//        cout<<cur_l<<","<<this->children.size()<<endl;
    }

    cell* get_next_children(int &iter){ // TODO optimize the efficient of this:: might not be important
        double u=1.0;
        double dis=(bounds[1]-bounds[0])/2;  // child should be divided by 2
        double l=1.0-dis*dim;
        unsigned int child_p_num=(1<<dim);
        while(iter<child_p_num){
            vector<double> lv(dim);
            unsigned iter_cp=iter;
            for (int j = 0; j < dim; ++j) {
                if(iter_cp%2==0){
                    lv[j]=bounds[j*2];
                }else{
                    lv[j]=(bounds[j*2]+bounds[j*2+1])/2.0;
                }
                iter_cp=iter_cp>>1;
            }
            double slv=sum(lv);

            if(u>slv && slv>l){
                vector<double> child_b(dim*2);
                for (int i = 0; i <dim ; ++i) {
                    child_b[i*2]=lv[i];
                    child_b[i*2+1]=lv[i]+dis;
                }
                cell_debug[cur_l+1]+=1;
                return new cell(child_b, cur_l+1, tar_l, 0, k, method);
            }else{
                ++iter;
            }
        }
        return nullptr;
    }

    void get_vertexes(){
        double sum_lb=0.0;
        for (int d = 0; d <dim; ++d) {
            sum_lb+=this->bounds[d*2];
        }

        double dis=this->bounds[1]-this->bounds[0];
        int num_ub=closeTo((1.0-sum_lb)/dis);

        assert(num_ub<dim);
        // select num_ub from dim
        std::vector<bool> v(dim);
        std::fill(v.begin(), v.begin() + num_ub, true);
        do {
            vector<double> one_cb;
            one_cb.reserve(dim);
            for (int i = 0; i < dim; ++i) {
                if (v[i]) {
                    one_cb.push_back(this->bounds[2*i+1]);
                }else{
                    one_cb.push_back(this->bounds[2*i]);
                }
            }
            this->vertexes.push_back(one_cb);
        } while (std::prev_permutation(v.begin(), v.end())); // forward next permutation
        vt_debug[cur_l]+=this->vertexes.size();
    }

    void CSA_insert(int p1, int p2, std::vector<std::vector<double>> &P){
        if(!(dmc[p1]>=k || dmc[p2]>=k)){
            int state=r_dominate(this->vertexes, P[p1], P[p2]);
            rtest_c+=1;
            if(state==DOMINATE){
                this->recursively_update_dmc(p2);
            } else if(state==DOMINATED){
                this->recursively_update_dmc(p1);
            }else{
                for(auto &child: this->children){
                    child->CSA_insert(p1, p2, P);
                }
            }
        }else{
            dmc_p_c[cur_l]+=1;
        }
    }

    void recursively_update_dmc(int i, int add=1){
        if(this->dmc[i]<k){
            this->dmc[i]+=add;
            for(auto &child: this->children){
                child->recursively_update_dmc(i, add);
            }
        }
    }

    void CSAp_insert(int p1, int p2, std::vector<std::vector<double>> &P){
        if( dmc[p1]>=k || dmc[p2]>=k ){
            dmc_p_c[cur_l]+=1;
            return;
        }
        if(rdo_graph[p1].find(p2)!=rdo_graph[p1].end() || rdo_graph[p2].find(p1)!=rdo_graph[p2].end()){
            dg_p_c[cur_l]+=1;
            return;
        }
        int state=r_dominate(vertexes, P[p1], P[p2]);
        rtest_c+=1;
        if(state==DOMINATE){
            this->CSAp_update(p1, p2);
        }else if(state==DOMINATED){
            this->CSAp_update(p2, p1);
        }else{
            for(auto &child:children){
                child->CSAp_insert(p1, p2, P);
            }
        }
    }

    inline void CSAp_update(int i, int j){
//        recursively_update_dmc(j);
        CSAp_update_graph(i, j);
        for(int id:rdo_graph[j]){
            if(rdo_graph[i].find(id)==rdo_graph[i].end()){ // a new option
//                recursively_update_dmc(id);
                CSAp_update_graph(i, id);
            }
        }
    }

    void CSAp_update_graph(int p1, int p2){
        if(dmc[p2]>=k || rdo_graph[p1].find(p2)!=rdo_graph[p1].end()){
            return;
        }
        if(dmc[p1]>=k){
            recursively_update_dmc(p2, k-dmc[p2]);
            return;
        }
        rdo_graph[p1].insert(p2);
        dmc[p2]+=1;
        for(auto &child: children){
            child->CSAp_update_graph(p1, p2);
        }
    }

    void MDA_insert(int id, std::vector<std::vector<double>> &P){
        double theta=-heap.front();
        double lb, ub;
        get_lb_ub(P[id], lb, ub);
        if(ub>theta){
            if(lb>theta){
                // heap push
                heap.push_back(-lb);
                std::push_heap(heap.begin(), heap.end());

                // heap pop
                std::pop_heap(heap.begin(), heap.end());
                heap.pop_back();
            }
            if(this->isLeaf()){
                s_rskyband.emplace_back(id, pair<double, double>(lb, ub));
            }else{
                for(auto &child: children){
                    child->MDA_insert(id, P);
                }
            }
        }
    }

    void MDA_superSet2RKS(std::vector<std::vector<double>> &P){
        double theta=-heap.front();
        vector<pair<int, pair<double, double>>> tmp; // TODO add to pseudo code
        for (auto &i:s_rskyband) {
            double ub=i.second.second;
            if(ub>theta){
                tmp.push_back(i);
            }
        }
        sort(tmp.begin(), tmp.end(), [](auto &a, auto &b){
           return a.second.second > b.second.second;
        });
        assert(tmp.size()>=k);
        vector<vector<double>> scores;

        for (int i = 0; i <k ; ++i) {
            this->rkskyband.push_back(tmp[i].first);
            this->rskyband_lb_MDA.emplace_back(tmp[i].second.first);
            vector<double> s;
            s.reserve(vertexes.size());
            for (auto &v:vertexes) {
                s.push_back(v*P[tmp[i].first]);
            }
            scores.push_back(s);
        }
#ifdef DEBUG
        for (int i = 0; i < k; ++i) {
            for (int j = i+1; j < k; ++j) {
                if(r_dominate(scores[i], scores[j])){
                    dHeat+=1;
                }
            }
        }
#endif
        for (int i = k; i <tmp.size() ; ++i) {
            int r_dominate_count=0;
            vector<double> s;
            for (int j = 0; j < rkskyband.size(); ++j) {
                if(tmp[i].second.first < rskyband_lb_MDA[j]){
                    if(s.empty()){
                        s.reserve(vertexes.size());
                        for (auto &v:vertexes) {
                            s.push_back(v*P[tmp[i].first]);
                        }
                    }
                    rtest_c+=1;
                    if(r_dominate(scores[j],s)){
                        r_dominate_count+=1;
                        if(r_dominate_count>=k){
                            break;
                        }
                    }
                }
            }
            if(r_dominate_count<k){
                this->dHeat+=r_dominate_count;
                this->rkskyband.push_back(tmp[i].first);
                this->rskyband_lb_MDA.emplace_back(tmp[i].second.first);
                if(s.empty()){
                    s.reserve(vertexes.size());
                    for (auto &v:vertexes) {
                        s.push_back(v*P[tmp[i].first]);
                    }
                }
                scores.push_back(s);
            }
        }
    }

    inline bool isLeaf() const{
//        return children.empty();
        return cur_l>=tar_l;
    }

    inline void get_lb_ub(std::vector<double> &p, double &lb, double &ub){
        lb=INFINITY;
        ub=0.0;
        double tmp;
        for(auto &v:vertexes){
            tmp=p*v;
            if(tmp>ub){
                ub=tmp;
            }
            if(tmp<lb){
                lb=tmp;
            }
        }
    }

    inline void get_scores(std::vector<double> &p, std::vector<double> &ret){
        ret.clear();
        ret.reserve(vertexes.size());
        for(auto &v:vertexes){
            ret.push_back(p*v);
        }
    }

    inline void get_lu_scores(std::vector<int> &ps, std::vector<std::vector<double>> &P,
                              std::vector<double> &lb, std::vector<double> &ub){
        lb=std::vector<double>(ps.size());
        ub=std::vector<double>(ps.size());
        for (int i = 0; i <ps.size() ; ++i) {
            this->get_lb_ub(P[ps[i]], lb[i], ub[i]);
        }
    }

    inline void get_all_scores(std::vector<std::vector<double>> &P,
                           std::vector<std::vector<double>> &scores){
        scores=std::vector<std::vector<double>>(P.size());
        for (int i = 0; i <P.size() ; ++i) {
            this->get_scores(P[i], scores[i]);
        }
    }

    void MDAp_insert(vector<int> &parent_sRSK, vector<vector<double>> &P, long cnt=0){
        assert(!parent_sRSK.empty());
        vector<double> lbs, ubs;
        this->get_lu_scores(parent_sRSK, P, lbs, ubs);
        double theta;
        for (double &lb:lbs) {
            theta=-heap.front();
            if(lb>theta){
                // heap push
                heap.push_back(-lb);
                std::push_heap(heap.begin(), heap.end());

                // heap pop
                std::pop_heap(heap.begin(), heap.end());
                heap.pop_back();
            }
        }
        theta=-heap.front();
        for (int i = 0; i <parent_sRSK.size(); ++i) {
            if(ubs[i]>theta){
                this->mdap_s_rksyband.push_back(parent_sRSK[i]);
                if(this->isLeaf()){ // a leaf node
                    this->s_rskyband.emplace_back(parent_sRSK[i], pair<double, double>(lbs[i], ubs[i]));
                }
            }
        }
        cnt+=this->mdap_s_rksyband.size();
        if(!isLeaf()){
            unsigned int child_p_num=(1<<dim);
            for (int i = 0; i < child_p_num; ++i) {
                cell* child=this->get_next_children(i);
                if(child!= nullptr){
                    child->MDAp_insert(this->mdap_s_rksyband, P, cnt);
                    delete (child);
                }
            }
//            this->get_children(0, method);
//            for(auto &child:children){
//                child->MDAp_insert(this->mdap_s_rksyband, P);
//                delete (child);
//            }
        }else{
            this->MDA_superSet2RKS(P);
            rsky_c+=this->rkskyband.size();
            uHeat=theta;
            rHeat=rkskyband.size();
#ifdef DEBUG
            cal_mm_am_uHeats(P);
            vector<double> heats={uHeat, MaxKMaxUHeat, AvgMaxUHeat, rHeat, dHeat};
            wheat.emplace_back(heats);
#endif
        }
        if(cnt>s_rsky_p_c){
            s_rsky_p_c=cnt;
        }
    }

    void cal_mm_am_uHeats(vector<vector<double>> &data){
        vector<double> mina;
        vector<double> maxa;
        mina.reserve(rkskyband.size());
        maxa.reserve(rkskyband.size());
        for (auto i:rkskyband) {
            double maxv=0;
            double minv=1e9;
            for(auto &v: vertexes){
                double tmp_score=v*data[i];
                if(maxv<tmp_score){
                    maxv=tmp_score;
                }
                if(minv>tmp_score){
                    minv=tmp_score;
                }
            }
            mina.push_back(minv);
            maxa.push_back(maxv);
        }
        sort(maxa.begin(), maxa.end(), greater<>());
        sort(mina.begin(), mina.end(), greater<>());
        MaxKMaxUHeat=maxa[k-1];
        AvgMaxUHeat=sum(maxa)/maxa.size();
    }

    void Baseline_insert(vector<vector<double>> &P){
        if(isLeaf()){
            assert(!P.empty());
            vector<vector<double>> s;
            this->get_all_scores(P, s);
            for (int i = 0; i <P.size() ; ++i) {
                int do_cnt=0;
                for (int j = 0; j < P.size(); ++j) {
                    if(i==j){
                        continue;
                    }
                    rtest_c+=1;
                    if(r_dominate(s[j], s[i])){
                        ++do_cnt;
                        if(do_cnt>=k){
                            break;
                        }
                    }
                }
                if(do_cnt<k){
                    this->rkskyband.push_back(i);
                }
            }
            long tmp=0;
            for (auto &v:s) {
                tmp+=v.size();
            }
            score_size=tmp>score_size?tmp:score_size;
        }else{
            unsigned int child_p_num=(1<<dim);
            for (int i = 0; i < child_p_num; ++i) {
                cell* child=this->get_next_children(i);
                if(child!= nullptr){
                    child->Baseline_insert(P);
                    delete (child);
                }
            }
        }
    }

//    void Baseline2_insert(vector<vector<double>> &P, Rtree *rtree_rt, unordered_map<long int, RtreeNode *> &ramTree){
//        if(isLeaf()){
//            long tmp=utk_rskyband(this->vertexes, this->dim, *rtree_rt, this->rkskyband, P, ramTree, this->k)*4;
//            score_size=tmp>score_size?tmp:score_size;
//        }else{
//            unsigned int child_p_num=(1<<dim);
//            for (int i = 0; i < child_p_num; ++i) {
//                cell* child=this->get_next_children(i);
//                if(child!= nullptr){
//                    child->Baseline2_insert(P, rtree_rt, ramTree);
//                    delete (child);
//                }
//            }
//        }
//    }
    void Baseline2_insert(vector<vector<double>> &P, Rtree *rtree_rt, unordered_map<long int, RtreeNode *> &ramTree){
        if(isLeaf()){
            long tmp=utk_rskyband(this->vertexes, this->dim, *rtree_rt, this->rkskyband, P, ramTree, this->k)*4;
            score_size=(tmp*2+s_rsky_p_c)>score_size?(tmp*2+s_rsky_p_c):score_size;
        }else{
//            // TODO delete me later
//            if(score_size!=0){
//                return;
//            }
//            // TODO delete me later
            unsigned int child_p_num=(1<<dim);
            for (int i = 0; i < child_p_num; ++i) {
//                // TODO delete me later
//                if(score_size!=0){
//                    return;
//                }
//                // TODO delete me later
                cell* child=this->get_next_children(i);
                if(child!= nullptr){
                    s_rsky_p_c+=child->dim*2+child->vertexes.size()*child->dim;
                    child->Baseline2_insert(P, rtree_rt, ramTree);
                    s_rsky_p_c-=child->dim*2+child->vertexes.size()*child->dim;
                    delete (child);
                }
            }
        }
    }


};


void get_all_leaves(cell &node, std::vector<cell*>& ret);

void cal_mem(cell &node);
#endif //HEATMAP_CELL_H

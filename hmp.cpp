
#include "hmp.h"
#include <numeric>
#include <rtree.h>
#include <chrono>


void Baseline(cell &root, std::vector<std::vector<double>> &P){
    for(auto &child:root.children){
        child->Baseline_insert(P);
    }
}

void Baseline2(cell &root, std::vector<std::vector<double>> &P){
    unordered_map<long int, RtreeNode *> ramTree;
    Rtree *rtree_rt= nullptr;
    build_rtree(rtree_rt, ramTree, P);
    cout<<"rtree build finish"<<endl;
    s_rsky_p_c+=root.dim*2+root.vertexes.size()*root.dim;

    for(auto &child:root.children){
        s_rsky_p_c+=child->dim*2+child->vertexes.size()*child->dim;
        child->Baseline2_insert(P, rtree_rt, ramTree);
    }
//    score_size+=ramTree.size()*1024;
}

vector<int> topk_single(int k,  std::vector<std::vector<double>> &P, vector<double>&w, double uHeat,
                 Rtree *rtree_rt, unordered_map<long int, RtreeNode *> &ramTree){
    RtreeNode* node;
    priority_queue<pair<double, int>> heap;
    double bound=uHeat;
    heap.emplace(INFINITY, rtree_rt->m_memory.m_rootPageID);
    int dim=w.size();
    vector<double> tmp_v(dim, 0.0);
    vector<int> topk;
    double tmp_score;
    double pageID;
    while(!heap.empty() || topk.size()<k){
        tmp_score=heap.top().first;
        pageID=heap.top().second;
        heap.pop();
        if (pageID >= MAXPAGEID){ // an option
            topk.push_back(pageID-MAXPAGEID);
        }else{ // an inter node
            node = ramTree[pageID];
            if (node->isLeaf()){ // this inter node contains leaves
                for (int i = 0; i < node->m_usedspace; i++){
                    tmp_score=P[node->m_entry[i]->m_id] * w;
                    if(tmp_score>=bound){
                        // the usual topk bound is set to 0
                        // we here can also set it as MaxMin_k
                        heap.emplace(tmp_score, node->m_entry[i]->m_id + MAXPAGEID);
                    }
                }
            }
            else{
                for (int i = 0; i < node->m_usedspace; i++){
                    auto &tmp=node->m_entry[i]->m_hc.getUpper();
                    for (int j = 0; j < dim; j++){
                        tmp_v[j]=tmp[j];
                    }
                    tmp_score=tmp_v*w;
                    if (tmp_score>=bound){
                        heap.emplace(tmp_score, node->m_entry[i]->m_id);
                    }
                }
            }
        }
    }
    assert(topk.size()==k);
    return topk;
}

vector<vector<double>> gen_weight_vector(int k, int dim, int num){
    vector<vector<double>> ws(num);
    srand(0); // a fixed random seed makes stable results

    for (int i=0;i<num;i++){ // generate user preference uniformly under the constraint \sum v_i=1,
        vector<double> v(dim,0.0);
        float res=1.0;
        for (int d=0;d<dim-1;d++){
            v[d] = res*(1.0-pow((double)random()/RAND_MAX,  1.0/(dim-d-1)));
            res-=v[d];
        }
        v[dim-1]=res;
        ws.push_back(v);
    }
    return ws;
}

double getMaxMin_k(cell &root, int k, vector<double> &w, vector<vector<double>> &P){
    vector<int> idx=vector<int>(P.size());
    std::iota(idx.begin(), idx.end(), 0);

    int dim=w.size();
    vector<double> cur_bounds(dim*2, 0);
    for (int i = 0; i < dim; ++i) {
        cur_bounds[dim*2+1]=1.0;
    }
    cell c=root; // empty init
    int cur_level=0;
    vector<int> new_idx=vector<int>();
    while(c.isLeaf()){ // time complexity d*h
        cur_level++;
        // get next level cell
        for (int i = 0; i <dim ; ++i) {  // update bound, time complexity d
            double mid=(cur_bounds[i*2]+cur_bounds[i*2+1])/2;
            if(w[i]>mid){
                cur_bounds[i*2]=mid;
            }else{
                cur_bounds[i*2+1]=mid;
            }
        }
        c=cell(cur_bounds, cur_level, root.tar_l, 0, k, mMDAp);

        vector<double> lbs, ubs;
        c.get_lu_scores(idx, P, lbs, ubs);
        double theta;
        vector<double> &heap=c.heap;
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
        if(c.isLeaf()){
//            for (int i = 0; i <idx.size(); ++i) {
//                if(ubs[i]>theta){
//                    c.s_rskyband.emplace_back(idx[i], pair<double, double>(lbs[i], ubs[i]));
//                }
//            }
//            c.MDA_superSet2RKS(P);
//            break;
//            # In fact, theta is MaxMin_k since the options with k maximum minimum utility must be in RKS
//            # why must be in RKS, since you can't find k options r-dominate each of them
//            # except these k-1 options with maximum minimum utility
//            # the other option's minimum score must be lower than k-th option's corresponding score
            return theta;
        }else{
            for (int i = 0; i <idx.size(); ++i) {
                if(ubs[i]>theta){
                    new_idx.push_back(idx[i]);
                }
            }
            idx=new_idx; // deep copy here
            new_idx=vector<int>();
            assert(!idx.empty());
        }
    }

}

void topk_multi(cell &root, int k,  std::vector<std::vector<double>> &P, int num,
                Rtree *rtree_rt, unordered_map<long int, RtreeNode *> &ramTree, bool uHeat){
    if(rtree_rt== nullptr){
        build_rtree(rtree_rt, ramTree, P);
    }
    if(P.empty()){
        cout<<"option dataset is empty!"<<endl;
        return;
    }
    int dim=P[0].size();
    vector<vector<double>> query_ws=gen_weight_vector(k, dim, num);
    vector<double> MaxMin_k;
    if(uHeat) {
        for (auto &w: query_ws) { // get MaxMin_k
            MaxMin_k.push_back(getMaxMin_k(root, k, w, P));
        }
    }else{
        MaxMin_k=vector<double>(num, 0.0);
    }
    auto ab = chrono::steady_clock::now();

    for(int i=0;i<num;i++){
        topk_single(k, P, query_ws[i], MaxMin_k[i], rtree_rt, ramTree);
    }

    auto ae = chrono::steady_clock::now();
    chrono::duration<double> elapsed_seconds= ae-ab;
    cout << "100 top-k query time cost: " << elapsed_seconds.count() << endl;
}

void CSA(cell &root, std::vector<std::vector<double>> &P){
    for (int i = 0; i <P.size() ; ++i) {
        for (int j = i+1; j <P.size() ; ++j) {
            root.CSA_insert(i, j, P);
//            for(auto &child:root.children){
//                child->CSA_insert(i, j,P);
//            }
        }
    }
    vector<cell *> all_leaves;
    get_all_leaves(root, all_leaves);
    for(auto &leaf: all_leaves){
        for (int i = 0; i < P.size(); ++i) {
            if(leaf->dmc[i]<leaf->k){
                leaf->rkskyband.push_back(i);
            }
        }
    }
}

void CSAp(cell &root, std::vector<std::vector<double>> &P){
    for (int i = 0; i <P.size() ; ++i) {
        for (int j = i+1; j <P.size() ; ++j) {
            root.CSAp_insert(i, j, P);
//            for(auto &child:root.children){
//                child->CSAp_insert(i,j, P);
//            }
        }
    }
    vector<cell *> all_leaves;
    get_all_leaves(root, all_leaves);
    for(auto &leaf: all_leaves){
        for (int i = 0; i < P.size(); ++i) {
            if(leaf->dmc[i]<leaf->k){
                leaf->rkskyband.push_back(i);
            }
        }
    }
}

void MDA(cell &root, std::vector<std::vector<double>> &P){
    for (int i = 0; i <P.size() ; ++i) {
//        root.MDA_insert(i, P);
        for(auto &child:root.children){
            child->MDA_insert(i, P);
        }
    }
    vector<cell *> all_leaves;
    get_all_leaves(root, all_leaves);
    for(auto &leaf: all_leaves){
        leaf->MDA_superSet2RKS(P);
    }
}

void MDAp(cell &root, std::vector<std::vector<double>> &P){
    vector<int> idx=vector<int>(P.size());
    std::iota(idx.begin(), idx.end(), 0);
    for(auto &child:root.children){
        child->MDAp_insert(idx, P);
    }
}




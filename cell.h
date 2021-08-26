//
// Created by 12859 on 2021/8/25.
//


#ifndef HEATMAP_CELL_H
#define HEATMAP_CELL_H

#include <vector>
#include <string>
#include "utils.h"
#include <unordered_set>

#define mCSA  1
#define mCSAp 2
#define mMDA  3
#define mHMDA 4

class cell{
public:
    std::vector<cell*> children;
    std::vector<int> rkskyband;
    int dim;
    int k;
    std::vector<double> bounds;
    std::vector<std::vector<double>> vertexes;
    std::vector<int> dmc; // "dominated" count array, dmc[i]=j means option i dominated by totally j options
    double uHeat; // the k_th largest minimum score
    double rHeat; // No. of RSK
    double dHeat; // No. of r-dominate relationships

    vector<unordered_set<int>> rdo_graph; // only used for CSA+, a "dominate" graph

    cell(std::vector<double> &b, int cur_level, int tar_level, int card, int K, int m=mHMDA){
        /*
         * \para b, bounds
         * \para cur_level, this cell would be in which level of the quad-tree
         * \para tar_level, the quad-tree will build how many level use
         * \para card, cardinality of the option dataset, used for init dmc size and so on
         */
        dim=b.size()/2;
        bounds=b;
        uHeat=0.0;
        rHeat=0.0;
        dHeat=0.0;
        this->k=K;
        dmc=std::vector<int>(card, 0);
        // generate this cell's vertexes
        this->get_vertexes();

        // begin to generate children
        if(cur_level<tar_level){
            if(dim>=3){
                std::vector<double> tmp(2*dim);
                for (int i = 0; i < dim; ++i) {
                    tmp[i*2]=b[i*2];
                    tmp[i*2+1]=b[i*2];
                }
                this->children.push_back(new cell(tmp, cur_level+1, tar_level, card, k, m));
            }
            for (int d = 0; d <dim ; ++d) {
                std::vector<double> tmp;
                for (int d2 = 0; d2 < dim; ++d2) {
                    double l=this->bounds[d*2], u=this->bounds[d*2+1];
                    if(d==d2){
                        tmp.push_back((l+u)/2.0);
                        tmp.push_back(u);
                    }else{
                        tmp.push_back(l);
                        tmp.push_back((l+u)/2.0);
                    }
                }
                this->children.push_back(new cell(tmp, cur_level+1, tar_level, card, k, m));
            }
        }
        if(m==mCSAp){
            this->rdo_graph=vector<unordered_set<int>>(card);
        }
    }

    void get_vertexes(){
        for (int d = 0; d <dim ; ++d) {
            std::vector<double> v(dim);
            for (int d2 = 0; d2 < dim; ++d2) {
                if(d==d2){
                    v[d2]=this->bounds[d*2+1]; // add upper bound
                }else{
                    v[d2]=this->bounds[d*2];   // add lower bound
                }
            }
            this->vertexes.push_back(v);
        }
    }

    void CSA_insert(int p1, int p2, std::vector<std::vector<double>> &P){
        if(!(dmc[p1]>=k || dmc[p2]>=k)){
            int state=r_dominate(this->vertexes, P[p1], P[p2]);
            if(state==DOMINATE){
                this->recursively_update_dmc(p2);
            } else if(state==DOMINATED){
                this->recursively_update_dmc(p1);
            }else{
                for(auto &child: this->children){
                    child->CSA_insert(p1, p2, P);
                }
            }
        }
    }

    void recursively_update_dmc(int i){
        if(this->dmc[i]<k){
            this->dmc[i]+=1;
            for(auto &child: this->children){
                child->recursively_update_dmc(i);
            }
        }
    }

    void CSAp_insert(int p1, int p2, std::vector<std::vector<double>> &P){
        if( dmc[p1]>k || dmc[p2]>=k ){
            return;
        }
        if(rdo_graph[p1].find(p2)!=rdo_graph[p1].end() || rdo_graph[p2].find(p1)!=rdo_graph[p2].end()){
            return;
        }
        int state=r_dominate(vertexes, P[p1], P[p2]);
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
        recursively_update_dmc(j);
        CSAp_update_graph(i, j);
        for(int id:rdo_graph[j]){
            if(rdo_graph[i].find(id)==rdo_graph[i].end()){ // a new option
                recursively_update_dmc(id);
                CSAp_update_graph(i, id);
            }
        }
    }

    void CSAp_update_graph(int p1, int p2){
        if(dmc[p2]>=k || rdo_graph[p1].find(p2)!=rdo_graph[p1].end()){
            return;
        }
        rdo_graph[p1].insert(p2);
        for(auto &child: children){
            child->CSAp_update_graph(p1, p2);
        }
    }


};


void get_all_leaves(cell &node, std::vector<cell*>& ret);

#endif //HEATMAP_CELL_H

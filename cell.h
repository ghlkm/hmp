//
// Created by 12859 on 2021/8/25.
//


#ifndef HEATMAP_CELL_H
#define HEATMAP_CELL_H

#include <vector>
#include <string>
#include "utils.h"


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

    cell(std::vector<double> &b, int cur_level, int tar_level, int card, int K, std::string &m= (std::string &) "MDA"){
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
                std::vector<double> tmp(dim);
                for (int i = 0; i < dim; ++i) {
                    tmp[i]=b[i*2];
                }
                this->children.push_back(new cell(tmp, cur_level+1, card, k, tar_level));
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
                this->children.push_back(new cell(tmp, cur_level+1, card, k, tar_level));
            }
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


};


void get_all_leaves(cell &node, std::vector<cell*>& ret){
    if(node.children.empty()){
        ret.push_back(&node);
    }else{
        for(auto &child: node.children){
            get_all_leaves(*child, ret);
        }
    }
}

#endif //HEATMAP_CELL_H

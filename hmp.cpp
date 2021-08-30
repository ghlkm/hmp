//
// Created by 12859 on 2021/8/25.
//

#include "hmp.h"
#include <numeric>

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


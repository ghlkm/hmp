//
// Created by 12859 on 2021/8/25.
//

#ifndef HEATMAP_UTILS_H
#define HEATMAP_UTILS_H

#include <vector>
#include "vector_operator.h"

#define DOMINATE 1
#define DOMINATED -1
#define NONDO 0

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

void helpmsg(const char* pgm);

const char* read(int a_argc, const char** a_argv, const char* a_param, const char* a_def);


#endif //HEATMAP_UTILS_H

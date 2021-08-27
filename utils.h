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

int r_dominate(std::vector<std::vector<double>> &vs, std::vector<double> &p1, std::vector<double> &p2);

void helpmsg(const char* pgm);

const char* read(int a_argc, const char** a_argv, const char* a_param, const char* a_def);

vector<vector<double>> read_options(const char* datafile, int dim);

inline bool v1_dominate_v2(vector<double>& v1, vector<double>& v2);

void kskyband_nortree(vector<int> &ret, vector<vector<double>> &data,int k);

bool r_dominate(std::vector<double> &s1, std::vector<double> &s2);

void kskyband_write(vector<vector<double>> &data, int k, string &filename, vector<vector<int>> &ret);

void kskyband_read (const string &filename, vector<vector<int>> &ret);
#endif //HEATMAP_UTILS_H

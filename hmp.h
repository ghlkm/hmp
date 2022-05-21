
#ifndef HEATMAP_HMP_H
#define HEATMAP_HMP_H
#include "cell.h"
#include <vector>


void Baseline(cell &root, std::vector<std::vector<double>> &P);

void Baseline2(cell &root, std::vector<std::vector<double>> &P);

void CSA(cell &root, std::vector<std::vector<double>> &P);

void MDA(cell &root, std::vector<std::vector<double>> &P);

void MDAp(cell &root, std::vector<std::vector<double>> &P);

void topk_multi(cell &root, int k,  std::vector<std::vector<double>> &P, int num,
                Rtree *rtree_rt, bool uHeat);

#endif //HEATMAP_HMP_H

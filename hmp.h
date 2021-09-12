//
// Created by 12859 on 2021/8/25.
//

#ifndef HEATMAP_HMP_H
#define HEATMAP_HMP_H
#include "cell.h"
#include <vector>


void Baseline(cell &root, std::vector<std::vector<double>> &P);

void CSA(cell &root, std::vector<std::vector<double>> &P);

void CSAp(cell &root, std::vector<std::vector<double>> &P);

void MDA(cell &root, std::vector<std::vector<double>> &P);

void MDAp(cell &root, std::vector<std::vector<double>> &P);

#endif //HEATMAP_HMP_H

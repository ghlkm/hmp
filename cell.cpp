//
// Created by 12859 on 2021/8/25.
//

#include "cell.h"


void get_all_leaves(cell &node, std::vector<cell*>& ret){
    if(node.children.empty()){
        ret.push_back(&node);
    }else{
        for(auto &child: node.children){
            get_all_leaves(*child, ret);
        }
    }
}
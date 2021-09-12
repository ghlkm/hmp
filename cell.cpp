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

long rsky_c=0; // CSA, CSA+, leaf
long dmc_c=0;  // CSA, CSA+, all
long rdo_g_c=0; // CSA+, dominate number
long s_rsky_c=0; // MDA, all
long s_rsky_p_c=0; // MDA+, max k
long rtest_c=0; // r-dominate test count

vector<long> dmc_p_c=vector<long>(10, 0);
vector<long> dg_p_c=vector<long>(10, 0);

void cal_mem(cell &node){
    if(node.method==mCSA){
        rsky_c+=node.rkskyband.size();
        dmc_c+=node.dmc.size();
    }else if(node.method==mCSAp){
        rsky_c+=node.rkskyband.size();
        dmc_c+=node.dmc.size();
        for(auto &i:node.rdo_graph){
            rdo_g_c+=i.size();
        }
    }else if(node.method==mMDA){
        rsky_c+=node.rkskyband.size();
        s_rsky_c+=node.s_rskyband.size();
    }else if(node.method==mMDAp){
//        mda_max+=node.mdap_s_rksyband.size();
//        s_rsky_p_c=s_rsky_p_c<mda_max?mda_max:s_rsky_p_c;
    }else{
        // pass
    }
    for(auto &child: node.children){
        cal_mem(*child);
    }
}

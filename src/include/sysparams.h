//
// Created by parthajit on 23/2/20.
//

#ifndef BPNET_BASIC_03_SYSPARAMS_H
#include <string>
using namespace std;
class sysparams{
public:
    int _from_size;
    int _to_size;
    int _exdeg; //Exists at least one vertex with degree.
    int _num_exdeg;
    int _total_count;
    int _num_cycles;
    string cifpymol;
    string corpymol;
    long res_from_size;
    long res_to_size;
    string is_overlap;
    string _outformat;
    string type;
    int cleaned_res;
    double wt_overlap_cutoff;
    string adj_file;
    string file_dir;
    string accn;
    int overlap_flag; 
    sysparams(){
            res_from_size = 0;
            res_to_size = 99999999;
            cifpymol = "TRUE";
            corpymol = "TRUE";
        _from_size = 3;
        _to_size = 999999;
	cleaned_res = 0;
        _exdeg = 2;
        _total_count = 0;
        _num_exdeg = 1;
        _num_cycles = 0;
        _outformat = "new";
        type = "BP";
        wt_overlap_cutoff = 0.0001;
        adj_file = "TRUE";
        is_overlap = "TRUE";
	overlap_flag = 1;
    }
};

#define BPNET_BASIC_03_SYSPARAMS_H

#endif //BPNET_BASIC_03_SYSPARAMS_H

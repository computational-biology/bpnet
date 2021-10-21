//
// Created by parthajit on 23/2/20.
//

#ifndef BPNET_BASIC_03_SYSPARAMS_H
#include <string>
using namespace std;
class sysparams{
      public:


	    char cifparam[512];
	    char accnparam[512];
	    char htparam[512];
	    char hdparam[512];
	    char hdvalparam[512];
	    char chainparam[512];
	    char chainvalparam[512];
	    char angparam[512];
	    char angvalparam[512];
	    char chparam[512];
	    char sgparam[512];
	    char evaltypeparam[512];
	    char corparam[50];
	    char nmrparam[50];
	    char nmrvalparam[50];
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
	    string ext;
	    int total_per_structure;
	    int overlap_flag; 
	    int overlap_method;
	    sysparams(){
		  strcpy(cifparam, "-dummyval");
		  strcpy(accnparam, "-dummyval");
		  strcpy(htparam, "-HT");
		  strcpy(hdparam, "-dummyval");
		  strcpy(hdvalparam, "-dummyval");
		  strcpy(chainparam, "-dummyval");
		  strcpy(chainvalparam, "-dummyval");
		  strcpy(angparam, "-dummyval");
		  strcpy(angvalparam, "-dummyval");
		  strcpy(chparam, "-dummyval");
		  strcpy(sgparam, "-dummyval");
		  strcpy(evaltypeparam, "-dummyval");
		  strcpy(corparam, "-dummyval");
		  strcpy(nmrparam, "-dummyval");
		  strcpy(nmrvalparam, "-dummyval");
		  total_per_structure = 0;
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
		  type = "OL";
		  wt_overlap_cutoff = 0.0001;
		  adj_file = "TRUE";
		  is_overlap = "TRUE";
		  overlap_flag = 1;
		  overlap_method = 0; // 0 for overlap value, 1 for c1p-c1p dist in case of overlap;
	    }
	    void print_params(FILE* fp)
	    {
		  fprintf(fp, "\n---------------------------P A R A M S    O P T E D ---------------------\n");
		        if(strcmp(htparam, "-dummyval") == 0){
			      fprintf(fp,"PARAM   HETATM NOT REQUESTED\n");
			}else{
			      fprintf(fp,"PARAM   HETATM REQUESTED\n");
			}
		        if(strcmp(hdparam, "-dummyval") == 0){
			      fprintf(fp,"PARAM   HYDROGEN BOND DIST(FOR BASE PAIR)    3.8A (DEFAULT)\n");
			}else{
			      fprintf(fp,"PARAM   HYDROGEN BOND DIST(FOR BASE PAIR)    %sA (REQUESTED)\n",hdvalparam);
			}
		        if(strcmp(angparam, "-dummyval") == 0){
			      fprintf(fp,"PARAM   ANGLE IN DEGREE(FOR BASE PAIR)    120 (DEFAULT)\n");
			}else{
			      fprintf(fp,"PARAM   ANGLE IN DEGREE(FOR BASE PAIR)    %s (REQUESTED)\n",angvalparam);
			}
		        if(strcmp(chparam, "-dummyval") == 0){
			      fprintf(fp,"PARAM   C-H...O/N MEDIATED BASE PAIR   REQUESTED\n");
			}else{
			      fprintf(fp,"PARAM   C-H...O/N MEDIATED BASE PAIR   NOT REQUESTED\n");
			}
		        if(strcmp(sgparam, "-dummyval") == 0){
			      fprintf(fp,"PARAM   SUGAR O2' MEDIATED BASE PAIR   REQUESTED\n");
			}else{
			      fprintf(fp,"PARAM   SUGAR O2' MEDIATED BASE PAIR   NOT REQUESTED\n");
			}
			if(overlap_flag == 1){
			      fprintf(fp,"PARAM   VERTEX CONNECTIVITY   OVERLAP BASED\n");
			}else{
			      fprintf(fp,"PARAM   VERTEX CONNECTIVITY   BASE-PAIR BASED\n");
			}
			if(overlap_flag == 1 && overlap_method ==0){
			      fprintf(fp,"PARAM   EDGE WEIGHT        BASE-BASE SURFACE OVERLAP (LARGER VALUE INDICATES CLOSENESS)\n");
			      fprintf(fp,"PARAM   WEIGHT CUT-OFF     %5.2lf\n", wt_overlap_cutoff);

			}else{
			      fprintf(fp,"PARAM   EDGE WEIGHT        C1'-C1' DISTANCE (SMALLER VALUE INDICATES CLOSENESS)\n");
			}
			fprintf(fp,"PARAM   NETSIZE %2d-%2d\n", _from_size, _to_size);
			fprintf(fp,"PARAM   REQUESTED AT LEAST %d VERTEX WITH DEGREE %d\n", _num_exdeg, _exdeg);
			if(_num_cycles == 0)
			      fprintf(fp,"PARAM   CYCLES NOT REQUESTED\n");
			else
			      fprintf(fp,"PARAM   CYCLES REQUESTED.  AT LEAST %d CYCLES.\n", _num_cycles);


	    }
};

#define BPNET_BASIC_03_SYSPARAMS_H

#endif //BPNET_BASIC_03_SYSPARAMS_H

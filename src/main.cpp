#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <time.h>
#include <assert.h>
#include <string.h>


#include <rnabp.h>
#include <stdio.h>
#include <stdlib.h>
#include <spgraph.h>
#include <djset.h>
#include <overlap.h>
#include <helix.h>
#include <iostream>





#include <stdlib.h>
#include <overlap.h>
#include <disjointset.h>
#include <graph.h>
#include <sec_seq.h>
#include <ntvariants.h>




#include <fstream>
#include <sysparams.h>

//#include <libpfind.h>

extern "C" void callbpfindc(char [],  char [], char [], char [], char [], char [], char [], char [], char [], char [], char [], char[], char[]);
using namespace std;
// using namespace boost;



void print_curr_date_time(FILE* fp){
    char s[200];
    time_t cal_time=time(NULL);
    struct tm* local_time = localtime(&cal_time);
    strftime(s,190,"%a, %b, %d, %Y %I::%M::%S %p\n",local_time);
    fprintf(fp,"%s",s);

}


















void tokenize(vector<string>& result, string input, string delim_list="\t "){
    vector<string> tmp_token;
    char s[400];// = input.c_str();
    strcpy(s,input.c_str());
    char delim[20];
    strcpy(delim,delim_list.c_str());
    char* p = strtok (s,delim);
    while (p != NULL) {
        //printf ("%s\n",p);
        tmp_token.push_back(p);
        p = strtok (NULL, delim_list.c_str());
    }



    //split(tmp_token, input, is_any_of(delim_list), token_compress_on);
    for(int i=0; i<tmp_token.size(); i++){
        if(tmp_token[i] != ""){
            result.push_back(tmp_token[i]);
        }
    }
}


void read_rob_file_and_tokenize(string accn, vector<vector<string> >& rob_file_token_list){

    std::ifstream infile;
    infile.open((accn+".rob").c_str(), ios::in);
    if (infile.is_open()) {
        std::string line;
        while (getline(infile, line)) {
            if(line.substr(0,4) == "OVLP"){
                vector<string> token;
                tokenize(token, line);
                rob_file_token_list.push_back(token);
            }
        }
        infile.close();
    }
}


void populate_edges_from_rob_format(Graph* g, vector<vector<string> >& rob_file_tokens,
                                    vector<vector<string> >& out_file_tokens, sysparams* syspar,
                                    ntvariants_t* ntvar){
    int num_res_out = (int)out_file_tokens.size();
    assert(g->get_num_vertex()  == num_res_out);
    for(int i=0; i<num_res_out; i++) {
        vector<string> &ith_line_tokens = out_file_tokens[i];
        int cif_from = atoi(ith_line_tokens[1].c_str());
	string inscode = ith_line_tokens[3];
	char ins;
	if(inscode.length() != 1){
		fprintf(stderr,"Error... ins-code having more than one characters\n");
		exit(EXIT_FAILURE);
	}
	ins = inscode.c_str()[0];
        string chain_name = ith_line_tokens[4];
        int tmp_corid = atoi(ith_line_tokens[0].c_str());
        string tmp_res_name = ith_line_tokens[2];
        g->add_cif_details(tmp_corid, cif_from, chain_name.c_str(), ins, tmp_res_name.c_str(), ntvar);
    }
        int num_res = (int)rob_file_tokens.size();
    for(int i=0; i<num_res; i++){
        vector<string>& ith_line_tokens = rob_file_tokens[i];
        assert((int)ith_line_tokens.size() >= 9);

        string tmp_cif = ith_line_tokens[2];
        int tmppos1 = tmp_cif.find(":");

        int cif_from = atoi((tmp_cif.substr(0, tmppos1)).c_str());
        int cif_to = atoi((tmp_cif.substr(tmppos1+1)).c_str());
        string chain_name = ith_line_tokens[4];
        tmppos1 = chain_name.find("-");
        string tmp_chain_from = chain_name.substr(0,tmppos1);
        string tmp_chain_to = chain_name.substr(tmppos1+1);
        string tmp_res_name = ith_line_tokens[3];
        tmppos1 = tmp_res_name.find(":");
        string tmp_final_res_from = tmp_res_name.substr(0,tmppos1);
        string tmp_final_res_to = tmp_res_name.substr(tmppos1+1);
        tmppos1 = ith_line_tokens[1].find(":");
        int v11 = atoi(ith_line_tokens[1].substr(0, tmppos1).c_str());
        int v12 = atoi(ith_line_tokens[1].substr(tmppos1+1).c_str());
        //printf("%d, %d, %s, %s\n",cif_from,tmp_chain.c_str(), tmp_final_res.c_str());
        //g->add_cif_details(v11, cif_from, tmp_chain_from.c_str(), tmp_final_res_from.c_str(), ntvar);
        //g->add_cif_details(v12, cif_to, tmp_chain_to.c_str(), tmp_final_res_to.c_str(), ntvar);


        double ovlpval = atof(ith_line_tokens[8].c_str());
        string bptype = ith_line_tokens[6];
        if(bptype == "--" && ovlpval < syspar->wt_overlap_cutoff) continue;




        int pos = ith_line_tokens[1].find(":");
        int v1 = atoi(ith_line_tokens[1].substr(0, pos).c_str());
        int v2 = atoi(ith_line_tokens[1].substr(pos+1).c_str());








        //string base_type_from = ith_line_tokens[3]+"-"+ith_line_tokens[5];
        string base_type_from = ith_line_tokens[5];
        int pos1 = ith_line_tokens[3].find(":");
        string frm_res = ith_line_tokens[3].substr(0,pos1);
        string to_res = ith_line_tokens[3].substr(pos1+1);
        //cout<<frm_res<<endl<<to_res<<endl;
        string pair_type = "";
        pos1 = ith_line_tokens[5].find(":");
        if(pos1 < 0 ){
            pair_type = ith_line_tokens[5];
        }else{
            string edge1 = ith_line_tokens[5].substr(0, 1);
            string edge2 = ith_line_tokens[5].substr(2, 1);
            string orient = ith_line_tokens[5].substr(3, 1);
            pair_type = edge2+":"+edge1+orient;
        }
        //exit(1);



        //cout<<"v1="<<v1<<", v2="<<v2<<", ovlpval="<<ovlpval<<"."<<endl;
        //string base_type_to = to_res+":"+frm_res+"-"+pair_type;
        string base_type_to = pair_type;
        //cout<< base_type_from<<","<<base_type_to<<endl;
        g->addUndirectedEdge(v1-1, v2-1);
        //string bp = ith_line_tokens[2]+":"+ith_line_tokens[7]+"-"+ith_line_tokens[10];
        //sDtring bp = ith_line_tokens[2]+":"+ith_line_tokens[7]+"-"+ith_line_tokens[10];
        //sDtring bp = ith_line_tokens[2]+":"+ith_line_tokens[7]+"-"+ith_line_tokens[10];
        //double eval = atof(ith_line_tokens[12].c_str());
        g->addBP_and_Evalue_from_rob(v1-1, v2-1, base_type_from, base_type_to, ovlpval);
    }
    //exit(1);

}





void read_out_file_and_tokenize(string accn, vector<vector<string> >& out_file_token_list){

    std::ifstream infile((accn+".out").c_str(), ios::in);
    if (infile.is_open()) {
        std::string line;
        while (getline(infile, line)) {
            if(line.at(0) == ' ' || isdigit(line.at(0))){
                vector<string> token;
                tokenize(token, line);
                out_file_token_list.push_back(token);
            }
        }
        infile.close();
    }
}

int get_cleaned_residue_size(string file_name){
    std::ifstream infile;
    infile.open(file_name.c_str(), ios::in);
    if (infile.is_open()) {
        std::string line;
        while (getline(infile, line)) {
            vector<string> token;
            tokenize(token, line);
            if(token.size()==6){
                if(token[1] == "Cleaned" && token[4] == "residues:"){
                    int n = std::atoi(token[5].c_str());
                    infile.close();
                    return n;
                }
            }

        }
        infile.close();
    }
    cerr<<"Error... Something is wrong in computing cleaned residue"<<endl;
    exit(11);
}

void populate_edges_new_out_format(Graph* g, vector<vector<string> >& out_file_tokens, ntvariants_t* ntvar){
    int num_res = (int)out_file_tokens.size();
    assert(g->get_num_vertex()  == num_res);
    for(int i=0; i<num_res; i++){
        vector<string>& ith_line_tokens = out_file_tokens[i];
        int cif_from = atoi(ith_line_tokens[1].c_str());
	string inscode = ith_line_tokens[3];
	char ins;
	if(inscode.length() != 1){
		fprintf(stderr,"Error... ins-code having more than one characters\n");
		exit(EXIT_FAILURE);
	}
	ins = inscode.c_str()[0];
        string chain_name = ith_line_tokens[4];
        int tmp_corid = atoi(ith_line_tokens[0].c_str());
        string tmp_res_name = ith_line_tokens[2];
        g->add_cif_details(tmp_corid, cif_from, chain_name.c_str(), ins, tmp_res_name.c_str(), ntvar);

        if(ith_line_tokens.size() > 5){
            if(ith_line_tokens.size() == 13){
                int v1 = atoi(ith_line_tokens[0].c_str());
                int v2 = atoi(ith_line_tokens[5].c_str());
                g->addUndirectedEdge(v1-1, v2-1);
                //string bp = ith_line_tokens[2]+":"+ith_line_tokens[7]+"-"+ith_line_tokens[10];
                string bp = ith_line_tokens[10];
                double eval = atof(ith_line_tokens[12].c_str());
                g->addBP_and_Evalue(v1-1, v2-1, bp, eval);



            }else if(ith_line_tokens.size() == 21){
                int v1 = atoi(ith_line_tokens[0].c_str());
                int v2 = atoi(ith_line_tokens[5].c_str());
                int v3 = atoi(ith_line_tokens[13].c_str());
                g->addUndirectedEdge(v1-1, v2-1);
                g->addUndirectedEdge(v1-1, v3-1);
                //string bp = ith_line_tokens[2]+":"+ith_line_tokens[7]+"-"+ith_line_tokens[10];
                string bp = ith_line_tokens[10];
                double eval = atof(ith_line_tokens[12].c_str());
                g->addBP_and_Evalue(v1-1, v2-1, bp, eval);
                //bp = ith_line_tokens[2]+":"+ith_line_tokens[15]+"-"+ith_line_tokens[18];
                bp = ith_line_tokens[18];
                eval = atof(ith_line_tokens[20].c_str());
                g->addBP_and_Evalue(v1-1, v3-1, bp, eval);
            }else if(ith_line_tokens.size() == 29){
                int v1 = atoi(ith_line_tokens[0].c_str());
                int v2 = atoi(ith_line_tokens[5].c_str());
                int v3 = atoi(ith_line_tokens[13].c_str());
                int v4 = atoi(ith_line_tokens[21].c_str());
                g->addUndirectedEdge(v1-1, v2-1);
                g->addUndirectedEdge(v1-1, v3-1);
                g->addUndirectedEdge(v1-1, v4-1);
                //string bp = ith_line_tokens[2]+":"+ith_line_tokens[7]+"-"+ith_line_tokens[10];
                string bp = ith_line_tokens[10];
                double eval = atof(ith_line_tokens[12].c_str());
                g->addBP_and_Evalue(v1-1, v2-1, bp, eval);
                //bp = ith_line_tokens[2]+":"+ith_line_tokens[15]+"-"+ith_line_tokens[18];
                bp = ith_line_tokens[18];
                eval = atof(ith_line_tokens[20].c_str());
                g->addBP_and_Evalue(v1-1, v3-1, bp, eval);
                //bp = ith_line_tokens[2]+":"+ith_line_tokens[23]+"-"+ith_line_tokens[26];
                bp = ith_line_tokens[26];
                eval = atof(ith_line_tokens[28].c_str());
                g->addBP_and_Evalue(v1-1, v4-1, bp, eval);
            }else if(ith_line_tokens.size() == 37){
                int v1 = atoi(ith_line_tokens[0].c_str());
                int v2 = atoi(ith_line_tokens[5].c_str());
                int v3 = atoi(ith_line_tokens[13].c_str());
                int v4 = atoi(ith_line_tokens[21].c_str());
                int v5 = atoi(ith_line_tokens[29].c_str());
                g->addUndirectedEdge(v1-1, v2-1);
                g->addUndirectedEdge(v1-1, v3-1);
                g->addUndirectedEdge(v1-1, v4-1);
                g->addUndirectedEdge(v1-1, v5-1);
                //string bp = ith_line_tokens[2]+":"+ith_line_tokens[7]+"-"+ith_line_tokens[10];
                string bp = ith_line_tokens[10];
                double eval = atof(ith_line_tokens[12].c_str());
                g->addBP_and_Evalue(v1-1, v2-1, bp, eval);
                //bp = ith_line_tokens[2]+":"+ith_line_tokens[15]+"-"+ith_line_tokens[18];
                bp = ith_line_tokens[18];
                eval = atof(ith_line_tokens[20].c_str());
                g->addBP_and_Evalue(v1-1, v3-1, bp, eval);
                //bp = ith_line_tokens[2]+":"+ith_line_tokens[23]+"-"+ith_line_tokens[26];
                bp = ith_line_tokens[26];
                eval = atof(ith_line_tokens[28].c_str());
                g->addBP_and_Evalue(v1-1, v4-1, bp, eval);
                //bp = ith_line_tokens[2]+":"+ith_line_tokens[31]+"-"+ith_line_tokens[34];
                bp = ith_line_tokens[34];
                eval = atof(ith_line_tokens[36].c_str());
                g->addBP_and_Evalue(v1-1, v5-1, bp, eval);
            }else{
                cerr<<"More than quadplet found!!!!!"<<endl;
                cerr<<"cor file serial: "<<ith_line_tokens[0]<<endl;
                exit(1);
            }
        }
    }
}

void populate_edges(Graph* g, vector<vector<string> >& out_file_tokens, ntvariants_t* ntvar){
    int num_res = (int)out_file_tokens.size();
    assert(g->get_num_vertex()  == num_res);
    for(int i=0; i<num_res; i++){
        vector<string>& ith_line_tokens = out_file_tokens[i];

        int cif_from = atoi(ith_line_tokens[1].c_str());
        string chain_name = ith_line_tokens[3];
        int tmp_corid = atoi(ith_line_tokens[0].c_str());
        string tmp_res_name = ith_line_tokens[2];
	char ins='?'; /* no provision for ins in old format out file. */
        g->add_cif_details(tmp_corid, cif_from, chain_name.c_str(), ins, tmp_res_name.c_str(), ntvar);

        if(ith_line_tokens.size() > 4){
            if(ith_line_tokens.size() == 11){
                int v1 = atoi(ith_line_tokens[0].c_str());
                int v2 = atoi(ith_line_tokens[4].c_str());
                g->addUndirectedEdge(v1-1, v2-1);
            }else if(ith_line_tokens.size() == 18){
                int v1 = atoi(ith_line_tokens[0].c_str());
                int v2 = atoi(ith_line_tokens[4].c_str());
                int v3 = atoi(ith_line_tokens[11].c_str());
                g->addUndirectedEdge(v1-1, v2-1);
                g->addUndirectedEdge(v1-1, v3-1);
            }else if(ith_line_tokens.size() == 25){
                int v1 = atoi(ith_line_tokens[0].c_str());
                int v2 = atoi(ith_line_tokens[4].c_str());
                int v3 = atoi(ith_line_tokens[11].c_str());
                int v4 = atoi(ith_line_tokens[18].c_str());
                g->addUndirectedEdge(v1-1, v2-1);
                g->addUndirectedEdge(v1-1, v3-1);
                g->addUndirectedEdge(v1-1, v4-1);
            }else if(ith_line_tokens.size() == 32){
                int v1 = atoi(ith_line_tokens[0].c_str());
                int v2 = atoi(ith_line_tokens[4].c_str());
                int v3 = atoi(ith_line_tokens[11].c_str());
                int v4 = atoi(ith_line_tokens[18].c_str());
                int v5 = atoi(ith_line_tokens[25].c_str());
                g->addUndirectedEdge(v1-1, v2-1);
                g->addUndirectedEdge(v1-1, v3-1);
                g->addUndirectedEdge(v1-1, v4-1);
                g->addUndirectedEdge(v1-1, v5-1);
            }else{
                cerr<<"More than quadplet found!!!!!"<<endl;
                cerr<<"cor file serial: "<<ith_line_tokens[0]<<endl;
                exit(1);
            }
        }
    }
}
void compute_base_pair_components(string accn, FILE* outputfile,
                                  sysparams* syspar, ntvariants_t* ntvar, ovlp_stat* stat){
    int num_cleaned_res = get_cleaned_residue_size(syspar->file_dir+accn+".out");
    if( syspar->res_from_size > num_cleaned_res || syspar->res_to_size < num_cleaned_res) return;
    cout<<"        NUMBER OF VERTERX IS: "<<num_cleaned_res<<endl;
    //if(run_bpfind == BPFINDRUN::TRUE){
    //    cout<<"        PAIR COMPUTATION STARTS"<<endl;
    //    string file_name = accn+".cif";
    //runbpfind((char*)file_name.c_str(), 1); // 1 for cif 0 for pdb
    //    cout<<"        PAIRCOMPUTATION ENDS SUCCESSFULLY"<<endl;
    //}
    cout<<"        GRAPH GENERATION FROM BASE-PAIR FILE STARTS"<<endl;
    vector<vector<string> > out_file_tokens;
    vector<vector<string> > out_file_tokens_for_rob;
    if(syspar->is_overlap == "TRUE"){
        syspar->type = "OL";
        read_rob_file_and_tokenize(syspar->file_dir+accn, out_file_tokens);
        read_out_file_and_tokenize(syspar->file_dir+accn, out_file_tokens_for_rob);
    }else{
        syspar->type = "BP";
        read_out_file_and_tokenize(syspar->file_dir+accn, out_file_tokens);
    }
    Graph g = Graph(num_cleaned_res, 1); // 1 means full graph with cif id etc.

    //printf("test1\n");
    /*
        This part of the code is needed because the bpfind versionm is old (does not generate ?. before pdbe visit of Bhattacharyya
        So, if .out file is supplied, it is sssumed to be of new version.
    */
    //if(run_bpfind == BPFINDRUN::TRUE){
    //    if(syspar->_outformat == "new"){
    //        populate_edges_new_out_format(&g, out_file_tokens);
    //    }else{
    //        populate_edges(&g, out_file_tokens);
    //    }

    //}else{
    if(syspar->is_overlap == "FALSE"){
        if(syspar->_outformat == "new"){
            //printf("test2\n");
            populate_edges_new_out_format(&g, out_file_tokens, ntvar);
        }else{
            populate_edges(&g, out_file_tokens, ntvar);
        }
    }else if(syspar->is_overlap == "TRUE"){
        populate_edges_from_rob_format(&g, out_file_tokens, out_file_tokens_for_rob, syspar, ntvar);
    }else{
        cerr<<"Error... out/rob supplied"<<endl;
        exit(1);
    }
    //}
    cout<<"        GRAPH GENERATION ENDS SUCCESSFULLY"<<endl;
    cout<<"        COMPONENT GENERATION STARTS"<<endl;
    g.computeConnectedComponents();
    cout<<"        COMPONENT GENERATION ENDS SUCCESSFULLY"<<endl;
    int num_components = g.num_components();
    cout<<"        FILE WRITING PROCESS STARTS"<<endl;

    string adj_file_name = syspar->file_dir+accn+".adj";
   /*if(syspar->overlap_flag == 1){
	 adj_file_name= syspar->file_dir+accn+"_ol.adj";
   }else{
	 adj_file_name= syspar->file_dir+accn+"_bp.adj";
   }*/
    FILE* adj_file;
    string pymol_file = syspar->file_dir+accn+".pml";
    if(syspar->adj_file == "TRUE"){
        adj_file = fopen(adj_file_name.c_str(),"w");
    }
    string edge_file_name = syspar->file_dir+accn+".edge";
//   if(syspar->overlap_flag == 1){
//	 edge_file_name= syspar->file_dir+accn+"_ol.edge";
//   }else{
//	 edge_file_name= syspar->file_dir+accn+"_bp.edge";
//   }
    FILE* edge_file = fopen(edge_file_name.c_str(),"w");
   fprintf(adj_file, "+====================================================================+\n");
fprintf(adj_file, "+                                                                    +\n");
fprintf(adj_file, "+                         B P N E T                                  +\n");
fprintf(adj_file, "+                       VERSION - 1-0-2                              +\n");
fprintf(adj_file, "+                                                                    +\n");
fprintf(adj_file, "+             A Software to Compute Multiplate                       +\n");
fprintf(adj_file, "+            Base-Pair Networsk in Nucleic Acids.                    +\n");
fprintf(adj_file, "+                                                                    +\n");
fprintf(adj_file, "+                                                                    +\n");
fprintf(adj_file, "+             BUG  REPORT: roy.parthajit@gmail.com                   +\n");
fprintf(adj_file, "+          TECH REPORT: dhananjay.bhattacharyya@saha.ac.in           +\n");
fprintf(adj_file, "+====================================================================+\n");
    fprintf(adj_file, "\n\n");
    fprintf(adj_file, "ABVR    CARD       Cardinality or size of the component.\n");
    fprintf(adj_file, "ABVR    SL         Serial No of the residue. This will match with _rna.pdb\n");
    fprintf(adj_file, "ABVR    EDG-TYPE   What type of edge it is.\n");
    fprintf(adj_file, "ABVR    WT         Weight of the edge. If TYP is OL, then Weight is Overlap, else c1'-c1' dist.\n");
    fprintf(adj_file, "ABVR    SP         Secondary Position in the chain, C for Coil, H for Helix etc.\n");
    fprintf(adj_file, "ABVR    TYP        Type of Network. OL for overlap based, BP for base-pair based.\n");
    fprintf(adj_file, "ABVR    RES        PDB/mmCIF residur number.\n");
    fprintf(adj_file, "ABVR    CHN        PDB/mmCIF chain number.\n");
    fprintf(adj_file, "ABVR    BS         Base Name.\n");
    fprintf(adj_file, "ABVR    MDL        Component Model number. It is the serial number of the first vertex.\n");
    fprintf(adj_file, "\n---------------------------P A R A M S    O P T E D ---------------------\n");
    if(syspar->adj_file == "TRUE"){
        fprintf(adj_file,"PARAM   NETSIZE %2d-%2d\n", syspar->_from_size, syspar->_to_size);
        fprintf(adj_file,"PARAM   DEGREE EXISTS %2d\n", syspar->_exdeg);
	if(syspar->overlap_flag == 1){
	      fprintf(adj_file,"PARAM   OVERLAP   REQUESTED\n");
	}else{
	      fprintf(adj_file,"PARAM   OVERLAP   NOT REQUESTED\n");
	}
    }

    fprintf(adj_file, "\n--------------------------- S U M M A R Y  R E P O R T ---------------------\n");
    fprintf(adj_file, "DATE    %s\n", OvlpGen::today().c_str());
    fprintf(adj_file,"ACCN    %s\n", syspar->accn.c_str());
    
    fprintf(adj_file, "TOTAL RESIDUE : %d\n", syspar->cleaned_res);
    fprintf(adj_file, "TOTAL COMPONENTS FOUND : %d\n", num_components);
    if(syspar->overlap_flag == TRUE){
	  fprintf(adj_file, "NO. OF BASE-PAIRS : %d\n", stat->cancnt + stat->noncancnt);
	  fprintf(adj_file, "NO. OF CANONICAL : %d\n", stat->cancnt);
	  fprintf(adj_file, "NO. OF NON-CAN   : %d\n", stat->noncancnt);
	  fprintf(adj_file, "NO. OF BI-FAR    : %d\n", stat->bfcnt);
	  fprintf(adj_file, "NO. OF ASTK : %d\n", stat->astkcnt);
	  fprintf(adj_file, "NO. OF OSTK : %d\n", stat->ostkcnt);
	  fprintf(adj_file, "NO. OF ADJA : %d\n", stat->adjacnt);
	  fprintf(adj_file, "NO. OF CLOS : %d\n", stat->closcnt);
	  fprintf(adj_file, "NO. OF CROS : %d\n", stat->croscnt);
	  fprintf(adj_file, "NO. OF PROX : %d\n", stat->proxcnt);
    }
    fprintf(edge_file,"REMARK QUERY Number of vertex %2d-%2d exists degree %d\n",syspar->_from_size,
            syspar->_to_size, syspar->_exdeg);
    sequence_t seq;
    sequence_create(&seq, (syspar->file_dir+accn+".dat").c_str() );
    cout<<"        SECONDARY SEQUENCE GENERATION ENDS SUCCESSFULLY\n";
    //exit(1);
    for(int i=0; i<num_components; i++) {
        vector<int>& comp = g.get_component(i);
        int comp_size = (int) comp.size();

        if(comp_size>=9 && syspar->is_overlap == "FALSE"){
            cout<<"Nine or more vertex compnent Found in base-pair network: "<<comp_size<<endl;
            //int x11;
            // cin>>x11;
        }
        if (comp_size >= syspar->_from_size && comp_size <= syspar->_to_size) {
            //cout<<"reaches here1"<<endl;
            g.fprint_component_adjacency(i, adj_file, edge_file, &seq, syspar);
            //cout<<"reaches here2"<<endl;
            fprintf(outputfile,"FOUND %4s %4d", accn.c_str(),comp_size);
            for(int j=0; j<comp_size ; j++){
                int cor_serial = comp.at(j) + 1 ; // This + 1 is because matrix index i is cor index i+1
                fprintf(outputfile," %6d (%2d) ", cor_serial, g.get_degree(cor_serial - 1)); // j implies cor id j+1
            }
            fprintf(outputfile,"\n");
        }
    }
    fprintf(stdout,"\n");
    if(syspar->adj_file == "TRUE"){
        fprintf(adj_file,"END\n");
        fclose(adj_file);
    }


    fprintf(edge_file,"END\n");
    fclose(edge_file);
    if(syspar->cifpymol == "TRUE"){
	  if(syspar->is_overlap == "TRUE"){
		g.gen_pymol_cif(pymol_file.c_str(), syspar);
	  }else{
		
		
		g.gen_pymol_basepair_cif(pymol_file.c_str(), syspar);
	  }
    }
    if(syspar->corpymol == "TRUE"){
	    g.gen_pymol_cor(pymol_file.c_str(), syspar);
    }
    g.fprint_component_summary(stdout);
    cout<<"        FILE WRITING PROCESS ENDS SUCCESSFULLY"<<endl;
    /*if(file_clean == FILE_CLEANUP::YES){
        cout<<"        SYSTEM CLEANUP STARTS"<<endl;
        string delete_file = accn+".out";
        remove(delete_file.c_str());
        delete_file = accn+".cor";
        remove(delete_file.c_str());
        delete_file = accn+".dat";
        remove(delete_file.c_str());
        delete_file = accn+".fasta";
        remove(delete_file.c_str());
        delete_file = accn+".hlx";
        remove(delete_file.c_str());
        delete_file = accn+".nup";
        remove(delete_file.c_str());
        remove("fort.22");
        remove("fort.23");
        cout<<"        SYSTEM CLEANUP ENDS"<<endl;
    }*/
    cout<<"ENDS ACCN: "<<accn<<endl;
    cout<<"----------------------------------------------"<<endl<<endl;
}

    void gen_pymol_basepair(struct nucbp* nbp, struct djset* set, 
		const char* pmlfile, sysparams* syspar){

	  char pymolline[2000];
	  int n = set->numset;
	  if(n == 0) return;
	  FILE* fp = fopen(pmlfile, "w");
	  assert(fp != NULL);
	  //fp = stdout;
	  fprintf(fp,"load %s.cif\n",syspar->accn.c_str());
	  int* visited = (int*) malloc (set->size * sizeof(int));
	  for(int i=0; i<set->size; ++i){
		visited[i] = 0;
	  }

	  int valid_comp = 0;
	  fprintf(fp, "select proteinall, polymer.protein\n");
	  fprintf(fp, "color olive,  proteinall\n");
	  fprintf(fp, "hide everything, proteinall\n");
	  fprintf(fp, "select waterall, solvent\n");
	  fprintf(fp, "hide everything, waterall\n");
	  fprintf(fp, "select nucleicall, polymer.nucleic\n");
	  fprintf(fp, "set_bond stick_radius, 0.50, nucleicall\n");
	  fprintf(fp, "set sphere_scale, 0.25, all\n");
	  fprintf(fp, "color violet,  nucleicall\n");
	  fprintf(fp, "set sphere_scale, 0.30, all\n");
	  fprintf(fp, "set cartoon_ring_mode, 0, nucleicall\n");
	  fprintf(fp, "set cartoon_ladder_mode, 0, nucleicall\n");
	  if(strcmp(syspar->chainvalparam, "-dummyval") != 0){
		fprintf(fp, "hide everything, nucleicall\n");
		fprintf(fp, "select ch, chain %s\n", syspar->chainvalparam);
		fprintf(fp, "show cartoon, ch\n");
	  }else{
		fprintf(fp, "show cartoon, nucleicall\n");
	  }
	  int color=0;
	  char colorbase[][20] = {"marine", "red", "green", "white", };
	  char colorline[][20] = {"white", "green", "red", "marine"};
	  int max_color=4;
	  for (int i = 0; i < set->size; ++i) {
		int vertex = i;
		if(visited[vertex] == 1) continue;
		int comp_size = djset_composize(set, vertex);
		for(int j=0; j<comp_size; ++j){
		      visited[vertex] = 1;
		      vertex = djset_next(set, vertex);
		}
		if(comp_size < syspar->_from_size || comp_size > syspar->_to_size) continue; 
		fprintf(fp, "select comp%d, ", vertex+1);
		
		for(int j=0; j<comp_size; ++j){
		      fprintf(fp, "(resi %d and chain %s) ", nbp[vertex].cifid, nbp[vertex].chain);
		      vertex = djset_next(set, vertex);
		}
		fprintf(fp,"\n");
		int index = color % max_color;
///		fprintf(fp, "color %s, comp%d\n", colorbase[index], vertex+1);
		if(comp_size == 2){
		      fprintf(fp, "color slate, comp%d\n", vertex+1);
		}else if(comp_size == 3){
		      fprintf(fp, "color green, comp%d\n", vertex+1);

		}else if(comp_size == 4){
		      fprintf(fp, "color red, comp%d\n", vertex+1);

		}else if(comp_size == 6){
		      fprintf(fp, "color white, comp%d\n", vertex+1);

		}else if(comp_size == 7){
		      fprintf(fp, "color yellow, comp%d\n", vertex+1);

		}else if(comp_size == 8){
		      fprintf(fp, "color limegreen, comp%d\n", vertex+1);

		}else{
		      fprintf(fp, "color brightorange, comp%d\n", vertex+1);
		}

		color ++;
		fprintf(fp, "set cartoon_ring_mode, 2,  comp%d\n", vertex+1);
		fprintf(fp, "set cartoon_ladder_mode, 1, comp%d\n", vertex+1);
		fprintf(fp, "show_as cartoon, comp%d\n", vertex+1);
		fprintf(fp, "show line, comp%d\n", vertex+1);
		fprintf(fp,"\n\n\n");
		vertex = i; //This is for line generation.
		for(int j=0; j<comp_size; ++j){
		      for(int k=0; k<nbp[vertex].numbp; ++k){
			    int othindex = nbp[vertex].oth_base_index[k];
			    fprintf(fp, "distance line%d,  (resi %d and chain %s and (name N3 or name C6)) , (resi %d and chain %s and (name N3 or name C6)), 15.0, 4\n", i+1, 
					nbp[vertex].cifid, nbp[vertex].chain,
					nbp[othindex].cifid, nbp[othindex].chain);
		      }
		      vertex = djset_next(set, vertex);
		}
		fprintf(fp, "set dash_width, 2.0, line%d\n", i+1);
		fprintf(fp, "set dash_gap, 0, line%d\n", i+1);
//		fprintf(fp, "set dash_color, %s, line%d\n", colorline[index], i+1);
		if(comp_size == 2){
		      fprintf(fp, "set dash_color, orange, line%d\n", i+1);
		}else if(comp_size == 3){
		      fprintf(fp, "set dash_color, red, line%d\n", i+1);

		}else if(comp_size == 4){
		      fprintf(fp, "set dash_color, green, line%d\n", i+1);

		}else if(comp_size == 6){
		      fprintf(fp, "set dash_color, marine, line%d\n", i+1);

		}else if(comp_size == 7){
		      fprintf(fp, "set dash_color, purple, line%d\n", i+1);

		}else if(comp_size == 8){
		      fprintf(fp, "set dash_color, wheat, line%d\n", i+1);

		}else{
		      fprintf(fp, "set dash_color, white, line%d\n", i+1);
		}

		fprintf(fp, "hide label, line%d\n", i+1);
		valid_comp ++;
	  }
	  free(visited);

	  fclose(fp);
    }
    void gen_pymol(struct nucbp* nbp, struct djset* set, 
		const char* pmlfile, sysparams* syspar, int iscif){
	    
	    char color[30][15] = {"red", "purple", "green", "blue", "orange", "white", "yellow", "gray", 
		                  "brightorange", "brown", "cyan", "lightblue", "lightteal", "marine",
	                          "oxygen", "raspberry", "sand", "smudge", "wheat", "sulfur", "slate", "violet"};
	    int numcolor = 22;
	    int n = set->numset;
	    if(n == 0) return;
	    FILE* fp = fopen(pmlfile, "w");
	    assert(fp != NULL);
	    //fp = stdout;
	    if(iscif == 1){
		  fprintf(fp,"load %s.cif\n",syspar->accn.c_str());
	    }else{
		  fprintf(fp,"load %s_rna.pdb\n",syspar->accn.c_str());
	    }
	    fprintf(fp, "\nselect proteinall, polymer.protein");
	    fprintf(fp, "\ncolor olive,  proteinall");
	    fprintf(fp, "\nhide everything, proteinall\n");
	    fprintf(fp, "\nselect nucleicall, polymer.nucleic\n");
	    if(strcmp(syspar->chainvalparam, "-dummyval") != 0){
		  fprintf(fp, "hide everything, nucleicall\n");
		  fprintf(fp, "select ch, chain %s\n", syspar->chainvalparam);
		  fprintf(fp, "show cartoon, ch\n");
	    }else{
		  fprintf(fp, "show cartoon, nucleicall\n");
	  }
	    int* visited = (int*) malloc (set->size * sizeof(int));
	    for(int i=0; i<set->size; ++i){
		  visited[i] = 0;
	    }

	    int valid_comp = 0;
	    for (int i = 0; i < set->size; ++i) {
		  int vertex = i;
		    if(visited[vertex] == 1) continue;
		    int comp_size = djset_composize(set, vertex);
		    for(int j=0; j<comp_size; ++j){
			  visited[vertex] = 1;
			  vertex = djset_next(set, vertex);
		    }
		    if(comp_size < syspar->_from_size || comp_size > syspar->_to_size) continue; 
		    fprintf(fp, "select comp%d, ", vertex+1);
		    for(int j=0; j<comp_size; ++j){
			  if(iscif == 1){
				fprintf(fp, "(resi %d and chain %s) ", nbp[vertex].cifid, nbp[vertex].chain);
			  }else{
				fprintf(fp, "resi %d ", vertex + 1);
			  }		
			  vertex = djset_next(set, vertex);
		    }
		    fprintf(fp,"\n");
		    fprintf(fp, "color %s, comp%d\n", color[valid_comp % numcolor], vertex+1);
		    fprintf(fp, "show  spheres, comp%d\n", vertex+1);
		    fprintf(fp,"\n");
		    valid_comp ++;
	    }
	    free(visited);
	    fprintf(fp,"sele sugbackbone, name P+OP1+OP2+O5*+C5*+C4*+O4*+C3*+O3*+C2*+O2*\n"); 
	    fprintf(fp,"hide spheres, sugbackbone\n");

	    fclose(fp);
    }
//    void gen_pymol_cor(const char* pymolpml, sysparams* syspar){
//	    char color[30][10] = {"red", "purple","green", "blue","orange", "white"};
//	    int numcolor = 6;
//	    int n = num_components();
//	    if(n == 0) return;
//	    FILE* fp = fopen(pymolpml, "w");
//	    assert(fp != NULL);
//	    //fp = stdout;
//	    fprintf(fp,"load %s_rna.pdb\n",syspar->accn.c_str());
//	    for (int i = 0; i < n; ++i) {
//		    vector<int> component = all_components.at(i);
//		    int comp_size = (int) component.size();
//		    if(comp_size < syspar->_from_size || comp_size > syspar->_to_size) continue; 
//		    if(comp_size == 1 && syspar->is_overlap == "TRUE"){
//			    fprintf(stderr, "Worning... In overlap mode an isolated vertex found\n");
//		    } 
//		    //if(comp_size == 1) continue;
//		    int v = component.at(0);
//		    fprintf(fp, "select comp%d, ",v+1);
//		    fprintf(fp, "resi %d", v+1);
//		    for(int j=1; j<comp_size; ++j){
//			    int v1 = component.at(j);
//			    fprintf(fp, "+%d", v1+1);
//		    }
//		    fprintf(fp,"\n");
//		    fprintf(fp, "color %s, comp%d\n", color[i%numcolor], v+1);
//		    fprintf(fp, "show  spheres, comp%d\n", v+1);
//		    fprintf(fp,"\n");
//	    }
//	    fprintf(fp,"show cartoon\n");
//	    fprintf(fp,"sele sugbackbone, name P+OP1+OP2+O5*+C5*+C4*+O4*+C3*+O3*+C2*+O2*\n"); 
//	    fprintf(fp,"hide spheres, sugbackbone\n");
//
//	    fclose(fp);
//    }


void show_help(){
    cout<<" This program is for finding the base-pair network In RNAs."<<endl;
    cout<<" The switches are:"<<endl;
    cout<<"     -netsize=5 (prints only network of size 5), Default is 3"<<endl;
    cout<<"     -netsize=5-10 (prints only network of size 5 to 10), Default is 3-30"<<endl;
    cout<<"     -exdeg=3  (Prints all networks with at least a vertex with degree 3"<<endl;
    cout<<"                (Default is 2)"<<endl;
    cout<<"     -cycles=2 (Computes the number of cycles, Main cycles), Default 0";
    cout<<"     -ovlpcutoff=20.0 for non-base pair overlap cutoff value"<<endl;
    cout<<" It works for multiple files, like, bpnet *.out -netsize=4 -exdeg=3"<<endl;
    cout<<"\n\n For detailed help, use --genhelp option. this will generate the help.md file in the current directory"<<endl;

}
void gen_helix(nucbp* nbp, sysparams* syspar, int ressize){

      string hlxfile = syspar->file_dir+syspar->accn+".hlx";

      FILE *fphlx;										/* output-file pointer */

      fphlx = fopen(hlxfile.c_str(), "w" );
      if ( fphlx == NULL ) {
	    fprintf ( stderr, "couldn't open file '%s'; %s\n",
			hlxfile.c_str(), strerror(errno) );
	    exit (EXIT_FAILURE);
      }


      struct pseudo_helix* all_pseudo_helix[1000];
      int psucount = 0;



      helix_pseudo_calc_all(all_pseudo_helix, &psucount, nbp, ressize);




      helix_pseudo_fprint_all(all_pseudo_helix, psucount, nbp, ressize, fphlx);




      struct helix* helix_array;
      int hlxcount = 0;


      helix_init_all(&helix_array, ressize);




      helix_calc_all(helix_array, &hlxcount, nbp, ressize);


      helix_print_all(helix_array, nbp, hlxcount, fphlx);



      helix_gen_pymol(all_pseudo_helix, psucount, helix_array, hlxcount, nbp, ressize, syspar);



      helix_free_all(helix_array);


      if( fclose(fphlx) == EOF ) {			/* close output file   */
	    fprintf ( stderr, "couldn't close file '%s'; %s\n",
			hlxfile.c_str(), strerror(errno) );
	    exit (EXIT_FAILURE);
      }

}
int main(int argc, char* argv[]) {

    sysparams syspar = sysparams();
    OvlpParameters::proximity_method = METHOD_OVERLAP;
    ntvariants_t ntvariants;
    ntvar_populate(&ntvariants);
    string* file_array  = new string[argc];
    int file_count      = 0;
    
    OvlpSurfaceDataFile * surfgen = NULL;
    OvlpAllSurfacePoints * all_surf_points = NULL;
    OvlpRNA_NucVatiants *nucVariants = NULL;
    OvlpRNA_NucVatiants *nucVariants_prox =  new OvlpRNA_NucVatiants();
    
    if(argc == 1){
        show_help();
        exit(1);
    }
    int overlapwtflag =0;
    cout<<"THE ENTIRE PROCESS STARTS"<<endl;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if(arg.substr(0,6) == "--help"){
            show_help();
            exit(1);
        }else if(arg.substr(0,7) == "-chain="){
	      strcpy(syspar.chainparam, "-ML");
	      strcpy(syspar.chainvalparam,arg.substr(7).c_str());
	}else if(arg.substr(0,8)=="-hbdist="){
	      strcpy(syspar.hdparam, "-HD");
	      strcpy(syspar.hdvalparam,arg.substr(8).c_str());
	}else if(arg.substr(0,8)=="-cutang="){
	      strcpy(syspar.angparam, "-VA");
	      strcpy(syspar.angvalparam, arg.substr(8).c_str());
	}else if(arg.substr(0,13)=="-sugmed=false"){
	      strcpy(syspar.sgparam, "-SG");
	}else if(arg.substr(0,12)=="-chmed=false"){
	      strcpy(syspar.chparam, "-CH");
	}else if(arg.substr(0,13)=="-hetatm=false"){
	      strcpy(syspar.htparam, "-dummyval");
	}else if(arg.substr(0, 11) =="-eval=false"){
	      strcpy(syspar.evaltypeparam, "-c1");
	}


	
	
	
	else if(arg.substr(0,9) == "-nettype="){
	      if(arg.substr(9,7) == "contact"){
		    syspar.overlap_flag = 1;
		    syspar.is_overlap = "TRUE";
		    syspar.type = "OL";
	      }else if(arg.substr(9,8) == "basepair"){
		    syspar.overlap_flag = 0;
		    syspar.is_overlap = "FALSE";
	      }else{    /* Exception Handling */ 
		    fprintf(stderr, "Error in flag -nettype, the value will be either 'contact' or 'basepair'\n");
		    exit(EXIT_FAILURE);
	      }
	}else if(arg.substr(0,10) == "-cifpymol="){
		if(arg.substr(10,4) == "true"){
			syspar.cifpymol = "TRUE";
		}else if(arg.substr(10,5) == "false"){
			syspar.cifpymol = "FALSE";
		}else{
			cerr<<"invalid value supplied to -cifpymol... supply true or false"<<endl;
			exit(EXIT_FAILURE);
		}
	}else if(arg.substr(0,10) == "-rnapymol="){
		if(arg.substr(10,4) == "true"){
			syspar.corpymol = "TRUE";
		}else if(arg.substr(10,5) == "false"){
			syspar.corpymol = "FALSE";
		}else{
			cerr<<"invalid value supplied to -corpymol... supply true or false"<<endl;
			exit(EXIT_FAILURE);
		}
	}else if(arg.substr(0,10) == "-wtcutoff="){
            syspar.wt_overlap_cutoff = atof(arg.substr(10).c_str());
	    overlapwtflag =1;
        }else if(arg.substr(0,8) == "-cycles=") {
            syspar._num_cycles = atoi(arg.substr(8).c_str());
        }else if(arg.substr(0,15) == "-wttype=c1p-c1p") {
	      strcpy(syspar.evaltypeparam, "-c1");

	    OvlpParameters::proximity_method = METHOD_C1PC1P;
            syspar.overlap_method = 1;
        }else if(arg.substr(0,7) == "-exdeg="){
            syspar._exdeg = atoi(arg.substr(7).c_str());
        }else if(arg.substr(0,10) == "-numexdeg="){
            syspar._num_exdeg = atoi(arg.substr(10).c_str());
        }else if(arg.substr(0,9) == "-netsize="){
            string subarg = arg.substr(9);
            int split = subarg.find("-");
            if(split < 0){
                syspar._from_size = syspar._to_size = atoi(subarg.c_str());
            }else{
                syspar._from_size = atoi(subarg.substr(0,split).c_str());
                syspar._to_size =  atoi(subarg.substr(split+1).c_str());
            }
            if((syspar._from_size==0 && syspar._to_size==0)||(syspar._from_size<0 || syspar._to_size<0)||
               (syspar._from_size>syspar._to_size)){
                cerr<<"Error in -netsize params..."<<endl;
                exit(19);
            }
        }else if (arg.substr(0, 1) == "-"){
            cerr<<"Error... Invalid switch "<<arg<<endl;
            exit(31);
        }else{
            file_array[file_count] = arg;
            file_count++;
        }
    }
    if(overlapwtflag == 1 &&  syspar.overlap_flag == 1 && syspar.overlap_method == 1){
		fprintf(stderr, "Error... C1'-C1' dist and weight cutoff cannot be opted simultaneously for overlap.\n");
		exit(EXIT_FAILURE);
    }
    if(overlapwtflag == 1 && syspar.overlap_flag == 0){
		fprintf(stderr, "Error... Weight cutoff cannot be opted for base pair networks.\n");
		exit(EXIT_FAILURE);
    }

//    FILE* fp = fopen("pairchain.net", "w");
//    fprintf(fp,"STARTS\n");
    
    char* nucdir = getenv("NUCLEIC_ACID_DIR");
	if(nucdir == NULL){
		fprintf(stderr, "Error... NUCLEIC_ACID_DIR not Defined.\n");
		exit(EXIT_FAILURE);
	}

	char nucfiledir[512];
	strcpy(nucfiledir,nucdir);
    
    if(syspar.overlap_flag == 1){
      OvlpAtomType base_type[] = {Carbon,Nitrogen,Oxygen};
      surfgen = new OvlpSurfaceDataFile(nucfiledir,base_type,3);
      all_surf_points = surfgen->generate_surface_points(3);
      nucVariants = new OvlpRNA_NucVatiants();
    }
    syspar.print_params(stdout);
    for(int i=0; i<file_count; i++){
	  int ressize;
        string file = file_array[i];
        int pos_sep = (int)file.find_last_of("/");
        int pos_dot = (int)file.find_last_of(".");
        if(pos_dot < 0){
            cerr<<"File extension not supplied. Please supply .pdb or .cif"<<endl;
            exit(1);
        }

	ovlp_stat stat;
	ovlp_stat_init(&stat);


        string accn = file.substr(pos_sep + 1, pos_dot - (pos_sep + 1));
        syspar.accn = accn;
        string ext = file.substr(pos_dot+1);
        syspar.ext = ext;
        syspar.file_dir = file.substr(0, pos_sep + 1);
        
        cout<<"STARTS ACCN: "<<accn<<endl;
	if(ext == "cif" || ext == "pdb" || ext == "out"){
	      if(ext == "out" && syspar.overlap_flag == 1){
		    fprintf(stderr, "Error... out mode works only for basepair network\n");
		    exit(EXIT_FAILURE);
	      }
	      if(ext == "out"){
		    fprintf(stderr, "Warning... out mode works if all the base pair related files \n");
		    fprintf(stderr, "           are present in the directory. Further you cannot \n");
		    fprintf(stderr, "          change the base pair or weight parameters. \n");
		    fprintf(stderr, "          You can only change network parameters. So, if\n");
		    fprintf(stderr, "          you are not sure what you want, run with cif/pdb extension.\n");
	      }
	      strcpy(syspar.accnparam,file.c_str());
	      if(ext == "cif" || ext == "pdb"){
		    if(ext == "cif"){
			  strcpy(syspar.cifparam, "-cif");
		    }
		    callbpfindc(syspar.cifparam, syspar.accnparam, syspar.htparam, 
				syspar.hdparam, syspar.hdvalparam, syspar.angparam, 
				syspar.angvalparam, syspar.chparam, syspar.sgparam, 
				syspar.corparam, syspar.evaltypeparam,
				syspar.chainparam, syspar.chainvalparam);
	      }

	      string file_path = syspar.file_dir+syspar.accn+".out";

	      ressize = get_cleaned_residue_size(file_path);
	      syspar.cleaned_res = ressize;
	      if(ressize == 0){
		    cout<<"\n        !!!!!!        NO NUCLEIC ACID FOUND FOR THIS STRUCTURE\n"<<endl;
		    cout<<"ENDS ACCN: "<<syspar.accn<<endl;
		    cout<<"----------------------------------------------"<<endl<<endl;
		    continue;
	      }
	      if(syspar.overlap_flag == 0){
		    int resinum;
		    //ovlp_residue_all_prox_comp(syspar.file_dir+syspar.accn, 4.0, nucVariants_prox, &resinum, &stat);
		    //overlap_gen_contact_map_for_base_pair(resinum, syspar.file_dir, syspar.accn, &stat);
		    struct graph g;
		    struct djset set;
		    struct nucbp* nbp;
		    graph_init(&g, ressize, UNDIRECTED);
		    nbp = (struct nucbp*) malloc ((ressize) * sizeof(struct nucbp));
		    string out = syspar.file_dir+syspar.accn+".out";
		    rnabp_scan_out(nbp, &ntvariants, ressize, out.c_str(), &g);
		    gen_helix(nbp, &syspar, ressize);
		    djset_init(&set, ressize);
		    graph_kruskal_component(&g, &set);
		    print_adjinfo(nbp, ressize, &set, &ntvariants, &syspar, &g, &stat);


		    int iscif = 0;
		    string pymolfile;
		    if(syspar.cifpymol == "TRUE"){
			  iscif = 1;
			  pymolfile = syspar.file_dir+syspar.accn+"_cif.pml";
			  if(syspar.is_overlap == "TRUE"){
				gen_pymol(nbp, &set, pymolfile.c_str(), &syspar, iscif);
			  }else{
				gen_pymol(nbp, &set, pymolfile.c_str(), &syspar, iscif);
//				gen_pymol_basepair(nbp, &set, pymolfile.c_str(), &syspar);

			  }

		    }
		    if(syspar.corpymol == "TRUE"){
			  iscif = 0;
			  pymolfile = syspar.file_dir+syspar.accn+"_rna.pml";
			  gen_pymol(nbp, &set, pymolfile.c_str(), &syspar, iscif);
		    }

		    graph_free(&g);
		    djset_free(&set);
		    free(nbp);
	      }else{
		    cout<<"        CONTACT COMPUTATION STARTS"<<endl;
		    int resinum; 
		    string rob_file_path = syspar.file_dir+syspar.accn+".rob";
		    cout<<"        OVERLAP COMPUTATION STARTS"<<endl;
		    OvlpParameters::pdb_accn = syspar.accn;
		    OvlpRNA_All_Residues * rna;
                    rna = ovlp_base_overlap_comp(syspar.file_dir+syspar.accn, syspar.wt_overlap_cutoff,
                    
                            surfgen,
                            all_surf_points,
                            nucVariants) ;
		    
		    ovlp_residue_all_prox_comp(syspar.file_dir+syspar.accn, 4.0, nucVariants_prox, &resinum, &stat);
		    overlap_gen_contact_map(resinum, syspar.file_dir, syspar.accn, &stat, rna, syspar.chainvalparam);
		    struct graph g;
		    struct djset set;
		    struct nucbp* nbp;
		    graph_init(&g, ressize, UNDIRECTED);
		    nbp = (struct nucbp*) malloc ((ressize) * sizeof(struct nucbp));
		    string out = syspar.file_dir+syspar.accn+".rob";
		    nucbp_scan_rob(nbp, &ntvariants, ressize, out.c_str(), &g);

		    // This part is for helix and pymol helix generation on contact base netwok mode.
		    struct nucbp* nbphlx;
		    nbphlx = (struct nucbp*) malloc ((ressize) * sizeof(struct nucbp));
		    struct graph ghlx;
		    graph_init(&ghlx, ressize, UNDIRECTED);
		    string outhlx = syspar.file_dir+syspar.accn+".out";
		    rnabp_scan_out(nbphlx, &ntvariants, ressize, outhlx.c_str(), &ghlx);
		    gen_helix(nbphlx, &syspar, ressize);
		    graph_free(&ghlx);
		    free(nbphlx);
		    // End of addition for helix generation.



		    djset_init(&set, ressize);
		    graph_kruskal_component(&g, &set);
		    print_adjinfo(nbp, ressize, &set, &ntvariants, &syspar, &g, &stat);

		    int iscif = 0;
		    string pymolfile;
		    if(syspar.cifpymol == "TRUE"){
			  iscif = 1;
			  pymolfile = syspar.file_dir+syspar.accn+"_cif.pml";
			  gen_pymol(nbp, &set, pymolfile.c_str(), &syspar, iscif);
		    }
		    if(syspar.corpymol == "TRUE"){
			  iscif = 0;
			  pymolfile = syspar.file_dir+syspar.accn+"_rna.pml";
			  gen_pymol(nbp, &set, pymolfile.c_str(), &syspar, iscif);
		    }

		    graph_free(&g);
		    djset_free(&set);
		    free(nbp);
		    delete rna;

	      }
        }else{
            cerr<<"Error..... invalid extension of the file "<<ext<<". Supply .cif or .pdb file"<<endl;
	      exit(EXIT_FAILURE);
	}
	FILE* summfp = stdout;
	fprintf(summfp, "---------------------------S U M M A R Y  R E P O R T ---------------------\n");
	fprintf(summfp, "TOTAL RESIDUE : %d\n", ressize);
	fprintf(summfp, "TOTAL NO. OF NETWORKS IN STRUCTURE : %d\n", syspar.total_per_structure);
	if(syspar.overlap_flag == TRUE){
	      fprintf(summfp, "NO. OF BASE-PAIRS : %d\n", stat.cancnt + stat.noncancnt);
	      fprintf(summfp, "NO. OF CANONICAL : %d\n", stat.cancnt);
	      fprintf(summfp, "NO. OF NON-CAN   : %d\n", stat.noncancnt);
	      fprintf(summfp, "NO. OF BI-FAR    : %d\n", stat.bfcnt);
	      fprintf(summfp, "NO. OF ASTK : %d\n", stat.astkcnt);
	      fprintf(summfp, "NO. OF OSTK : %d\n", stat.ostkcnt);
	      fprintf(summfp, "NO. OF ADJA : %d\n", stat.adjacnt);
	      fprintf(summfp, "NO. OF CLOS : %d\n", stat.closcnt);
	      fprintf(summfp, "NO. OF CROS : %d\n", stat.croscnt);
	      fprintf(summfp, "NO. OF PROX : %d\n", stat.proxcnt);
	}else{
	      fprintf(summfp, "CONTACT NETWORK NOT REQUESTED\n");
	}
	cout<<"ENDS ACCN: "<<accn<<endl;
	cout<<"----------------------------------------------"<<endl<<endl;
    }
    if(syspar.overlap_flag == 1){
      delete surfgen;
      delete all_surf_points;
      delete nucVariants;
    }
    printf("Total %d components found on your query\n",syspar._total_count);
    delete [] file_array;
    cout<<"THE ENTIRE PROCESS COMPLETES"<<endl;
    return 0;
}

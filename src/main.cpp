#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <time.h>
#include <assert.h>
#include <string.h>







#include <stdlib.h>
#include <overlap.h>
#include <disjointset.h>
#include <graph.h>
#include <sec_seq.h>
#include <ntvariants.h>




#include <fstream>
#include <sysparams.h>

//#include <libpfind.h>

extern "C" void callbpfindc(char [],  char [], char [], char [], char [], char [], char [], char [], char [], char []);
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
    Graph g = Graph(num_cleaned_res, 1); // 0 means full graph with cif id etc.

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

    fprintf(adj_file, "\n---------------------------S U M M A R Y  R E P O R T ---------------------\n");
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
	    g.gen_pymol_cif(pymol_file.c_str(), syspar);
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



void show_help(){
    cout<<" This program is for finding the base-pair network In RNAs."<<endl;
    cout<<" It assumes a .out or .rob file as its input."<<endl;
    cout<<" In case of .rob file it expects .out file as well."<<endl;
    cout<<" In both the cases, it assumes the .dat file."<<endl;
    cout<<" The switches are:"<<endl;
    cout<<"     -netsize=5 (prints only network of size 5), Default is 3"<<endl;
    cout<<"     -netsize=5-10 (prints only network of size 5 to 10), Default is 3-30"<<endl;
    cout<<"     -exdeg=3  (Prints all networks with at least a vertex with degree 3"<<endl;
    cout<<"                (Default is 2)"<<endl;
    cout<<"     -cycles=2 (Computes the number of cycles, Main cycles), Default 0";
    cout<<"     -outformat=new/old (default new)"<<endl;
    cout<<"     -ovlpcutoff=20.0 for non-base pair overlap cutoff value"<<endl;
    cout<<"     -adj=true for creation of adjacency matrix. Default false"<<endl;
    cout<<"\n\n It generates a .edge file for edge lists, .adj(optional) file for adjacency"<<endl;
    cout<<" matrix for every file. It also generates a pairchain.net for all files."<<endl;
    cout<<" It works for multiple files, like, bpnet *.out -netsize=4 -exdeg=3"<<endl;

}
int main(int argc, char* argv[]) {

    sysparams syspar = sysparams();
    ntvariants_t ntvariants;
    ntvar_populate(&ntvariants);
    string* file_array  = new string[argc];
    int file_count      = 0;
    char cifparam[512] = "-dummyval";
    char accnparam[512] = "-dummyval";
    char htparam[512] = "-HT";
    char hdparam[512] = "-dummyval";
    char hdvalparam[512] = "-dummyval";
    char angparam[512] = "-dummyval";
    char angvalparam[512] = "-dummyval";
    char chparam[512] = "-dummyval";
    char sgparam[512] = "-dummyval";
    char corparam[50] = "-dummyval";
    
    OvlpSurfaceDataFile * surfgen = NULL;
    OvlpAllSurfacePoints * all_surf_points = NULL;
    OvlpRNA_NucVatiants *nucVariants = NULL;
    OvlpRNA_NucVatiants *nucVariants_prox =  new OvlpRNA_NucVatiants();
    
    if(argc == 1){
        show_help();
        exit(1);
    }
    cout<<"THE ENTIRE PROCESS STARTS"<<endl;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if(arg.substr(0,5) == "-help"){
            show_help();
            exit(1);
        }else if(arg.substr(0,8)=="-hbdist="){
	      strcpy(hdparam, "-HD");
	      strcpy(hdvalparam,arg.substr(8).c_str());
	}else if(arg.substr(0,8)=="-cutang="){
	      strcpy(angparam, "-VA");
	      strcpy(angvalparam, arg.substr(8).c_str());
	}else if(arg.substr(0,13)=="-sugmed=false"){
	      strcpy(sgparam, "-SG");
	}else if(arg.substr(0,12)=="-chmed=false"){
	      strcpy(chparam, "-CH");
	}else if(arg.substr(0,13)=="-hetatm=false"){
	      strcpy(chparam, "-dummyval");
	}


	
	
	
	else if(arg.substr(0,9) == "-overlap="){
	      if(arg.substr(9,4) == "true"){
		    syspar.overlap_flag = 1;
		    syspar.is_overlap = "TRUE";
	      }else if(arg.substr(9,5) == "false"){
		    syspar.overlap_flag = 0;
		    syspar.is_overlap = "FALSE";
	      }else{    /* Exception Handling */ 
		    fprintf(stderr, "Error in flag -overlap, the value will be either true or false\n");
		    exit(EXIT_FAILURE);
	      }
	}else if(arg.substr(0,5) == "-adj="){
            if(arg.substr(5)== "true"){
                syspar.adj_file = "TRUE";
            }else if(arg.substr(5)== "false"){
                syspar.adj_file = "FALSE";
            }else{
                cerr<<"invalid value supplied to -adj... supply true or false"<<endl;
                exit(EXIT_FAILURE);
            }
        }else if(arg.substr(0,10) == "-cifpymol="){
		if(arg.substr(10,4) == "true"){
			syspar.cifpymol = "TRUE";
			syspar.corpymol = "FALSE";
		}else if(arg.substr(10,5) == "false"){
			syspar.cifpymol = "FALSE";
		}else{
			cerr<<"invalid value supplied to -cifpymol... supply true or false"<<endl;
			exit(EXIT_FAILURE);
		}
	}else if(arg.substr(0,10) == "-rnapymol="){
		if(arg.substr(10,4) == "true"){
			syspar.corpymol = "TRUE";
			syspar.cifpymol = "FALSE";
		}else if(arg.substr(10,5) == "false"){
			syspar.corpymol = "FALSE";
		}else{
			cerr<<"invalid value supplied to -corpymol... supply true or false"<<endl;
			exit(EXIT_FAILURE);
		}
	}/*else if(arg.substr(0,8) == "-numres="){
            string subarg = arg.substr(8);
            int split = subarg.find("-");
            if(split < 0){
                cerr<<"Error in -numres params..."<<endl;
            }else{
                syspar.res_from_size = atoi(subarg.substr(0,split).c_str());
                syspar.res_to_size =  atoi(subarg.substr(split+1).c_str());
            }
            if((syspar.res_from_size == 0 && syspar.res_to_size == 0)||(syspar.res_from_size<0 || syspar.res_to_size<0)|| (syspar.res_from_size > syspar.res_to_size)){
                cerr<<"Error in -numres params..."<<endl;
                exit(19);
            }
        }*/else if(arg.substr(0,12) == "-ovlpcutoff="){
            cerr<<"Worning.... Overlap cutoff works only for overlap and is ignored for base-pair networks"<<endl;
            syspar.wt_overlap_cutoff = atof(arg.substr(12).c_str());
        }else if(arg.substr(0,8) == "-cycles=") {
            syspar._num_cycles = atoi(arg.substr(8).c_str());
        }else if(arg.substr(0,7) == "-exdeg="){
            syspar._exdeg = atoi(arg.substr(7).c_str());
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
    cout<<" THE QUERY: From ="<<syspar._from_size<<", to="<<syspar._to_size<<" with minsize="<<syspar._exdeg<<endl;
    FILE* fp = fopen("pairchain.net", "w");
    fprintf(fp,"STARTS\n");
    
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
    for(int i=0; i<file_count; i++){
	  int ressize;
        string file = file_array[i];
        int pos_sep = (int)file.find_last_of("/");
        int pos_dot = (int)file.find_last_of(".");
        if(pos_dot < 0){
            cerr<<"File extension not supplied. Please supply .pdb or .cif, out or .rob file"<<endl;
            exit(1);
        }

	ovlp_stat stat;
	ovlp_stat_init(&stat);


        string accn = file.substr(pos_sep + 1, pos_dot - (pos_sep + 1));
        syspar.accn = accn;
        string ext = file.substr(pos_dot+1);
        syspar.file_dir = file.substr(0, pos_sep + 1);
        
        cout<<"STARTS ACCN: "<<accn<<endl;
        /*if(ext == "out"){
	      string file_path = syspar.file_dir+syspar.accn+".out";
	      FILE* fp = fopen(file_path.c_str(),"r");
	      if(fp == NULL){     
		    fprintf(stderr, "Error... out format chosen, but directory does not contain .out or .dat file\n");
		    exit(EXIT_FAILURE);
	      }
	      fclose(fp);

	      file_path = syspar.file_dir+syspar.accn+".dat";
	      fp = fopen(file_path.c_str(),"r");
	      if(fp == NULL){     
		    fprintf(stderr, "Error... out format chosen, but directory does not contain .out or .dat file\n");
		    exit(EXIT_FAILURE);
	      }
	      fclose(fp);
            compute_base_pair_components(accn, fp, &syspar, &ntvariants);
        }
        else */ 
	if(ext == "cif" || ext == "pdb"){
	      strcpy(accnparam,file.c_str());
	      /*if(syspar.overlap_flag == 0){
		    strcpy(corparam,"-NOCOR");
	      }*/
	      if(ext == "cif" || ext == "pdb"){
		    if(ext == "cif"){
			  strcpy(cifparam, "-cif");
		    }
//		    printf("Starting BPFIND\n"); // bhatta
		    callbpfindc(cifparam, accnparam, htparam, 
				hdparam, hdvalparam, angparam, 
				angvalparam, chparam, sgparam, 
				corparam );
//		    printf("Finished BPFIND\n"); // bhatta
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
		    compute_base_pair_components(accn, fp, &syspar, &ntvariants, &stat);
	      }else{
		    /*syspar.overlap_flag = 0;
		    syspar.is_overlap = "FALSE";
		    compute_base_pair_components(accn, fp, &syspar, &ntvariants);
		    syspar.overlap_flag = 1;
		    syspar.is_overlap = "TRUE";*/
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
		    
		    //call_overlap(rob_file_path); 
		    
		    ovlp_residue_all_prox_comp(syspar.file_dir+syspar.accn, 4.0, nucVariants_prox, &resinum, &stat);
		    overlap_gen_contact_map(resinum, syspar.file_dir, syspar.accn, &stat, rna);
		    //rna->gen_dist_map();
		    compute_base_pair_components(accn, fp, &syspar, &ntvariants, &stat);
		    delete rna;

	      }

            // Compute_base_pair_components(accn, fp, &syspar);
        }else{
            cerr<<"Error..... invalid extension of the file "<<ext<<". Supply .cif or .pdb file"<<endl;
	      exit(EXIT_FAILURE);
	}
	/*
	return 0;
	if(ext == "out"){
	      string file_path = syspar.file_dir+syspar.accn+".out";
	      FILE* fp = fopen(file_path.c_str(),"r");
	      if(fp == NULL){     
		    fprintf(stderr, "Error... out format chosen, but directory does not contain .out or .dat file\n");
		    exit(EXIT_FAILURE);
	      }
	      fclose(fp);

	      file_path = syspar.file_dir+syspar.accn+".dat";
	      fp = fopen(file_path.c_str(),"r");
	      if(fp == NULL){     
		    fprintf(stderr, "Error... out format chosen, but directory does not contain .out or .dat file\n");
		    exit(EXIT_FAILURE);
	      }
	      fclose(fp);
            compute_base_pair_components(accn, fp, &syspar, &ntvariants);
        }else if(syspar.overlap_flag == 1){
	      string rob_file_path = syspar.file_dir+syspar.accn+".rob";
	      
	      syspar.is_overlap = "TRUE";
	      compute_base_pair_components(accn, fp, &syspar, &ntvariants);
	}else if(ext == "rob"){
	      string file_path = syspar.file_dir+syspar.accn+".out";
	      FILE* fp = fopen(file_path.c_str(),"r");
	      if(fp == NULL){     
		    fprintf(stderr, "Error... Lazy overlap chosen, but directory does not contain .out, .dat or .rob file\n");
		    exit(EXIT_FAILURE);
	      }
	      fclose(fp);

	      file_path = syspar.file_dir+syspar.accn+".dat";
	      fp = fopen(file_path.c_str(),"r");
	      if(fp == NULL){     
		    fprintf(stderr, "Error... Lazy overlap chosen, but directory does not contain .out, .dat or .rob file\n");
		    exit(EXIT_FAILURE);
	      }
	      fclose(fp);
	      file_path = syspar.file_dir+syspar.accn+".rob";
	      fp = fopen(file_path.c_str(),"r");
	      if(fp == NULL){     
		    fprintf(stderr, "Error... Lazy overlap chosen, but directory does not contain .out, .dat or .rob file\n");
		    exit(EXIT_FAILURE);
	      }
	      fclose(fp);
	      syspar.is_overlap = "TRUE";
	      compute_base_pair_components(accn, fp, &syspar, &ntvariants);
        }else{
            cerr<<"Error..... invalid extension of the file "<<ext<<". Supply .out or .rob file"<<endl;
            exit(1);
        }
    */
	FILE* summfp = stdout;
	fprintf(summfp, "---------------------------S U M M A R Y  R E P O R T ---------------------\n");
	fprintf(summfp, "TOTAL RESIDUE : %d\n", ressize);
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
	      fprintf(summfp, "Overlap not requested");


	}
    }
    if(syspar.overlap_flag == 1){
      delete surfgen;
      delete all_surf_points;
      delete nucVariants;
    }
    fprintf(fp,"ENDS");
    fclose(fp);
    printf("Total %d components found on your query\n",syspar._total_count);
    delete [] file_array;
    cout<<"THE ENTIRE PROCESS COMPLETES"<<endl;
    return 0;
}

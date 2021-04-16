#ifndef __BPNET_RNABP_H__
#define __BPNET_RNABP_H__
#include <stdio.h>
#include <stdlib.h>

#include <spgraph.h>
#include <ntvariants.h>
#include <djset.h>
#include <sysparams.h>
#include <sec_seq.h>
#include <overlap.h>



#define MAX_EDGE (7)
struct nucbp{
      int cifid;
      char resname[4];
      char resclass;
      char chain[4];
      char ins;
      char name[MAX_EDGE][5];  //W:WC
      char type[MAX_EDGE][3];  // BP, TP, BF
      double eval[MAX_EDGE];

      int numbp;
};

void nucbp_scan_rob(struct nucbp *self, ntvariants_t* ntvar, int numres, const char *robfile, struct graph* g)
{
      FILE* fp = fopen(robfile, "r");
      if(fp == NULL){    /* Exception Handling */ 
	    fprintf(stderr, "Error in %s at line %d... Unable to open file %s\n", __FILE__, __LINE__, robfile);
	    exit(EXIT_FAILURE);
      }
      char line[256];
      while(fgets(line, 256, fp) != NULL){
	    if(strncmp(line,"NUMROW", 6) == 0) break;

      }

      char sep[] = "\t :-\n";
      char *token;
      token = strtok(line, sep); //skip NUMROW word.
      token = strtok(NULL, sep); //Consider number of row.
      int numrow = atoi(token);
      for(int i=0; i<numres; ++i){
	    self[i].numbp = 0;
      }

      for(int k=0; k<numrow; ++k){	    
	    if(fgets(line, 256, fp) == NULL){
		  fprintf(stderr, "Error... Abnormal end of file encountered\n");
		  exit(EXIT_FAILURE);
	    }
	    token = strtok(line, sep);

	    token = strtok(NULL, sep);
	    int i = atoi(token);
	    i --;
	    self[i].numbp ++;

	    token = strtok(NULL, sep);
	    int j = atoi(token);
	    j --;
	    self[j].numbp++;


	    token = strtok(NULL, sep);
	    self[i].ins = token[0];

	    token = strtok(NULL, sep);
	    self[i].cifid = atoi(token);


	    token = strtok(NULL, sep);
	    self[j].cifid = atoi(token);

	    token = strtok(NULL, sep);
	    self[j].ins = token[0];

	    token = strtok(NULL, sep);
	    strcpy(self[i].resname, token);
	    self[i].resclass = ntvar_getchar(ntvar, self[i].resname);

	    token = strtok(NULL, sep);
	    strcpy(self[j].resname, token);
	    self[j].resclass = ntvar_getchar(ntvar, self[j].resname);




	    token = strtok(NULL, sep);
	    strcpy(self[i].chain, token);

	    token = strtok(NULL, sep);
	    strcpy(self[j].chain, token);

	    // reading base pair name, : is a valid char here.
	    token = strtok(NULL, "\t \n");
	    char bpname[6];
	    strcpy(bpname, token);
	    int pos = graph_set_edge(g, i, j);
	    strcpy(self[i].name[pos], token);

	    token = strtok(NULL, "\t  \n");
	    strcpy(self[i].type[pos], token);

	    //skipping :
	    token = strtok(NULL, "\t  \n");

	    // Reading overlap value 
	    token = strtok(NULL, "\t  \n");
	    self[i].eval[pos] = atof(token);

	    graph_set_wt(g, i, pos , self[i].eval[pos]);

	    int pos2 = graph_set_edge(g, j, i);
	    if(bpname[1] == ':') { // If base pair, then reverse the name W:HT  to H:WT
		  char ch = bpname[2];
		  bpname[2] = bpname[0];
		  bpname[0] = ch;
		  strcpy(self[j].name[pos2], bpname);
	    }else{
		  strcpy(self[j].name[pos2], bpname);
	    }
	    strcpy(self[j].type[pos2], self[i].type[pos]);
	    self[j].eval[pos2] = self[i].eval[pos];
	    graph_set_wt(g, j, pos2, self[i].eval[pos]);

      }
      fclose(fp);

}
void rnabp_scan_out(struct nucbp *self, ntvariants_t* ntvar, int numres, const char *outfile, struct graph* g) 
{
      FILE* fp = fopen(outfile, "r");
      if(fp == NULL){    /* Exception Handling */ 
	    fprintf(stderr, "Error in %s at line %d... Unable to open file %s\n", __FILE__, __LINE__, outfile);
	    exit(EXIT_FAILURE);
      }
      char line[256];
      while(fgets(line, 256, fp) != NULL){
	    if(line[0] == '#') continue;
	    break; /* pass the lines */
      }

      int i;
      char sep[] = "\t \n";
      char *token;
      while(line != NULL){
	    if(strlen(line)<=6) break;

	    token = strtok(line, sep);
	    int i = atoi(token);
	    i --;
	    self[i].numbp = 0;

	    token = strtok(NULL, sep);
	    self[i].cifid = atoi(token);

	    token = strtok(NULL, sep);
	    strcpy(self[i].resname, token);
	    self[i].resclass= ntvar_getchar(ntvar, self[i].resname);

	    token = strtok(NULL, sep);
	    if(strlen(token) != 1){
		  fprintf(stderr,"Error... ins-code having more than one characters\n");
		  exit(EXIT_FAILURE);
	    }
	    self[i].ins = token[0];

	    token = strtok(NULL, sep);
	    strcpy(self[i].chain, token);

	    /* For base pairs */
	    while((token = strtok(NULL, sep)) != NULL){

		  self[i].numbp++;
		  int j = atoi(token);
		  j --;
		  int pos = graph_set_edge(g, i, j);


		  token = strtok(NULL, sep);

		  token = strtok(NULL, sep);

		  token = strtok(NULL, sep);

		  token = strtok(NULL, sep);

		  token = strtok(NULL, sep);
		  strcpy(self[i].name[pos], token);

		  token = strtok(NULL, sep);
		  strcpy(self[i].type[pos], token);

		  token = strtok(NULL, sep);
		  self[i].eval[pos] = atof(token);
		  graph_set_wt(g, i, pos, self[i].eval[pos]);

	    }
	    fgets(line, 256, fp);

      }
      fclose(fp);
}


int get_numres(char* outfile)
{
      FILE* fp = fopen(outfile, "r");
      assert(fp != NULL);
      if(fp == NULL){    /* Exception Handling */ 
	    fprintf(stderr, "Error in %s at line %d... Unable to open file %s\n", __FILE__, __LINE__, outfile);
	    exit(EXIT_FAILURE);
      }
      char line[256];
      int flag = 0;
      int nres = -1;
      while(fgets(line, 256, fp) != NULL){
	    if(line[0] != '#') break;
	    if(strncmp(line, "#HEADER   Cleaned number", 24) != 0) continue;
	    nres = atoi(line+38);
	    flag = 1;
	    break;
      }
      if(flag == 0){
	    fprintf(stderr, "Error.... Something wrong in out file. No cleaned residue entry found.\n");
	    fclose(fp);
	    exit(EXIT_FAILURE);
      }
      return nres;
}


void print_compo_edge_info(struct nucbp* self, struct graph* g,
	    struct djset* set,  sysparams* syspar, sequence_t* seq, int comp_size, int vertex, FILE* file_adj)
{
      int v1=vertex;
      int v2=vertex;
      int* array = (int*) malloc (comp_size * sizeof(int));
      djset_sort(set, vertex, array);
//      for(int i=0; i<comp_size; ++i){
//	    fprintf(file_adj, "%d, ",array[i]);
//      }
      fprintf(file_adj, "\n");
      fprintf(file_adj, "HEAD    SL    SL DEG  EDG-TYPE      WT  SP  TYP  RES    CHN  RES    CHN  BS  BS\n");
      fprintf(file_adj, "        --    -- ---  --------      --  --  ---  ---    ---  ---    ---  --  --\n");

      for(int i=0; i<comp_size; ++i){
	    v1 = array[i];
	    int deg = g->deg[v1];
	    for(int j=0; j<deg; ++j){
		  v2 = graph_edge_at(g, v1, j);
		  int index = j;

		  char ith_ins[5];
		  ith_ins[1]= self[v1].ins;
		  char jth_ins[5];
		  jth_ins[1]= self[v2].ins;

		  ith_ins[0] = '/';
		  jth_ins[0] = '/';

		  if(ith_ins[1] == '?'){
			ith_ins[1] = ' ';
			ith_ins[0] = ' ';
		  }
		  if(jth_ins[1] == '?'){
			jth_ins[1] = ' ';
			jth_ins[0] = ' ';
		  }

		  ith_ins[2] = '\0';
		  jth_ins[2] = '\0';


		  string bpname = self[v1].name[index];
		  double enerval = graph_get_wt(g, v1, index); 
		  int deg = graph_deg(g, v2);
		  fprintf(file_adj, "EDGE %5d %5d  %2d  %c:%c-%4s %8.2lf %c:%c %s %5d%2s %3s %5d%2s %3s %3s "
			      "%3s\n",
			      v1+1, 
			      v2+1,
			      deg,
			      self[v1].resclass,
			      self[v2].resclass,
			      bpname.c_str(),
			      enerval,
			      seq->sec_seq_str[v1],
			      seq->sec_seq_str[v2],
			      syspar->type.c_str(),
			      self[v1].cifid,
			      ith_ins,
			      self[v1].chain,
			      self[v2].cifid,
			      jth_ins,
			      self[v2].chain,
			      self[v1].resname,
			      self[v2].resname);
	    }
      }
      free(array);
      fprintf(file_adj,"#\n");
}

void print_intro(FILE* adj_file)
{
      fprintf(adj_file, "        +====================================================================+\n");
      fprintf(adj_file, "        +                                                                    +\n");
      fprintf(adj_file, "        +                         B P N E T                                  +\n");
      fprintf(adj_file, "        +                       VERSION - 1-0-2                              +\n");
      fprintf(adj_file, "        +                                                                    +\n");
      fprintf(adj_file, "        +             A Software to Compute Multiplate                       +\n");
      fprintf(adj_file, "        +            Base-Pair Networsk in Nucleic Acids.                    +\n");
      fprintf(adj_file, "        +                                                                    +\n");
      fprintf(adj_file, "        +                                                                    +\n");
      fprintf(adj_file, "        +             BUG  REPORT: roy.parthajit@gmail.com                   +\n");
      fprintf(adj_file, "        +          TECH REPORT: dhananjay.bhattacharyya@saha.ac.in           +\n");
      fprintf(adj_file, "        +====================================================================+\n");
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
}

void print_summary(sysparams* syspar, ovlp_stat* stat, int num_comp, FILE* adj_file)
{
      fprintf(adj_file, "\n---------------------------S U M M A R Y  R E P O R T ---------------------\n");
      fprintf(adj_file, "DATE    %s\n", OvlpGen::today().c_str());
      fprintf(adj_file,"ACCN    %s\n", syspar->accn.c_str());

      fprintf(adj_file, "TOTAL RESIDUE : %d\n", syspar->cleaned_res);
//      fprintf(adj_file, "TOTAL COMPONENTS FOUND : %d\n", syspar->total_per_structure);
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
}

//void print_adjacency(struct rnabp* self, 
//	    struct graph* g,
//	    sysparams* syspar, 
//	    FILE* file_adj,
//	    int num_comp,
//	    int vertex)
//{
//                fprintf(file_adj,"HEAD  TYPE  MDL RES DEG");
//                for(int i=0; i<comp_size; ++i){
//                    fprintf(file_adj, "%6d ", cor_res_no[i]);
//                }
//                fprintf(file_adj,"\n");
//            for(int i=0; i<comp_size; ++i){
//		  fprintf(file_adj, "ADJC   %s %5d  ", syspar->type.c_str(), vertex+1);
//		  fprintf(file_adj, "%c   %1d ",res_class_name[cor_res_no[i] - 1], comp_deg[i]);
//		  int val;
//		  for(int j=0; j<comp_size; ++j){
//			val = graph_edge_index(g, )
//			fprintf(file_adj, "%6d ", );
//		  }
//		  fprintf(file_adj, "\n");
//		  fprintf(file_adj,"#\n");
//	    }
//}
//void print_weight(struct rnabp* self, struct graph* g, int comp_size, sysparams* syspar, FILE* adj_file)
//{
//
//                for(int i=0; i<comp_size; ++i){
//                    fprintf(file_adj,"WEGT   %s %5d        ", syspar->type.c_str(), cor_res_no[0]);
//                    for(int j=0; j<comp_size; ++j){
//                        fprintf(file_adj, "%6.2lf ", adjenergy[cor_res_no[i]-1][cor_res_no[j]-1]);
//                    }
//                    fprintf(file_adj, "\n");
//                }
//                fprintf(file_adj,"TER\n");
//}

void print_adjinfo(struct nucbp* self, int nres, struct djset* set, 
	    ntvariants_t* ntvar, sysparams* syspar, struct graph* g, ovlp_stat* stat)
{

      FILE* adj_file = fopen((syspar->file_dir+syspar->accn+".adj").c_str(), "w");
      print_intro(adj_file);
      syspar->print_params(adj_file);
      int num_comp = set->numset;
      print_summary(syspar, stat,  num_comp, adj_file);
      sequence_t seq;
      sequence_create(&seq, (syspar->file_dir+syspar->accn+".dat").c_str() );
      cout<<"        SECONDARY SEQUENCE GENERATION ENDS SUCCESSFULLY\n";
      int* check = (int*) malloc (nres * sizeof(int));
      for(int i=0; i<nres; ++i){
	    check[i] = 0;
      }

      syspar->total_per_structure = 0; 
      for(int i=0; i<nres; ++i){
	    if(check[i] == 1) continue;

	    int comp_size = djset_composize(set, i);
	    int cycles    = djset_cycles(set, i);
	    int pos = i;
	    int exdeg = 0;
	    int num_exdeg = 0;
	    for(int j=0; j<comp_size; ++j){
		  check[pos] = 1;
		  if(syspar->_exdeg == graph_deg(g, pos)){
			exdeg = 1;
			num_exdeg ++;
		  }
		  pos = djset_next(set, pos);
	    }
	    printf("exdeg=%d, numexdeg=%d, compsize=%d, cicles=%d\n", exdeg, num_exdeg, comp_size, cycles);
	    if(exdeg == 0) continue;
	    if(num_exdeg < syspar->_num_exdeg) continue;
	    if(comp_size < syspar->_from_size || comp_size > syspar->_to_size) continue;
	    if(cycles < syspar->_num_cycles) continue;
	    char comp_name[10];
	    graph_compo_isomorph_name(g, set, i, comp_name);
	    fprintf(adj_file,"COMP   %s %5d  %5s\n", syspar->accn.c_str(), i+1, comp_name);
	    fprintf(adj_file,"CARD   %5d\n",comp_size);
	    print_compo_edge_info(self, g, set,  syspar, &seq, comp_size, i, adj_file);
	    syspar->_total_count ++;
	    syspar->total_per_structure ++; 
      }
      fprintf(adj_file,"END\n");
      free(check);
      fclose(adj_file);
}




#endif// __BPNET_RNABP_H__

//
// Created by parthajit on 23/2/20.
//

#ifndef BPNET_BASIC_03_GRAPH_H
#include <vector>
#include <string>
#include <assert.h>
#include <iostream>
#include <stdio.h>
#include <disjointset.h>
#include <sysparams.h>
#include <sec_seq.h>
#include <ntvariants.h>


using namespace std;

void isompophic_name(int **A, int *deg, int size, char *name){
	strcpy(name, "##");
	if(size == 1){
		strcpy(name,"1G");
	}else if(size == 2){
		strcpy(name,"B1");
	}else if(size == 3){
		int d1 = 0;
		int d2 = 0;
		for(int i=0; i<size; ++i){
			if(deg[i] == 1){
				d1 ++;
			}else if(deg[i] == 2){
				d2 ++;
			}else{
				fprintf(stderr, "Error... Impossible degree incountered\n");
				exit(EXIT_FAILURE);
			}
		}
		if(d1 == 2 && d2 == 1){
			strcpy(name,"T1");
		}else if(d1 == 0 && d2 == 3){
			strcpy(name,"T2");
		}else{
			fprintf(stderr, "Error... Invalid degree sum incountered\n");
			exit(EXIT_FAILURE);
		}
	}else if(size == 4){
		int d1 = 0;
		int d2 = 0;
		int d3 = 0;
		for(int i=0; i<size; ++i){
			if(deg[i] == 1){
				d1 ++;
			}else if(deg[i] == 2){
				d2 ++;
			}else if(deg[i] == 3){
				d3 ++;
			}else{
				fprintf(stderr, "Error... Impossible degree incountered\n");
				exit(EXIT_FAILURE);
			}
		}
		if(d1 == 2 && d2 == 2 && d3 == 0){
			strcpy(name,"Q1");
		}else if(d1 == 3 && d2 == 0 && d3 == 1){
			strcpy(name,"Q2");
		}else if(d1 == 1 && d2 == 2 && d3 == 1){
			strcpy(name,"Q3");
		}else if(d1 == 0 && d2 == 4 && d3 == 0){
			strcpy(name,"Q4");
		}else if(d1 == 0 && d2 == 2 && d3 == 2){
			strcpy(name,"Q5");
		}else if(d1 == 0 && d2 == 0 && d3 == 4){
			strcpy(name,"Q6");
		}else{
			fprintf(stderr, "Error... Invalid degree sum incountered\n");
			exit(EXIT_FAILURE);
		}
	}else if(size == 5){
		int d1 = 0;
		int d2 = 0;
		int d3 = 0;
		int d4 = 0;
		strcpy(name, "P?");
		for(int i=0; i<size; ++i){
			if(deg[i] == 1){
				d1 ++;
			}else if(deg[i] == 2){
				d2 ++;
			}else if(deg[i] == 3){
				d3 ++;
			}else if(deg[i] == 4){
				d4 ++;
			}else{
				fprintf(stderr, "Error... Impossible degree incountered\n");
				exit(EXIT_FAILURE);
			}
		}
		if(d1 == 2 && d2 == 3 && d3 == 0 && d4 == 0){
			strcpy(name,"P1");
		}else if(d1 == 3 && d2 == 1 && d3 == 1 && d4 == 0){
			strcpy(name,"P2");
		}else if(d1 == 1 && d2 == 3 && d3 == 1 && d4 == 0){ /* P3 or P7 case */
			strcpy(name,"P@");
			for(int i=0; i<size; ++i){
				if(deg[i] == 1){
					for(int j=0; j<size; j++){
						if(A[i][j] == 1){ 
							/* check if leaf node is connected to 
							 * degree 2 or degree 3 */
							if(deg[j] == 2){
								strcpy(name, "P3");
							}else{
								strcpy(name, "P7");
							}
						}
					}
					break;
				}
			}
		}else if(d1 == 1 && d2 == 1 && d3 == 3 && d4 == 0){
			strcpy(name,"P4");
		}else if(d1 == 2 && d2 == 1 && d3 == 2 && d4 == 0){
			strcpy(name,"P5");
		}else if(d1 == 0 && d2 == 3 && d3 == 2 && d4 == 0){
			strcpy(name,"P6");
		}else if(d1 == 0 && d2 == 5 && d3 == 0 && d4 == 0){
			strcpy(name,"P8"); /* P7 done earlier */
		}else if(d1 == 0 && d2 == 2 && d3 == 2 && d4 == 1){
			strcpy(name,"P9");
		}else if(d1 == 1 && d2 == 2 && d3 == 1 && d4 == 1){
			strcpy(name,"P10");
		}else if(d1 == 2 && d2 == 2 && d3 == 0 && d4 == 1){
			strcpy(name,"P11");
		}else if(d1 == 4 && d2 == 0 && d3 == 0 && d4 == 1){
			strcpy(name,"P12");
		}else{
			fprintf(stderr, "Error... Invalid degree sum incountered\n");
			exit(EXIT_FAILURE);
		}
	}else if(size == 6){
		strcpy(name, "6#");
		int d1 = 0;
		int d2 = 0;
		int d3 = 0;
		int d4 = 0;
		int d5 = 0;
		for(int i=0; i<size; ++i){
			if(deg[i] == 1){
				d1 ++;
			}else if(deg[i] == 2){
				d2 ++;
			}else if(deg[i] == 3){
				d3 ++;
			}else if(deg[i] == 4){
				d4 ++;
			}else if(deg[i] == 5){
				d5 ++;
			}else{
				fprintf(stderr, "Error... Impossible degree incountered\n");
				exit(EXIT_FAILURE);
			}
		}
		if(d1 == 2 && d2 == 4 && d3+d4+d5 == 0){
			strcpy(name,"S112");
		}else if( d1 == 4 && d3 == 2 && d2+d4+d5 == 0){
			strcpy(name, "S109");
		}else if(d1 == 3 && d2 == 2 && d3 == 1 && d4+d5 == 0){
			strcpy(name,"S111?");
			int twodeg=0;
			for(int i=0; i<size; ++i){
				if(deg[i] == 3){
					for(int j=0; j<size; ++j){
						if(A[i][j] == 1){
							if(deg[j] == 2){
								twodeg ++;
							}
						}
					}
					break;
				}
			}
			if(twodeg == 2){ /* Two Twodeg vertex connected to 3 deg vertex is graph 110 
					    else graph 111*/
				strcpy(name,"S110");
			}else{
				strcpy(name,"S111");
			}
		}else if(d1 == 4 && d3 == 2 && d2+d4+d5 ==0){
			strcpy(name,"S109");
		}
	}else if(size == 7){
	    strcpy(name, "7#");
		int d1 = 0;
		int d2 = 0;
		int d3 = 0;
		int d4 = 0;
		int d5 = 0;
		int d6 = 0;
		for(int i=0; i<size; ++i){
			if(deg[i] == 1){
				d1 ++;
			}else if(deg[i] == 2){
				d2 ++;
			}else if(deg[i] == 3){
				d3 ++;
			}else if(deg[i] == 4){
				d4 ++;
			}else if(deg[i] == 5){
				d5 ++;
			}else if(deg[i] == 6){
				d6 ++;
			}else{
				fprintf(stderr, "Error... Impossible degree incountered\n");
				exit(EXIT_FAILURE);
			}
		}
		if(d1 == 2 && d2 == 5){
			strcpy(name, "7L");
		}
	}else if(size == 8){
	    strcpy(name, "8#");
		int d1 = 0;
		int d2 = 0;
		int d3 = 0;
		int d4 = 0;
		int d5 = 0;
		int d6 = 0;
		int d7 = 0;
		for(int i=0; i<size; ++i){
			if(deg[i] == 1){
				d1 ++;
			}else if(deg[i] == 2){
				d2 ++;
			}else if(deg[i] == 3){
				d3 ++;
			}else if(deg[i] == 4){
				d4 ++;
			}else if(deg[i] == 5){
				d5 ++;
			}else if(deg[i] == 6){
				d6 ++;
			}else if(deg[i] == 7){
				d6 ++;
			}else{
				fprintf(stderr, "Error... Impossible degree incountered\n");
				exit(EXIT_FAILURE);
			}
		}
		if(d1 == 2 && d2 == 6){
			strcpy(name, "8L");
		}
	}else if(size == 9){
	    strcpy(name, "9#");
		int d1 = 0;
		int d2 = 0;
		int dn = 0;
		for(int i=0; i<size; ++i){
			if(deg[i] == 1){
				d1 ++;
			}else if(deg[i] == 2){
				d2 ++;
			}else{
				dn ++;
				//fprintf(stderr, "Error... Impossible degree incountered\n");
				//exit(EXIT_FAILURE);
			}
		}
		if(d1 == 2 && d2 == 7){
			strcpy(name,"9L");
		}
	}else if(size == 10){
	    strcpy(name, "10#");
		int d1 = 0;
		int d2 = 0;
		int dn = 0;
		for(int i=0; i<size; ++i){
			if(deg[i] == 1){
				d1 ++;
			}else if(deg[i] == 2){
				d2 ++;
			}else{
				dn ++;
//				fprintf(stderr, "Error... Impossible degree incountered\n");
//				exit(EXIT_FAILURE);
			}
		}
		if(d1 == 2 && d2 == 8){
			strcpy(name,"10L");
		}
	}else if(size == 11){
	    strcpy(name, "11#");
		int d1 = 0;
		int d2 = 0;
		int dn = 0;
		for(int i=0; i<size; ++i){
			if(deg[i] == 1){
				d1 ++;
			}else if(deg[i] == 2){
				d2 ++;
			}else{
				dn ++;
//				fprintf(stderr, "Error... Impossible degree incountered\n");
//				exit(EXIT_FAILURE);
			}
		}
		if(d1 == 2 && d2 == 9){
			strcpy(name,"11L");
		}
	}else if(size == 12){
	    strcpy(name, "12#");
		int d1 = 0;
		int d2 = 0;
		int dn = 0;
		for(int i=0; i<size; ++i){
			if(deg[i] == 1){
				d1 ++;
			}else if(deg[i] == 2){
				d2 ++;
			}else{
				dn ++;
//				fprintf(stderr, "Error... Impossible degree incountered\n");
//				exit(EXIT_FAILURE);
			}
		}
		if(d1 == 2 && d2 == 10){
			strcpy(name,"12L");
		}
	}else if(size > 12){
		sprintf(name,"%dG",size);
	}
}

#define EMPTY (-1)
#define MAX_DEG (6)
class Graph{
    int V;    // No. of vertices

    // Pointer to an array containing adjacency lists
    int* Adj;
    int* Degree;
    string * adjbp;
    double* adjenergy;
    int* cifid;
    char** chain;
    char* ins_code;
    char** res_act_name;
    char*  res_class_name;  /* A or a ,U or u for modified bases etc */
    int graph_type; /*1=full graph with cifid, residue name etc or 0=graph for conputation*/
    vector<vector<int> > all_components;
    char* valid_comp;
    vector<string> component_name;
    void DFS_visit(int v1, bool visited[], vector<int>& component){
        visited[v1] = true;
        component.push_back(v1);

        // Recur for all the vertices
        // adjacent to this vertex
        for(int i = 0; i < V; i++) {
	      if(get_edge_index(v1, i) >= 0){
                if (!visited[i]) {
                    DFS_visit(i, visited, component);
                }
            }

        }
    }
public:
    int get_edge_index(int i, int j){
	  int adjval;
	  for(int k=0; i<MAX_DEG; ++k){
		adjval = Adj[i*V+k];
		if(adjval == EMPTY) return -1;
		if(Adj[i*V+k] == j) return k;// k is returned for reference to other access like weight;
	  }
	  return -1;
    }
    int set_edge(int i, int j){
	  int pos;
	  for(int k=0; k<MAX_DEG; ++k){
		pos = i*V + k;
		if(Adj[pos] == -1){
		      Adj[pos] = j;
		      return k;
		}
	  }
	  {    /* Exception Handling */ 
		fprintf(stderr, "Falat error in func %s on file %s at line %d... MAX DEGREE overflows.\n", __func__, __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	  }
	  return -1;
    }
    inline double get_weight(int i, int index){
	  return adjenergy[i*V + index];
    }
    inline void set_wetght(int i, int index, double wt){
	  adjenergy[i*V + index] = wt;
	  return;
    }


    inline string get_bp(int i, int index){
	  return adjbp[i*V + index];
    }
    inline void set_bp(int i, int index, string bpname){
	  adjbp[i*V + index] = bpname;
	  return;
    }




    void fprint_component_summary(FILE* fp){
    	int s = (int)all_components.size();
		fprintf(fp, "        SIZE:");
		int cnt = 0;
		for(int i=0 ;i< s; ++i){
			if(valid_comp[i] == 'T'){
				cnt++;
				fprintf(fp, "%d   ",(int)all_components.at(i).size());
			}
		}
		fprintf(fp, "\n");
		fprintf(fp, "        TOTAL COMPONENTS:%d\n",cnt);
    }
    void gen_pymol_basepair_cif(const char* pmlfile, sysparams* syspar){

	  
	  int n = num_components();
	  if(n == 0) return;
	  FILE* fp = fopen(pmlfile, "w");
	  if(fp == NULL){    /* Exception Handling */ 
		fprintf(stderr, "Error in function %s. (File: %s, Line %d)... cannot open file for writing pml script.\n", __func__, __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	  }
	  //fp = stdout;
	  fprintf(fp,"load %s.cif\n",syspar->accn.c_str());
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
	  fprintf(fp, "show cartoon, nucleicall\n");
	  fprintf(fp, "set cartoon_ring_mode, 0, nucleicall\n");
	  fprintf(fp, "set cartoon_ladder_mode, 1, nucleicall\n");
	  fprintf(fp, "show_as cartoon, cartoon\n");

	  for (int i = 0; i < n; ++i) {
		vector<int> component = all_components.at(i);
		int comp_size = (int) component.size();
		if(comp_size < syspar->_from_size || comp_size > syspar->_to_size) continue; 
		if(comp_size == 1 && syspar->is_overlap == "TRUE"){
		      fprintf(stderr, "Worning... In overlap mode an isolated vertex found\n");
		} 
		//if(comp_size == 1) continue;
		int v = component.at(0);
		fprintf(fp, "select comp%d, ",v+1);
		fprintf(fp, "(resi %d and chain %s) ", this->cifid[v], this->chain[v] );
		for(int j=1; j<comp_size; ++j){
		      int v1 = component.at(j);
		      fprintf(fp, "(resi %d and chain %s) ", this->cifid[v1], this->chain[v1] );
		}
		fprintf(fp,"\n");
		fprintf(fp, "util.cbak comp%d\n", v+1);
		fprintf(fp,"\n");
	  }
	  fclose(fp);
    }
    void gen_pymol_cif(const char* pmlfile, sysparams* syspar){
	    char color[30][10] = {"red", "purple","green", "blue","orange", "white"};
	    int numcolor = 6;
	    int n = num_components();
	    if(n == 0) return;
	    FILE* fp = fopen(pmlfile, "w");
	    assert(fp != NULL);
	    //fp = stdout;
	    fprintf(fp,"load %s.cif\n",syspar->accn.c_str());
	    for (int i = 0; i < n; ++i) {
		    vector<int> component = all_components.at(i);
		    int comp_size = (int) component.size();
		    if(comp_size < syspar->_from_size || comp_size > syspar->_to_size) continue; 
		    if(comp_size == 1 && syspar->is_overlap == "TRUE"){
			    fprintf(stderr, "Worning... In overlap mode an isolated vertex found\n");
		    } 
		    //if(comp_size == 1) continue;
		    int v = component.at(0);
		    fprintf(fp, "select comp%d, ",v+1);
		    fprintf(fp, "(resi %d and chain %s) ", this->cifid[v], this->chain[v] );
		    for(int j=1; j<comp_size; ++j){
			    int v1 = component.at(j);
			    fprintf(fp, "(resi %d and chain %s) ", this->cifid[v1], this->chain[v1] );
		    }
		    fprintf(fp,"\n");
		    fprintf(fp, "color %s, comp%d\n", color[i%numcolor], v+1);
		    fprintf(fp, "show  spheres, comp%d\n", v+1);
		    fprintf(fp,"\n");
	    }
	    fprintf(fp,"sele sugbackbone, name P+OP1+OP2+O5*+C5*+C4*+O4*+C3*+O3*+C2*+O2*\n"); 
	    fprintf(fp,"hide spheres, sugbackbone\n");

	    fclose(fp);
    }
    void gen_pymol_cor(const char* pymolpml, sysparams* syspar){
	    char color[30][10] = {"red", "purple","green", "blue","orange", "white"};
	    int numcolor = 6;
	    int n = num_components();
	    if(n == 0) return;
	    FILE* fp = fopen(pymolpml, "w");
	    assert(fp != NULL);
	    //fp = stdout;
	    fprintf(fp,"load %s_rna.pdb\n",syspar->accn.c_str());
	    for (int i = 0; i < n; ++i) {
		    vector<int> component = all_components.at(i);
		    int comp_size = (int) component.size();
		    if(comp_size < syspar->_from_size || comp_size > syspar->_to_size) continue; 
		    if(comp_size == 1 && syspar->is_overlap == "TRUE"){
			    fprintf(stderr, "Worning... In overlap mode an isolated vertex found\n");
		    } 
		    //if(comp_size == 1) continue;
		    int v = component.at(0);
		    fprintf(fp, "select comp%d, ",v+1);
		    fprintf(fp, "resi %d", v+1);
		    for(int j=1; j<comp_size; ++j){
			    int v1 = component.at(j);
			    fprintf(fp, "+%d", v1+1);
		    }
		    fprintf(fp,"\n");
		    fprintf(fp, "color %s, comp%d\n", color[i%numcolor], v+1);
		    fprintf(fp, "show  spheres, comp%d\n", v+1);
		    fprintf(fp,"\n");
	    }
	    fprintf(fp,"show cartoon\n");
	    fprintf(fp,"sele sugbackbone, name P+OP1+OP2+O5*+C5*+C4*+O4*+C3*+O3*+C2*+O2*\n"); 
	    fprintf(fp,"hide spheres, sugbackbone\n");

	    fclose(fp);
    }

    explicit Graph(int V, int graph_type){
        this->graph_type = graph_type;
        this->V = V;
        if(graph_type == 0){
            Adj = new int[V * MAX_DEG];
            Degree = new int[V];
            for (int k = 0; k < V ; ++k) {
                Degree[k] = 0;
            }
            for(int i = 0; i <V; ++i) {
                for(int j=0; j<MAX_DEG; ++j){
		      Adj[i*V + j] = -1;
		}
            }
            adjbp = NULL;
            adjenergy = NULL;
            cifid = NULL;
            chain = NULL;
            ins_code = NULL;
            res_act_name = NULL;
            res_class_name = NULL;
        }else if(graph_type == 1){
            Adj = new int[V* MAX_DEG];
            adjbp = new string[V * MAX_DEG];
            adjenergy = new double[V * MAX_DEG];
            Degree = new int[V];
            cifid = new int[V];
            chain = new char*[V];
            ins_code = new char[V];
            res_act_name = new char*[V];
            res_class_name = new char[V];

            for(int i = 0; i < V; ++i) {
                chain[i] = new char[5];
                res_act_name[i] = new char[5];
            }
            for (int k = 0; k < V ; ++k) {
                Degree[k] = 0;
                cifid[k] = -1;
                chain[k][0] = '\0';
                res_act_name[k][0] = '\0';
                res_class_name[k] = '#';

            }
            for(int i=0; i<V; i++){
                for(int j=0; j<MAX_DEG; j++){
                    Adj[i*V + j] = -1;
		    adjbp[i*V + j] = "";
		    adjenergy[i*V + j] = 0.0;
                }
            }
        }
    }
    ~Graph(){
        // cout<<"Destructor invoked"<<endl;
        delete [] Adj;
        delete [] Degree;
        if(graph_type == 1){
            delete [] adjbp;
            delete [] adjenergy;
            delete [] cifid;
            delete [] chain;
            delete [] res_act_name;
            delete [] res_class_name;
            delete [] ins_code;
            free(valid_comp);
        }

        //cout<<"      MEMORY CLEANUP DONE SUCCESSFULLY"<<endl;
    }

    void addDirectedEdge(int v1, int v2){
	  if( ! (v1>=0 && v1<V && v2>=0 && v2 <V)){    /* Exception Handling */ 
		fprintf(stderr, "Error in %s at line %d... vertex overflow/underflow.\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	  }
        set_edge(v1, v2);
        Degree[v1] = Degree[v1] + 1;
    }
    void add_cif_details(int corid, int cif, const char* cif_chain, char ins, const char* res_name, ntvariants_t*
	ntvar){
        cifid[corid - 1] = cif;
        ins_code[corid-1] = ins;
        strcpy(chain[corid - 1],cif_chain);
        strcpy(this->res_act_name[corid - 1],res_name);
        char base_name = '#';
        if(normal_guanine(res_name) == true){
            base_name = 'G'; /* Upper case */
        }else if(normal_cytosine(res_name) == TRUE){
            base_name = 'C';
        }else if(normal_adenine(res_name) == TRUE){
            base_name = 'A';
        }else if(normal_uracil(res_name) == TRUE){
            base_name = 'U';
        }else if(ntvar_is_guavar(ntvar, res_name) == TRUE){
            base_name = 'g'; /* lower case */
        }else if(ntvar_is_cytvar(ntvar, res_name) == TRUE){
            base_name = 'c';
        }else if(ntvar_is_adevar(ntvar, res_name) == TRUE){
            base_name = 'a';
        }else if(ntvar_is_uravar(ntvar, res_name) == TRUE){
            base_name = 'u';
        }else{
            fprintf(stderr, "Error... Unkonwn residue encountered (%s).\n", res_name);
            //exit(EXIT_FAILURE);
        }
        res_class_name[corid -1] = base_name;
        //printf("%d, %d, %s, %s\n",corid-1,cif, chain[corid-1],res_act_name[corid-1]);

    }

    int count_cycles(){
        Disjoint_Set set = Disjoint_Set(this->V);
        int count = 0;
        for (int i = 0; i < V; ++i) {
            for (int j = i+1; j < V; ++j) {
                if(get_edge_index(i, j) >=0 ){ // if there is edge.
                    int leader_i = set.Find(i);
                    int leader_j = set.Find(j);
                    //cout<<"i="<<leader_i<<", "<<leader_j<<endl;
                    if(leader_i == leader_j){
                        count ++;
                    }
                    set.Union(leader_i, leader_j);
                }
            }
        }
        return count;
    }
    void addBP_and_Evalue(int v1, int v2, string bp, double eval){
        assert(v1>=0 && v1<V && v2>=0 && v2 <V);
	int index = get_edge_index(v1,v2);
        if(get_bp(v1, index) == ""){
	    set_bp(v1, index, bp);
	    set_wetght(v1, index, eval);
        }

//adjbp[v2][v1] = bp;
        //adjenergy[v2][v1] = eval;
    }
    void addBP_and_Evalue_from_rob(int v1, int v2, string bpfrom, string bpto, double eval){
        assert(v1>=0 && v1<V && v2>=0 && v2 <V);
	int index = get_edge_index(v1,v2);
        if(get_bp(v1, index) == ""){
	    set_bp(v1, index, bpfrom);
	    set_wetght(v1, index, eval);
        }
	index = get_edge_index(v2,v1);
        if(get_bp(v2, index) == ""){
	    set_bp(v2, index, bpto);
	    set_wetght(v2, index, eval);
        }
    }
    void addUndirectedEdge(int v1, int v2){
        assert(v1>=0 && v1<V && v2>=0 && v2 <V);
	if(get_edge_index(v1, v2) == -1){
	      set_edge(v1, v2);
	      Degree[v1] = Degree[v1] + 1;
	}
	if(get_edge_index(v2, v1) == -1){
	      set_edge(v2, v1);
	      Degree[v2] = Degree[v2] + 1;
	}
    }
    int get_degree(int v1){
        return Degree[v1];  // i=0 returns cor id 1;
    }
    void computeConnectedComponents(){
        bool *visited = new bool[V];
        for(int i = 0; i < V;i++){
            visited[i] = false;
        }


        for (int v=0; v<V; v++) {
            vector<int> component;
            if (! visited[v] ) {
                DFS_visit(v, visited, component);
                all_components.push_back(component);
            }

        }
        delete [] visited;
        int s = (int)all_components.size();
        valid_comp = (char*) malloc( s * sizeof(char));
        for(int i=0; i<s; ++i){
        	valid_comp[i] = 'F';
        }
    }

    void compute_iso_class(){
        int n = num_components();
        for (int i = 0; i < n; ++i) {
            vector<int> component = all_components.at(i);
            int comp_size = (int) component.size();


        }
    }
    void fprint_component_adjacency(int component_index, FILE* file_adj, FILE* file_edge_list, sequence_t*
    seq, sysparams* syspar){
        vector<int> component = all_components.at(component_index);
        int comp_size = (int) component.size();

        if(!(comp_size <= syspar->_to_size && comp_size >= syspar->_from_size)) return;
        //assert(comp_size<30);
        int** comp_adj;
	int*  comp_deg;
        comp_adj = new int*[comp_size];
        comp_deg = new int[comp_size];
        for (int i = 0; i < comp_size; ++i) {
            comp_adj[i] = new int[comp_size];
        }
        int* cor_res_no = new int[comp_size];
        for(int i=0; i<comp_size; ++i){
            int cor_serial = component.at(i) + 1;
            cor_res_no[i] = cor_serial;
        }
        int mindeg = 0;
        Graph gcomp = Graph(comp_size, 0);  // 0 for Partial graph

        for(int i=0; i<comp_size; ++i){
            int serial_i = component.at(i);
	    int degree = this->get_degree(serial_i);
	    comp_deg[i] = degree;
            if(degree >= syspar->_exdeg){
                mindeg ++;
            }
            for(int j=0; j<comp_size; ++j){
                int serial_j = component.at(j);
		if(get_edge_index(serial_i, serial_j) >=0){
		      comp_adj[i][j] = 1;
		}else{
		      comp_adj[i][j] = 0;
		}
                if(comp_adj[i][j] == 1){
                    gcomp.addUndirectedEdge(i, j);  /* Creating comonent graph for cycle detection*/
                }
            }
        }
        int count = gcomp.count_cycles();
        if(count >0 && syspar->is_overlap == "FALSE"){
            cout<<count<<"                  Cycles found!! at COMP "<<cor_res_no[0]<<" with component size "
																				""<<comp_size<<endl;
        }

	fprintf(file_adj,"\n\n---------------- C O M P O N E N T      I N F O      S T A R T S -----------------\n\n");
        if(mindeg >= syspar->_num_exdeg && count >= syspar->_num_cycles){
            this->valid_comp[component_index] = 'T';
            syspar->_total_count++;
            if(syspar->adj_file == "TRUE"){
		    char comp_name[10];
				isompophic_name(comp_adj, comp_deg, comp_size, comp_name);
                fprintf(file_adj,"COMP   %s %5d  %5s\n", syspar->accn.c_str(), cor_res_no[0], comp_name);
                fprintf(file_adj,"CARD   %5d\n",comp_size);
            }

            fprintf(file_edge_list,"COMP   %s %5d\n", syspar->accn.c_str(), cor_res_no[0]);
            if(mindeg>2){
                if(syspar->adj_file == "TRUE"){
                    fprintf(file_adj,"INFO    DEGREE %2d COUNT %2d\n",syspar->_exdeg, mindeg);
                }
                fprintf(file_edge_list,"INFO   DEGREE %2d COUNT %2d\n",syspar->_exdeg, mindeg);
            }
            if(count>0){
                if(syspar->adj_file == "TRUE"){
                    fprintf(file_adj,"INFO   CYCLE %5d\n",count);
                }
                fprintf(file_edge_list,"INFO   CYCLE %5d\n",count);
            }
            fprintf(file_edge_list,"META CARD %5d\n",comp_size);
	    fprintf(file_adj, "HEAD    SL    SL  EDG-TYPE      WT  SP  TYP  RES    CHN  RES    CHN  BS  BS\n");
	    fprintf(file_adj, "        --    --  --------      --  --  ---  ---    ---  ---    ---  --  --\n");
            for(int i=0; i<comp_size; ++i){
                for(int j=0; j<comp_size; ++j){
//                    if(syspar->adj_file == "TRUE"){
//                        fprintf(file_adj, "%6d ", comp_adj[i][j]);
//                    }
                    if(j>=i && comp_adj[i][j] != 0){
                    	char ith_ins[5];
                    	ith_ins[1]= ins_code[cor_res_no[i]-1];
						char jth_ins[5];
						jth_ins[1]= ins_code[cor_res_no[j]-1];

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


						string bpname = get_bp(cor_res_no[i]-1, cor_res_no[j]-1);
						double enerval = get_weight(cor_res_no[i]-1, cor_res_no[j]-1);
                        fprintf(file_adj, "EDGE %5d %5d  %c:%c-%4s %8.2lf %c:%c %s %5d%2s %3s %5d%2s %3s %3s "
												"%3s\n",
                                cor_res_no[i], cor_res_no[j],
                                res_class_name[cor_res_no[i]-1],
                                res_class_name[cor_res_no[j]-1],
                                bpname.c_str(),
                                enerval,
				seq->sec_seq_str[cor_res_no[i] -1 ],
                                seq->sec_seq_str[cor_res_no[j] -1],
                                syspar->type.c_str(),
                                cifid[cor_res_no[i] -1],
								ith_ins,
								chain[cor_res_no[i] -1],
                                cifid[cor_res_no[j] -1],
								jth_ins,
								chain[cor_res_no[j] -1],
                                res_act_name[cor_res_no[i] -1],
								res_act_name[cor_res_no[j] -1]);
                    }
                }
//                if(syspar->adj_file == "TRUE"){
//		    fprintf(file_adj, "\n");
//                }

            }
	    fprintf(file_adj,"#\n");
            if(syspar->adj_file == "TRUE"){
                fprintf(file_adj,"HEAD  TYPE  MDL RES DEG");
                for(int i=0; i<comp_size; ++i){
                    fprintf(file_adj, "%6d ", cor_res_no[i]);
                }
                fprintf(file_adj,"\n");
            }
            for(int i=0; i<comp_size; ++i){
                if(syspar->adj_file == "TRUE") {
                    fprintf(file_adj, "ADJC   %s %5d  ", syspar->type.c_str(), cor_res_no[0]);
                }
		if(syspar->adj_file == "TRUE"){
			fprintf(file_adj, "%c   %1d ",res_class_name[cor_res_no[i] - 1], comp_deg[i]);
		}
                for(int j=0; j<comp_size; ++j){
                    if(syspar->adj_file == "TRUE"){
                        fprintf(file_adj, "%6d ", comp_adj[i][j]);
                    }
		}
                if(syspar->adj_file == "TRUE"){
		    fprintf(file_adj, "\n");
                }
	    }
            if(syspar->adj_file == "TRUE"){
                fprintf(file_adj,"#\n");
                for(int i=0; i<comp_size; ++i){
                    fprintf(file_adj,"WEGT   %s %5d        ", syspar->type.c_str(), cor_res_no[0]);
                    for(int j=0; j<comp_size; ++j){
                        fprintf(file_adj, "%6.2lf ", adjenergy[(cor_res_no[i]-1)*V + (cor_res_no[j]-1)]);
                    }
                    fprintf(file_adj, "\n");
                }
                fprintf(file_adj,"TER\n");
                fprintf(file_edge_list,"TER\n");
            }



            delete[] cor_res_no;
            for (int i = 0; i < comp_size; ++i) {
                delete [] comp_adj[i];
            }
            delete [] comp_adj;
	    delete [] comp_deg;
        }

    }


    int num_components(){
        return (int)all_components.size();
    }
    vector<int>& get_component(int i){
        return all_components[i];
    }
    int get_num_vertex(){
        return V;
    }


};
#define BPNET_BASIC_03_GRAPH_H

#endif //BPNET_BASIC_03_GRAPH_H

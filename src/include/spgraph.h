//
// Created by parthajit on 23/2/20.
//

#ifndef BPNET_BASIC_03_SPGRAPH_H
#define BPNET_BASIC_03_SPGRAPH_H
#include <string>
#include <assert.h>
#include <iostream>
#include <stdio.h>
#include "djset.h"
#include <sysparams.h>
#include <sec_seq.h>
#include <ntvariants.h>

#define EMPTY (-1)
#define MAX_DEGS (10)
enum gtype {DIRECTED =0, UNDIRECTED=1};






struct graph{
      int v;
      int* adj;
      int* deg;
      double* w;
      double memsize;
      enum gtype type;
};
void graph_free(struct graph* self)
{
      free(self->adj);   
      free(self->deg); 
      free(self->w);   
}

void graph_init(struct graph* self, int vnum, enum gtype type)
{
      int sz = MAX_DEGS * vnum;
      self->memsize = sz;
      self->v   = vnum;
      self->adj =    (int*) malloc(sz * sizeof(int));
      self->w   = (double*) malloc(sz  * sizeof(double));
      self->deg =    (int*) malloc(vnum  * sizeof(int));
      self->type = type;

      for (int k = 0; k < vnum ; ++k) {
	    self->deg[k] = 0;
      }

//      for(int i = 0; i < sz; ++i) {
//		  self->adj[i] = -1;
//      }
}

int graph_edge_at(struct graph* self, int vertex, int index){
      return self->adj[MAX_DEGS * vertex + index];
}


int graph_edge_index(struct graph* self, int i, int j)
{
      int val;
      int start = MAX_DEGS * i;
      int end   = start + self->deg[i];
      for(int k=start; k<end; ++k){
	    val = self->adj[k];
	    if(val == j) return k - start;// Returned for reference to other access like weight;
      }
      return -1;
}

int graph_set_edge(struct graph* self, int i, int j){
      if(self->deg[i] >= MAX_DEGS){    /* Exception Handling */ 
	    fprintf(stderr, "Falat error in func %s on file %s at line %d... MAX DEGREE overflows.\n",
			__func__, __FILE__, __LINE__);
	    exit(EXIT_FAILURE);
      }
      int index = MAX_DEGS * i + self->deg[i];
      self->adj[index] = j;
      self->deg[i] += 1;
      return (self->deg[i] - 1);
}

double graph_get_wt(struct graph* self, int i, int index)
{
      return self->w[i*MAX_DEGS + index];
}

void graph_set_wt(struct graph* self, int i, int index, double wt)
{
      self->w[i*MAX_DEGS + index] = wt;
      return;
}

int graph_deg(struct graph* self, int vertex)
{
      return self->deg[vertex];
}

void graph_kruskal_component(struct graph* self, struct djset* set)
{
      int vertex;
      for(int i=0; i<self->v; ++i){
	    int index = MAX_DEGS * i;
	    for(int j=0; j<self->deg[i]; ++j){
		  vertex = self->adj[index + j];
		  if(vertex > i){
			djset_union(set, i, vertex);
		  }
	    }
      }
}

void graph_compo_isomorph_name(struct graph* self, struct djset* set, int vertex, char* name){
      int size = djset_composize(set, vertex);


	strcpy(name, "##");
	if(size == 1){
		strcpy(name,"1G");
	}else if(size == 2){
		strcpy(name,"B1");
	}else if(size == 3){
		int d1 = 0;
		int d2 = 0;
		int v = vertex;
		for(int i=0; i<size; ++i){
			if(self->deg[v] == 1){
				d1 ++;
			}else if(self->deg[v] == 2){
				d2 ++;
			}else{
				fprintf(stderr, "Error... Impossible degree incountered\n");
				exit(EXIT_FAILURE);
			}
			v = djset_next(set, v);
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
		int v = vertex;
		for(int i=0; i<size; ++i){
			if(self->deg[v] == 1){
				d1 ++;
			}else if(self->deg[v] == 2){
				d2 ++;
			}else if(self->deg[v] == 3){
				d3 ++;
			}else{
				fprintf(stderr, "Error... Impossible degree incountered\n");
				exit(EXIT_FAILURE);
			}
			v = djset_next(set, v);
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
		int v = vertex;
		for(int i=0; i<size; ++i){
			if(self->deg[v] == 1){
				d1 ++;
			}else if(self->deg[v] == 2){
				d2 ++;
			}else if(self->deg[v] == 3){
				d3 ++;
			}else if(self->deg[v] == 4){
				d4 ++;
			}else{
				fprintf(stderr, "Error... Impossible degree incountered\n");
				exit(EXIT_FAILURE);
			}
			v = djset_next(set, v);
		}
		if(d1 == 2 && d2 == 3 && d3 == 0 && d4 == 0){
			strcpy(name,"P1");
		}else if(d1 == 3 && d2 == 1 && d3 == 1 && d4 == 0){
			strcpy(name,"P2");
		}else if(d1 == 1 && d2 == 3 && d3 == 1 && d4 == 0){ /* P3 or P7 case */
			strcpy(name,"P@");
			int u = vertex;
			for(int i=0; i<size; ++i){
				if(self->deg[u] == 1){
				      int v = djset_next(set, u); 
					for(int j=1; j<size; j++){
						if(graph_edge_index(self, u, v) >= 0){ 
							/* check if leaf node is connected to 
							 * degree 2 or degree 3 */
							if(self->deg[v] == 2){
								strcpy(name, "P3");
							}else{
								strcpy(name, "P7");
							}
						}
						v = djset_next(set, v);
					}
					break;
				}
				u = djset_next(set, u);
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
		int v = vertex;
		for(int i=0; i<size; ++i){
			if(self->deg[v] == 1){
				d1 ++;
			}else if(self->deg[v] == 2){
				d2 ++;
			}else if(self->deg[v] == 3){
				d3 ++;
			}else if(self->deg[v] == 4){
				d4 ++;
			}else if(self->deg[v] == 5){
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
			int u = vertex;
			for(int i=0; i<size; ++i){
				if(self->deg[u] == 3){
				      int v = djset_next(set, u);
					for(int j=1; j<size; ++j){
					      if(graph_edge_index(self, u, v) >= 0){
							if(self->deg[v] == 2){
								twodeg ++;
							}
						}
						v = djset_next(set, v);
					}
					break;
				}
				u = djset_next(set, u);
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
		int v = vertex;
		for(int i=0; i<size; ++i){
			if(self->deg[v] == 1){
				d1 ++;
			}else if(self->deg[v] == 2){
				d2 ++;
			}else if(self->deg[v] == 3){
				d3 ++;
			}else if(self->deg[v] == 4){
				d4 ++;
			}else if(self->deg[v] == 5){
				d5 ++;
			}else if(self->deg[v] == 6){
				d6 ++;
			}else{
				fprintf(stderr, "Error... Impossible degree incountered\n");
				exit(EXIT_FAILURE);
			}
			v = djset_next(set, v);
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
		int v = vertex;
		for(int i=0; i<size; ++i){
			if(self->deg[v] == 1){
				d1 ++;
			}else if(self->deg[v] == 2){
				d2 ++;
			}else if(self->deg[v] == 3){
				d3 ++;
			}else if(self->deg[v] == 4){
				d4 ++;
			}else if(self->deg[v] == 5){
				d5 ++;
			}else if(self->deg[v] == 6){
				d6 ++;
			}else if(self->deg[v] == 7){
				d6 ++;
			}else{
				fprintf(stderr, "Error... Impossible degree incountered\n");
				exit(EXIT_FAILURE);
			}
			v = djset_next(set, v);
		}
		if(d1 == 2 && d2 == 6){
			strcpy(name, "8L");
		}
	}else if(size == 9){
	    strcpy(name, "9#");
		int d1 = 0;
		int d2 = 0;
		int dn = 0;
		int v = vertex;
		for(int i=0; i<size; ++i){
			if(self->deg[v] == 1){
				d1 ++;
			}else if(self->deg[v] == 2){
				d2 ++;
			}else{
				dn ++;
				//fprintf(stderr, "Error... Impossible degree incountered\n");
				//exit(EXIT_FAILURE);
			}
			v = djset_next(set, v);
		}
		if(d1 == 2 && d2 == 7){
			strcpy(name,"9L");
		}
	}else if(size == 10){
	    strcpy(name, "10#");
		int d1 = 0;
		int d2 = 0;
		int dn = 0;
		int v = djset_next(set, v);
		for(int i=0; i<size; ++i){
			if(self->deg[v] == 1){
				d1 ++;
			}else if(self->deg[v] == 2){
				d2 ++;
			}else{
				dn ++;
//				fprintf(stderr, "Error... Impossible degree incountered\n");
//				exit(EXIT_FAILURE);
			}
			v = djset_next(set, v);
		}
		if(d1 == 2 && d2 == 8){
			strcpy(name,"10L");
		}
	}else if(size == 11){
	    strcpy(name, "11#");
		int d1 = 0;
		int d2 = 0;
		int dn = 0;
		int v = vertex;
		for(int i=0; i<size; ++i){
			if(self->deg[v] == 1){
				d1 ++;
			}else if(self->deg[v] == 2){
				d2 ++;
			}else{
				dn ++;
//				fprintf(stderr, "Error... Impossible degree incountered\n");
//				exit(EXIT_FAILURE);
			}
			v = djset_next(set, v);
		}
		if(d1 == 2 && d2 == 9){
			strcpy(name,"11L");
		}
	}else if(size == 12){
	    strcpy(name, "12#");
		int d1 = 0;
		int d2 = 0;
		int dn = 0;
		int v = vertex;
		for(int i=0; i<size; ++i){
			if(self->deg[v] == 1){
				d1 ++;
			}else if(self->deg[v] == 2){
				d2 ++;
			}else{
				dn ++;
//				fprintf(stderr, "Error... Impossible degree incountered\n");
//				exit(EXIT_FAILURE);
			}
			v = djset_next(set, v);
		}
		if(d1 == 2 && d2 == 10){
			strcpy(name,"12L");
		}
	}else if(size > 12){
		sprintf(name,"%dG",size);
	}
}
#endif //BPNET_BASIC_03_SPGRAPH_H

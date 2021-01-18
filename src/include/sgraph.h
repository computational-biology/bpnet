//
// Created by parthajit on 7/10/20.
//

#ifndef BPNET_BASIC_03_SGRAPH_H
#define BPNET_BASIC_03_SGRAPH_H

#include<stdio.h>
#include<stdlib.h>
#define MAXDEG 6
#define TRUE 1
#define FALSE 0

typedef struct{
	size_t** adj;
	int v;
	double* weight;
	int isdir;
	int iswtd;
}spgraph_t;

void spgraph_init(spgraph_t* g, int v, int is_directed, int is_weighted){
	g->v   = v;
	g->isdir = is_directed;
	g->iswtd = is_weighted;
	adj = (size_t**) malloc(MAXDEG * sizeof(size_t));
	for(int i=0; i<MAXDEG; ++i){
		g->adj[i] = (int*) malloc(v * sizeof(int));
	}
	for(int i=0; i<MAXDEG; ++i){
		g->wt[i]  = (double*) malloc(v * sizeof(double));
	}
	g->directed = is_directed;
}

void spgraph_setedge(int v1, int v2, double weight){
	;
}

void spgraph_free(spgraph_t* g){
	for(int i=0; i<MAXDEG; ++i){
		free(g->adj[i]);
	}
	for(int i=0; i<MAXDEG; ++i){
		free(g->wt[i]);
	}
}





#endif //BPNET_BASIC_03_SGRAPH_H

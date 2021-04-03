//
// Created by parthajit on 23/2/20.
//

#ifndef BPNET_BASIC_03_DJSET_H
#define BPNET_BASIC_03_DJSET_H

#include <stdlib.h>
#include <stdio.h>


struct djset{
      int* leader;
      int* rank;
      int* link;
      int* cardinal;
      int* cycles;
      int size;
      int numset;

};

int intcmp (const void* a, const void* b) {
      return ( *(int*)a - *(int*)b );
}

void djset_sort(struct djset* self, int elem, int* array){
      int start = elem;
      static int flag =0;
      int i = 0;
      do{
	    array[i] = start;

	    i++;
	    start = self->link[start];
      }while(start != elem);
      qsort(array, i, sizeof(int), intcmp);
}
void djset_init(struct djset* self, int size)
{
      self->size   = size;
      self->numset = size;

      self->leader     = (int*) malloc (size * sizeof(int));
      self->rank       = (int*) malloc (size * sizeof(int));
      self->link       = (int*) malloc (size * sizeof(int));
      self->cardinal   = (int*) malloc (size * sizeof(int));
      self->cycles   = (int*) malloc (size * sizeof(int));

      for (int i = 0; i < size; ++i) {
	    self->leader[i] = i;
      }
      for (int i = 0; i < size; ++i) {
	    self->link[i] = i;
      }
      for (int i = 0; i < size; ++i) {
	    self->rank[i] = 0;
      }
      for (int i = 0; i < size; ++i) {
	    self->cardinal[i] = 1;
      }
      for (int i = 0; i < size; ++i) {
	    self->cycles[i] = 0;
      }
}
void djset_free(struct djset* self)
{
      free(self->link);
      free(self->rank);
      free(self->leader);
      free(self->cardinal);
      free(self->cycles);
}

int djset_find(struct djset* self, int v)
{
      while(self->leader[v] != v){
	    v = self->leader[v];
      }
      return v;
}

void djset_link(struct djset* self, int x, int y)
{
      self->numset --;
      int tmp = self->link[y];     
      self->link[y] = self->link[x];
      self->link[x] = tmp;
      if(self->rank[x] > self->rank[y]){
	    self->leader[y] = x;
	    self->cardinal[x] += self->cardinal[y];
	    self->cycles[x]   += self->cycles[y];
      }else{
	    self->leader[x] = y;
	    self->cardinal[y] += self->cardinal[x];
	    self->cycles[y]   += self->cycles[x];
	    if(self->rank[x] == self->rank[y]){
		  self->rank[y] ++;
	    }
      }
}
void djset_union(struct djset* self, int u, int v)
{
      int lx = djset_find(self, u);
      int ly = djset_find(self, v);
      if(lx != ly){
	    djset_link(self, lx, ly);
      }else{
	    self->cycles[lx] ++;
      }
}
int djset_cycles(struct djset* self, int vertex){
      int leader = djset_find(self, vertex);
      return self->cycles[leader];
}
int djset_composize(struct djset* self, int vertex){
      int leader = djset_find(self, vertex);
      return self->cardinal[leader];
}
int djset_next(struct djset* self, int vertex){
      return self->link[vertex];
}

#endif //BPNET_BASIC_03_DJSET_H

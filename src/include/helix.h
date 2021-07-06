/*
 * =====================================================================================
 *
 *       Filename:  helix.h
 *
 *    Description:  This header file manipulates hlx file. 
 *
 *        Version:  1.0
 *        Created:  27/06/21 01:34:03 PM IST
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  PARTHAJIT ROY (PR), roy.parthajit@gmail.com
 *   Organization:  The University of Burdwan
 *
 * =====================================================================================
 */



#ifndef  __helix_H__
#define  __helix_H__

#include <stdlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>

#include <rnabp.h>
#include <spgraph.h>
#include <stack.h>
#include <sysparams.h>


#define TRUR 1
#define FALSE 0


#define HELIX_MAX 30
#define HAIRPIN_LOOP_MAX 6

struct helix{
      int i[HELIX_MAX];
      int j[HELIX_MAX];
      int size;
      int is_hairpin;
      int pinsize;
      char loop_res[HAIRPIN_LOOP_MAX];
      char loopname[64];
};				/* ----------  end of struct helix  ---------- */

struct pseudo_helix{
      struct helix hlx;
      struct pseudo_helix* branch[6];
      int count;
};


struct pymol_param{
      char hlxcol[20];
      char loopcol[20];
      char psuhlxcol[20];
      char nonhlxcol[20];
      char psujunccol[20];
      int psudocount;
      int loopcount;
      int helixcount;
};

void helix_fprint(struct helix* self, struct nucbp* nbp, FILE* fp){
      if(self->is_hairpin == TRUE){
	    fprintf(fp, "TYPE  STEM-LOOP  with (%s) loop structure.\n", self->loopname);
      }else{
	    fprintf(fp, "TYPE  PURE HELIX structure.\n");
      }
      for(int k=0; k<self->size; ++k){
	    if(self->is_hairpin == TRUE){
		  fprintf(fp, "STEM  ");
	    }else{
		  fprintf(fp, "HELX  ");
	    }
	    fprintf(fp, "%c %c  ", nbp[self->i[k]].resclass, nbp[self->j[k]].resclass);
	    fprintf(fp, "%5d %-4s", nbp[self->i[k]].cifid,nbp[self->i[k]].chain) ;
	    fprintf(fp, "%5d %-4s", nbp[self->j[k]].cifid,nbp[self->j[k]].chain) ;
	    fprintf(fp, "\n");
      }
      if(self->is_hairpin == TRUE){
	    fprintf(fp, "LOOP  ");
	    for(int i=0; i<self->pinsize; ++i){
		  fprintf(fp, "%c",self->loop_res[i]);
	    }
	    fprintf(fp, "\n");
      }
      fprintf(fp, "#\n");
}

void helix_fprint1(struct helix* self, struct nucbp* nbp, FILE* fp){
      fprintf(fp, "%4s.%-6d ", nbp[self->i[0]].chain,nbp[self->i[0]].cifid) ;
      for(int k=0; k<self->size; ++k){
	    fprintf(fp, "%c", nbp[self->i[k]].resclass);
      }
      fprintf(fp, " %4s.%-6d    helix size=%d", nbp[self->i[self->size-1]].chain,nbp[self->i[self->size-1]].cifid, self->size) ;
      if(self->is_hairpin == TRUE){
	    fprintf(fp, ", %s, ", self->loopname);
	    for(int i=0; i<self->pinsize; ++i){
		  fprintf(fp, "%c",self->loop_res[i]);
	    }
      }
      fprintf(fp, "\n");

      fprintf(fp, "            ");
      for(int k=0; k<self->size; ++k){
	    fprintf(fp, "|");
      }
      if(self->is_hairpin == TRUE){
		  fprintf(fp, " )");
//	    for(int k=0; k<self->pinsize/2; ++k){
//		  fprintf(fp, ":");
//	    }
//	    if(self->pinsize % 2 == 1){
//		  fprintf(fp, ")");
//	    }
      }
      fprintf(fp, "\n");
      fprintf(fp, "%4s.%-6d ", nbp[self->j[0]].chain,nbp[self->j[0]].cifid) ;
      for(int k=0; k<self->size; ++k){
	    fprintf(fp, "%c", nbp[self->j[k]].resclass);
      }
      fprintf(fp, " %4s.%-6d", nbp[self->j[self->size-1]].chain,nbp[self->j[self->size-1]].cifid) ;
      fprintf(fp, "\n\n");

}

int is_paired(struct nucbp* nbp, int i, int j){
      for(int k=0; k<nbp[i].numbp; ++k){
	    if(nbp[i].oth_base_index[k] == j) return 1;
      }
      return 0;
}

int is_next_resi(struct nucbp* nbp, int i, int nexti){
      if(nexti-i > 1) return FALSE;
      if(strcmp(nbp[i].chain, nbp[nexti].chain) != 0) return FALSE;
      else if(nbp[i].cifid < nbp[nexti].cifid) return TRUE;
      else if(nbp[i].cifid == nbp[nexti].cifid && nbp[i].ins <  nbp[nexti].ins) return TRUE;
      else return FALSE;
}

int is_prev_resi(struct nucbp* nbp, int i, int previ){
      if(i-previ >1) return FALSE;
      if(strcmp(nbp[i].chain, nbp[previ].chain) != 0) return FALSE;
      else if(nbp[i].cifid >  nbp[previ].cifid) return TRUE;
      else if(nbp[i].cifid == nbp[previ].cifid && nbp[i].ins > nbp[previ].ins) return TRUE;
      else return FALSE;
}




/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  helix_comp
 *  Description:  This function computes pure helix from
 *  		  rnabp data base.  
 * =====================================================================================
 */

void helix_calc(struct helix* self, struct nucbp* nbp, struct graph* g, int i, int j){
      /* The base case */
      if(graph_edge_index(g, i, j) != -1) return;
      if(is_paired(nbp, i,j) == FALSE) return;
      else{
	    graph_set_edge(g, i, j);
	    graph_set_edge(g, j, i);
	    if(self->size == HELIX_MAX){    /* Exception Handling */ 
		  fprintf(stderr, "Error in function %s. (File: %s, Line %d)... Helix Overflow Encountered.\n", __func__, __FILE__, __LINE__);
		  exit(EXIT_FAILURE);
	    }	       
	    self->i[self->size]=i;
	    self->j[self->size]=j;
	    self->size ++;
	    if(is_next_resi(nbp, i, i+1) == TRUE && is_prev_resi(nbp, j, j-1) == TRUE){
		  helix_calc(self, nbp, g, i+1, j-1);
	    }
      }
      return;
}		/* -----  end of function helix_calc  ----- */


struct pseudo_helix* pseudo_helix_getnode(){
      struct pseudo_helix* node;
      node = (struct pseudo_helix*) malloc (sizeof(struct pseudo_helix));
      node->hlx.size = 0;
      node->count = 0;
      for(int i=0; i<6; ++i){
	    node->branch[i] = NULL;
      }
      return node;
}
void pseudo_helix_free(struct pseudo_helix* self){
      for(int i=0; i<self->count; ++i){
	    pseudo_helix_free(self->branch[i]);
      }
      free(self);
      return;
}
struct pseudo_helix* pseudo_helix_calc(struct nucbp* nbp, struct graph* g, int i, int j){
      struct pseudo_helix * self = pseudo_helix_getnode();
      struct helix* h = &self->hlx;
      helix_calc(h, nbp, g, i, j);
      if(self->hlx.size < 3){
	    self->hlx.size = 0;
	    free(self);
	    return NULL;
      }
      
      int  i_last = h->i[h->size - 1];
      int  j_last = h->j[h->size - 1];
      if(is_next_resi(nbp, i_last, i_last+1) == TRUE){
	    for(int k=0; k<nbp[i_last + 1].numbp; ++k){
		  struct pseudo_helix* tmphlx = pseudo_helix_calc(nbp, g, i_last+1, nbp[i_last +1].oth_base_index[k]);
		  if(tmphlx != NULL){
			self->branch[self->count] = tmphlx;
			self->count ++;
		  }
	    }
      }

      if(is_prev_resi(nbp, j_last, j_last - 1) == TRUE){
	    for(int k=0; k<nbp[j_last - 1].numbp; ++k){
		  struct pseudo_helix* tmphlx = pseudo_helix_calc(nbp, g, nbp[j_last - 1].oth_base_index[k], j_last - 1);
		  if(tmphlx != NULL){
			self->branch[self->count] = tmphlx;
			self->count ++;
		  }
	    }
      }
      return self;
}

void helix_pseudo_fprint(struct pseudo_helix* self, struct nucbp* nbp, struct stack* stk, FILE* fp){
      stack_push(stk, self);
      if(self->count == 2){
	    fprintf(fp, "More than one branch found...\n");
      }
      if(self->count == 0){
	    struct pseudo_helix* tmphlx = (pseudo_helix*) stk->elem[0];
//            struct helix* f_elem = &tmphlx->hlx;
//
//	    for(int i=0; i<=stk->top; ++i){
//		  tmphlx = (pseudo_helix*) stk->elem[i];
//		  fprintf(fp, "%5d~", tmphlx->hlx.i[0] + 1) ;
//		  for(int k=0; k<tmphlx->hlx.size; ++k){
//			fprintf(fp, " ");
//		  }
//	    }
//	    fprintf(fp,"\n");
	    for(int i=0; i<=stk->top; ++i){
		  tmphlx = (pseudo_helix*) stk->elem[i];
		  for(int k=0; k<tmphlx->hlx.size; ++k){
			fprintf(fp,"PHLX  ");
			fprintf(fp, "%5d:%-5d  ", tmphlx->hlx.i[k]+1, tmphlx->hlx.j[k]+1);
			fprintf(fp, "%c=%c", nbp[tmphlx->hlx.i[k]].resclass, nbp[tmphlx->hlx.j[k]].resclass);
			fprintf(fp, "  %5d %-4s", nbp[tmphlx->hlx.i[k]].cifid,nbp[tmphlx->hlx.i[k]].chain);
			fprintf(fp, "  %5d %-4s", nbp[tmphlx->hlx.j[k]].cifid,nbp[tmphlx->hlx.j[k]].chain);
			if(i<stk->top && k == tmphlx->hlx.size -1){
			      fprintf(fp," <---Break\n");
			}else{
			      fprintf(fp,"\n");
			}
		  }
	    }
	    fprintf(fp,"#\n\n");

	    
//	    for(int i=0; i<=stk->top; ++i){
//		  tmphlx = (pseudo_helix*) stk->elem[i];
//		  for(int k=0; k<tmphlx->hlx.size; ++k){
//			fprintf(fp, "|");
//		  }
//		  fprintf(fp,"~");
//	    }
//	    fprintf(fp,"\n");
//
//	    
//	    for(int i=0; i<=stk->top; ++i){
//		  tmphlx = (pseudo_helix*) stk->elem[i];
//		  for(int k=0; k<tmphlx->hlx.size; ++k){
//			fprintf(fp, "%c", nbp[tmphlx->hlx.j[k]].resclass);
//		  }
//		  fprintf(fp,"~");
//	    }
      }

	    
     

      
      for(int k=0; k<self->count; ++k){
	    helix_pseudo_fprint(self->branch[k], nbp, stk, fp);
      }
      
      
      stack_pop(stk);

      
}


void helix_pseudo_calc_all(struct pseudo_helix* all_pseudo_helix[], int* counter, struct nucbp* nbp, int nres)
{

      int cnt = 0;
      struct graph g;
      graph_init(&g, nres, UNDIRECTED); 
      for(int i=0; i<nres; ++i){
	    for(int k=0; k<nbp[i].numbp; ++k){
		  int j = nbp[i].oth_base_index[k];
		  if(graph_edge_index(&g, i, j) != -1) continue;

		  struct pseudo_helix* tmp = pseudo_helix_calc(nbp, &g, i, j);
		  if(tmp != NULL){
			if(tmp->count == 0) {
			      free(tmp);
			}else{
			      all_pseudo_helix[cnt] = tmp;
			      cnt ++;
			}

		  }
	    }
      }
      *counter = cnt;
      graph_free(&g);
}

void helix_pseudo_fprint_all(struct pseudo_helix* all_pseudo_helix[], int counter, struct nucbp* nbp, int nres, FILE* fp){
      struct stack stk;
      stack_init(&stk, nres/10);

      
      for(int i=0; i<counter; ++i){
	    helix_pseudo_fprint(all_pseudo_helix[i], nbp, &stk, fp);
      }

      
      stack_free(&stk);
}
void helix_init_all(struct helix** self, int maxsize){
      *self = (struct helix*) malloc ( maxsize * sizeof(struct helix) );
      if ( *self == NULL ) {
	    fprintf ( stderr, "\nDynamic memory allocation failed\n" );
	    exit (EXIT_FAILURE);
      }
}

void get_loop_name(char* seq, int len, char* name){
      if(len == 4){
	    if(toupper(seq[0]) == 'G' 
			&&
	       (toupper(seq[2]) == 'G' || toupper(seq[2]) == 'A') 
	                &&
	        toupper(seq[3]) == 'A'){
		  strcpy(name, "GNRA Tetra");
	    }else if(toupper(seq[0]) == 'U' 
			&&
	             toupper(seq[2]) == 'C'
	                &&
	             toupper(seq[3]) == 'G'){
		  strcpy(name, "UNCG Tetra");
	    }else{
		  strcpy(name,"Tetra");
	    }
      }else{
	    name[0] = '\0';
      }
}

int helix_check_hairpin(struct helix* self, struct nucbp* nbp){
      const int loop_size = HAIRPIN_LOOP_MAX;
      self->is_hairpin = FALSE;
      self->pinsize = 0;
      char* ich = nbp[self->i[self->size-1]].chain;
      char* jch = nbp[self->j[self->size-1]].chain;

      /*
      int icif = nbp[self->i[self->size-1]].cifid;
      int jcif = nbp[self->j[self->size-1]].cifid;

      int iins = nbp[self->i[self->size-1]].ins;
      int jins = nbp[self->j[self->size-1]].ins;
      */
      if(strcmp(ich, jch) != 0) return FALSE;

      int pin_size = abs(self->i[self->size-1] - self->j[self->size-1]) - 1;
      if(pin_size <= loop_size){
	    self->is_hairpin = TRUE;
	    self->pinsize = pin_size;
	    int index = self->i[self->size - 1] + 1; // Starting of loop.
	    if(index > self->j[self->size-1]){
		  index = self->j[self->size - 1] + 1; // Starting of Loop.
	    }
	    for(int i=0; i<self->pinsize; ++i){
		  self->loop_res[i] = nbp[index+i].resclass;
	    }
	    get_loop_name(self->loop_res, self->pinsize, self->loopname);

	    
	    return TRUE;
      }
      return FALSE;
}

void helix_calc_all(struct helix* self, int* hlxcount, struct nucbp* nbp, int nres){
      struct graph g;
      *hlxcount = 0;
      graph_init(&g, nres, UNDIRECTED); 

      
      for(int i=0; i<nres; ++i){
	    
	    
	    for(int k=0; k<nbp[i].numbp; ++k){
		  int j = nbp[i].oth_base_index[k];
		  if(graph_edge_index(&g, i, j) != -1) continue;
		  self[*hlxcount].size = 0;
		  helix_calc(self+ *hlxcount, nbp, &g, i, j);
		  helix_check_hairpin(self + *hlxcount, nbp);

		  
		  if(self[*hlxcount].size >= 3){
			*hlxcount = *hlxcount + 1;
		  }
	    }
      }

      
      graph_free(&g);
}
void helix_print_all(struct helix* self, struct nucbp* nbp, int count, FILE* fp){
      for(int i=0; i<count; ++i){
	    helix_fprint(self+i, nbp, fp);
      }
}

//void helix_pymol_gen_all(struct helix* self,  struct pymol_param* pyparam, struct nucbp* nbp, int count, FILE* fp){
//      int hlxcount = 0;
//      int loopcnt = 0;
//      for(int i=0; i<count; ++i){
//	    void helix_pymol_gen(self+i, pyparam, nbp, &hlxcount, &loopcount, fp){
//	    helix_fprint(self+i, nbp, fp);
//      }
//}
//
void helix_free_all(struct helix* self){
      free ( self );
}


void pymol_param_init(struct pymol_param* self){
      strcpy(self->hlxcol, "green");
      strcpy(self->loopcol, "white");
      strcpy(self->psuhlxcol, "red");
      strcpy(self->nonhlxcol, "blue");
      strcpy(self->psujunccol, "magenta");
      self->psudocount = 0;
      self->loopcount = 0;
      self->helixcount = 0;
}

void helix_gen_psudo_pml(struct pseudo_helix* self, struct nucbp* nbp, struct stack* stk, struct pymol_param* pyparam, FILE* fp){
      stack_push(stk, self);
      if(self->count == 0){
//	    fprintf(fp, "select psuhelix%d, ", pyparam->psudocount);
	    struct pseudo_helix* tmphlx = (pseudo_helix*) stk->elem[0];
	    char psubreak[2048];
	    char psuline[512];
	    strcpy(psubreak, "");
	    int junccount = 0;
	    for(int i=0; i<=stk->top; ++i){
		  tmphlx = (pseudo_helix*) stk->elem[i];
		  for(int k=0; k<tmphlx->hlx.size; ++k){
//			fprintf(fp, "(resi %d and chain %s) ", nbp[tmphlx->hlx.i[k]].cifid,nbp[tmphlx->hlx.i[k]].chain);
//			fprintf(fp, "(resi %d and chain %s) ", nbp[tmphlx->hlx.j[k]].cifid,nbp[tmphlx->hlx.j[k]].chain);
			if(i<stk->top && k == tmphlx->hlx.size -1){
			      sprintf(psuline, "select psubreak%d, ", junccount);
			      strcat(psubreak, psuline);
			      sprintf(psuline, "(resi %d and chain %s) ", nbp[tmphlx->hlx.i[k]].cifid,nbp[tmphlx->hlx.i[k]].chain);
			      strcat(psubreak, psuline);
			      sprintf(psuline, "(resi %d and chain %s) ", nbp[tmphlx->hlx.j[k]].cifid,nbp[tmphlx->hlx.j[k]].chain);
			      strcat(psubreak, psuline);
			      sprintf(psuline, "\ncolor %s, psubreak%d\n", pyparam->psujunccol, junccount);
			      strcat(psubreak, psuline);
			      sprintf(psuline, "show_as surface, psubreak%d\n\n", junccount);
			      strcat(psubreak, psuline);
			      junccount ++;
			}
		  }
	    }
	    fprintf(fp,"\n");

//	    fprintf(fp, "color %s, psuhelix%d\n", pyparam->psuhlxcol, pyparam->psudocount);
//	    fprintf(fp, "show spheres, psuhelix%d\n", pyparam->psudocount );

	    fprintf(fp,"\n");
	    fprintf(fp, "%s\n", psubreak);

	    fprintf(fp,"\n");
	    pyparam->psudocount = pyparam->psudocount + 1;
      }
      for(int k=0; k<self->count; ++k){
	    helix_gen_psudo_pml(self->branch[k], nbp, stk, pyparam, fp);
      }
      stack_pop(stk);
}

void helix_pymol_gen(struct helix* self, struct pymol_param* pyparam, struct nucbp* nbp, FILE* fp){
//      if(self->is_hairpin == TRUE){
//	    fprintf(fp, "select stemloop%d, ", pyparam->loopcount);
//      }else{
//	    fprintf(fp, "select helix%d, ", pyparam->helixcount);
//      }
      fprintf(fp, "select helix%d, ", pyparam->helixcount);
      for(int k=0; k<self->size; ++k){
	    fprintf(fp, "(resi %d and chain %s) ", nbp[self->i[k]].cifid, nbp[self->i[k]].chain) ;
	    fprintf(fp, "(resi %d and chain %s) ", nbp[self->j[k]].cifid, nbp[self->j[k]].chain) ;
      }
      fprintf(fp, "\ncolor %s, helix%d\n", pyparam->hlxcol, pyparam->helixcount);
      fprintf(fp, "set cartoon_ring_mode, 3, helix%d\n", pyparam->helixcount);
      fprintf(fp, "set cartoon_ladder_mode, 1, helix%d\n", pyparam->helixcount);
      fprintf(fp, "show_as cartoon, helix%d\n\n", pyparam->helixcount);
      pyparam->helixcount = pyparam->helixcount + 1;
      if(self->is_hairpin == TRUE){
	    fprintf(fp, "select loop%d, ", pyparam->loopcount);
	    int small = (self->i[self->size-1]< self->j[self->size-1]) ? self->i[self->size-1] : self->j[self->size-1];
	    for(int i=1; i<=self->pinsize; ++i){
		  fprintf(fp, "(resi %d and chain %s) ", nbp[small+i].cifid, nbp[small+i].chain);
	    }
	    fprintf(fp, "\ncolor %s, loop%d\n", pyparam->loopcol, pyparam->loopcount);
	    fprintf(fp, "show_as cartoon, loop%d\n\n", pyparam->loopcount);
	    fprintf(fp, "show stick, loop%d\n\n", pyparam->loopcount);
	    pyparam->loopcount = pyparam->loopcount + 1; 
      }
//      if(self->is_hairpin == TRUE){
//	    fprintf(fp, "\ncolor %s, stemloop%d\n", pyparam->loopcol, pyparam->loopcount);
//	    fprintf(fp, "show sphere, stemloop%d\n\n", pyparam->loopcount);
//	    pyparam->loopcount = pyparam->loopcount + 1; 
//      }else{
//	    fprintf(fp, "\ncolor %s, helix%d\n", pyparam->loopcol, pyparam->helixcount);
//	    fprintf(fp, "show sphere, helix%d\n\n", pyparam->helixcount);
//	    pyparam->helixcount = pyparam->helixcount + 1;
//
//      }
      fprintf(fp, "\n\n");
}
void helix_gen_pymol(struct pseudo_helix* psuhlx[], int psuhlxcount, struct helix* hlx, int helixcount, struct nucbp* nbp, sysparams* syspar){
      struct pymol_param pymolparam;
      pymol_param_init(&pymolparam);
      string inputfile_name  = syspar->accn+"."+syspar->ext;
      string outputfile_name = syspar->file_dir+syspar->accn+"_helix.pml";
      FILE* fp = fopen(outputfile_name.c_str(), "w");
      fprintf(fp, "load %s\n", inputfile_name.c_str());

      fprintf(fp, "\nselect proteinall, polymer.protein");
      fprintf(fp, "\ncolor olive,  proteinall");
      fprintf(fp, "\nset cartoon_transparency, 0.5, proteinall");
      fprintf(fp, "\nshow cartoon, proteinall\n");

      fprintf(fp, "\nselect nucleicall, polymer.nucleic");
      fprintf(fp, "\nset_bond stick_radius, 0.50, nucleicall");
      fprintf(fp, "\ncolor yellow,  nucleicall");
      fprintf(fp, "\nset sphere_scale, 0.50, nucleicall");
//      fprintf(fp, "\nshow sticks, nucleicall");
//      fprintf(fp, "\nshow spheres, nucleicall");
//      fprintf(fp, "\nset cartoon_ladder_mode,0");
//      fprintf(fp, "\nset cartoon_ring_finder,0");
      fprintf(fp, "\nshow cartoon, nucleicall\n");

      for(int i=0; i<helixcount; ++i){
	    helix_pymol_gen(hlx+i, &pymolparam, nbp, fp);

      }

      struct stack stk;
      stack_init(&stk, psuhlxcount+10);
      for(int i=0; i<psuhlxcount; ++i){
	    helix_gen_psudo_pml(psuhlx[i], nbp, &stk, &pymolparam,  fp);
      }
      stack_free(&stk);
      fclose(fp);
}


#endif   /* ----- #ifndef __helix_H__  ----- */

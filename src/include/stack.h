/*
 * =====================================================================================
 *
 *       Filename:  stack.h
 *
 *    Description:  This library is for generic stack that stores
 *                  void pointers. The users have to cast it accordingly.  
 *
 *        Version:  1.0
 *        Created:  02/07/21 08:05:47 PM IST
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  PARTHAJIT ROY (PR), roy.parthajit@gmail.com
 *   Organization:  The University of Burdwan
 *
 * =====================================================================================
 */



struct stack{
      void** elem;
      int top;
      int smax;
};
void stack_init(struct stack* self, int size){
      self->smax = size;
      self->elem = (void**) malloc (size * sizeof(void*));
      self->top  = -1; 
}

void stack_push(struct stack* self, void* val){
      if(self->top == self->smax-1){    /* Exception Handling */ 
	    fprintf(stderr, "Error in function %s. (File: %s, Line %d)... Stack overflows.\n", __func__, __FILE__, __LINE__);
	    exit(EXIT_FAILURE);
      }
      self->elem[++self->top] = val;
}
void* stack_pop(struct stack* self){
      return self->elem[self->top--];
}
void stack_free(struct stack* self){
      free(self->elem);
}

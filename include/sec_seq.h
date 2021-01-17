//
// Created by parthajit on 29/2/20.
//

#ifndef BPNET_BASIC_03_SEC_SEQ_H
#define BPNET_BASIC_03_SEC_SEQ_H

#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>

typedef struct{
    char* sec_seq_str;
    int size;
}sequence_t;

void sequence_create(sequence_t* seq, const char* dat_file_name){
    FILE* fp = fopen(dat_file_name, "r");
    assert(fp != NULL);
    char tmpseq[20000];
    seq->size = 0;
    char line[100];
    int len;
    cout<<"One\n";
    while(fgets(line, sizeof(line), fp) != NULL) {
    cout<<"coming inside="<<line<<"END"<<endl;
        if(line[0] == '>') continue;
        len = sprintf(tmpseq+seq->size,"%s", line);
        cout<<"Line="<<line<<", len="<<len<<"."<<endl;
        seq->size += len - 1; // To discard \n
    }
    cout<<"Two\n";
    /* Allocate the memory of the struture and copy */

    seq->sec_seq_str = (char*) malloc((seq->size+2)* sizeof(char));
    //cout<<"SEq size = "<<seq->size<<endl;
    sprintf(seq->sec_seq_str,"%s",tmpseq);
    seq->sec_seq_str[seq->size] = '\0';
    //printf("The val=%s.\n",seq->sec_seq_str);
    //printf("SEQ size = %d\n", seq->size);
    assert(fp != NULL);
    fclose(fp);
}
void sequence_destroy(sequence_t* seq){
    free(seq->sec_seq_str);
}

#endif //BPNET_BASIC_03_SEC_SEQ_H

//
// Created by parthajit on 26/9/20.
//

#ifndef BPNET_BASIC_03_NTVARIANTS_H
#define BPNET_BASIC_03_NTVARIANTS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#define TRUE 1
#define FALSE 0
typedef struct{
	char AdeVariants[200][5];
	char GuaVariants[200][5];
	char CytVariants[200][5];
	char UraVariants[200][5];
	int AdeVarCnt;
	int GuaVarCnt;
	int CytVarCnt;
	int UraVarCnt;
}ntvariants_t;

int ntvar_is_adevar(ntvariants_t* ntvar, const char* base){
    for(int i=0; i<ntvar->AdeVarCnt; ++i){
        if(strcmp(ntvar->AdeVariants[i], base) == 0){
            return TRUE;
        }
    }
    return FALSE;
}

int ntvar_is_guavar(ntvariants_t* ntvar, const char* base){
    for(int i=0; i<ntvar->GuaVarCnt; ++i){
        if(strcmp(ntvar->GuaVariants[i], base) == 0){
            return TRUE;
        }
    }
    return FALSE;
}
int ntvar_is_cytvar(ntvariants_t* ntvar, const char* base){
    for(int i=0; i<ntvar->CytVarCnt; ++i){
        if(strcmp(ntvar->CytVariants[i], base) == 0){
            return TRUE;
        }
    }
    return FALSE;
}
int ntvar_is_uravar(ntvariants_t* ntvar, const char* base){
    for(int i=0; i<ntvar->UraVarCnt; ++i){
        if(strcmp(ntvar->UraVariants[i], base) == 0){
            return TRUE;
        }
    }
    return FALSE;
}




void ntvar_populate(ntvariants_t* ntvar){
	char* nucdir = getenv("NUCLEIC_ACID_DIR");
	if(nucdir == NULL){
		fprintf(stderr, "Error... NUCLEIC_ACID_DIR not Defined.\n");
		exit(EXIT_FAILURE);
	}

	char filename[512];
	strcpy(filename,nucdir);
	strcat(filename,"/AdeVariants.name");
	FILE * fp = fopen(filename, "r");
	assert(fp != NULL);
	char line[20];
	ntvar->AdeVarCnt = 0;
	while(fgets(line, 20, fp) != NULL){
		int len = (int)strlen(line);
		line[len-1] = '\0'; // To stop newline char
		strcpy(ntvar->AdeVariants[ntvar->AdeVarCnt], line);
		ntvar->AdeVarCnt++;
	}
	fclose(fp);

	strcpy(filename,nucdir);
	strcat(filename,"/GuaVariants.name");
	fp = fopen(filename, "r");
	assert(fp != NULL);
	ntvar->GuaVarCnt = 0;
	while(fgets(line, 20, fp) != NULL){
		int len = (int)strlen(line);
		line[len-1] = '\0'; // To stop newline char
		strcpy(ntvar->GuaVariants[ntvar->GuaVarCnt], line);
		ntvar->GuaVarCnt++;
	}
	fclose(fp);

	strcpy(filename,nucdir);
	strcat(filename,"/CytVariants.name");
	fp = fopen(filename, "r");
	assert(fp != NULL);
	ntvar->CytVarCnt = 0;
	while(fgets(line, 20, fp) != NULL){
		int len = (int)strlen(line);
		line[len-1] = '\0'; // To stop newline char
		strcpy(ntvar->CytVariants[ntvar->CytVarCnt], line);
		ntvar->CytVarCnt++;
	}
	fclose(fp);

	strcpy(filename,nucdir);
	strcat(filename,"/UraVariants.name");
	fp = fopen(filename, "r");
	assert(fp != NULL);
	ntvar->UraVarCnt = 0;
	while(fgets(line, 20, fp) != NULL){
		int len = (int)strlen(line);
		line[len-1] = '\0'; // To stop newline char
		strcpy(ntvar->UraVariants[ntvar->UraVarCnt], line);
		ntvar->UraVarCnt++;
	}
	strcpy(ntvar->UraVariants[ntvar->UraVarCnt], "PSU");
	ntvar->UraVarCnt++;
	fclose(fp);
}

int normal_adenine(const char* base){
    if(strcmp(base, "A") == 0 || strcmp(base, "DA") == 0){
        return TRUE;
    }else{
        return FALSE;
    }
}

int normal_guanine(const char* base){
    if(strcmp(base, "G") == 0 || strcmp(base, "DG") == 0){
        return TRUE;
    }else{
        return FALSE;
    }
}

int normal_cytosine(const char* base){
    if(strcmp(base, "C") == 0 || strcmp(base, "DC") == 0){
        return TRUE;
    }else{
        return FALSE;
    }
}

int normal_uracil(const char* base){
    if(strcmp(base, "U") == 0 || strcmp(base, "DT") == 0){
        return TRUE;
    }else{
        return FALSE;
    }
}

char ntvar_getchar(ntvariants_t * ntvar, char* res_name)
{
      char base_name;

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
            exit(EXIT_FAILURE);
        }
	return base_name;

}


#endif //BPNET_BASIC_03_NTVARIANTS_H

#include "utils.h"
#include <stdio.h>
#include <stdlib.h>

double **read_matrix(char *source, int n){

    FILE *file = fopen(source, "r");

    if(file == NULL){
        return NULL;
    }

    double **arr = malloc(n * sizeof(double *));
    if(arr == NULL){
        fclose(file);
        return NULL;
    }

    for(int i=0; i<n; i++){
        //Creo memoria para cada renglon y la almaceno en arr
        arr[i] = malloc(n * sizeof(double));

        for(int j=0; j<n; j++){

            if(fscanf(file, "%lf", &arr[i][j]) != 1){ //Seria lo mismo si pongo *(arr+i)+j
                printf("Error al leer la matriz, revisa el .txt\n");
                break;
            }
        }
    }

    fclose(file);

    return arr;


}

double **create_matrix(int n){

    double **L = (double **)malloc(n * sizeof(double *));
    
    if(L == NULL){
        printf("Error, no se creó la matriz\n");
        return NULL;
    }
    
    for(int i=0; i<n; i++){

        L[i] = calloc(n, sizeof(double));

        if(L[i] == NULL){
        printf("Error, no se creó la matriz\n");
        for(int m=0; m<i; m++) free(L[m]);
        free(L);
        return NULL;
        }

    }

    return L;

}

void free_matrix(double **M, int n) {
    for (int i = 0; i < n; i++) free(M[i]);
    free(M);
}

void display_matrix(double **M, int n){

    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){

            printf("%.2lf ", M[i][j]);

        }
        printf("\n");
    }
}





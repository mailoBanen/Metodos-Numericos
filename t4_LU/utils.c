#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MAT_TOLERANCE 1e-9

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
    printf("\n");
}

void multiply(double **A, double **B, double **output, int n){

    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            
            double valor = 0;
            for(int k=0; k<n; k++){

                valor += A[i][k] * B[k][j];

            }

            output[i][j] = valor;

        }
    }
}

int is_equal(double **A, double **B, int n){


    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){

            if(fabs(A[i][j] - B[i][j]) > MAT_TOLERANCE){

                printf("Las matrices no son iguales.\n");
                printf("A[%d][%d] = %f  B[%d][%d] = %f\n", i, j, A[i][j], i, j, B[i][j]);
                return 0;
            }
        }
    }

    return 1;
}


#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MAT_TOLERANCE 1e-5

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

double *read_vec(char *source, int n){

    FILE *file = fopen(source, "r");

    if(file == NULL){
        return NULL;
    }

    //Creo un puntero de tamaño n, y sino retorna nulo
    double *arr = (double *)malloc(n * sizeof(double));
    if(arr == NULL){
        fclose(file);
        return NULL;
    }

    for(int i=0; i<n; i++){

            //Escaneo el archivo y voy posicionando los valores en el arreglo
            if(fscanf(file, "%lf", &arr[i]) != 1){ //Seria lo mismo si pongo *(arr+i)+j
                printf("Error al leer la matriz, revisa el .txt\n");
                break;
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

double *create_vec(int n){

    double *V = (double *)calloc(n, sizeof(double));
    if(V == NULL){
        printf("Error creando vector de tamaño %d\n", n);
        return NULL;
    } 
    return V;
}

void fill_vec(double *vec, double num, int n){
    for(int i=0; i<n; i++){
        vec[i] = num;
    }
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

void display_vec(double *b, int n){

    for(int i=0; i<n; i++){
            printf("%.3lf\n", b[i]);
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

void matXvec(double **A, double *x, double *vout, int n){

    for(int i=0; i<n; i++){
        double suma=0;
        for(int j=0; j<n; j++){

            suma += A[i][j] * x[j];

        }

        vout[i] = suma;
    }
}




int is_equal_m(double **A, double **B, int n){

    double error_max=0;

    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){

            if(fabs(A[i][j] - B[i][j]) > error_max){
                error_max = fabs(A[i][j] - B[i][j]);
            }

            if(fabs(A[i][j] - B[i][j]) > MAT_TOLERANCE){

                printf("Las matrices no son iguales.\n");
                printf("A[%d][%d] = %f  B[%d][%d] = %f\n", i, j, A[i][j], i, j, B[i][j]);
                return 0;
            }
        }
    }

    printf("Error: %e\n", error_max);

    return 1;
}

int is_equal_v(double *x, double *y, int n){

    double error_max = 0;

    for(int i=0; i<n; i++){

        if(fabs(x[i] - y[i]) > error_max){
            error_max = fabs(x[i] - y[i]);
        }

        if(fabs(x[i] - y[i]) > MAT_TOLERANCE){

            printf("Vectores no iguales: x[%d]=%f , y[%d]=%f\n", i, x[i], i, y[i]);
            return 0;
        }
    }

    printf("Son iguales\ndiff: %e\n", error_max);
    return 1;
}


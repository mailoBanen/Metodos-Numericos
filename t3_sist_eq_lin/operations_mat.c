#include "operations_mat.h"
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

void display_matrix(double **M, int n){

    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){

            printf("%.2lf ", M[i][j]);

        }
        printf("\n");
    }
}

void display_vec(double *b, int n){

    for(int i=0; i<n; i++){
            printf("%.3lf\n", b[i]);
    }
}

void free_matrix(double **M, int n){

    for(int i=0; i<n; i++){

        free(M[i]);
    }

    free(M);
}

void free_vec(double *b){
    free(b);
}


//Operaciones
double *diagonal(double **D, double *b, int n){

    double *x = (double *)malloc(n * sizeof(double));

    for(int i=0; i<n; i++){
        x[i] = b[i] / D[i][i];
    }

    return x;
}

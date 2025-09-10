#include <stdio.h>
#include <stdlib.h>

double **read_matrix(char *source, int n);
void display_matrix(double **M, int n);
void free_matrix(double **M, int n);

int main(){

    double **matrix;
    int n=15; 

    matrix = read_matrix("sistemasDeEqLineales/D.txt", n);
    display_matrix(matrix, n);
    free_matrix(matrix, n);

    return 0;
}


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

void display_matrix(double **M, int n){

    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){

            printf("%.2lf ", M[i][j]);

        }
        printf("\n");
    }
}

void free_matrix(double **M, int n){

    for(int i=0; i<n; i++){

        free(M[i]);
    }

    free(M);
}






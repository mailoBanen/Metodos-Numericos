#include "operations_mat.h"
#include <stdio.h>
#include <stdlib.h>

int main(){

    double **matrix, *vec, *x;

    //Tamaño de mi sistema
    int n=15; 

    //Leo de mis archivos la matriz y el vector
    matrix = read_matrix("sistemasDeEqLineales/A.txt", n);
    vec = read_vec("sistemasDeEqLineales/b.txt", n);
    
    //Resuelvo la matriz diagonal
    x = gauss_pivot(matrix, vec, n);

    //display_matrix(matrix, n);
    free_matrix(matrix, n);

    display_vec(x, n);

    free_vec(x);


    return 0;
}
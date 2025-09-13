#include <stdio.h>
#include <stdlib.h>
#include "operations_mat.h"

int main(){

    double **matrix, *vec, *x;
    //Tamaño de mi sistema
    int n=15; 

    //Leo de mis archivos la matriz y el vector
    matrix = read_matrix("sistemasDeEqLineales/L.txt", n);
    vec = read_vec("sistemasDeEqLineales/b_l.txt", n);
    
    //Resuelvo la matriz diagonal
    x = lower(matrix, vec, n);

    //display_matrix(matrix, n);
    free_matrix(matrix, n);

    display_vec(x, n);

    free_vec(x);
    free_vec(vec);

    return 0;
}
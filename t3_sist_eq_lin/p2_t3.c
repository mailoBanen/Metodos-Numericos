#include <stdio.h>
#include <stdlib.h>
#include "operations_mat.h"


int main(){

    double **matrix, *vec, *x;
    //Tamaño de mi sistema
    int n=15; 

    //Leo de mis archivos la matriz y el vector
    matrix = read_matrix("sistemasDeEqLineales/U.txt", n);
    vec = read_vec("sistemasDeEqLineales/b_u.txt", n);
    
    //Resuelvo la matriz diagonal
    x = upper(matrix, vec, n);

    //display_matrix(matrix, n);
    free_matrix(matrix, n);

    display_vec(x, n);

    free_vec(x);

    return 0;
}
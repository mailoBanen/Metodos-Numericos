#include <stdio.h>
#include <stdlib.h>
#include "operations_mat.h"

int main(){


    double **matrix, *vec, *x;
    int n=15; 

    matrix = read_matrix("sistemasDeEqLineales/D.txt", n);
    vec = read_vec("sistemasDeEqLineales/b_d.txt", n);
    x = diagonal(matrix, vec, n);

    //display_matrix(matrix, n);
    free_matrix(matrix, n);

    display_vec(x, n);


    free_vec(x);
    free_vec(vec);

    return 0;
}
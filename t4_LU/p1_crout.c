#include "utils.h"
#include "decomposition.h"
#include <stdio.h>
#include <stdlib.h>

int main(){

    int n=50;
    double **A =  read_matrix("sistemasTarea4/A.txt", n);
    double *b = read_vec("sistemasTarea4/b.txt", n);
    double **L, **U;

    //Creo las matrices donde descompondré
    L = create_matrix(n);
    U = create_matrix(n);

    //Factorizo LU
    crouth(A, L, U, n);

    /*Creo el vector solucion*/
    double *x = create_vec(n);

    /*Resuelvo LU con factorizacion Crouth*/
    solve_crouth(L, U, x, b, n);
    
    /*Creo el vector resultado para comparar Ax = r = b*/
    double *r = create_vec(n);
    matXvec(A, x, r, n);

    /*Comparo si r y b son iguales*/
    is_equal_v(r, b, n);

    free_matrix(A, n);
    free_matrix(L, n);
    free_matrix(U, n);
    free(x);
    free(r);

    return 0;
}
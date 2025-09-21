#include "utils.h"
#include "decomposition.h"
#include <stdio.h>
#include <stdlib.h>


int main(){

    int n=50;
    double **A =  read_matrix("sistemasTarea4/SPD.txt", n);
    double *b = read_vec("sistemasTarea4/b_spd.txt", n);
     
    double **L, **LT;

    /*Creo mis matrices L y LT*/
    L = create_matrix(n);
    LT = create_matrix(n);

    /*Factorizo con cholesky*/
    cholesky(A, L, LT, n);

    double *x = create_vec(n);

    /*Resuelvo el sistema con cholesky*/
    solve_LU(L, LT, x, b, n);

    double *r = create_vec(n);
    matXvec(A, x, r, n);

    /*Confirmo que la solucion sea igual al vector b*/
    is_equal_v(b, r, n);

    free_matrix(A, n);
    free_matrix(L, n);
    free_matrix(LT, n);
    free(r);
    free(x);


    return 0;
    
}
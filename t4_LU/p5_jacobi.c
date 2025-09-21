#include "utils.h"
#include "decomposition.h"
#include <stdio.h>
#include <stdlib.h>


int main(){

    int n=125;
    double **A =  read_matrix("sistemasTarea4/A_big.txt", n);
    double *b = read_vec("sistemasTarea4/b_big.txt", n);

    /*Vector solucion creado*/
    double *x = calloc(n, sizeof(double));

    /*Aplico jacobi y guardo en x*/
    jacobi(A, x, b, n, 1000, 1e-5);

    double *comp = create_vec(n);

    matXvec(A, x, comp, n);
    is_equal_v(comp, b, n);


    free_matrix(A, n);
    free(x);
    free(comp);
    free(b);


    return 0;
    
}
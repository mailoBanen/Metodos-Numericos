#include "utils.h"
#include "decomposition.h"
#include <stdio.h>
#include <stdlib.h>


int main(){

    int n=50;
    double **A =  read_matrix("sistemasTarea4/SPD.txt", n);
     
    
    double **L, **U;

    L = create_matrix(n);
    U = create_matrix(n);

    cholesky(A, L, U, n);
    double **R = create_matrix(n);
    multiply(L, U, R, n);

    if(is_equal(A, R, n)){
        printf("Son iguales\n");
    }

    free_matrix(A, n);
    free_matrix(L, n);
    free_matrix(U, n);
    free_matrix(R, n);

    return 0;
    
}
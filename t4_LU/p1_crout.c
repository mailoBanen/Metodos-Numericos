#include "utils.h"
#include "decomposition.h"
#include <stdio.h>
#include <stdlib.h>


int main(){

    int n=50;
    double **A =  read_matrix("sistemasTarea4/sistemasTarea4/A.txt", n);
    double **L, **U;

    L = create_matrix(n);
    U = create_matrix(n);

    crouth(A, L, U, n);
    //display_matrix(L, n);
    display_matrix(U, n);

    free_matrix(A, n);
    free_matrix(L, n);
    free_matrix(U, n);

    return 0;
    
}
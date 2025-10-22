#include <stdio.h>
#include <stdlib.h>
#include "utils.h"
#include "methods.h"


int main(){

    int rows, cols;
    double **A = read_matrix("archivos_tarea7/A1.txt", &rows, &cols);
    double *b = read_vec("archivos_tarea7/A1.txt", &rows);
    double *x = create_vec(rows);

    display_matrix(A, rows);
    solveQR(A, x, b, rows);
    display_vec(x, rows);
    err_solution(A, x, b, rows);


    


    return 0;
}
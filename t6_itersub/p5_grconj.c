#include <stdlib.h>
#include <stdio.h>
#include "utils.h"
#include "funciones.h"


int main(){

    int rows, cols;
    int n;

    double **A  = read_matrix("archivos_tarea6/A2.txt", &rows, &cols);
    double *b = read_vec("archivos_tarea6/b2.txt", &n);

    double *x_0 = create_vec(rows);

    gradient_conj(A, x_0, b, rows, 1000, 1e-10);


    //display_vec(x_0, rows);
    err_solution(A, x_0, b, rows);


    return 0;
}
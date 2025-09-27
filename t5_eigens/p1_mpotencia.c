#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "eigens.h"
#include "utils.h"



int main(){

    int rows, cols;
    double **A = read_matrix("archivos_tarea5/Eigen_5x5.txt", &rows, &cols);
    display_matrix(A, rows);

    /*Matriz para guardar todos los valores calculados*/
    double **eigens = create_matrix(rows+1);

    /*Cuantos quiero calcular*/
    int k=4;

    /*Vector inicial*/
    double *v_0 = create_vec(rows);
    fill_vec(v_0, 1/sqrt(rows), rows);

    k_power_method(k, A, v_0, eigens, rows, 1000, 1e-12);

    display_matrix(eigens, rows+1);

    
    free_matrix(A, rows);
    free_matrix(eigens, rows+1);
    free(v_0);

    return 0;
}
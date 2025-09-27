#include <stdio.h>
#include <stdlib.h>
#include "eigens.h"
#include "utils.h"

int main(){

    int rows, cols;
    double **A = read_matrix("archivos_tarea5/Eigen_3x3.txt", &rows, &cols);
    //display_matrix(A, rows);

    /*Creo la matriz inversa */
    double **A_inv = create_matrix(rows);
    inverse_cholesky(A, A_inv, rows);
    display_matrix(A, rows);
  

    /*Valores iniciales*/
    double *x = create_vec(rows);
    x[0] = 1;

    double lambda = 1;


    k_biggest_evals(2, A, x, &lambda, rows, 1000, 1e-6);

    /*Aplico metodo de la potencia inversa*/
    // inverse_power_method(A_inv, x, &lambda, rows, 100, 1e-12);

    // printf("El eigen minimo es: %lf \n", lambda);
    // display_vec(x, rows);

    // /*Confirmo que sea un valor y vector valido*/
    // double *r1=create_vec(rows);
    // double *r2=create_vec(rows);
    // matXvec(A, x, r1, rows);
    // numXvec(lambda, x, r2, rows);

    // is_equal_v(r1, r2, rows);

    free_matrix(A, rows);
    free_matrix(A_inv, rows);
    free(x);

    return 0;
}
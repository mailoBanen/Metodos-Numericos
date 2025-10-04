#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "utils.h"
#include "funciones.h"



int main(){

    int rows, cols;
    double **A = read_matrix("archivos_tarea6/Eigen_125x125.txt", &rows, &cols);
    int k=7;
    
    double *v_0 = create_vec(rows);
    fill_vec(v_0, 1.0/sqrt(rows), rows);

    double **eigens = create_matrix(rows+1);
    k_power_method(k, A, v_0, eigens, rows, 10000, 1e-10);
    
    /*Error de los valores obtenidos*/
    for(int i=0; i<k; i++){

        error_evec(A, eigens[i], eigens[rows][i], rows);

    }

    //display_matrix(eigens, rows+1);
    
   

    return 0;
}
#include <stdlib.h>
#include <stdio.h>
#include "utils.h"
#include "funciones.h"



int main(){

    int rows, cols;
    double **A = read_matrix("archivos_tarea6/Eigen_125x125.txt", &rows, &cols);
    // display_matrix(A, rows);
    
    /*Las columnas de Phi son el numero de eigenvectores a calcular (colPhi)*/
    int rowPhi = rows, colPhi = 7;
    double **Phi = create_phi(rowPhi, colPhi);
    // printf("Phi incial:\n");
    // display_mat_rect(Phi, rowPhi, colPhi);


    /*Vector que va a guardar los eigenvalores*/
    double *lambdas = create_vec(rows);

    
    /*Phi va a tener como columnas los eigenvectores de A
    Si A es de nxn y quiero k vectores necesito que 
    Phi nxk para tener k vectores*/
    iter_subspace(A, rows, cols, Phi, rowPhi, colPhi, lambdas, 1000, 1e-5);

    display_vec(lambdas, colPhi);
    //display_mat_rect(Phi, rowPhi, colPhi);

    double *r1 = create_vec(rows);

    /*Error de los valores obtenidos*/
    for(int i=0; i<colPhi; i++){

        get_column(Phi, r1, rows, i);
        error_evec(A, r1, lambdas[i], rows);

    }
    

    return 0;
}
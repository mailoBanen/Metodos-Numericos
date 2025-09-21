#include "utils.h"
#include "decomposition.h"
#include <stdio.h>
#include <stdlib.h>

int main(){

    int n=4;
    double phi_0 = 0, phi_n = 100, K=1, Q=5, L=1;

    /*Creo los vectores de la matriz A*/
    double **vec_A = malloc(2 * sizeof(double *));
    double *main_diag_A = create_vec(n);
    double *sup_diag_A = create_vec(n);
    fill_vec(main_diag_A, 2, n);
    fill_vec(sup_diag_A, -1, n);

    vec_A[0] = main_diag_A;
    vec_A[1] = sup_diag_A;


    /*Muestro vectores de A*/
    /* for(int i=0; i<2; i++){
        for(int j=0; j<n; j++){
            printf("%f ", vec_A[i][j]);
        }
        printf("\n");
    } */

    /*Creo los vectores de la matriz L*/
    double **vec_L = malloc(2 * sizeof(double *));
    double *main_diag_L = create_vec(n);
    double *sup_diag_L = create_vec(n);
    vec_L[0] = main_diag_L;
    vec_L[1] = sup_diag_L;

    /*Muestro los vectores de L*/
    /* printf("\n");
    for(int i=0; i<2; i++){
        for(int j=0; j<n; j++){
            printf("%f ", vec_L[i][j]);
        }
        printf("\n");
    } */

    /*Ya tengo A, y tengo creado L, ahora hay que factorizar
    pero en el formato vector*/
    cholesky_tridiag_vec(vec_A, vec_L, n);

    /* printf("\n");
    for(int i=0; i<2; i++){
        for(int j=0; j<n; j++){
            printf("%f ", vec_L[i][j]);
        }
        printf("\n");
    } */

    /*Aqui ya tengo la factorizacion de A*/
    double delta_x  = L/(n+1), dx2 = (delta_x*delta_x);

    /*Creo todas las variables necesarias para solucionar el sistema*/
    double *b = create_vec(n);
    b[0]= Q/K*(dx2) + phi_0;
    b[n-1] = Q/K*(dx2) + phi_n;

    /*Inicializo el vector b*/
    for(int i=1; i<n-1; i++){
        b[i] = Q/K*dx2;
    }

    double *x = create_vec(n);

    /*Ya habiendo factorizado por cholesky, resuelvo el sistema
    con (LL^t)x=b*/
    solve_cholesky_tridiag_vec(vec_L, b, x, n);


    /*Comparo los resultados obtenidos*/
    double *r = create_vec(n);

    /*Multiplico la matriz Ax = r = b*/
    for(int i=0; i<n; i++){
        double suma=0;
        for(int j=0; j<n; j++){

            double val=0;
            if(i==j){
                val=2;
            }else if(j == i+1 || j==i-1){
                val = -1;
            }else{
                val=0;
            }

            suma += val*x[j];

        }
        r[i] = suma;
    }
    
    /*Verifico si la multiplicacion es igual al resultado.*/
    is_equal_v(r, b, n);

    return 0;
}
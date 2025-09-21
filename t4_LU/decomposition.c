#include "decomposition.h"
#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define ZERO_VAL 1e-12


void lower(double **L, double *x,  double *b, int n){

    x[0] = b[0]/L[0][0]; 

    double suma=0;

    for(int i=1; i<n; i++){
        suma = 0.0L;
        for(int j=0; j<i; j++){

            suma += L[i][j] * x[j];

        }

        //En caso de haber un numero practicamente zero
        //en la diagonal
        if(fabsl(L[i][i]) < ZERO_VAL){
            return;
        }

        x[i] = (b[i] - suma)/L[i][i];
    }
}

void upper(double **U, double *x, double *b, int n){

    x[n-1] = b[n-1]/U[n-1][n-1]; 

    double suma=0;

    for(int i=n-2; i>=0; i--){
        suma = 0.0L;
        for(int j=i+1; j<=n-1; j++){

            suma += U[i][j] * x[j];

        }

        //En caso de haber un numero practicamente zero
        //en la diagonal
        if(fabsl(U[i][i]) < ZERO_VAL){
            return;
        }

        x[i] = (b[i] - suma)/U[i][i];
    }
}



void crouth(double **A, double **L, double **U, int n){

    /*Inicializo la diagonal principal con 1*/
    for(int i=0; i<n; i++){
        U[i][i] = 1.0;
    }

    //Recorro la diagonal
    for(int i=0; i<n; i++){

        //Cada diagonal inicio desde el elemento de la diagonal, hacia abajo
        //con los valores de "l"
        for(int k=i; k<n; k++){/*k es el iterador para recorrer la columna*/

            double suma=0;

            for(int m=0; m<i; m++){

                suma += L[k][m] * U[m][i];

            }

            L[k][i] = A[k][i] - suma;

        }

        //Despues en esa mismo diagonal, calculo todos los valores de u
        for(int k=i+1; k<n; k++){/*k es el iterador para recorrer la fila
                                    exceptuando el elemento de la diagonal*/

            double suma=0;

            for(int m=0; m<i; m++){

                suma += L[i][m] * U[m][k];

            }

            U[i][k] = (A[i][k] - suma) / L[i][i];

        }
    }
}

void solve_crouth(double **L, double **U, double *x, double *b, int n){

    double *y = create_vec(n);
    lower(L, y, b, n);

    upper(U, x, y, n);
    free(y);
}

void solve_LU(double **L, double **U, double *x, double *b, int n){

    double *y = create_vec(n);
    lower(L, y, b, n);

    upper(U, x, y, n);
    free(y);
}



void cholesky(double **A, double **L, double **U, int n){

    L[0][0] = sqrt(A[0][0]);
    U[0][0] = sqrt(A[0][0]);

    for(int i=1; i<n; i++){

        L[i][0] = A[i][0]/L[0][0];
        U[0][i] = L[i][0];

    }

    if(isnan(L[0][0])){

        printf("No es posible usar cholesky\n");
        exit(2);
    }  


    for(int j=1; j<n-1; j++){

        double suma = 0;

        for (int k = 0; k < j; k++)
        {
            suma += L[j][k] * L[j][k];
        }

        L[j][j] = sqrt(A[j][j] - suma);
        U[j][j] = L[j][j];

        for(int i=j+1; i<n; i++){
            
            double suma = 0;

            for(int k=0; k<i; k++){

                suma += L[i][k] * L[j][k];

            }


            L[i][j] = (A[i][j] - suma) / L[j][j];

            if(isnan(L[i][i])){

                printf("No es posible usar cholesky\n");
                exit(2);
            }

            U[j][i] = L[i][j];

        }
    }

    double suma = 0;
    for(int k=0; k<n-1; k++){

        suma += (L[n-1][k] * L[n-1][k]);

    }

    L[n-1][n-1] = sqrt(A[n-1][n-1] - suma);
    U[n-1][n-1] = L[n-1][n-1];

}

void cholesky_tridiag(double **A, double **L, double **LT, int n){

    L[0][0] = sqrt(A[0][0]);
    LT[0][0] = L[0][0];

    int j = 0;

    for(int i=1; i<n; i++){

        L[i][j] = A[i][j]/L[j][j];
        LT[j][i] = L[i][j];

        L[i][i] = sqrt(A[i][i] - (L[i][i-1] * L[i][i-1]) );
        LT[i][i] = L[i][i];

        j++;
    }
}

void cholesky_tridiag_vec(double **A, double **L, int n){

    /*La diagonal principal es en 0, la otra en 1*/

    L[0][0] = sqrt(A[0][0]);

    for(int i=1; i<n; i++){

        L[1][i] = A[1][i]/L[0][i-1];

        L[0][i] = sqrt(A[0][i] - (L[1][i] * L[1][i]));
    }
}

void solve_cholesky_tridiag_vec(double **L, const double *b, double *x, int n) {
    double *y = (double*)malloc(sizeof(double) * n);
    if (y==NULL) return;

    // Hacia adelante Ly = b
    y[0] = b[0] / L[0][0];
    for (int i = 1; i < n; ++i) {
        y[i] = (b[i] - L[1][i] * y[i-1]) / L[0][i];
    }

    //Hacia atrás L^Tx = Ux = y
    x[n-1] = y[n-1] / L[0][n-1];
    for (int i = n-2; i >= 0; --i) {
        x[i] = (y[i] - L[1][i+1] * x[i+1]) / L[0][i];
    }

    free(y);
}




void jacobi(double **A, double *x, double *b, int n, int n_iters, double tolerance){

    if(n<=0 || n_iters<=0){
        exit(2);
    }

    //Creo el vector siguiente del metodo
    double *x_next = (double *)malloc(n * sizeof(double)); 
    if(x_next == NULL) exit(1);

    /*Creo copia de cada uno de los vectores para no perder referencia
    y solo hacer swapeo de vectores*/
    double *x_curr = x;
    double *cx_next = x_next;

    /*Bandera de convergencia y contador de iteraciones*/
    int k=0, converged = 0;

    while(k < n_iters && converged == 0){

        /*El error se inicializa en cero*/
        double error = 0;

        /*Calculo la primera iteracion del vector solución*/
        for(int i=0; i<n; i++){

            if(fabs(A[i][i]) < ZERO_VAL){

                printf("Error: Cero en la diagonal\n");
                exit(2);
            }
            
            double suma=0;

            for(int j=0; j<n; j++){

                if(j == i) continue;

                suma += A[i][j]*x_curr[j];

            }

            cx_next[i] = (b[i] - suma)/A[i][i];

            double curr_err_abs = fabs(cx_next[i] - x_curr[i]);
            
            error = (curr_err_abs > error) ? curr_err_abs : error;

        }


        //printf("absoluto error: %e\n", error);
        /*Si encuentro una solucion valida, la regreso*/
        if(error < tolerance){

            printf("Solucion aproximada k: %d\n", k);
            if(x != cx_next){
                memcpy(x, cx_next, n*sizeof(double));
            }
            free(x_next);
            converged=1;
            return;
        }

        k++;

        /*Si no encontré una solucion, mi vector siguiente, ahora es el actual
        y busco la siguiente iteracion*/
        double *temp = x_curr;
        x_curr = cx_next;
        cx_next = temp;

    }

    if(!converged){
        printf("Se terminaron las iteraciones y no convergio\n");
        if(x != x_curr){
            memcpy(x, x_curr, n*sizeof(double));
        }
        free(x_next);
    }
}

void gauss_seidel(double **A, double *x, double *b, int n, int n_iters, double tolerance){

    double *x_old = (double *)malloc(n * sizeof(double));
    if(x_old == NULL){
        printf("Error: memoria insuficiente\n");
        exit(1);
    }

    int k=0;

    while(k < n_iters){
        
        /*El error se inicializa en cero*/
        double error = 0;

        memcpy(x_old, x, n*sizeof(double));

        /*Calculo la primera iteracion del vector solución*/
        for(int i=0; i<n; i++){

            double xi_old = x[i];

            if(fabs(A[i][i]) < ZERO_VAL){

                printf("Error: Cero en la diagonal\n");
                exit(2);
            }
            
            double suma=0;

            /*Vectores recientemente actualizados*/
            for(int j=0; j<i; j++){

                suma += A[i][j]*x[j];

            }

            /*Vectores no calculados aun*/
            for(int j=i+1; j<n; j++){

                suma+= A[i][j]*x_old[j];
            }


            x[i] = (b[i] - suma)/A[i][i];

            double err_abs = fabs(x[i] - xi_old);

            error = (err_abs > error) ? err_abs : error;

        }


        /*Si encuentro una solucion valida, la regreso*/
        if(error < tolerance){

            printf("Solucion aproximada n: %d\n", k);
            break;
        }

        k++;

    }

    if(k == n_iters){
    printf("Se terminaron las iteraciones\n");
    }

    free(x_old);
    
}

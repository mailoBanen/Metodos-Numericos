#include "decomposition.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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
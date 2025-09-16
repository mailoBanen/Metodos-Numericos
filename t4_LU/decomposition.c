#include "decomposition.h"
#include <stdio.h>
#include <stdlib.h>

void crouth(double **A, double **L, double **U, int n){

    /*Inicializo la diagonal principal con 1*/
    for(int i=0; i<n; i++){
        U[i][i] = 1.0;
    }

    //Recorro la diagonal
    for(int i=0; i<n; i++){

        //Cada diagonal inicio desde el elemento de la diagonal, hacia abajo
        //con los valores de l
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
#include "utils.h"
#include "decomposition.h"
#include <stdio.h>
#include <stdlib.h>

// Función para crear una matriz tridiagonal n x n
double** crear_matriz_tridiagonal(int n, double diag_principal, double diag_superior, double diag_inferior) {
    // Reservar memoria para el array de punteros (filas)
    double** matriz = (double**)malloc(n * sizeof(double*));
    if (matriz == NULL) {
        printf("Error: No se pudo asignar memoria para las filas\n");
        return NULL;
    }
    
    // Reservar memoria para cada fila
    for (int i = 0; i < n; i++) {
        matriz[i] = (double*)malloc(n * sizeof(double));
        if (matriz[i] == NULL) {
            printf("Error: No se pudo asignar memoria para la fila %d\n", i);
            // Liberar memoria ya asignada antes de salir
            for (int j = 0; j < i; j++) {
                free(matriz[j]);
            }
            free(matriz);
            return NULL;
        }
        
        // Inicializar toda la fila a cero
        for (int j = 0; j < n; j++) {
            matriz[i][j] = 0.0;
        }
    }
    
    // Llenar las tres diagonales
    for (int i = 0; i < n; i++) {
        // Diagonal principal
        matriz[i][i] = diag_principal;
        
        // Diagonal superior (excepto para la última fila)
        if (i < n - 1) {
            matriz[i][i + 1] = diag_superior;
        }
        
        // Diagonal inferior (excepto para la primera fila)
        if (i > 0) {
            matriz[i][i - 1] = diag_inferior;
        }
    }
    
    return matriz;
}


int main(){

    int n=4;
    //double **A =  read_matrix("sistemasTarea4/SPD.txt", n);
    //double **A = create_matrix(n);
    
    double **A = crear_matriz_tridiagonal(n, 2.0, 1.0, 1.0);
    display_matrix(A, n);

    double **L, **LT;

    L = create_matrix(n);
    LT = create_matrix(n);

    cholesky_tridiag(A, L, LT, n);
    
    double **R = create_matrix(n);
    multiply(L, LT, R, n);

    display_matrix(L, n);
    display_matrix(LT, n);

    if(is_equal(A, R, n)){
        printf("Son iguales\n");
    }

    free_matrix(A, n);
    free_matrix(L, n);
    free_matrix(LT, n);
    free_matrix(R, n);

    return 0;
    
}
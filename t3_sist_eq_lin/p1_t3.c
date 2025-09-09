#include <stdio.h>
#include <stdlib.h>

//La funcion recibe una matriz 2D, un vector, y el tamaño de la matriz.
double *diagonal(double **D, double *b, int size);


int main(){

    //Creo el vector b
    double *b = (double *)malloc(3 * sizeof(double));

    //Creo la matriz D
    double **D = (double **)malloc(3 * sizeof(double *));
    for (int i = 0; i < 3; i++) {
        D[i] = (double *)malloc(3 * sizeof(double));
    }

    //Lleno la diagonal de la matriz con valores conocidos
    //El vector b con solo "2"'s
    for(int i=0; i<3; i++){

        D[i][i] = i+1;
        b[i] = 2;
    }

    double *x = diagonal(D, b, 3);

    //Imprimo los valores para verificar valores
    for (size_t i = 0; i < 3; i++)
    {
        printf("%f\n", x[i]);
    }

    //libero x
    free(x);

    //Libero los renglones de D
    for (int i = 0; i < 3; i++) {
        free(D[i]);
    }

    //Libero D
    free(D);


    return 0;
}


double *diagonal(double **D, double *b, int n){

    double *x = (double *)malloc(n * sizeof(double));

    for(int i=0; i<n; i++){
        x[i] = b[i] / D[i][i];
    }

    return x;
}
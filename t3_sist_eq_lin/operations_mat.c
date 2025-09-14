#include "operations_mat.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MATRIX_TOLERANCE 1e-14


double **read_matrix(char *source, int n){

    FILE *file = fopen(source, "r");

    if(file == NULL){
        return NULL;
    }

    double **arr = malloc(n * sizeof(double *));
    if(arr == NULL){
        fclose(file);
        return NULL;
    }

    for(int i=0; i<n; i++){
        //Creo memoria para cada renglon y la almaceno en arr
        arr[i] = malloc(n * sizeof(double));

        for(int j=0; j<n; j++){

            if(fscanf(file, "%lf", &arr[i][j]) != 1){ //Seria lo mismo si pongo *(arr+i)+j
                printf("Error al leer la matriz, revisa el .txt\n");
                break;
            }
        }
    }

    fclose(file);

    return arr;

}

double *read_vec(char *source, int n){

    FILE *file = fopen(source, "r");

    if(file == NULL){
        return NULL;
    }

    //Creo un puntero de tamaño n, y sino retorna nulo
    double *arr = (double *)malloc(n * sizeof(double));
    if(arr == NULL){
        fclose(file);
        return NULL;
    }

    for(int i=0; i<n; i++){

            //Escaneo el archivo y voy posicionando los valores en el arreglo
            if(fscanf(file, "%lf", &arr[i]) != 1){ //Seria lo mismo si pongo *(arr+i)+j
                printf("Error al leer la matriz, revisa el .txt\n");
                break;
            }
        
    }

    fclose(file);

    return arr;
}

void display_matrix(double **M, int n){

    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){

            printf("%.2lf ", M[i][j]);

        }
        printf("\n");
    }
}

void display_vec(double *b, int n){

    for(int i=0; i<n; i++){
            printf("%.3lf\n", b[i]);
    }
}

void free_matrix(double **M, int n){

    for(int i=0; i<n; i++){

        free(M[i]);
    }

    free(M);
}

void free_vec(double *b){
    free(b);
}

void sum_arrs(double *arr1, double *arr2, double *result, double  n){

    for(int i=0; i<n; i++){

        result[i] = arr1[i] + arr2[i];

    }

}

void mult_arr(double *arr1, double k, double *result, int n){

    for(int i=0; i<n; i++){

        result[i] = arr1[i] * k;

    }
}

void swap_row(double **A, double *b, int i, int in, int n){

    double *temp = malloc((n+1) * sizeof(double));

    for(int j=0; j<n; j++){

        temp[j] = A[in][j];
        A[in][j] = A[i][j];
        A[i][j] = temp[j];

    }

    temp[n] = b[in];
    b[in] = b[i];
    b[i] = temp[n];

    free(temp);

}

void swap_col(double **A, int *perm, int j, int jn, int n){

    double temp=0;

    for(int i = 0; i < n; i++) {
    temp = A[i][j];
    A[i][j] = A[i][jn];
    A[i][jn] = temp;

    }

    int temp_perm; 
    temp_perm = perm[j];
    perm[j] = perm[jn];
    perm[jn] = temp_perm;
}

void find_max(double **A, double *b, int i, int j, int *perm, int n){

    //Inicio valores para maximo valor, y su respectiva posicion
    double max_val = fabs(A[i][j]);
    int max_row = i;
    int max_col = j;

    //Bandera en caso de que hubo un maximo en un renglon debajo
    int row_flag = 0;

    // Busco primero en los renglones debajo
    for (int in = i; in < n; in++) {
        if (fabs(A[in][j]) > max_val) {
            max_val = fabs(A[in][j]);
            max_row = in;
            row_flag = 1;
        }
    }

    //Si encontre un maximo en otro renglon, la bandera se activo
    //e intercambio renglones.
    if (row_flag) {
        swap_row(A, b, i, max_row, n);
        return;
    }

    // Si no encontre un mejor elemento en los renglones, reviso toda la submatriz.
    max_val = fabs(A[i][j]);
    for (int in = i; in < n; in++) {
        for (int jn = j; jn < n; jn++) {
            if (fabs(A[in][jn]) > max_val) {
                max_val = fabs(A[in][jn]);
                max_row = in;
                max_col = jn;
            }
        }
    }

    //Si se encontro un maximo en otra columna, cambio
    if (max_col != j) {
        swap_col(A, perm, j, max_col, n);
    }

    //Tambien cambio el renglon.
    if (max_row != i) {
        swap_row(A, b, i, max_row, n);
    }
}




//Operaciones
double *diagonal(double **D, double *b, int n){

    double *x = (double *)malloc(n * sizeof(double));

    for(int i=0; i<n; i++){
        x[i] = b[i] / D[i][i];
    }

    return x;
}

double *lower(double **L, double *b, int n){

    double *x = (double *)malloc(n * sizeof(double));
    x[0] = b[0]/L[0][0]; 

    double suma=0;

    for(int i=0; i<n; i++){
        suma = 0.0L;
        for(int j=0; j<i; j++){

            suma += L[i][j] * x[j];

        }

        //En caso de haber un numero practicamente zero
        //en la diagonal
        if(fabsl(L[i][i]) < MATRIX_TOLERANCE){
            free(x);
            return NULL;
        }

        x[i] = (b[i] - suma)/L[i][i];
    }

    return x;
}

double *upper(double **U, double *b, int n){

    double *x = (double *)malloc(n * sizeof(double));
    x[n-1] = b[n-1]/U[n-1][n-1]; 

    double suma=0;

    for(int i=n-2; i>=0; i--){
        suma = 0.0L;
        for(int j=i+1; j<=n-1; j++){

            suma += U[i][j] * x[j];

        }

        //En caso de haber un numero practicamente zero
        //en la diagonal
        if(fabsl(U[i][i]) < MATRIX_TOLERANCE){
            free(x);
            return NULL;
        }

        x[i] = (b[i] - suma)/U[i][i];
    }

    return x;
}

double *gauss_elim(double **A, double *b, int n){

    //constante escalamiento
    double m;

    double *aux = malloc(n * sizeof(double));
    if(aux == NULL) return NULL;
    double *x;

    //Recorro toda la diagonal con i
    for(int i=0; i<n; i++){
        //Con j recorro los renglones de la matriz
        for(int j=i+1; j<n; j++){

            //Para cada renglon, hay un m diferente en funcion del pivote 
            //y del valor 
            m = (A[j][i]/A[i][i]);

            //Multiplico el arreglo por una constante y guardo el resultado en un arreglo auxiliar
            mult_arr(A[i], -m, aux, n);

            //Despues sumo estos dos arreglos y los guardo en el que voy a modificar
            sum_arrs(A[j], aux, A[j],n);
            b[j] = b[j] - m*b[i] ;

        }
    }

    free(aux);

    x = upper(A, b, n);

    return x;

}

double *gauss_pivot(double **A, double *b, int n){

    //constante escalamiento
    double m;
    //Matriz de permutaciones
    int *perm = malloc(n * sizeof(int));
    //Inicio vector de permutacion.
    for(int i=0; i<n; i++)
    {    perm[i] = i;}
    
    //Matriz auxiliar para hacer operaciones con renglones.
    double *aux = malloc(n * sizeof(double));
    if(aux == NULL) return NULL;

    //Vector x de solucion
    double *x;

    //Recorro toda la diagonal con i
    for(int i=0; i<n; i++){
        //Con j recorro los renglones de la matriz
        for(int j=i+1; j<n; j++){

            //Si hay un cero en la diagonal, realizar pivoteo
            if(fabs(A[i][i]) < MATRIX_TOLERANCE){
                
                //Se modifica la misma matriz y se guardan las permutaciones
                find_max(A, b, i, i, perm, n);
            }

            //Para cada renglon, hay un m diferente en funcion del pivote 
            //y del valor 
            m = (A[j][i]/A[i][i]);

            //Multiplico el arreglo por una constante y guardo el resultado en un arreglo auxiliar
            mult_arr(A[i], -m, aux, n);

            //Despues sumo estos dos arreglos y los guardo en el que voy a modificar
            sum_arrs(A[j], aux, A[j],n);
            b[j] = b[j] - m*b[i] ;

        }
    }

    free(aux);

    //Vector respuesta ya ordenado
    double* x_ordered = (double*)malloc(n * sizeof(double));

    x = upper(A, b, n);

    for (int i = 0; i < n; i++) {
        x_ordered[perm[i]] = x[i];
    }

    free(x);
    free(perm);

    return x_ordered;
}


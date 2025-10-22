#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "methods.h"
#include "utils.h"
#define ZERO_VAL 1e-16

// //Funcion original
// double f(double x){

//     //printf("%Lf\n", x);
//     return (1e6L*powl(e, x)) + ((435e3L)/(x))*(powl(e, x) - 1) - 1564000.0L;
// }

// //Derivada de la funcion
// double df(double x){

//     return (((pow(e,x))*(1e6L*powl(x,2) + 435e3L*x - 435e3L)) + (435000.0L))/(powl(x,2));
// }


/************Tarea 2 metodos numericos*****************/

double secant(double (*func)(double), double x0, double x1, int n_iters, double tol){

    /*Los puntos no necesitan estar alrededor de la raiz*/
    double x_next;

    for(int i=0; i<n_iters; i++){

        if(fabsl(func(x1) - func(x0)) < 1e-16){

            printf("Error: Div casi por cero\n");
            return -1;

        }

            x_next = x1 - func(x1)*(x1-x0)/(func(x1) - func(x0));

        if(fabsl(x_next-x1)<tol){

            //printf("%Lf\n", x_next);
            return x_next;
            
        }else{

            x0 = x1;
            x1 = x_next;
            
        }
    }
    printf("Iteraciones terminadas.\n");
     return -1;

}

double newton(double (*func)(double), double (*dfunc)(double), double x0, int n_iters, double tol){

    double x_next;
    int i=0;

    x_next = x0 - (func(x0))/(dfunc(x0));

    while(fabsl(x_next-x0)>tol && n_iters>i){

        x0 = x_next;
        x_next = x0 - ((func(x0))/(dfunc(x0)));
        i++;

    }
    return x0;

}

double false_position(double (*func)(double), double a, double b, int n_iters, double tol){

    double c, fc, fa, fb;

    //Inicializo las variables a usar.
    c = b - func(b)*(b-a)/(func(b)-func(a));
    fa = func(a);
    fb = func(b);
    fc = func(c);

    //revisar que tengan signos opuestos
    if((fa*fb)>0){
        return -1;
    }

    for(int i=0; i<n_iters;i++){

        //Calculo el nuevo punto y su evaluacion
        c = b - func(b)*(b-a)/(func(b)-func(a));
        fc = func(c);

        //Verifico si el absoluto es cercano a cero por debajo de la tolerancia
        if(fabsl(fc) < tol){

            printf("Raiz encontrada en %d iters\n", i);
            return c;
        }

        //Si con el punto "a" hay cambio de signo, recorto b, sino reviso al contrario.
        if((fa * fc) < 0){

            fb = fc;
            b = c;

        }else{

            fa = fc;
            a = c;
        }
    }
    printf("Iteraciones agotadas.");
    return -1.0;

}

double bisection(double (*func)(double), double a, double b, int n_iters, double tol){

    double c, fc, fa, fb;

    //Inicializo las variables a usar.
    c = (a+b)/(2);
    fa = func(a);
    fb = func(b);
    fc = func(c);

    //revisar que tengan signos opuestos
    if((fa*fb)>0){
        return -1;
    }

    for(int i=0; i<n_iters;i++){

        //Calculo el nuevo punto y su evaluacion
        c = (a+b)/2;
        fc = func(c);

        //Verifico si el absoluto es cercano a cero por debajo de la tolerancia
        if(fabsl(fc) < tol){

            printf("Raiz encontrada en %d iters\n", i);
            return c;
        }

        //Si con el punto "a" hay cambio de signo, recorto b, sino reviso al contrario.
        if((fa * fc) < 0){

            fb = fc;
            b = c;

        }else{

            fa = fc;
            a = c;
        }
    }
    printf("Iteraciones agotadas.");
    return -1.0;

}

double fix_point(double (*func)(double), double p0_, double n_iters, double tol){

    /*En func va alguna formulacion diferente de f(x)*/

    double p, p0 = p0_;

    for(int i=0; i<n_iters; i++){

        //Cambia dependiento de la funcion 
        p = func(p0);

        if(fabsl(p-p0) < tol){

            printf("raiz encontrada en %d iters\n",i);
            return p;
        
        }else{

            p0 = p;
        }

        //printf("%e\n",p);

    }

    printf("no se encontró raiz\n");

    return -1;

}

/****************************************************** */


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


/*************Tarea 3 sistemas ecuaciones********************* */

void diagonal(double **D, double *b, double *x,  int n){

    /*x es el vector respuesta de salida*/

    for(int i=0; i<n; i++){
        x[i] = b[i] / D[i][i];
    }

}

void lower(double **L, double *x, double *b, int n){

    x[0] = b[0]/L[0][0]; 

    double suma=0;

    for(int i=0; i<n; i++){
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

void gauss_elim(double **A, double *b, double *x, int n){

    //constante escalamiento
    double m;

    double *aux = malloc(n * sizeof(double));
    if(aux == NULL) return;

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

    upper(A, x, b, n);

}

void gauss_pivot(double **A, double *b, double *x, int n){

    //constante escalamiento
    double m;
    //Matriz de permutaciones
    int *perm = malloc(n * sizeof(int));
    //Inicio vector de permutacion.
    for(int i=0; i<n; i++)
    {    perm[i] = i;}
    
    //Matriz auxiliar para hacer operaciones con renglones.
    double *aux = malloc(n * sizeof(double));
    if(aux == NULL) return;


    //Recorro toda la diagonal con i
    for(int i=0; i<n; i++){
        //Con j recorro los renglones de la matriz
        for(int j=i+1; j<n; j++){

            //Si hay un cero en la diagonal, realizar pivoteo
            if(fabs(A[i][i]) < ZERO_VAL){
                
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
    if(x_ordered == NULL){

        printf("No se aloco memoria correctamente para x_ordered\n");
        free(perm);
        return;
    }

    upper(A, x, b, n);

    for (int i = 0; i < n; i++) {
        x_ordered[perm[i]] = x[i];
    }

    for(int i=0; i<n; i++) x[i] = x_ordered[i];

    free(perm);
    free(x_ordered);

}

/************************************************************/


/***************Tareas 4 factorizaciones, Jacobi y GS*************/

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

    /*Manualmente los primeros valores*/
    L[0][0] = sqrt(A[0][0]);
    U[0][0] = sqrt(A[0][0]);

    /*Calculo la primera columna*/
    for(int i=1; i<n; i++){

        L[i][0] = A[i][0]/L[0][0];
        U[0][i] = L[i][0];

    }

    if(isnan(L[0][0])){

        printf("No es posible usar cholesky\n");
        exit(2);
    }  

    /*Empiezo a partir del segundo elemento*/
    for(int j=1; j<n-1; j++){

        double suma = 0;

        /*Primero los valores de la diagonal*/
        for (int k = 0; k < j; k++)
        {
            suma += L[j][k] * L[j][k];
        }

        L[j][j] = sqrt(A[j][j] - suma);
        U[j][j] = L[j][j];

        /*Despues los debajo de la diagonal*/
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

    /*El ultimo se calcula manualmente*/
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

    /*A tambien está en este formato*/
    /*La diagonal principal es en 0, la otra en 1*/

    L[0][0] = sqrt(A[0][0]);

    for(int i=1; i<n; i++){

        L[1][i] = A[1][i]/L[0][i-1];

        L[0][i] = sqrt(A[0][i] - (L[1][i] * L[1][i]));
    }
}

void solve_cholesky_tridiag_vec(double **L, const double *b, double *x, int n){
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

/***************************************************** */


/*****************Tarea5 Potencia y potencia inversa*************/

void k_power_method(int k, double **A, double *v_0, double **eigens, int n, int n_iters, double tol){

    /*Inicio con 0 valores calculados*/
    int n_evals = 0;

    /*Inicio mi contador de iteraciones*/
    int iter=0;

    /*Vector auxiliar siguiente*/
    double *v_n = create_vec(n);

    /*Mi vector inicial lo inicio como unitario*/
    fill_vec(v_0, 1/sqrt(n), n);

    /*Valor que contiene el error entre calculo
    de los eigenvalores*/
    double diff=tol+1;
    double lambda_next;
    double vec_error;


    /*Mientras tenga menos e-valores de los que me pidieron*/
    while(n_evals < k){

        /*Inicio el vector auxiliar en cada calculo*/
        fill_vec(v_0, 1/sqrt(n), n);
        double lambda = 1;

        /*Quito las proyecciones desde el principio*/
        //quitar_proj(v_0, eigens, n_evals, n);

        iter=0;
        diff = tol+1;
        vec_error = tol+1;

        /*Mientras el error sea muy grande y no se haya
        llegado a la iteracion maxima*/
        while((diff > tol) && iter < n_iters){

            /*Quito la proyeccion en v_0*/
            quitar_proj(v_0, eigens, n_evals, n);

            /*Multiplico A*v_0 = v_n*/
            matXvec(A, v_0, v_n, n);

            quitar_proj(v_n, eigens, n_evals, n);

            /*Calculo lambda*/
            double num;
            double den;
            dot_product(v_n, v_n, &num, n);
            dot_product(v_0, v_n, &den, n);
            lambda_next = num/den;

            /*Normalizo v_n*/
            normalize_l2(v_n, n);
            quitar_proj(v_n, eigens, n_evals, n);

            /*Calculo el error*/
            diff = fabs(lambda_next - (lambda));
            diff_vecs(v_0, v_n, &vec_error, n);

            /*El siguiente ahora es el actual*/
            for (int i = 0; i < n; i++){

                v_0[i] = v_n[i];
            }

            iter++;
        }

        /*Si ya convergi a un valor salgo del ciclo
        y guardo cada valor en su respectiva posicion
        De 0 a n-1 filas contienen e-vectores
        En la posicion n estan los e-valores*/
        for(int j=0; j<n; j++){

            /*Guardo el eigenvector*/
            eigens[n_evals][j] = v_0[j];

        }

        /*Guardo el eigenvalor*/
        eigens[n][n_evals] = lambda_next;

        /*Aumento el contador porque ya calcule uno mas*/
        n_evals++;

    }
}

void k_inverse_power_method(int k, double **A, double *v_0, double **eigens, int n, int n_iters, double tol){

    /*Calculo la matriz inversa*/
    double **A_inv = create_matrix(n);
    inverse_cholesky(A, A_inv, n);
    printf("inversa terminada\n");


    /*Inicio con 0 valores calculados*/
    int n_evals = 0;

    /*Inicio mi contador de iteraciones*/
    int iter=0;

    /*Vector auxiliar siguiente*/
    double *v_n = create_vec(n);

    /*Mi vector inicial lo inicio como unitario*/
    fill_vec(v_0, 1/sqrt(n), n);

    /*Valor que contiene el error entre calculo
    de los eigenvalores*/
    double diff=tol+1;

    double lambda_next;


    /*Mientras tenga menos e-valores de los que me pidieron*/
    while(n_evals < k){

        /*Inicio el vector auxiliar en cada calculo*/
        fill_vec(v_0, 1/sqrt(n), n);
        double lambda = 1;

        /*Quito las proyecciones desde el principio*/
        //quitar_proj(v_0, eigens, n_evals, n);

        iter=0;
        diff = tol+1;

        /*Mientras el error sea muy grande y no se haya
        llegado a la iteracion maxima*/
        while(diff > tol && iter < n_iters){

            /*Quito la proyeccion en v_0*/
            quitar_proj(v_0, eigens, n_evals, n);

            /*Multiplico A*v_0 = v_n*/
            matXvec(A_inv, v_0, v_n, n);

            quitar_proj(v_n, eigens, n_evals, n);

            /*Calculo lambda*/
            double num;
            double den;
            dot_product(v_0, v_n, &num, n);
            dot_product(v_n, v_n, &den, n);
            lambda_next = num/den;

            /*Normalizo v_n*/
            normalize_l2(v_n, n);
            quitar_proj(v_n, eigens, n_evals, n);

            /*Calculo el error*/
            diff = fabs(lambda_next - (lambda));

            /*El siguiente ahora es el actual*/
            for (int i = 0; i < n; i++){

                v_0[i] = v_n[i];
            }

            iter++;
        }

        /*Si ya convergi a un valor salgo del ciclo
        y guardo cada valor en su respectiva posicion
        De 0 a n-1 filas contienen e-vectores
        En la posicion n estan los e-valores*/
        for(int j=0; j<n; j++){

            /*Guardo el eigenvector*/
            eigens[n_evals][j] = v_0[j];

        }

        /*Guardo el eigenvalor*/
        eigens[n][n_evals] = lambda_next;

        /*Aumento el contador porque ya calcule uno mas*/
        n_evals++;

    }

    free_matrix(A_inv, n);

}


void k_inverse_power_method2(int k, double **A, double *v_0, double **eigens, int n, int n_iters, double tol){

    /*Calculo la matriz inversa*/
    double **L = create_matrix(n);
    double **LT = create_matrix(n);
    // printf("inversa terminada\n");

    cholesky(A, L, LT, n);


    /*Inicio con 0 valores calculados*/
    int n_evals = 0;

    /*Inicio mi contador de iteraciones*/
    int iter=0;

    /*Vector auxiliar siguiente*/
    double *v_n = create_vec(n);

    /*Mi vector inicial lo inicio como unitario*/
    fill_vec(v_0, 1/sqrt(n), n);

    /*Valor que contiene el error entre calculo
    de los eigenvalores*/
    double diff=tol+1;

    double lambda_next;


    /*Mientras tenga menos e-valores de los que me pidieron*/
    while(n_evals < k){

        /*Inicio el vector auxiliar en cada calculo*/
        fill_vec(v_0, 1/sqrt(n), n);
        double lambda = 1;

        /*Quito las proyecciones desde el principio*/
        //quitar_proj(v_0, eigens, n_evals, n);

        iter=0;
        diff = tol+1;

        /*Mientras el error sea muy grande y no se haya
        llegado a la iteracion maxima*/
        while(diff > tol && iter < n_iters){

            /*Quito la proyeccion en v_0*/
            quitar_proj(v_0, eigens, n_evals, n);

            /*Resuelvo A*v_n = v_0*/
            solve_LU(L, LT, v_n, v_0, n);

            quitar_proj(v_n, eigens, n_evals, n);

            /*Calculo lambda*/
            double num;
            double den;
            dot_product(v_0, v_n, &num, n);
            dot_product(v_n, v_n, &den, n);
            lambda_next = num/den;

            /*Normalizo v_n*/
            normalize_l2(v_n, n);
            quitar_proj(v_n, eigens, n_evals, n);

            /*Calculo el error*/
            diff = fabs(lambda_next - (lambda));

            /*El siguiente ahora es el actual*/
            for (int i = 0; i < n; i++){

                v_0[i] = v_n[i];
            }

            iter++;
        }

        /*Si ya convergi a un valor salgo del ciclo
        y guardo cada valor en su respectiva posicion
        De 0 a n-1 filas contienen e-vectores
        En la posicion n estan los e-valores*/
        for(int j=0; j<n; j++){

            /*Guardo el eigenvector*/
            eigens[n_evals][j] = v_0[j];

        }

        /*Guardo el eigenvalor*/
        eigens[n][n_evals] = lambda_next;

        /*Aumento el contador porque ya calcule uno mas*/
        n_evals++;

    }

    free_matrix(L, n);
    free_matrix(LT, n);

}




void inverse_cholesky(double **A, double **A_inv, int n){

    double **L, **LT;

    L = create_matrix(n);
    LT = create_matrix(n);

    cholesky(A, L, LT, n);

    /*Creo mi vector identidad auxiliar*/
    double *col_identity = create_vec(n);
    double *col_inverse = create_vec(n);

    for(int i=0; i<n; i++){

        /*Inicio el vector b para la primera columnda de la inversa */
        col_identity[i] = 1;

        /*Resuelvo el sistema factorizado con b = i-esima columna de I*/
        solve_LU(L, LT, col_inverse, col_identity, n);

        /*Ahora que tengo mi columna de la inversa, la copio*/
        for(int j=0; j<n; j++){

            A_inv[j][i] = col_inverse[j];

        }

        col_identity[i] = 0;
    }

    free(col_identity);
    free(col_inverse);
    free_matrix(L, n);
    free_matrix(LT, n);
}

/*****************************************************************/



/***********Tarea 6 iteracion subespacio y gradiente conjugado******/

void power_method_mod(double **A, double **eigens, int i_0, int n_evals, int n, int n_iters){

    /*n_evals, me dice cuantos eigenvectores ya tengo, estos estan guardados como las filas de eigens.
    v_0 va a ser mi vectori inicial y eigens tendra los otros vectores a restar*/

    /*Creo una copia del vector en i_0 para trabajar con la copia*/
    double *v_0 = create_vec(n);
    for(int i = 0; i<n; i++) v_0[i] = eigens[i_0][i];

    /*Inicio mi contador de iteraciones*/
    int iter=0;

    /*Vector auxiliar siguiente*/
    double *v_n = create_vec(n);

    /*Valor que contiene el error entre calculo
    de los eigenvalores*/
    double diff=1;
    //double lambda_next;
    double vec_error;

    /************Solo quiero calcular un eigenvector*********************/

    double lambda = 1;

    iter=0;
    diff =1;
    vec_error =1;

    /*Mientras no se hagan las iteraciones sugeridas, seguir.
    Solo se hacen unas cuantas iteraciones, no muchisimas*/
    while(iter < n_iters){

        /*Quito las proyecciones en v_0, de los otros eigenvectores, v_0
        es el eigenvalor con indice i_0, por eso no le quito las contribuciones de ese
        sino seria quitarse a si mismo*/
        quitar_proj_except(v_0, eigens, i_0, n_evals, n);

        /*Multiplico A*v_0 = v_n*/
        matXvec(A, v_0, v_n, n);

        /*A vn tambien le quito las contribuciones de los demas
        menos del que quiero modificar*/
        quitar_proj_except(v_n, eigens, i_0, n_evals, n);

        /*Calculo lambda*/
        // double num;
        // double den;
        // dot_product(v_n, v_n, &num, n);
        // dot_product(v_0, v_n, &den, n);
        // lambda_next = num/den;

        /*Normalizo v_n y le quito la respectiva proyeccion*/
        normalize_l2(v_n, n);

        /*Calculo el error*/
        // diff = fabs(lambda_next - (lambda));
        diff_vecs(v_0, v_n, &vec_error, n);

        /*El siguiente ahora es el actual*/
        for (int i = 0; i < n; i++){

            v_0[i] = v_n[i];
        }

        iter++;

    }

    /*Aunque no haya convergido guardo el vector y el valor
    y guardo cada valor en su respectiva posicion
    De 0 a n-1 filas contienen e-vectores*/
    for(int j=0; j<n; j++){

        /*Guardo el eigenvector*/
        eigens[i_0][j] = v_0[j];

    }

    /***************Solo calculo un eigenvector********************* 
      porque el el que voy a sustituir dentro de eigens
    */


    free(v_0);
    free(v_n);
}

void jacobi_rotations(double **A, int n, double **e_vecs, int max_iters, double tol){

    if(n<2){
        printf("Tamaño invalido de la matriz.\n");
    }

    int n_iters = 0;
    double max_val = fabs(A[1][0]);
    int max_i, max_j;
    double curr_val;
    double *v_aux = create_vec(n);


    /*Mientras haya iteraciones disponibles
    y el mayor valor fuera de la diagonal sea muy grande*/
    while(n_iters < max_iters && max_val > tol){

        max_val = 0;
        max_i = 1;  
        max_j = 0;


        /*Encuentro el valor maximo fuera de la diagonal y su posicion en la matriz
        como es simetrica solo tengo que revisar hasta antes de llegar a la diagonal*/
        for(int i=1; i<n; i++){
            for(int j=0; j<i; j++){

                curr_val = fabs(A[i][j]);
                if(curr_val > max_val){

                    max_val = curr_val;
                    max_i = i; max_j = j;

                }

            }

        }

        if(max_val <= tol) break;

        /*Calculo el angulo de rotacion*/
        double alpha = A[max_i][max_i] - A[max_j][max_j];
        double beta = 2.0 * A[max_i][max_j];

        double theta = (1.0L/2.0L)*atan2(beta, alpha);


        /*Calculo las trigonometricas del angulo*/
        double c = cos(theta);
        double s = sin(theta);

        /*Guardo valores antes de modificar la matriz*/
        double temp_aii = A[max_i][max_i];
        double temp_ajj = A[max_j][max_j];
        double temp_aij = A[max_i][max_j];

        /*Calculo los elementos nuevos en la diagonal*/
        A[max_i][max_i] = (temp_aii*c*c) + (temp_ajj*s*s) + (2*temp_aij*c*s);
        A[max_j][max_j] = (temp_aii*s*s) + (temp_ajj*c*c) - (2*temp_aij*c*s);

        /*Los elementos pivote se hacer cero porque en esto nos basamos para el calculo de theta*/
        A[max_i][max_j] = 0.0;
        A[max_j][max_i] = 0.0;

        /*Modfico los valores que no son el pivote ni la diagonal*/
        for(int k=0; k<n; k++){

            if(k != max_i && k != max_j){

                double temp_aik = A[max_i][k];
                double temp_ajk = A[max_j][k];


                A[max_i][k] = (temp_aik*c) + (temp_ajk*s);
                A[k][max_i] = A[max_i][k];

                A[max_j][k] = (temp_ajk*c) - (temp_aik*s);
                A[k][max_j] = A[max_j][k]; 

            }

        }

        /*Con jacobi voy a tener los eigenvectores como filas de la matriz identidad
        no como columnas*/

        for(int k=0; k<n; k++){

            double temp_eki = e_vecs[k][max_i];
            double temp_ekj = e_vecs[k][max_j];


            e_vecs[k][max_i] = temp_eki*c + temp_ekj*s;
            e_vecs[k][max_j] = -temp_eki*s + temp_ekj*c;

        }

        n_iters++;

    }

    free(v_aux);
}

void gradient_conj(double **A, double *x, double *b, int n, int max_iters, double tol){

    double *r = create_vec(n);
    double *p = create_vec(n);
    double *Ap = create_vec(n);

    /* r = b - A*x */
    matXvec(A, x, r, n);
    substract_vec(b, r, r, n); 

    /* p = r */
    for (int i = 0; i < n; ++i) p[i] = r[i];

    /* rr_prod = r*r para la norma l2
    rr < tol^2 para que termine */
    double rr_prod;
    dot_product(r, r, &rr_prod, n); 

    double tol2 = tol * tol;

    int iter = 0;

    while (rr_prod > tol2 && iter < max_iters) {
        /*************Inicia el algoritmo******************* */
        /* Ap = A * p */
        matXvec(A, p, Ap, n);

        double denom;
        dot_product(p, Ap, &denom, n);

        double alpha = rr_prod / denom;

        /* x = x + alpha * p */
        numXvec(alpha, p, Ap, n);  
        add_vecs(Ap, x, x, n);        

        /*********Aqui ya tengo actualizado x********/

        /* r = r - alpha * A*p */
        matXvec(A, p, Ap, n);       
        numXvec(alpha, Ap, Ap, n);   
        substract_vec(r, Ap, r, n);  
        /*r contiene el nuevo residuo*/

        double rr_next;
        dot_product(r, r, &rr_next, n);

        double beta = rr_next / rr_prod;

        /* p = r + beta * p */
        numXvec(beta, p, Ap, n);   
        add_vecs(r, Ap, p, n);  

        rr_prod = rr_next;
        iter++;

    }

    free(r);
    free(p);
    free(Ap);

}

void iter_subspace(double **A, int rowA, int colA, double **Phi, int rowPhi, int colPhi, 
                                                double *lambdas, int n_iters, double tol){


        /*Los eigenvectores van a terminar en Phi, pero de derecha a izquierda
        Tambie los eigenvalores, porque en la matriz Bmm salen en la diagonal 
        de abajo para arribe*/

        /*B_mm cuadrada va a contener los eigenvalores*/
        double **Bmm = create_matrix(colPhi);

        /*Phi_mn * A_nn = Maux_mn*/
        double **M_aux = create_rect(colPhi, colA);

        /*Auxiliar para Phi transpuesta*/
        double **PhiT = create_rect(colPhi, rowPhi);

        /*Matriz para guardar los eigenvectores de Bmm
        tiene que ser unitaria al principio*/
        double **Phi_v = create_identity(colPhi);
        double **Phi_next = create_rect(rowPhi, colPhi);



        for(int iters=0; iters<n_iters; iters++){

            /****************Aqui es donde empieza el algoritmo****************/
            if(iters % 100 == 0) printf("iteracion %d\n", iters);

            /*PhiT * A = Maux -> Maux * Phi = Bmm*/
            /*Genero mi matriz Bmm*/
            matXmat_rect(PhiT, colPhi, rowPhi, A, rowA, colA, M_aux);
            matXmat_rect(M_aux, colPhi, colA, Phi, rowPhi, colPhi, Bmm);


            /*Reinicio Phi_v con la identidad, que va a guardar los eigenvecs*/
            for(int i=0; i<colPhi; i++){
                for(int j=0; j<colPhi; j++){

                    if(i == j) Phi_v[i][j] = 1.0;
                    else Phi_v[i][j] = 0.0;

                }
            }

            /*Bmm se modifica y contiene los eigenvalores pero solo quiero los eigenvectores
            Los valores de Bmm no me interesan, Bmm contiene en su diagonal los eigenvalores
            Phi_v tiene que ser la identidad inicialmente para que guarde los valores correctos*/
            jacobi_rotations(Bmm, colPhi, Phi_v, 1000, 1e-8);

            // transpose(Phi_v, colPhi);

            /*Genero siguiente iteracion Phi * Phi_v = Phi_next*/
            matXmat_rect(Phi, rowPhi, colPhi, Phi_v, colPhi, colPhi, Phi_next);

            /**************Aqui ya tengo calculado Phi_next***************************/

            /*Despues de calcular Phi_next ahora quiero aplicar el metodo de la potencia
            a cada una de las columnas de Phi next*/

            /*Guardo el nuevo Phi (Phi_next), en la transpuesta para tener los eigenvectores
             como filas y poder trabajar con las funciones hechas
            La informacion que habia en la transpuesta ya no la vamos a usar asi que nos sirve*/
            transpose_rect(Phi_next, rowPhi, colPhi, PhiT, colPhi, rowPhi);

            /*A cada columna de Phi (fila de PhiT) le hago el metodo de potencia modificado.
            Solo unas pocas iteraciones*/
            for(int i=0; i<colPhi; i++){

                power_method_mod(A, PhiT, i, colPhi, rowPhi, 3);

            }

            /*Ortogonalizo PhiT, sus filas quedan ortonormales*/
            ortogonalize_mat(PhiT, colPhi, rowPhi);

            /*Guardo las modificaciones en Phi para la siguiente iteracion*/
            transpose_rect(PhiT, colPhi, rowPhi, Phi, rowPhi, colPhi);

            /*Aqui tengo PhiT con la informacion nueva
            Phi con la informacion nueva
            Y eso es todo lo que importa para la siguiente iteracion*/

            /*Aplico potencia a cada columna de Phi_next, pocas iters
            A * colPhi
            El eigen vector obtenido, va a ser la nueva columna de Phi
            por lo que tengo que copiar los valores en colPhi y despues meterla en Phi*/

        }
        
        for(int i=0; i<colPhi; i++){

            lambdas[i] = Bmm[i][i];

        }

        // printf("Bmm final:\n");
        // display_matrix(Bmm, colPhi);

        // printf("Phi final:\n");
        // display_mat_rect(Phi, rowPhi, colPhi);

        free_matrix(Phi_v, colPhi);
        free_matrix(Bmm, colPhi);
        free_matrix(M_aux, colPhi);
        free_matrix(PhiT, colPhi);
        free_matrix(Phi_next, rowPhi);


}


void inv_iter_subspace(double **A, int rowA, int colA, double **Phi, int rowPhi, int colPhi, 
                                                    double *lambdas, int n_iters, double tol){

    /*B_mm cuadrada va a contener los eigenvalores*/
    double **Bmm = create_matrix(colPhi);

    /*Phi_mn * A_nn = Maux_mn*/
    double **M_aux = create_rect(colPhi, colA);

    /*Auxiliar para Phi transpuesta*/
    double **PhiT = create_phi(colPhi, rowPhi);

    /*Matriz para guardar los eigenvectores de Bmm
    tiene que ser unitaria al principio*/
    double **Phi_v = create_identity(colPhi);
    double **Phi_next = create_rect(rowPhi, colPhi);

    double **A_inv = create_matrix(rowA);
    inverse_cholesky(A, A_inv, rowA);
    // printf("Matriz inversa de A:\n");
    // display_matrix(A_inv, rowA);

    // printf("Matriz A:\n");
    // display_matrix(A, rowA);


    for(int iters=0; iters<n_iters; iters++){

        /************************Aqui es donde empieza el algoritmo**************************** */

        /*PhiT * A = Maux -> Maux * Phi = Bmm*/
        /*Genero mi matriz Bmm*/
        matXmat_rect(PhiT, colPhi, rowPhi, A_inv, rowA, colA, M_aux);
        matXmat_rect(M_aux, colPhi, colA, Phi, rowPhi, colPhi, Bmm);


        /*Reinicio Phi_v con la identidad, que va a guardar los eigenvecs*/
        for(int i=0; i<colPhi; i++){
            for(int j=0; j<colPhi; j++){

                if(i == j) Phi_v[i][j] = 1.0;
                else Phi_v[i][j] = 0.0;

            }
        }

        /*Bmm se modifica y contiene los eigenvalores pero solo quiero los eigenvectores
        Los valores de Bmm no me interesan, Bmm contiene en su diagonal los eigenvalores
        Phi_v tiene que ser la identidad inicialmente para que guarde los valores correctos*/
        jacobi_rotations(Bmm, colPhi, Phi_v, 1000, 1e-12);

        /*No es necesario corregir la inversa con jacobi, solo hasta el final, porque no usamos
        los eigenvalores para hacer calculos, solo los vectores*/
        //for(int i=0; i<colPhi; i++) Bmm[i][i] = 1.0/Bmm[i][i];

        /*Genero siguiente iteracion Phi * Phi_v = Phi_next*/
        matXmat_rect(Phi, rowPhi, colPhi, Phi_v, colPhi, colPhi, Phi_next);

        /**************Aqui ya tengo calculado Phi_next***************************/

        /*Despues de calcular Phi_next ahora quiero aplicar el metodo de la potencia
        a cada una de las columnas de Phi next*/

        /*Guardo el nuevo Phi (Phi_next), en la transpuesta para tener los eigenvectores como filas y poder trabajar con las funciones hechas
        La informacion que habia en la transpuesta ya no la vamos a usar asi que nos sirve*/
        transpose_rect(Phi_next, rowPhi, colPhi, PhiT, colPhi, rowPhi);

        /*A cada columna de Phi (fila de PhiT) le hago el metodo de potencia modificado, unas pocas iteraciones*/
        for(int i=0; i<colPhi; i++){

            power_method_mod(A_inv, PhiT, i, colPhi, rowPhi, 3);

        }

        /*Ortogonalizo PhiT, sus filas quedan ortonormales*/
        ortogonalize_mat(PhiT, colPhi, rowPhi);

        /*Guardo las modificaciones en Phi para la siguiente iteracion*/
        transpose_rect(PhiT, colPhi, rowPhi, Phi, rowPhi, colPhi);

        /*Aqui tengo PhiT con la informacion nueva
        Phi con la informacion nueva
        Y eso es todo lo que importa para la siguiente iteracion*/

        /*Aplico potencia a cada columna de Phi_next, pocas iters
        A * colPhi
        El eigen vector obtenido, va a ser la nueva columna de Phi
        por lo que tengo que copiar los valores en colPhi y despues meterla en Phi*/

    }
    
    /*Obtenemos los eigenvalores correctos*/
    for(int i=0; i<colPhi; i++){

        lambdas[i] = 1.0/Bmm[i][i];

    }

    // printf("Bmm final:\n");
    // display_matrix(Bmm, colPhi);

    // printf("Phi final:\n");
    // display_mat_rect(Phi, rowPhi, colPhi);

    free_matrix(Phi_v, colPhi);
    free_matrix(Bmm, colPhi);
    free_matrix(M_aux, colPhi);
    free_matrix(PhiT, colPhi);
    free_matrix(Phi_next, rowPhi);
    free_matrix(A_inv, rowA);


}


/********************************************************************** */


/*********************Tarea 7 QR******************************/


void QR(double **A, double **Q, double **R, int n){
    /*Toma como entrada la matriz A, que factoriza Q y R*/

    /*Primero hay que ortogonalizar A*/

    double *a = create_vec(n);
    double *a_prev = create_vec(n);
    double dot_p = 0;

    for(int i=0; i<n; i++){

        /*Tomo la i-esima columna de A*/
        get_column(A, a, n, i);

        for(int j=0; j<i; j++){
            /*Quito las contribuciones anteriores*/
            get_column(Q, a_prev, n, j);

            /*Producto punto entre el actual y uno de los anteriores*/
            dot_product(a, a_prev, &dot_p, n);

            R[j][i] = dot_p;

            /*Multiplico el producto punto por el previo*/
            numXvec(dot_p, a_prev, a_prev, n);

            /*Le quito la contribucion al actual*/
            substract_vec(a, a_prev, a, n);


        }

        dot_product(a, a, &dot_p, n);
        R[i][i] = sqrt(dot_p);
        normalize_l2(a, n);
        set_column(Q, a, n, i);

    }

    free(a);
    free(a_prev);

}

void transpose_inplace(double **Q, double **QT, int n){

    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){

            QT[j][i] = Q[i][j];

        }

    }

}



void solveQR(double **A, double *x, double *b, int n){


    double **Q = create_matrix(n);
    double **R = create_matrix(n);
    double **QT = create_matrix(n);
    double *y = create_vec(n);

    QR(A, Q, R, n);
    transpose_inplace(Q, QT, n);                                                                                                                                                                                             
    matXvec(QT, b, y, n);
    upper(R, x, y, n);
    display_matrix(Q, n);
    display_matrix(QT, n);

    free_matrix(Q, n);
    free_matrix(R, n);
    free_matrix(QT, n);
    free(y);

}


























#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "eigens.h"
#include "utils.h"
#define ZERO_VAL 1e-10

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

void solve_LU(double **L, double **U, double *x, double *b, int n){

    double *y = create_vec(n);
    lower(L, y, b, n);

    upper(U, x, y, n);
    free(y);
}




void normalize_inf(double *v, int n){

    /*Busco el valor mas grande del vector*/
    double max = fabs(v[0]);

    for(int i=1; i<n; i++){

        /*Si encuentro un valor mas grande lo intercambio*/
        if(fabs(v[i]) > max){

            max = fabs(v[i]);

        }
    }

    /*Normalizo todo el vector con el maximo*/
    for(int i=0; i<n; i++){
        v[i] = v[i]/max;
    }
}

void normalize_l2(double *v, int n){

    double norm2 = 0.0;
    for (int i = 0; i < n; ++i) norm2 += v[i] * v[i];

    norm2 = sqrt(norm2);

    if (norm2 <= 1e-16) return; /* vector casi nulo */
    for (int i = 0; i < n; ++i) v[i] /= norm2;

}

void diff_vecs(double *u, double *v, double *diff_vec, int n){

    double max = fabs(u[0] - v[0]);

    for(int i=1; i<n; i++){

        double val = fabs(u[i] - v[i]);

        if(val > max){
            max = val;
        }
    }

    (*diff_vec)=max;
}

void quitar_proj(double *v, double **e_vecs, int n_vecs, int n){

    /*Si no e calculado ningun eigenvec, no hago nada*/
    if(n_vecs == 0){
        return;
    }

    /*Variables auxiliares*/
    double dot_prod = 0;
    double *v_aux = create_vec(n);

    /*Quito la contribucion de los n_vecs ya calculados*/
    for(int i=0; i<n_vecs; i++){

        /*Producto punto entre el vector y el i-eigenvec
        lo guardo en la variable dot_prod*/
        dot_product(v, e_vecs[i], &dot_prod, n);

        /*Multiplico el dot_prod por el i-eigenvec
        y lo guardo en el v_aux*/
        numXvec(dot_prod, e_vecs[i], v_aux, n);

        /*Al vector original le resto la contribucion
        del i-eigenvec*/
        for(int j=0; j<n; j++){

            v[j] = v[j] - v_aux[j];

        }

    }

    free(v_aux);
}




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




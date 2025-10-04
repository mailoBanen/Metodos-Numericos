#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "funciones.h"
#include "utils.h"
#define ZERO_VAL 1e-14


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
            jacobi_rotations(Bmm, colPhi, Phi_v, 1000, 1e-12);

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
        col_identity[i] = 1.0;

        /*Resuelvo el sistema factorizado con b = i-esima columna de I*/
        solve_LU(L, LT, col_inverse, col_identity, n);

        /*Ahora que tengo mi columna de la inversa, la copio*/
        for(int j=0; j<n; j++){

            A_inv[j][i] = col_inverse[j];

        }

        col_identity[i] = 0.0;
    }

    free(col_identity);
    free(col_inverse);
    free_matrix(L, n);
    free_matrix(LT, n);
}




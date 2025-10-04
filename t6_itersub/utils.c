#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define ZERO_VAL 1e-10

double **read_matrix(char *source, int *rows, int *cols){

    FILE *file = fopen(source, "r");

    if(file == NULL){
        return NULL;
    }

    fscanf(file, "%d", rows);
    fscanf(file, "%d", cols);

    double **arr = malloc((*rows) * sizeof(double *));
    if(arr == NULL){
        fclose(file);
        return NULL;
    }

    for(int i=0; i<(*rows); i++){
        //Creo memoria para cada renglon y la almaceno en arr
        arr[i] = malloc((*cols) * sizeof(double));

        for(int j=0; j<(*cols); j++){

            if(fscanf(file, "%lf", &arr[i][j]) != 1){ //Seria lo mismo si pongo *(arr+i)+j
                printf("Error al leer la matriz, revisa el .txt\n");
                break;
            }
        }
    }

    fclose(file);

    return arr;
}

double *read_vec(char *source, int *n){

    FILE *file = fopen(source, "r");

    if(file == NULL){
        return NULL;
    }

    fscanf(file, "%d", n);

    //Creo un puntero de tamaño n, y sino retorna nulo
    double *arr = (double *)malloc((*n) * sizeof(double));
    if(arr == NULL){
        fclose(file);
        return NULL;
    }

    for(int i=0; i<(*n); i++){

            //Escaneo el archivo y voy posicionando los valores en el arreglo
            if(fscanf(file, "%lf", &arr[i]) != 1){ //Seria lo mismo si pongo *(arr+i)+j
                printf("Error al leer la matriz, revisa el .txt\n");
                break;
            }
        
    }

    fclose(file);

    return arr;
}




double **create_matrix(int n){

    double **L = (double **)malloc(n * sizeof(double *));
    
    if(L == NULL){
        printf("Error, no se creó la matriz\n");
        return NULL;
    }
    
    for(int i=0; i<n; i++){

        L[i] = calloc(n, sizeof(double));

        if(L[i] == NULL){
        printf("Error, no se creó la matriz\n");
        for(int m=0; m<i; m++) free(L[m]);
        free(L);
        return NULL;
        }

    }

    return L;

}

double *create_vec(int n){

    double *V = (double *)calloc(n, sizeof(double));
    if(V == NULL){
        printf("Error creando vector de tamaño %d\n", n);
        return NULL;
    } 
    return V;
}

double **create_identity(int n){

    double **I = create_matrix(n);

    for(int i=0; i<n; i++){

        I[i][i] = 1.0L;
        
    }

    return I;
}

double max_arg_vec(double *v, int n){
    /*De un vector v tamaño n, retornar el numero mas grande*/

    double max = fabs(v[0]);
    double curr_val;

    for (int i = 1; i < n; i++)
    {
        
        curr_val = fabs(v[i]);
        if(curr_val > max){

            max = curr_val;

        }
        
    }

    return max;
    
}

double **create_phi(int n, int m){
    /*Con esto debo de crear una matriz de tamaño nxm que tenga columnas ortonormales
    simplemente una matriz identidad cortada*/

    if(n>=m){

        double **Phi = (double **)malloc(n * sizeof(double *));

        if(Phi==NULL){
            printf("Error, no se creó la matriz Phi\n"); 
            return NULL;
        } 

        for(int i=0; i<n; i++){

            Phi[i] = calloc(m, sizeof(double));

            if(Phi[i] == NULL){

                printf("Error, no se creó la matriz Phi\n"); 

                for(int j=0; j<i; j++){
                    free(Phi[j]);
                }
                free(Phi);
                return NULL;

            }

        }

        for(int i=0; i<m; i++){
            Phi[i][i] = 1.0;
        }

        return Phi;

    }else{


        double **Phi = (double **)malloc(n * sizeof(double *));

        if(Phi==NULL){
            printf("Error, no se creó la matriz Phi\n"); 
            return NULL;
        } 

        for(int i=0; i<n; i++){

            Phi[i] = calloc(m, sizeof(double));

            if(Phi[i] == NULL){

                printf("Error, no se creó la matriz Phi\n"); 

                for(int j=0; j<i; j++){
                    free(Phi[j]);
                }
                free(Phi);
                return NULL;

            }

        }

        for(int i=0; i<n; i++){
            Phi[i][i] = 1.0;
        }

        return Phi;

    }

}

double **create_rect(int n, int m){

    if(n>=m){

        double **Phi = (double **)malloc(n * sizeof(double *));

        if(Phi==NULL){
            printf("Error, no se creó la matriz Phi\n"); 
            return NULL;
        } 

        for(int i=0; i<n; i++){

            Phi[i] = calloc(m, sizeof(double));

            if(Phi[i] == NULL){

                printf("Error, no se creó la matriz Phi\n"); 

                for(int j=0; j<i; j++){
                    free(Phi[j]);
                }
                free(Phi);
                return NULL;

            }

        }

        return Phi;

    }else{


        double **Phi = (double **)malloc(n * sizeof(double *));

        if(Phi==NULL){
            printf("Error, no se creó la matriz Phi\n"); 
            return NULL;
        } 

        for(int i=0; i<n; i++){

            Phi[i] = calloc(m, sizeof(double));

            if(Phi[i] == NULL){

                printf("Error, no se creó la matriz Phi\n"); 

                for(int j=0; j<i; j++){
                    free(Phi[j]);
                }
                free(Phi);
                return NULL;

            }

        }

        return Phi;

    }

}

void fill_vec(double *vec, double num, int n){
    for(int i=0; i<n; i++){
        vec[i] = num;
    }
}





void free_matrix(double **M, int n) {
    for (int i = 0; i < n; i++) free(M[i]);
    free(M);
}

void display_matrix(double **M, int n){

    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){

            printf("%.4lf ", M[i][j]);

        }
        printf("\n");
    }
    printf("\n");
}

void display_vec(double *b, int n){

    for(int i=0; i<n; i++){
            printf("%e\n", b[i]);
    }
    printf("\n");
}

void display_mat_rect(double **Phi, int n, int m){


    for(int i=0; i<n; i++){
        for(int j=0; j<m; j++){

            printf("%.4lf ", Phi[i][j]);

        }
        printf("\n"); 
    }
    printf("\n");

}




void get_column(double **A, double *column, int n, int j_col){

    for(int i=0; i<n; i++){

        column[i] = A[i][j_col];

    }

}

void set_column(double **A, double *column, int n, int j_col){

    for(int i=0; i<n; i++){

        A[i][j_col] = column[i];

    }

}




void transpose_rect(double **A, int a_rows, int a_cols, double **AT, int at_rows, int at_cols){

    if(a_rows != at_cols || a_cols != at_rows){
        printf("No son dimensiones similares.\n");
        return;
    }

    for(int i=0; i<a_rows; i++){
        for(int j=0; j<a_cols; j++){

            AT[j][i] = A[i][j];


        }
    }

}

void multiply(double **A, double **B, double **output, int n){

    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            
            double valor = 0;
            for(int k=0; k<n; k++){

                valor += A[i][k] * B[k][j];

            }

            output[i][j] = valor;

        }
    }
}

void matXvec(double **A, double *x, double *vout, int n){

    for(int i=0; i<n; i++){
        double suma=0;
        for(int j=0; j<n; j++){

            suma += A[i][j] * x[j];

        }

        vout[i] = suma;
    }
}

void dot_product(double *x, double *y, double *n_out, int n){

    double suma=0;

    for(int i=0; i<n; i++){

        suma += x[i] * y[i];

    }

    *n_out = suma;

}

void numXvec(double num, double *x, double *x_out, int n){
    
    for(int i=0; i<n; i++){
        x_out[i] = x[i]*num; 
    }

}

void substract_vec(double *u, double *v, double *v_out, int n){

    for(int i=0; i<n; i++){

        v_out[i] = u[i] - v[i];

    }

}

void add_vecs(double *u, double *v, double *v_out, int n){

    for (int i = 0; i < n; i++)
    {
        v_out[i] = u[i] + v[i];
    }

}

void matXmat_rect(double **A, int rowA, int colA, 
                           double **B, int rowB, int colB, double **C){
    
    //Verificar la condición columnas de A == filas de B
    if (colA != rowB) {
        printf("Error: El numero de columnas de A debe ser igual al numero de filas de B.\n");
        return;
    }

    //El resultado C tendrá dimensiones (filA x colB)

    int filC = rowA;
    int colC = colB;
    
    //C[i][j] = sum (A[i][k] * B[k][j]) para k desde 0 hasta colA-1 (o filB-1)
    
    for (int i = 0; i < rowA; i++) {      
        for (int j = 0; j < colB; j++) {   
            C[i][j] = 0.0;                 

            for (int k = 0; k < colA; k++) {
                
                C[i][j] += A[i][k] * B[k][j];

            }
        }
    }
}




void outter_product(double *x, double *y, double **A_out, int n){

    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){

            A_out[i][j] = x[i]*y[j];

        }
    }
}

void numXmat(double num, double **A, double **A_out, int n){

    for(int i=0; i<n; i++){
        double *Ai = A[i];
        double *Ai_out = A_out[i];
        for(int j=0; j<n; j++){

            Ai_out[j] = Ai[j]*num;

        }
    }
}

void normalize_l2(double *v, int n){

    double norm2 = 0.0;
    for (int i = 0; i < n; ++i) norm2 += v[i] * v[i];

    norm2 = sqrt(norm2);

    if (norm2 <= 1e-16) return; /* vector casi nulo */
    for (int i = 0; i < n; ++i) v[i] /= norm2;

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

void diff_mat(double **A, double **B, double *diff_mat, int n){

    double max = fabs(A[0][0] - B[0][0]);

    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){

            double curr_val = fabs(A[i][j] - B[i][j]);
            if(curr_val > max){

                max = curr_val;

            }


        }
        
    }
    

    (*diff_mat) = max;

}




int is_equal_m(double **A, double **B, int n){

    double error_max=0;

    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){

            if(fabs(A[i][j] - B[i][j]) > error_max){
                error_max = fabs(A[i][j] - B[i][j]);
            }

            if(fabs(A[i][j] - B[i][j]) > ZERO_VAL){

                printf("Las matrices no son iguales.\n");
                printf("A[%d][%d] = %f  B[%d][%d] = %f\n", i, j, A[i][j], i, j, B[i][j]);
                return 0;
            }
        }
    }

    printf("Error: %e\n", error_max);

    return 1;
}

int is_equal_v(double *x, double *y, int n){

    double error_max = 0;

    for(int i=0; i<n; i++){
        double diff = fabs(x[i] - y[i]);
        
        if(diff > error_max){
            error_max = diff;
        }

        if(diff > ZERO_VAL){
            printf("Vectores no iguales - Primera diferencia en índice %d: %e\n", i, diff);
            printf("Error máximo encontrado: %e\n", error_max);
            return 0;
        }
    }

    printf("Son iguales - Error máximo: %e\n", error_max);
    return 1;
}

int is_identity(double **A, int n){

    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){

            if(i==j){

                if(fabs(A[i][j]-1) > ZERO_VAL){

                    printf("No es la matriz identidad\n");
                    return -1;

                }
            }else{

                if(fabs(A[i][j]) > ZERO_VAL){

                    return -1;
                }

            }
        }
    }

    return 1;
}

int error_evec(double **A, double *v, double lambda, int n){

    /*Vectores de resultados*/
    double *r1 = create_vec(n);
    double *r2 = create_vec(n);

    /*Calculos Av = lambda*v*/
    matXvec(A, v, r1, n);
    numXvec(lambda, v, r2, n);

    /*comparo vectores para encontrar los resultados maximos*/
    double max_err = fabs(r1[0] - r2[0]);
    double curr_err;

    for(int i=1; i<n; i++){

        curr_err = fabs(r1[i] - r2[i]);
        if(curr_err > max_err){

            max_err = curr_err;
        }

    }

    printf("La diferencia maxima es: %e\n", max_err);

    free(r1);
    free(r2);

    return 0;
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

void ortogonalize_mat(double **eigens, int rows, int cols){

    double *v_aux = create_vec(cols);
    double dot_prod = 0;

    /*Quiero procesar cada fila en orden*/
    for (int i = 0; i < rows; i++) {
        
        /* Hacer la fila i ortogonal a todas las filas anteriores 0 a i-1 */
        for (int j = 0; j < i; j++) {
            
            /* Calcular proyección de A[i] sobre A[j] */
            dot_product(eigens[i], eigens[j], &dot_prod, cols);
            
            /* Restar la proyección: A[i] = A[i] - proj */
            numXvec(dot_prod, eigens[j], v_aux, cols);

            for (int k = 0; k < cols; k++) {
                eigens[i][k] = eigens[i][k] - v_aux[k];
            }
        }
        
        // Normalizar la fila i después de hacerla ortogonal
        normalize_l2(eigens[i], cols);
    }

    free(v_aux);
}

void quitar_proj_except(double *v, double **e_vecs, int except_i, int n_vecs, int n){

    /*A un vector le quiero quitar todas las proyecciones, excepto una*/
    /*Si no e calculado ningun eigenvec, no hago nada*/
    if(n_vecs == 0){
        return;
    }

    /*Variables auxiliares*/
    double dot_prod = 0;
    double *v_aux = create_vec(n);

    /*Quito la contribucion de los n_vecs ya calculados*/
    for(int i=0; i<n_vecs; i++){

        /*A mi vector v le quito todas las proyecciones, menos la del vector con
        el indice except_i, porque v está en el idice except_i */
        if(i != except_i){

            /*Producto punto entre el vector y el i-eigenvec
            lo guardo en la variable dot_prod*/
            dot_product(v, e_vecs[i], &dot_prod, n);

            /*Al vector original le resto la contribucion
            del i-eigenvec*/
            for(int j=0; j<n; j++){

                v[j] = v[j] - (e_vecs[i][j] * dot_prod);

            }

        }

    }
    free(v_aux);
}

void err_solution(double **A, double *x, double *b, int n){

    double *r = create_vec(n);
    double diff_ = 0;

    matXvec(A, x, r, n);

    diff_vecs(r, b, &diff_, n);

    printf("El error maximo del sistema es: %e\n", diff_);

    free(r);

}







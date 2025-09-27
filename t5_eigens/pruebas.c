#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// Definición de constantes
#define N 3          // Dimensión de la matriz (ej. 3x3)
#define MAX_ITER 1000 // Número máximo de iteraciones
#define TOLERANCIA 1e-6 // Criterio de convergencia

// -------------------------------------------------------------------
// 1. FUNCIONES AUXILIARES DE VECTORES Y MATRICES
// -------------------------------------------------------------------

// Inicializa un vector con valores aleatorios
void inicializar_vector(double *v) {
    for (int i = 0; i < N; i++) {
        v[i] = (double)rand() / RAND_MAX; // Valor entre 0 y 1
    }
}

// Normaliza un vector (lo convierte a longitud unitaria)
double norma(double *v) {
    double sum = 0.0;
    for (int i = 0; i < N; i++) {
        sum += v[i] * v[i];
    }
    return sqrt(sum);
}

void normalizar(double *v) {
    double norm = norma(v);
    if (norm < 1e-12) return; // Evitar división por cero
    for (int i = 0; i < N; i++) {
        v[i] /= norm;
    }
}

// Producto punto (dot product)
double producto_punto(double *v1, double *v2) {
    double sum = 0.0;
    for (int i = 0; i < N; i++) {
        sum += v1[i] * v2[i];
    }
    return sum;
}

// Multiplicación matriz-vector: z = A * w
void mult_matriz_vector(double A[N][N], double *w, double *z) {
    for (int i = 0; i < N; i++) {
        z[i] = 0.0;
        for (int j = 0; j < N; j++) {
            z[i] += A[i][j] * w[j];
        }
    }
}

// -------------------------------------------------------------------
// 2. FUNCIÓN DE DEFLACIÓN (Gram-Schmidt Modificado)
// -------------------------------------------------------------------

// -------------------------------------------------------------------
// 2. FUNCIÓN DE DEFLACIÓN (Gram-Schmidt Modificado)
//    (Mantenemos la misma función, solo cambiamos dónde se llama)
// -------------------------------------------------------------------

// Aplica el proceso de Gram-Schmidt para eliminar la proyección
// del vector w sobre los eigenvectores ya encontrados (V).
// En el método de potencia, w es el vector resultante (z).
void deflacion_gram_schmidt(double *w, double **V, int num_encontrados) {
    for (int k = 0; k < num_encontrados; k++) {
        // Coeficiente de proyección: alpha = dot(w, v_k)
        double alpha = producto_punto(w, V[k]);
        
        // Restar la proyección: w = w - alpha * v_k
        for (int i = 0; i < N; i++) {
            w[i] -= alpha * V[k][i];
        }
    }
    // NOTA: No normalizamos aquí, la normalización se hace después en el bucle principal.
}

// -------------------------------------------------------------------
// 3. MÉTODO DE LA POTENCIA (Eigenvalores más Grandes) - CORREGIDO
// -------------------------------------------------------------------

void metodo_potencia_mas_grandes(double A[N][N], int k_eigen, double *eigenvalores, double **eigenvectores) {
    
    double *w_actual = (double*)malloc(N * sizeof(double));
    double *z = (double*)malloc(N * sizeof(double)); // z = A * w
    
    for (int k = 0; k < k_eigen; k++) {
        inicializar_vector(w_actual);
        normalizar(w_actual);
        
        for (int iter = 0; iter < MAX_ITER; iter++) {
            // Guardar el w_actual anterior para la condición de parada
            double *w_anterior = (double*)malloc(N * sizeof(double));
            for(int i=0; i<N; i++) w_anterior[i] = w_actual[i];
            
            // 1. POTENCIA: z = A * w_actual
            mult_matriz_vector(A, w_actual, z);
            
            // 2. DEFLACIÓN (Gram-Schmidt): ¡CORRECCIÓN CLAVE!
            //    Asegura que el vector 'z' esté en el subespacio ortogonal.
            deflacion_gram_schmidt(z, eigenvectores, k); 
            
            // 3. NORMALIZACIÓN: w_actual = z / ||z||
            for(int i=0; i<N; i++) w_actual[i] = z[i];
            normalizar(w_actual);

            // 4. CRITERIO DE PARADA
            double diff = 0.0;
            for(int i=0; i<N; i++) diff += pow(w_actual[i] - w_anterior[i], 2);
            free(w_anterior);

            if (sqrt(diff) < TOLERANCIA) {
                // 5. CÁLCULO DEL EIGENVALOR (Cociente de Rayleigh)
                mult_matriz_vector(A, w_actual, z);
                eigenvalores[k] = producto_punto(w_actual, z);
                
                // Almacenar el eigenvector encontrado
                for(int i=0; i<N; i++) eigenvectores[k][i] = w_actual[i];
                break;
            }
        }
        // Si el bucle termina por MAX_ITER, el código simplemente guarda la última aproximación
    }
    
    free(w_actual);
    free(z);
}

// -------------------------------------------------------------------
// 4. MÉTODO DE LA POTENCIA INVERSA (Eigenvalores más Chicos)
// -------------------------------------------------------------------

// Esta función requiere una matriz inversa (A_inversa)
// En la práctica, se utiliza una factorización LU para resolver el sistema A*z = w
// Aquí, para simplicidad del ejemplo, asumiremos que A_inversa ya fue calculada.
void metodo_potencia_inversa_mas_chicos(double A_inv[N][N], int k_eigen, double *eigenvalores, double **eigenvectores) {

    // (El proceso es idéntico al de la potencia, solo que se usa A_inv)
    double *w_actual = (double*)malloc(N * sizeof(double));
    double *z = (double*)malloc(N * sizeof(double));
    
    for (int k = 0; k < k_eigen; k++) {
        inicializar_vector(w_actual);
        normalizar(w_actual);
        
        for (int iter = 0; iter < MAX_ITER; iter++) {
            deflacion_gram_schmidt(w_actual, eigenvectores, k);
            
            double *w_anterior = (double*)malloc(N * sizeof(double));
            for(int i=0; i<N; i++) w_anterior[i] = w_actual[i];
            
            // Multiplicación por la inversa: z = A_inv * w
            mult_matriz_vector(A_inv, w_actual, z);
            
            // Normalización
            for(int i=0; i<N; i++) w_actual[i] = z[i];
            normalizar(w_actual);
            
            // Criterio de Parada
            double diff = 0.0;
            for(int i=0; i<N; i++) diff += pow(w_actual[i] - w_anterior[i], 2);
            free(w_anterior);

            if (sqrt(diff) < TOLERANCIA) {
                // El eigenvalor de A es el recíproco del eigenvalor de A_inv
                // lambda_inv = dot(w_actual, A_inv * w_actual)
                mult_matriz_vector(A_inv, w_actual, z);
                double lambda_inv = producto_punto(w_actual, z);
                
                eigenvalores[k] = 1.0 / lambda_inv; // El valor propio de la matriz A
                
                // Almacenar el eigenvector
                for(int i=0; i<N; i++) eigenvectores[k][i] = w_actual[i];
                break;
            }
        }
    }
    
    free(w_actual);
    free(z);
}

// -------------------------------------------------------------------
// 5. FUNCIÓN PRINCIPAL (MAIN) y EJEMPLO
// -------------------------------------------------------------------

int main() {
    srand(time(NULL)); // Inicializar la semilla aleatoria
    
    // Matriz de ejemplo (simétrica para garantizar ortogonalidad de eigenvectores)
    // Eigenvalores: 4, 1, -1. Eigenvectores (sin normalizar): (1, 0, 1), (0, 1, 0), (1, 0, -1)
    double A[N][N] = {
        {1, 0, 3},
        {0, 1, 0},
        {3, 0, 1}
    };

    // Matriz inversa A^-1 (usada para la potencia inversa)
    // Eigenvalores: 1/4, 1, -1.
    double A_inv[N][N] = {
        {-0.25, 0.0, 0.75},
        { 0.0,  1.0, 0.0},
        { 0.75, 0.0, -0.25}
    };
    
    int k_to_find = 2; // Queremos encontrar los 2 más grandes y los 2 más chicos
    
    // Asignación de memoria para resultados
    double *e_val_grande = (double*)malloc(k_to_find * sizeof(double));
    double **e_vec_grande = (double**)malloc(k_to_find * sizeof(double*));
    double *e_val_chico = (double*)malloc(k_to_find * sizeof(double));
    double **e_vec_chico = (double**)malloc(k_to_find * sizeof(double*));

    for(int i=0; i<k_to_find; i++) {
        e_vec_grande[i] = (double*)malloc(N * sizeof(double));
        e_vec_chico[i] = (double*)malloc(N * sizeof(double));
    }
    
    // ----------------------------------------------------------------
    // EJECUTAR MÉTODOS
    // ----------------------------------------------------------------
    
    // I. Más Grandes (Método de la Potencia)
    metodo_potencia_mas_grandes(A, k_to_find, e_val_grande, e_vec_grande);

    // II. Más Chicos (Método de la Potencia Inversa)
    metodo_potencia_inversa_mas_chicos(A_inv, k_to_find, e_val_chico, e_vec_chico);
    
    // ----------------------------------------------------------------
    // IMPRIMIR RESULTADOS
    // ----------------------------------------------------------------
    
    printf("--- Eigenvalores MAS GRANDES (Metodo de la Potencia con Deflacion) ---\n");
    for (int k = 0; k < k_to_find; k++) {
        printf("Lambda_%d: %.4f, Vector: (", k + 1, e_val_grande[k]);
        for (int i = 0; i < N; i++) {
            printf("%.4f%s", e_vec_grande[k][i], (i == N - 1) ? "" : ", ");
        }
        printf(")\n");
    }
    
    printf("\n--- Eigenvalores MAS CHICOS (Metodo de la Potencia Inversa con Deflacion) ---\n");
    for (int k = 0; k < k_to_find; k++) {
        printf("Lambda_%d: %.4f, Vector: (", k + 1, e_val_chico[k]);
        for (int i = 0; i < N; i++) {
            printf("%.4f%s", e_vec_chico[k][i], (i == N - 1) ? "" : ", ");
        }
        printf(")\n");
    }

    // Liberar memoria
    for(int i=0; i<k_to_find; i++) {
        free(e_vec_grande[i]);
        free(e_vec_chico[i]);
    }
    free(e_vec_grande);
    free(e_val_grande);
    free(e_vec_chico);
    free(e_val_chico);
    
    return 0;
}
#ifndef OPERATIONS_MAT_H
#define OPERATIONS_MAT_H


double **read_matrix(char *source, int n);
double *read_vec(char *source, int n);

void display_matrix(double **M, int n);
void display_vec(double *b, int n);

void free_matrix(double **M, int n);
void free_vec(double *b);

void sum_arrs(double *arr1, double *arr2, double *result, double  n);
void mult_arr(double *arr1, double k, double *result, int n);


//Operaciones
double *diagonal(double **D, double *b, int n);
double *lower(double **L, double *b, int n);
double *upper(double **U, double *b, int n);
double *gauss_elim(double **A, double *b, int n);

#endif
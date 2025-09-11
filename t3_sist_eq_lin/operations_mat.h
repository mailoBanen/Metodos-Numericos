#ifndef OPERATIONS_MAT_H
#define OPERATIONS_MAT_H


double **read_matrix(char *source, int n);
double *read_vec(char *source, int n);

void display_matrix(double **M, int n);
void display_vec(double *b, int n);

void free_matrix(double **M, int n);
void free_vec(double *b);


//Operaciones
double *diagonal(double **D, double *b, int size);

#endif
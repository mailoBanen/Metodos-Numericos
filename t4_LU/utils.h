#ifndef UTILS_H
#define UTILS_H


double **read_matrix(char *source, int n);
double *read_vec(char *source, int n);



double **create_matrix(int n);
double *create_vec(int n);
void fill_vec(double *vec, double num, int n);


void free_matrix(double **M, int n);
void display_matrix(double **M, int n);
void display_vec(double *b, int n);

void multiply(double **A, double **B, double **output, int n);
void matXvec(double **A, double *x, double *vout, int n);

int is_equal_m(double **A, double **B, int n);
int is_equal_v(double *x, double *y, int n);

#endif
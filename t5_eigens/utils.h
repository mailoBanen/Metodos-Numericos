#ifndef UTILS_H
#define UTILS_H


double **read_matrix(char *source, int *rows, int *cols);
double *read_vec(char *source, int n);



double **create_matrix(int n);
double *create_vec(int n);
void fill_vec(double *vec, double num, int n);


void free_matrix(double **M, int n);
void display_matrix(double **M, int n);
void display_vec(double *b, int n);

void multiply(double **A, double **B, double **output, int n);
void matXvec(double **A, double *x, double *vout, int n);
void dot_product(double *x, double *y, double *n_out, int n);
void numXvec(double num, double *x, double *x_out, int n);
void outter_product(double *x, double *y, double **A_out, int n);
void numXmat(double num, double **A, double **A_out, int n);

int is_equal_m(double **A, double **B, int n);
int is_equal_v(double *x, double *y, int n);
int is_identity(double **A, int n);

#endif
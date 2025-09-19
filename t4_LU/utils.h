#ifndef UTILS_H
#define UTILS_H


double **read_matrix(char *source, int n);
double **create_matrix(int n);
void free_matrix(double **M, int n);
void display_matrix(double **M, int n);
void multiply(double **A, double **B, double **output, int n);
int is_equal(double **A, double **B, int n);

#endif
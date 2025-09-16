#ifndef UTILS_H
#define UTILS_H


double **read_matrix(char *source, int n);
double **create_matrix(int n);
void free_matrix(double **M, int n);
void display_matrix(double **M, int n);

#endif
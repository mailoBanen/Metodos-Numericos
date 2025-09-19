#ifndef DECOMPOSITION_H
#define DECOMPOSITION_H

void crouth(double **A, double **L, double **U, int n);

void cholesky(double **A, double **L, double **U, int n);

void cholesky_tridiag(double **A, double **L, double **LT, int n);

#endif 
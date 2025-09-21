#ifndef DECOMPOSITION_H
#define DECOMPOSITION_H


void lower(double **L, double *x, double *b, int n);
void upper(double **U, double *x, double *b, int n);



void crouth(double **A, double **L, double **U, int n);
void solve_crouth(double **L, double **U, double *x, double *b, int n);
void solve_LU(double **L, double **U, double *x, double *b, int n);


void cholesky(double **A, double **L, double **U, int n);
void cholesky_tridiag(double **A, double **L, double **LT, int n);
void cholesky_tridiag_vec(double **A, double **L, int n);
void solve_cholesky_tridiag_vec(double **L, const double *b, double *x, int n);


void jacobi(double **A, double *x, double *b, int n, int n_iters, double tolerance);
void gauss_seidel(double **A, double *x, double *b, int n, int n_iters, double tolerance);


#endif 
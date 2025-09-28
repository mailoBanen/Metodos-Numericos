#ifndef EIGENS_H
#define EIGENS_H

void lower(double **L, double *x,  double *b, int n);
void upper(double **U, double *x, double *b, int n);

void cholesky(double **A, double **L, double **U, int n);
void solve_LU(double **L, double **U, double *x, double *b, int n);
void normalize_inf(double *v, int n);
void quitar_proj(double *v, double **e_vecs, int n_vecs, int n);

void k_power_method(int k, double **A, double *v_0, double **eigens, int n, int n_iters, double tol);
void k_inverse_power_method(int k, double **A, double *v_0, double **eigens, int n, int n_iters, double tol);
void k_biggest_evals(int k, double **A, double *x, double *lambda, int n, int n_iters, double tol);


void inverse_cholesky(double **A, double **A_inv, int n);


#endif
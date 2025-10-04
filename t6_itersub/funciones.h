#ifndef FUNCIONES_H
#define FUNCIONES_H


void power_method_mod(double **A, double **eigens,int i_0, int n_evals, int n, int n_iters);
void jacobi_rotations(double **A, int n, double **e_vecs, int max_iters, double tol);
// void gradient_conj(double **A, double *x_0, double *b, int n, int max_iters, double tol);
void gradient_conj(double **A, double *x, double *b, int n, int max_iters, double tol);
void iter_subspace(double **A, int rowA, int colA, double **Phi, int rowPhi, int colPhi, double *lambdas, int n_iters, double tol);
void inv_iter_subspace(double **A, int rowA, int colA, double **Phi, int rowPhi, int colPhi, double *lambdas, int n_iters, double tol);


void lower(double **L, double *x,  double *b, int n);
void upper(double **U, double *x, double *b, int n);

void cholesky(double **A, double **L, double **U, int n);
void solve_LU(double **L, double **U, double *x, double *b, int n);


void k_power_method(int k, double **A, double *v_0, double **eigens, int n, int n_iters, double tol);
void k_inverse_power_method(int k, double **A, double *v_0, double **eigens, int n, int n_iters, double tol);
void inverse_cholesky(double **A, double **A_inv, int n);



#endif
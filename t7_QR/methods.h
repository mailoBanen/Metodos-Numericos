#ifndef METHODS_H
#define METHODS_H


#define PI 3.14159265358979323846
#define E 2.71828182845904523536
#define GOLDEN_RATIO 1.618033988749895


double secant(double (*func)(double), double x0, double x1, int n_iters, double tol);
double newton(double (*func)(double), double (*dfunc)(double), double x0, int n_iters, double tol);
double false_position(double (*func)(double), double a, double b, int n_iters, double tol);
double bisection(double (*func)(double), double a, double b, int n_iters, double tol);
double fix_point(double (*func)(double), double p0_, double n_iters, double tol);



void sum_arrs(double *arr1, double *arr2, double *result, double  n);
void mult_arr(double *arr1, double k, double *result, int n);
void swap_row(double **A, double *b, int i, int in, int n);
void swap_col(double **A, int *perm, int j, int jn, int n);
void find_max(double **A, double *b, int i, int j, int *perm, int n);



void diagonal(double **D, double *b, double *x, int n);
void lower(double **L, double *x, double *b, int n);
void upper(double **U, double *x, double *b, int n);
void gauss_elim(double **A, double *b, double *x, int n);
void gauss_pivot(double **A, double *b, double *x, int n);




void crouth(double **A, double **L, double **U, int n);
void solve_crouth(double **L, double **U, double *x, double *b, int n);
void solve_LU(double **L, double **U, double *x, double *b, int n);
void cholesky(double **A, double **L, double **U, int n);
void cholesky_tridiag(double **A, double **L, double **LT, int n);
void cholesky_tridiag_vec(double **A, double **L, int n);
void solve_cholesky_tridiag_vec(double **L, const double *b, double *x, int n);
void jacobi(double **A, double *x, double *b, int n, int n_iters, double tolerance);
void gauss_seidel(double **A, double *x, double *b, int n, int n_iters, double tolerance);



void k_power_method(int k, double **A, double *v_0, double **eigens, int n, int n_iters, double tol);
void k_inverse_power_method(int k, double **A, double *v_0, double **eigens, int n, int n_iters, double tol);
void k_inverse_power_method2(int k, double **A, double *v_0, double **eigens, int n, int n_iters, double tol);
void inverse_cholesky(double **A, double **A_inv, int n);


void power_method_mod(double **A, double **eigens, int i_0, int n_evals, int n, int n_iters);
void jacobi_rotations(double **A, int n, double **e_vecs, int max_iters, double tol);
void gradient_conj(double **A, double *x, double *b, int n, int max_iters, double tol);
void iter_subspace(double **A, int rowA, int colA, double **Phi, int rowPhi, int colPhi, 
                                                double *lambdas, int n_iters, double tol);

void inv_iter_subspace(double **A, int rowA, int colA, double **Phi, int rowPhi, int colPhi, 
                                                    double *lambdas, int n_iters, double tol);
                                                    







void QR(double **A, double **Q, double **R, int n);
void transpose_inplace(double **Q, double **QT, int n);
void solveQR(double **A, double *x, double *b, int n);








#endif
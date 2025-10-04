#ifndef UTILS_H
#define UTILS_H


double **read_matrix(char *source, int *rows, int *cols);
double *read_vec(char *source, int *n);



double **create_matrix(int n);
double *create_vec(int n);
double **create_identity(int n);
double max_arg_vec(double *v, int n);
double **create_phi(int n, int m);
double **create_rect(int n, int m);


void fill_vec(double *vec, double num, int n);


void free_matrix(double **M, int n);
void display_matrix(double **M, int n);
void display_vec(double *b, int n);
void display_mat_rect(double **Phi, int n, int m);
void get_column(double **A, double *column, int n, int j_col);
void set_column(double **A, double *column, int n, int j_col);
void transpose_rect(double **A, int arows, int acols, double **AT, int atrows, int atcols);



void multiply(double **A, double **B, double **output, int n);
void matXvec(double **A, double *x, double *vout, int n);
void dot_product(double *x, double *y, double *n_out, int n);
void numXvec(double num, double *x, double *x_out, int n);
void substract_vec(double *u, double *v, double *v_out, int n);
void add_vecs(double *u, double *v, double *v_out, int n);
void matXmat_rect(double **MatrizA, int filA, int colA, 
                           double **MatrizB, int filB, int colB, double **C);

void outter_product(double *x, double *y, double **A_out, int n);
void numXmat(double num, double **A, double **A_out, int n);
void normalize_l2(double *v, int n);
void normalize_inf(double *v, int n);
void diff_vecs(double *u, double *v, double *diff_vec, int n);
void diff_mat(double **A, double **B, double *diff_mat, int n);




int is_equal_m(double **A, double **B, int n);
int is_equal_v(double *x, double *y, int n);
int is_identity(double **A, int n);
int error_evec(double **A, double *v, double lambda, int n);
void quitar_proj(double *v, double **e_vecs, int n_vecs, int n);
void ortogonalize_mat(double **eigens, int rows, int cols);
void quitar_proj_except(double *v, double **e_vecs, int except_i, int n_vecs, int n);
void err_solution(double **A, double *x, double *b, int n);


#endif
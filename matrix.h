#define EPS 1e-9

enum { STR_SIZE = 256 };

double **create_matrix(int n);
void free_matrix(double **matrix, int n);
double *create_vector(int n);
void free_vector(double *vector);
double **copy_matrix(double **matrix, int n);
double **create_np1_matrix(double **matrix, double *vector, int n);
double **create_nx2n_matrix(double **matrix, int n);
double **inverse_matrix(double **matrix, int n);
void fill_matrix(double **matrix, double *vector, int n, FILE *fp);
int find_row(double **matrix, int n, int row);
void swap_rows(double **r_1, double **r_2);
int find_max_element(double **matrix, int n, int row);
void swap_elements(double *elm_1, double *elm_2);
void swap_columns(double **matrix, int n, int clm_1, int clm_2);
void strsub(double **matrix, int n, int row_1, int row_2, double k);
void gauss_direct(double **matrix, int num_row, int num_clm, int row);
void gauss_reverse(double **matrix, int num_row, int num_clm, int row);
double determinant(double **matrix, int n);
double *gauss(double **matrix, double *vector, int n);
double *modified_gauss(double **matrix, double *vector, int n);
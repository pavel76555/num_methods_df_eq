#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "matrix.h"

double **create_matrix(int n)
{
    double **matrix = calloc(n, sizeof(*matrix));
    for (int i = 0; i < n; i++) {
        matrix[i] = calloc(n, sizeof(double));
    }

    return matrix;
}

double *create_vector(int n)
{
    double *vector = calloc(n, sizeof(*vector));
    return vector;
}

void free_matrix(double **matrix, int n)
{
    for (int i = 0; i < n; i++) {
        free(matrix[i]);
    }

    free(matrix);
}

void free_vector(double *vector)
{
    free(vector);
}

double **copy_matrix(double **matrix, int n)
{
    double **matrix_copy = create_matrix(n);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            matrix_copy[i][j] = matrix[i][j];
        }
    }

    return matrix_copy;
}

int find_row(double **matrix, int n, int row)
{
    for (int i = row + 1; i < n; i++) {
        if (fabs(matrix[i][row]) > EPS) {
            return i;
        }
    }

    return row;
}

void swap_rows(double **r_1, double **r_2)
{
    double *tmp = *r_1;
    *r_1 = *r_2;
    *r_2 = tmp;
}

void swap_elements(double *elm_1, double *elm_2)
{
    double tmp = *elm_1;
    *elm_1 = *elm_2;
    *elm_2 = tmp;
}

void swap_columns(double **matrix, int n, int clm_1, int clm_2)
{
    for (int i = 0; i < n; i++) {
        swap_elements(&matrix[i][clm_1], &matrix[i][clm_2]);
    }
}

void strsub(double **matrix, int n, int row_1, int row_2, double k)
{
    for (int i = 0; i < n; i++) {
        matrix[row_2][i] -= matrix[row_1][i] * k;
    }
}

int find_max_element(double **matrix, int n, int row)
{
    int i_max = row;
    double max = fabs(matrix[row][i_max]);

    for (int i = row + 1; i < n; i++) {
        if (fabs(matrix[row][i]) > max) {
            i_max = i;
            max = fabs(matrix[row][i]);
        }
    }

    return i_max;
}

void fill_matrix(double **matrix, double *vector, int n, FILE *fp)
{
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            fscanf(fp, "%lf", &matrix[i][j]);
        }

        fscanf(fp, "%lf", &vector[i]);
    }
}

double **create_np1_matrix(double **matrix, double *vector, int n)
{
    double **np1_matrix = copy_matrix(matrix, n);

    for (int i = 0; i < n; i++) {
        np1_matrix[i] = realloc(np1_matrix[i], (n + 1) * sizeof(double));
        np1_matrix[i][n] = vector[i];
    }

    return np1_matrix;
}

double **create_nx2n_matrix(double **matrix, int n)
{
    double **nx2n_matrix = copy_matrix(matrix, n);

    for (int i = 0; i < n; i++) {
        nx2n_matrix[i] = realloc(nx2n_matrix[i], (2 * n) * sizeof(double));
    }

    for (int i = 0; i < n; i++) {
        for (int j = n; j < 2 * n; j++) {
            if (i + n == j) {
                (nx2n_matrix[i][j] = 1);
            } else {
                (nx2n_matrix[i][j] = 0);
            }
        }
    }

    return nx2n_matrix;
}

double **inverse_matrix(double **matrix, int n)
{
    double **nx2n_matrix = create_nx2n_matrix(matrix, n);

    for (int i = 0; i < n; i++) {
        gauss_direct(nx2n_matrix, n, 2 * n, i);
    }

    for (int i = 0; i < n; i++) {
        gauss_reverse(nx2n_matrix, n, 2 * n, i);
    }

    double **inv_matrix = create_matrix(n);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            inv_matrix[i][j] = nx2n_matrix[i][j + n];
        }
    }

    free_matrix(nx2n_matrix, n);

    return inv_matrix;
}

void gauss_direct(double **matrix, int num_row, int num_clm, int row)
{
    if (matrix[row][row] == 0) {
        int c_row = find_row(matrix, num_row, row);
        swap_rows(&matrix[row], &matrix[c_row]);
    }

    double res = matrix[row][row];

    for (int i = row; i < num_clm; i++) {
        matrix[row][i] /= res;
    }

    for (int i = row + 1; i < num_row; i++) {
        strsub(matrix, num_clm, row, i, matrix[i][row]);
    }
}

void gauss_reverse(double **matrix, int num_row, int num_clm, int row)
{
    for (int i = row + 1; i < num_row; i++) {
        strsub(matrix, num_clm, i, row, matrix[row][i]);
    }
}

double determinant(double **matrix, int n)
{
    double det = 1;
    double **matrix_copy = copy_matrix(matrix, n);

    for (int i = 0; i < n; i++) {
        int k;
        if (fabs(matrix_copy[i][i]) < EPS) {
            k = find_row(matrix_copy, n, i);

            swap_rows(&matrix_copy[i], &matrix_copy[k]);
            det *= -1;
        }

        for (int j = i + 1; j < n; j++) {
            if (fabs(matrix_copy[i][i]) < EPS) {
                det = 0;
                free_matrix(matrix_copy, n);
                return det;

            }

            strsub(matrix_copy, n, i, j, matrix_copy[j][i] / matrix_copy[i][i]);
        }
    }

    for (int i = 0; i < n; i++) {
        det *= matrix_copy[i][i];
    }

    free_matrix(matrix_copy, n);

    return det;
}

double *gauss(double **matrix, double *vector, int n)
{
    double **np1_matrix = create_np1_matrix(matrix, vector, n);

    for (int i = 0; i < n; i++) {
        gauss_direct(np1_matrix, n, n + 1, i);
    }

    for (int i = 0; i < n; i++) {
        gauss_reverse(np1_matrix, n, n + 1, i);
    }

    double *ans = create_vector(n);

    for (int i = 0; i < n; i++) {
        ans[i] = np1_matrix[i][n];
    }

    free_matrix(np1_matrix, n);

    return ans;
}

double *modified_gauss(double **matrix, double *vector, int n)
{
    double **np1_matrix = create_np1_matrix(matrix, vector, n);

    int *order = calloc(n, sizeof(*order));

    for (int i = 0; i < n; i++) {
        order[i] = i;
    }

    for (int k = 0; k < n; k++) {
        if (np1_matrix[k][k] == 0) {
            int c_row = find_row(matrix, n, k);
            swap_rows(&np1_matrix[k], &np1_matrix[c_row]);
        }

        int i_max = find_max_element(np1_matrix, n, k);
        swap_columns(np1_matrix, n, k, i_max);


        int tmp = order[k];
        order[k] = order[i_max];
        order[i_max] = tmp;

        double res = np1_matrix[k][k];

        for (int i = k; i < n + 1; i++) {
            np1_matrix[k][i] /= res;
        }

        for (int i = k + 1; i < n; i++) {
            strsub(np1_matrix, n + 1, k, i, np1_matrix[i][k]);
        }
    }

    for (int i = 0; i < n; i++) {
        gauss_reverse(np1_matrix, n, n + 1, i);
    }

    double *ans = create_vector(n);

    for (int i = 0; i < n; i++) {
        ans[i] = np1_matrix[order[i]][n];
    }

    free_matrix(np1_matrix, n);
    free(order);

    return ans;
}

double norm(double **matrix, int n)
{
    double max = 0;
    for (int i = 0; i < n; i++) {
        max += matrix[i][0];
    }

    double cur_max;

    for (int j = 1; j < n; j++) {
        cur_max = 0;
        for (int i = 0; i < n; i++) {
            cur_max += matrix[i][j];
        }

        if (cur_max > max) {
            max = cur_max;
        }
    }

    return max;
}

double condition_num(double **matrix, int n)
{
    double **inv_matrix = inverse_matrix(matrix, n);
    double ans = norm(matrix, n) * norm(inv_matrix, n);
    free_matrix(inv_matrix, n);
    return ans;
}

int main(int argc, char const *argv[])
{
    if (argc < 2) {
        printf("No file!\n");
        return 0;
    }

    char str[STR_SIZE] = {};
    snprintf(str, STR_SIZE, "tests/%s", argv[1]);
    FILE *fp = fopen(str, "r");
    int n;
    fscanf(fp, "%d", &n);

    double **matrix = create_matrix(n);
    double *vector = create_vector(n);

    fill_matrix(matrix, vector, n, fp);
    fclose(fp);

    double det = determinant(matrix, n);
    printf("Определитель: %-14g\n\n", det);
    if (det != 0) {

        printf("Число обусловленности: %g\n\n", condition_num(matrix, n));
        printf("Обратная матрица:\n");

        double **inv_matrix = inverse_matrix(matrix, n);

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                printf("%-14g ", inv_matrix[i][j]);
            }
            printf("\n");
        }

        printf("\n");
        free_matrix(inv_matrix, n);

        printf("Решение системы:\n");

        double *ans = gauss(matrix, vector, n);

        if (ans != NULL) {
            for (int i = 0; i < n; i++) {
                printf("x%d = %g\n", i + 1, ans[i]);
            }
        }

        printf("\n");
        free(ans);

        printf("Решение системы мод. методом Гаусса: \n");

        ans = modified_gauss(matrix, vector, n);

        if (ans != NULL) {
            for (int i = 0; i < n; i++) {
                printf("x%d = %g\n", i + 1, ans[i]);
            }
        }

        printf("\n");

        free(ans);
    }

    free_matrix(matrix, n);
    free_vector(vector);

    return 0;
}

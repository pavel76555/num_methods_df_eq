#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define EPS 1e-9
enum { STR_SIZE = 256 };

double **create_matrix(int n)
{
    double **matrix = calloc(n, sizeof(*matrix));
    for (int i = 0; i < n; i++) {
        matrix[i] = calloc(n, sizeof(double));
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            matrix[i][j] = 0.;
        }
    }

    return matrix;
}

double *create_vector(int n)
{
    double *vector = calloc(n, sizeof(*vector));
    for (int i = 0; i < n; i++) {
        vector[i] = 0.;
    }

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

void fill_matrix(double **matrix, double *vector, int n, FILE *fp)
{
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            fscanf(fp, "%lf", &matrix[i][j]);
        }

        fscanf(fp, "%lf", &vector[i]);
    }
}

double **mul_matr_matr(double **matrix_1, double **matrix_2, int n)
{
    double **mul = create_matrix(n);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {
                mul[i][j] += matrix_1[i][k] * matrix_2[k][j];
            }
        }
    }

    return mul;
}

double *mul_matr_vect(double **matrix, double *vector, int n)
{
    double *mul = create_vector(n);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            mul[i] += matrix[i][j] * vector[j];
        }
    }

    return mul;
}

double **transpose_matrix(double **matrix, int n)
{
    double **trp = create_matrix(n);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            trp[i][j] = matrix[j][i];
        }
    }

    return trp;
}

double difference(double *vector_1, double *vector_2, int n)
{
    double res = 0;

    for (int i = 0; i < n; ++i) {
        res += (vector_1[i] - vector_2[i]) * (vector_1[i] - vector_2[i]);
    }

    return sqrt(res);
}

void relaxation(double **matrix, double *vector, int n, double w)
{
    int cnt = 0;
    double **trp_matrix = transpose_matrix(matrix, n);
    double **mul_matrix = mul_matr_matr(trp_matrix, matrix, n);
    double *vector_change = mul_matr_vect(trp_matrix, vector, n);
    double *vector_ans = create_vector(n);
    double *vector_prev = create_vector(n);

    for (int i = 0; i < n; i++) {
        vector_ans[i] = 0;
        vector_prev[i] = 1;
    }

    while (difference(vector_ans, vector_prev, n) >= EPS) {
        cnt++;

        for (int i = 0; i < n; i++) {
            vector_prev[i] = vector_ans[i];
        }

        for (int i = 0; i < n; i++) {
            double sum = 0;

            for (int j = 0; j < i; j++) {
                sum += (mul_matrix[i][j] * vector_ans[j]);
            }

            for (int j = i; j < n; j++) {
                sum += (mul_matrix[i][j] * vector_prev[j]);
            }

            vector_ans[i] = w * (vector_change[i] - sum) / mul_matrix[i][i] + vector_prev[i];
        }
    }

    // printf("Число итераций: %d\n\n", cnt);
    // printf("Значение w: %g\n\n", w);

    // for (int i = 0; i < n; i++) {
    //     printf("x%d = %lf\n", i + 1, vector_ans[i]);
    // }

    printf("%.2lf %d\n", w, cnt);

    free_matrix(trp_matrix, n);
    free_matrix(mul_matrix, n);
    free_vector(vector_change);
    free_vector(vector_prev);
    free_vector(vector_ans);
}

int
main(int argc, char const **argv)
{
    if (argc < 2) {
        printf("No file!\n");
        return 0;
    }

    double w;
    // scanf("%lf", &w);
    // if (w <= 0. || w >= 2.) {
    //     printf("Incorrect w\n");
    //     return 0;
    // }

    char str[STR_SIZE] = {};
    snprintf(str, STR_SIZE, "tests/%s", argv[1]);
    FILE *fp = fopen(str, "r");
    int n;
    fscanf(fp, "%d", &n);

    double **matrix = create_matrix(n);
    double *vector = create_vector(n);

    fill_matrix(matrix, vector, n, fp);
    fclose(fp);

    for (w = 0.05; w < 2.0; w += 0.05) {
        relaxation(matrix, vector, n, w);
    }

    // relaxation(matrix, vector, n, w);

    free_matrix(matrix, n);
    free_vector(vector);
    return 0;
}

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double p(double x) {
    return -1;
}

double q(double x) {
    return 0;
}

double f(double x) {
    return 0;
}

double
find_a(double x_0, double h, int i)
{
    return (h * h * q(x_0 + h * i) - h * p(x_0 + h * i) + 1);
}

double
find_b(double x_0, double h, int i)
{
    return (h * p(x_0 + h * i) - 2);
}

double
find_c(double x_0, double h, int i)
{
    return 1;
}

double
find_d(double x_0, double h, int i)
{
    return (-(h * h * f(x_0 + i * h)));
}


double*
create_vector(int size)
{
    double *vector = (double*)malloc(size * sizeof(double));
    return vector;
}

void
destroy_vector(double *vector)
{
    free(vector);
}

void
fin_dif(double (*p)(double), double (*q)(double), double (*f)(double),
        double sigma_1, double gamma_1, double delta_1,
        double sigma_2, double gamma_2, double delta_2,
        double x_0, double x_n, int n)
{
    double h = (x_n - x_0) / n;

    double *alpha = create_vector(n + 1);
    double *beta = create_vector(n + 1);

    alpha[0] = -gamma_1 / (h * sigma_1 - gamma_1);
    beta[0] = (h * delta_1) /(h * sigma_1 - gamma_1);

    for (int i = 1; i < n; i++) {
        alpha[i] = -find_c(x_0, h, i) / (find_a(x_0, h, i) * alpha[i - 1] + find_b(x_0, h, i));
        beta[i] = (find_d(x_0, h, i) - find_a(x_0, h, i) * beta[i - 1]) / (find_a(x_0, h, i) * alpha[i - 1] + find_b(x_0, h, i));
    }

    double a_n = -gamma_2;
    double b_n = h * sigma_2 + gamma_2;
    double d_n = h *delta_2;
    beta[n] = (d_n - a_n * beta[n - 1]) / (a_n * alpha[n - 1] + b_n);

    double *y = create_vector(n + 1);
    y[n] = beta[n];

    for (int i = n - 1; i >= 0; i--) {
        y[i] = alpha[i] * y[i + 1] + beta[i];
        // printf("x = %lf\ty = %lf\n", x_0 + h * i, y[i]);
        printf("%lf %lf\n", x_0 + h * i, y[i]);
    }

    destroy_vector(alpha);
    destroy_vector(beta);
    destroy_vector(y);

    // double h = (x_n - x_0) / n;
    // double x = x_0 + h;
    // double *alpha = create_vector(n);
    // double *beta = create_vector(n);

    // alpha[0] = - gamma_1 / (h * sigma_1 - gamma_1);
    // beta[0] = delta_1 / (sigma_1 - gamma_1 / h);

    // for (int i = 1; i < n; i++) {
    //     double k_1 = 1 / (h * h) + p(x) / (2 * h);
    //     double k_2 = 1 / (h * h) - p(x) / (2 * h);
    //     double k_3 = -2 / (h * h) + q(x);

    //     alpha[i] = -k_1 / (k_2 * alpha[i - 1] + k_3);
    //     beta[i] = (f(x) - k_2 * beta[i - 1]) / (k_3 + k_2 * alpha[i - 1]);
    //     x += h;
    // }
    
    // double y = (gamma_2 * beta[n - 1] + h * delta_2) / (gamma_2 * (1 - alpha[n - 1]) + h * sigma_2);

    // for (int i = n - 1; i >= 0; i--) {
    //     // printf("x = %lf\ty = %lf\n", x, y);
    //     printf("%lf;%lf\n", x, y);
    //     x -= h;
    //     y = alpha[i] * y + beta[i];
    // }
    
    // // printf("x = %lf\ty = %lf\n", x, y);
    // printf("%lf;%lf\n", x, y);

    // destroy_vector(alpha);
    // destroy_vector(beta);
}

int main(int argc, char *argv[])
{
    fin_dif(p, q, f, 1, 0, -1, -1, 1, 2, 0, 1, 10);
    printf("\n\n");
    fin_dif(p, q, f, 1, 0, -1, -1, 1, 2, 0, 1, 100);
    
    return 0; 
}

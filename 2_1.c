#include <math.h>
#include <stdio.h>

double f(double x, double y)
{
    return exp(-sin(x)) - y * cos(x);
}

double
f1(double x, double u, double v)
{
    return (cos(x + 1.1 * v) + u);
}

double
f2(double x, double u, double v)
{
    return ((-1) * v * v + 2.1 * u + 1.1);
}

void RK2(double (*function)(double, double), double a, double b, double n, double x_0, double y_0)
{
    double x, y;
    double h = (b - a) / n;
    double f_prev, x_prev, y_prev = y_0;

    for (int i = 0; i < n; i++) {
        x = x_0 + h * (i + 1);
        x_prev = x_0 + h * i;
        f_prev = function(x_prev, y_prev);
        y = y_prev + (f_prev + function(x, y_prev + f_prev * h)) * h / 2;
        y_prev = y;

        printf("%lf %lf\n", x, y);
    }
}

void RK4(double (*function)(double, double), double a, double b, double n, double x_0, double y_0)
{
    double x, y;
    double h = (b - a) / n;
    double x_prev, y_prev = y_0;
    double k_1, k_2, k_3, k_4;

    for (int i = 0; i < n; i++) {
        x = x_0 + h * (i + 1);
        x_prev = x_0 + h * i;

        k_1 = h * function(x_prev, y_prev);
        k_2 = h * function(x_prev + h / 2, y_prev + k_1 / 2);
        k_3 = h * function(x_prev + h / 2, y_prev + k_2 / 2);
        k_4 = h * function(x_prev + h, y_prev + k_3);

        y = y_prev + (k_1 + 2 * k_2 + 2 * k_3 + k_4) / 6;
        y_prev = y;

        printf("%lf %lf\n", x, y);

    }
}

void RK2_s(double (*function_1)(double, double, double), double (*function_2)(double, double, double),
    double a, double b, double n, double x_0, double y1_0, double y2_0)
{
    double h = (b - a) / n;
    double x, u, v, f_u, f_v;
    double x_prev, u_prev = y1_0, v_prev = y2_0;

    for (int i = 0; i < n; i++) {
        x = x_0 + h * (i + 1);
        x_prev = x_0 + h * i;
        f_u = u_prev + h * function_1(x_prev, u_prev, v_prev);
        f_v = v_prev + h * function_2(x_prev, u_prev, v_prev);
        u = u_prev + h / 2 * (function_1(x_prev, u_prev, v_prev) + function_1(x, f_u, f_v));
        v = v_prev + h / 2 * (function_2(x_prev, u_prev, v_prev) + function_2(x, f_u, f_v));
        u_prev = u;
        v_prev = v;

        printf("%lf %lf %lf\n", x, v);
    }
}

void
RK4_s(double (*function_1)(double, double, double),
    double (*function_2)(double, double, double),
    double a, double b, double n, double x_0, double y1_0, double y2_0)
{
    double h = (b - a) / n;
    double x, u, v, f_u, f_v;
    double x_prev, u_prev = y1_0, v_prev = y2_0;
    double k_1, k_2, k_3, k_4, m_1, m_2, m_3, m_4;

    for (int i = 0; i < n; i++) {
        x = x_0 + h * (i + 1);
        x_prev = x_0 + h * i;

        k_1 = h * function_1(x_prev, u_prev, v_prev),
            m_1 = h * function_2(x_prev, u_prev, v_prev);

        k_2 = h * function_1(x_prev + h / 2, u_prev + k_1 / 2, v_prev + m_1 / 2);
        m_2 = h * function_2(x_prev + h / 2, u_prev + k_1 / 2, v_prev + m_1 / 2);

        k_3 = h * function_1(x_prev + h / 2, u_prev + k_2 / 2, v_prev + m_2 / 2);
        m_3 = h * function_2(x_prev + h / 2, u_prev + k_2 / 2, v_prev + m_2 / 2);

        k_4 = h * function_1(x_prev + h, u_prev + k_3, v_prev + m_3);
        m_4 = h * function_2(x_prev + h, u_prev + k_3, v_prev + m_3);

        u = u_prev + (k_1 + 2 * k_2 + 2 * k_3 + k_4) / 6;
        v = v_prev + (m_1 + 2 * m_2 + 2 * m_3 + m_4) / 6;

        u_prev = u;
        v_prev = v;

        printf("%lf %lf %lf\n", x, u, v);
    }
}

int main(int argc, const char *argv[])
{
    double a = 0;
    double b = 3;
    int n = 100;

    // Алгоритм Рунге-Кутта второго порядка точности для одного уравнения
    printf("\nДля одного уравнения второй порядок точности дает точки:\n");
    RK2(f, a, b, n, 0, 10);
    printf("\n\n");

    // Алгоритм Рунге-Кутта четвертого порядка точности для одного уравнения
    printf("\nДля одного уравнения четвёртый порядок точности дает точки:\n");
    RK4(f, a, b, n, 0, 0);
    printf("\n\n");


    // Алгоритм Рунге-Кутта второго порядка точности для системы из 2 уравнений
    printf("\nДля системы из двух уравнений второй порядок точности дает точки:\n");
    RK2_s(f1, f2, a, b, n, 0, 0.25, 1);
    printf("\n\n");

    // //Алгоритм Рунге-Кутта четвертого порядка точности для системы из 2 уравнений
    printf("\nДля системы из двух уравнений четвертый порядок точности дает точки:\n");
    RK4_s(f1, f2, a, b, n, 0, 0.25, 1);
    printf("\n\n");

    return 0;
}

#include "Main_Header.h"

int main()
{
    SetConsoleCP(1251);
    SetConsoleOutputCP(1251);

    const int N = 10;
    const int m = 1;
    double x_vals[N] = { 1.20, 1.82, 3.31, 7.24, 8.92, 9.12, 10.97, 14.45, 22.91, 36.31 };
    double v_vals[N] = { 630, 890, 2350, 3810, 4630, 4820, 5230, 7500, 11800, 19600 };

    // Вычисление сумм для нормальных уравнений
    vector<double> sumX(2 * m + 1, 0.0);
    for (int k = 0; k <= 2 * m; ++k)
    {
        for (int i = 0; i < N; ++i)
        {
            sumX[k] += pow(x_vals[i], k);
        }
    }

    // Формирование матрицы A и вектора b
    vector<vector<double>> A(m + 1, vector<double>(m + 1));
    vector<double> b_right(m + 1, 0.0);
    for (int i = 0; i <= m; ++i)
    {
        for (int j = 0; j <= m; ++j)
        {
            if (i == 0 && j == 0) A[i][j] = N;
            else A[i][j] = sumX[i + j];
        }

        for (int t = 0; t < N; ++t)
        {
            b_right[i] += v_vals[t] * pow(x_vals[t], i);
        }
    }

    // Решаем систему нормальных уравнений
    vector<double> coeffs = GaussMethod(A, b_right);
    double a = coeffs[0];
    double b = coeffs[1];  // коэффициент b

    // Вычисление остаточной дисперсии и sigma
    double s2 = 0.0;
    for (int i = 0; i < N; ++i)
    {
        double v_pred = a + b * x_vals[i];
        s2 += pow(v_vals[i] - v_pred, 2);
    }
    s2 /= (N - m - 1);
    double sigma = sqrt(s2);

    // Вывод результатов
    cout << fixed << setprecision(3);
    cout << "Линейная аппроксимация: v = a + b*x\n";
    cout << "a = " << a << " км/с\n";
    cout << "b = " << b << " км/(с*Мпк)\n";
    cout << "Остаточная дисперсия S^2: " << s2 << "\n";
    cout << "Среднеквадратичное отклонение: " << sigma << "\n\n";

    cout << "Таблица результатов:\n";
    cout << "-----------------------------------------------------------------\n";
    cout << "|   x (Мпк) |   v (км/с) |  Предсказанное v (км/с) | Отклонение |\n";
    cout << "-----------------------------------------------------------------\n";
    for (int i = 0; i < N; ++i) {
        double v_pred = a + b * x_vals[i];
        double delta = v_vals[i] - v_pred;
        cout << "| " << setw(9) << x_vals[i] << " | " << setw(10) << v_vals[i] << " | " << setw(23) << v_pred << " | " << setw(10) << delta << " |\n";
    }
    cout << "-----------------------------------------------------------------\n";

    return 0;
}
#include "Main_Header.h"

// Решение СЛАУ методом Гаусса для произвольной размерности
vector<double> GaussMethod(const vector<vector<double>>& A, const vector<double>& b)
{
    int n = A.size();
    vector<vector<double>> M(n, vector<double>(n + 1));
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j) M[i][j] = A[i][j];
        M[i][n] = b[i];
    }
    // Прямой ход
    for (int k = 0; k < n; ++k)
    {
        // выбор опорного элемента
        int piv = k;
        for (int i = k + 1; i < n; ++i)
        {
            if (fabs(M[i][k]) > fabs(M[piv][k]))
                piv = i;
        }
        swap(M[k], M[piv]);
        // вычитание строк
        for (int i = k + 1; i < n; ++i)
        {
            double factor = M[i][k] / M[k][k];
            for (int j = k; j <= n; ++j)
            {
                M[i][j] -= factor * M[k][j];
            }
        }
    }
    // Обратный ход
    vector<double> x(n);
    for (int i = n - 1; i >= 0; --i)
    {
        double sum = M[i][n];
        for (int j = i + 1; j < n; ++j)
        {
            sum -= M[i][j] * x[j];
        }
        x[i] = sum / M[i][i];
    }
    return x;
}
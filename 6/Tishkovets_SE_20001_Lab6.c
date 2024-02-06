#define _USE_MATH_DEFINES
#pragma warning(disable:4996)
#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <locale.h>
#include <math.h>
#include <stdlib.h>
#define MEMORY_ERROR 2 //код ошибки выделения памяти

//динамически выделяет память для матрицы размера n x m
//n - длина столбца (кол-во строк); m - длина строки (кол-во столбцов)
double** createMatrix(int n, int m) {
    double** matrix = (double**)malloc(n * sizeof(double*));
    if (matrix == NULL) { //проверка на ошибку с выделением памяти
        return MEMORY_ERROR;
    }
    for (int i = 0; i < n; i++) {
        matrix[i] = (double*)malloc(m * sizeof(double));
    }
    return matrix;
}

//считывает из файла матрицу в двумерный массив
double** readMatrix(FILE* file, double** A, int n, int m) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            fscanf(file, "%lf", &A[i][j]);
        }
    }
    return A;
}

//очищает память, выделенную для матрицы 
void freeMatrix(double** matrix, int n) {
    for (int i = 0; i < n; i++) {
        free(matrix[i]);
    }
    free(matrix);
}

//печатает матрицу
void print_matrix(double** matrix, int n, int m) {
    printf("\nМатрица:\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            if (j == n) printf("| "); //для отделения вектора свободных чисел
            printf("%15.10lf ", matrix[i][j]); //ширина поля 5, 1 цифра после запятой, выравнивание по правому краю
        }
        printf("\n");
    }
    printf("\n");
}

//умножает матрицу на вектор столбец
double* multiplyMatrixByColumn(double** matrix, double* column, int n) {
    double* res = createMatrix(n, 1);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            res[i] += matrix[i][j] * column[j];
        }
    }
    return res;
}

//умножает один вектор столбец (транспонированный) на другой
double multiplyColumnByColumn(double* x1, double* x2, int n) {
    double res = 0.0;
    for (int i = 0; i < n; i++) {
        res += x1[i] * x2[i];
    }
    return res;
}

//считает евклидову норму вектора
double Norm(double* x, int n) {
    double norm = 0.0;
    for (int i = 0; i < n; i++) { //находит сумму квадратов
        norm += x[i] * x[i];
    }
    norm = sqrt(norm); //извлекает корень
    return norm;
}

//реализует степенной метод с оптимальным сдвигом
double power_method_shift(double** A, double* l_real, double eps, int* iter, int n) {
    double* x_prev = (double*)malloc(n * sizeof(double*)); //начальное приближение
    double* x = (double*)malloc(n * sizeof(double*)); //следующее приближение
    for (int i = 0; i < n; i++) {
        x_prev[i] = 1.0; //идентификация начального приближения
    }

    double shift = (l_real[1] + l_real[n - 1]) / 2; //оптимальный сдвиг

    double** B = createMatrix(n, n); //новая матрица, уже со сдвигом
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i != j) B[i][j] = A[i][j];
            else B[i][j] = A[i][j] - shift; // B = A - shift * E
        }
    }

    double lambda = 0; //лямбда
    double lambda_prev = 10; //лямбда на предыдущей итерации
    int Iter = 0; //число итераций

    while (fabs(lambda - lambda_prev) > eps) {
        x = multiplyMatrixByColumn(B, x_prev, n); //следующее приближение

        double temp1 = multiplyColumnByColumn(x, x_prev, n); //числитель
       // double temp2 = multiplyColumnByColumn(x_prev, x_prev, n); //знаменатель

        lambda_prev = lambda;
        lambda = temp1; //обновляет лямбду lambda_k+1 = (x_k+1, x_k) / (x_k, x_k)
            // так как векторы нормируются, temp2 = 1

        double norm = Norm(x, n); //находит евклидову норму вектора
        for (int i = 0; i < n; i++) {
            x[i] /= norm; //нормирует вектор
        }

        for (int i = 0; i < n; i++) { //обновляет приближение
            x_prev[i] = x[i];
        }

        Iter++;
    }

    *iter = Iter;

    return lambda + shift;
}

int main(void) {
    setlocale(LC_CTYPE, "Russian");
    int n = 10;
    int Iter = 0;
    double eps = 0;

    /* ИССЛЕДОВАНИЕ ЗАВИСИМОСТИ АБСОЛЮТНОЙ ПОГРЕШНОСТИ И ЧИСЛА ИТЕРАЦИЙ ОТ ТОЧНОСТИ */

    FILE* fp; //файл с 3 матрицами
    fopen_s(&fp, "file1.txt", "r"); //открывает файл из MatLab
    if (!fp) { //проверка успешного открытия файла
        printf("Ошибка при открытии файла\n");
        exit(-1); //выходит из программы с кодом -1
    }

    FILE* f1; //файл для записи
    fopen_s(&f1, "norm6.csv", "w"); //открывает файл для записи норм
    if (!f1) { //проверка успешного открытия файла
        printf("Ошибка при открытии файла\n");
        exit(-1); //выходит из программы с кодом -1
    }

    double lambda_real[10] = { 100, 10, 9, 8, 7, 6, 5, 4, 3, 2 }; //столбец настоящих сч

    for (int s = 0; s < 3; s++) {
        double** A = createMatrix(n, n);
        A = readMatrix(fp, A, n, n);
        print_matrix(A, n, n);

        for (int d = 1; d <= 10; d++) {
            eps = pow(10, -d);

            //вызов степенного метода со сдвигом
            double lambda = power_method_shift(A, lambda_real, eps, &Iter, n);

            double error = fabs(lambda - lambda_real[0]);

            printf("Максимальное сч: %.15lf\n", lambda);

            /* заполняет таблицу Excel - ';' для перехода к следующей ячейке */
            fprintf(f1, "%.10lf;%d;%.15lf\n", eps, Iter, error);
        }
        fprintf(f1, "\n");
    }

    fclose(f1);
    fclose(fp);

    return 0;
}
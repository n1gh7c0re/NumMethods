#define _USE_MATH_DEFINES
#pragma warning(disable:4996)
#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <locale.h>
#include <math.h>
#include <stdlib.h>
#include <windows.h>
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

//очищает память, выделенную для матрицы 
void freeMatrix(double** matrix, int n) {
    for (int i = 0; i < n; i++) {
        free(matrix[i]);
    }
    free(matrix);
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

//печатает матрицу
void print_matrix(double** matrix, int n, int m) {
    printf("\nМатрица:\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            if (j == n) printf("| "); //для отделения вектора свободных чисел
            printf("%15.8lf ", matrix[i][j]); //ширина поля 5, 1 цифра после запятой, выравнивание по правому краю
        }
        printf("\n");
    }
    printf("\n");
}

//перемножает две квадратные матрицы одной размерности
double** multiplyMatrixByMatrix(double** A, double** B, int n) {
    double** C = createMatrix(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            C[i][j] = 0;
            for (int k = 0; k < n; k++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return C;
}

//реализует метод Леверье для отыскания собственных чисел
double leverrier(double* p, double** A, int n) {
    LARGE_INTEGER freq, start, end;
    QueryPerformanceFrequency(&freq);
    QueryPerformanceCounter(&start); //запоминает время начала выполнения метода

    double* S = malloc(n * sizeof(double)); //массив следов матрицы
    for (int i = 0; i < n; i++) {
        S[i] = 0.0; //индексация
        p[i] = 0.0; //индексация
    }

    double** B = createMatrix(n, n); //вспомогательная матрица для нахождения степеней матрицы А
    for (int d = 0; d < n; d++) { //копирует элементы матрицы А в матрицу В
        for (int l = 0; l < n; l++) {
            B[d][l] = A[d][l];
        }
    }

    for (int k = 0; k < n; k++) { //находит след матрицы А степени k = 1, 2, ..., n
        for (int j = 0; j < n; j++) {
            S[k] += B[j][j];
        }
        B = multiplyMatrixByMatrix(B, A, n);
    }

    //получает коэффициенты полинома p(k) = 1/k * [ S(k) - p(1) * S(k-1) - ... - p(k-1) * S(1) ]
    for (int i = 0; i < n; i++) {
        p[i] = S[i] / (i + 1);
        for (int j = 0; j < i; j++) {
            p[i] -= p[j] * S[i - j - 1] / (i + 1);
        }
    }

    QueryPerformanceCounter(&end); //запоминает время конца выполнения метода
    //находит время в секундах на исполнение метода
    double t = (double)(end.QuadPart - start.QuadPart) / (double)freq.QuadPart;
    return t;
}

int main(void) {
    setlocale(LC_CTYPE, "Russian");
    int n = 10; //размерность СЛАУ
    int a = 1;
    int c = 10; //границы диапазона значений элементов матрицы

    /* ИССЛЕДОВАНЕ ЗАВИСИМОСТИ ВРЕМЕНИ ВЫПОЛНЕНИЯ МЕТОДА ОТ РАЗМЕРНОСТИ МАТРИЦЫ */

    //FILE* f2; //файл для записи
    //fopen_s(&f2, "time5.csv", "w"); //открывает файл для записи времени
    //if (!f2) { //проверка успешного открытия файла
    //    printf("Ошибка при открытии файла\n");
    //    exit(-1); //выходит из программы с кодом -1
    //}

    //for (int d = 0; d < 6; d++) {
    //    double** A = createMatrix(n, n);
    //    for (int i = 0; i < n; i++) { //заполняет матрицу случайными числами
    //        for (int j = 0; j < n; j++) {
    //            A[i][j] = (double)rand() / RAND_MAX * (c - a) + a;
    //        }
    //    }
    //    double* p = malloc(n * sizeof(double)); //массив промежуточных коэффициентов полинома

    //    double t = leverrier(p, A, n); //вызывает метод Леверье

    //    /* заполняет таблицу Excel - ';' для перехода к следующей ячейке */
    //    fprintf(f2, "%lf;%d\n", t, n);
    //    printf("%lf\n", t);

    //    freeMatrix(A, n);
    //    free(p);
    //    n *= 2;
    //}
    //fclose(f2);
    
    /* ИССЛЕДОВАНИЕ СИММЕТРИЧНЫХ МАТРИЦ */

    //FILE* f1; //файл для записи
    //fopen_s(&f1, "polynomial1.txt", "w"); //открывает файл для записи коэффициентов симметричной матрицы
    //if (!f1) { //проверка успешного открытия файла
    //    printf("Ошибка при открытии файла\n");
    //    exit(-1); //выходит из программы с кодом -1
    //}

    //FILE* fp; //файл с матрицами
    //fopen_s(&fp, "file1.txt", "r"); //открывает файл из MatLab
    //if (!fp) { //проверка успешного открытия файла
    //    printf("Ошибка при открытии файла\n");
    //    exit(-1); //выходит из программы с кодом -1
    //}

    //for (int d = 0; d < 3; d++) {
    //    double** A = createMatrix(n, n);
    //    A = readMatrix(fp, A, n, n);
    //    print_matrix(A, n, n);

    //    double* p = malloc(n * sizeof(double)); //массив промежуточных коэффициентов полинома

    //    double t = leverrier(p, A, n); //вызывает метод Леверье

    //    double* p_true = malloc((n + 1) * sizeof(double)); //массив окончательных коэффициентов полинома
    //    p_true[0] = 1; //согласно формуле первый коээффициент равен 1


    //    for (int i = 0; i < n; i++) {
    //        p_true[i + 1] = -p[i]; //меняет знак
    //    }

    //    for (int i = 0; i < n + 1; i++) { //записывает в текстовый файл коэффициенты
    //        fprintf(f1, "%.15lf ", p_true[i]);
    //    }
    //    fprintf(f1, "\n");

    //    freeMatrix(A, n);
    //    free(p);
    //    free(p_true);
    //}
    //fclose(f1);
    //fclose(fp);

    /* ИССЛЕДОВАНИЕ НЕ СИММЕТРИЧНЫХ МАТРИЦ */

    FILE* f3; //файл для записи
    fopen_s(&f3, "polynomial2.txt", "w"); //открывает файл для записи коэффициентов несимметричной матрицы
    if (!f3) { //проверка успешного открытия файла
        printf("Ошибка при открытии файла\n");
        exit(-1); //выходит из программы с кодом -1
    }

    FILE* fp2; //файл с матрицами
    fopen_s(&fp2, "file2.txt", "r"); //открывает файл из MatLab
    if (!fp2) { //проверка успешного открытия файла
        printf("Ошибка при открытии файла\n");
        exit(-1); //выходит из программы с кодом -1
    }

    for (int d = 0; d < 3; d++) {
        double** A = createMatrix(n, n);
        A = readMatrix(fp2, A, n, n);
        print_matrix(A, n, n);

        double* p = malloc(n * sizeof(double)); //массив промежуточных коэффициентов полинома

        double t = leverrier(p, A, n); //вызывает метод Леверье

        double* p_true = malloc((n + 1) * sizeof(double)); //массив окончательных коэффициентов полинома
        p_true[0] = 1; //согласно формуле первый коээффициент равен 1


        for (int i = 0; i < n; i++) {
            p_true[i + 1] = -p[i]; //меняет знак
        }

        for (int i = 0; i < n + 1; i++) { //записывает в текстовый файл коэффициенты
            fprintf(f3, "%.15lf ", p_true[i]);
        }
        fprintf(f3, "\n");

        freeMatrix(A, n);
        free(p);
        free(p_true);
    }
    fclose(f3);
    fclose(fp2);

    return 0;
}
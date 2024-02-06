#define _USE_MATH_DEFINES
#pragma warning(disable:4996)
#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <locale.h>
#include <math.h>
#include <stdlib.h>
#include <windows.h>
#include <time.h>
#define MEMORY_ERROR 2 //код ошибки выделения памяти
#define DET_ERROR -2   //код нарушения условия применимости метода

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

//считает евклидову норму вектора
double Norm(double* x, int n) {
    double norm = 0.0;
    for (int i = 0; i < n; i++) { //находит сумму квадратов
        norm += x[i] * x[i];
    }
    norm = sqrt(norm); //извлекает корень
    return norm;
}

//считает евклидову норму разности двух векторов
double NormDiff(double* x1, double* x2, int n) {
    double norm = 0.0;
    for (int i = 0; i < n; i++) { //находит сумму квадратов
        norm += (x1[i] - x2[i]) * (x1[i] - x2[i]);
    }
    norm = sqrt(norm); //извлекает корень
    return norm;
}

//транспонирует матрицу А размерности n
double** transpose_matrix(double** A, int n) {
    double** B = createMatrix(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            B[j][i] = A[i][j];
        }
    }
    return B;
}

/* получает обратную матрицу с помощью метода Гаусса - например:
   1 2 3 | 1 0 0     1 0 0 | 2 8 9
   4 5 6 | 0 1 0 --> 0 1 0 | 3 6 5
   7 8 9 | 0 0 1     0 0 1 | 1 7 4 */
double** inverse_matrix(double** A, int n) {
    double r, a;
    double** matrix = createMatrix(n, 2 * n);
    for (int i = 0; i < n; i++) { //копирует элементы исходной матрицы
        for (int j = 0; j < n; j++) {
            matrix[i][j] = A[i][j];
        }
    }

    for (int i = 0; i < n; i++) {
        //получает единичную матрицу справа от исходной
        for (int j = n; j < 2 * n; j++) {
            if (i == (j - n)) matrix[i][j] = 1;
            else matrix[i][j] = 0;
        }
    }

    //прямой ход метода Гаусса
    for (int i = 0; i < n; i++) { //идет по столбцам
        for (int j = 0; j < n; j++) { //по строкам
            if (i != j) { //получает нули над и под диагональю, т.к. j начинается с 0
                r = matrix[j][i] / matrix[i][i];
                for (int k = 0; k < 2 * n; k++) {
                    matrix[j][k] -= r * matrix[i][k];
                }
            }
        }
    }

    //получает единичную матрицу слева
    for (int i = 0; i < n; i++) {
        a = matrix[i][i];
        for (int j = 0; j < 2 * n; j++) {
            matrix[i][j] /= a;
        }
    }

    double** A_1 = createMatrix(n, n); //обратная матрица
    for (int i = 0; i < n; i++) { //получает обратную квадратную матрицу размерности n
        for (int j = 0; j < n; j++) {
            A_1[i][j] = matrix[i][n + j];
        }
    }
    freeMatrix(matrix, n);

    return A_1;
}

//реализует модифицированный метод Грама-Шмидта
double gramSchmidt(double** A, double** Q, double** R, int n) {
    LARGE_INTEGER freq, start, end;
    QueryPerformanceFrequency(&freq);
    QueryPerformanceCounter(&start); //запоминает время начала выполнения метода

    for (int j = 0; j < n; j++) {
        for (int s = 0; s < n; s++) { //копирует j-й столбец матрицы A в j-й столбец матрицы Q
            Q[s][j] = A[s][j];
        }

        double* Q_j = (double*)malloc(n * sizeof(double*));
        for (int k = 0; k < n; k++) { //копирует j-й столбец матрицы Q
            Q_j[k] = Q[k][j];
        }

        for (int i = 0; i < j; i++) {
            double* Q_i = (double*)malloc(n * sizeof(double*));
            for (int z = 0; z < n; z++) { //копирует i-й столбец матрицы Q
                Q_i[z] = Q[z][i];
            }

            R[i][j] = multiplyColumnByColumn(Q_i, Q_j, n);
            for (int t = 0; t < n; t++) {
                Q_j[t] = Q_j[t] - R[i][j] * Q_i[t];
            }
            free(Q_i);
        }

        R[j][j] = Norm(Q_j, n);
        if (R[j][j] == 0) return DET_ERROR; //столбцы матрицы А не являются линейно независимыми

        for (int l = 0; l < n; l++) {
            Q_j[l] = Q_j[l] / R[j][j];
        }

        for (int h = 0; h < n; h++) { //заполняем j-й уже новый столбец матрицы Q
            Q[h][j] = Q_j[h];
        }
        free(Q_j);
    }

    QueryPerformanceCounter(&end); //запоминает время конца выполнения метода
    //находит время в секундах на исполнение метода
    double t = (double)(end.QuadPart - start.QuadPart) / (double)freq.QuadPart;
    return t;
}

int main(void) {
    setlocale(LC_CTYPE, "Russian");
    int n = 10;
    int m = n + 1;
    int a = -100, b = 100; //границы диапазона значений элементов матрицы

    /* ИССЛЕДОВАНЕ ЗАВИСИМОСТИ ВРЕМЕНИ ВЫПОЛНЕНИЯ МЕТОДА ОТ РАЗМЕРНОСТИ МАТРИЦЫ */

    //FILE* f2; //файл для записи
    //fopen_s(&f2, "time3.csv", "w"); //открывает файл для записи норм
    //if (!f2) { //проверка успешного открытия файла
    //    printf("Ошибка при открытии файла\n");
    //    exit(-1); //выходит из программы с кодом -1
    //}

    //for (int d = 0; d < 8; d++) {
    //    double** matrix = createMatrix(n, m);
    //    for (int i = 0; i < n; i++) { //заполняет матрицу случайными числами
    //        for (int j = 0; j < m; j++) {
    //            matrix[i][j] = (double)rand() / RAND_MAX * (b - a) + a;
    //        }
    //    }

    //    double* b_real = createMatrix(n, 1); //столбец свободных членов исходной матрицы
    //    for (int z = 0; z < n; z++) { //копирует все элементы
    //        b_real[z] = matrix[z][n];
    //    }

    //    double** Q = createMatrix(n, n); //ортогональная матрица
    //    double** R = createMatrix(n, n); //матрица перехода от ортонормированного базиса к исходному
    //    for (int i = 0; i < n; i++) { //заполняет матрицы R и Q нулями
    //        for (int j = 0; j < n; j++) {
    //            R[i][j] = 0.0;
    //            Q[i][j] = 0.0;
    //        }
    //    }

    //    double t = gramSchmidt(matrix, Q, R, n);
    //    if (t != DET_ERROR) {
    //        double** QT = transpose_matrix(Q, n); //транспонированная матрица 
    //        double** R_1 = inverse_matrix(R, n); //обратная матрица

    //        double* x = (double*)malloc(n * sizeof(double*)); //вектор решений
    //        x = multiplyMatrixByColumn(multiplyMatrixByMatrix(R_1, QT, n), b_real, n);

    //        /* заполняет таблицу Excel - ';' для перехода к следующей ячейке */
    //        fprintf(f2, "%lf;%d\n", t, n);
    //    }
    //    else printf("Определитель матрицы равен 0"); //проверка условий применимости метода

    //    n *= 2;
    //    m = n + 1;
    //}
    //fclose(f2); //закрывает файл

    /* ИССЛЕДОВАНИЕ ЗАВИСИМОСТИ АБСОЛЮТНОЙ ПОГРЕШНОСТИ И НЕВЯЗКИ ОТ ЧИСЛА ОБУСЛОВЛЕННОСТИ */

    FILE* fp; //файл с матрицами
    fopen_s(&fp, "file.txt", "r"); //открывает файл из MatLab
    if (!fp) { //проверка успешного открытия файла
        printf("Ошибка при открытии файла\n");
        exit(-1); //выходит из программы с кодом -1
    }

    FILE* f1; //файл для записи
    fopen_s(&f1, "norm3.csv", "w"); //открывает файл для записи норм
    if (!f1) { //проверка успешного открытия файла
        printf("Ошибка при открытии файла\n");
        exit(-1); //выходит из программы с кодом -1
    }

    double x_real[10] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 }; //столбец настоящих решений

    for (int d = 0; d < 8; d++) {
        double** A = createMatrix(n, m);
        A = readMatrix(fp, A, n, m);

        double* b_real = createMatrix(n, 1); //столбец свободных членов исходной матрицы
        for (int z = 0; z < n; z++) { //копирует все элементы
            b_real[z] = A[z][n];
        }

        double cond = pow(10, d); //число обусловленности
        printf("Исходная:\n");
        print_matrix(A, n, m);

        double** Q = createMatrix(n, n); //ортогональная матрица
        double** R = createMatrix(n, n); //верхняя треугольная матрица
        for (int i = 0; i < n; i++) { //заполняет матрицы R и Q нулями
            for (int j = 0; j < n; j++) {
                R[i][j] = 0.0;
                Q[i][j] = 0.0;
            }
        }

        double t = gramSchmidt(A, Q, R, n);
        if (t != DET_ERROR) {
            double** QT = transpose_matrix(Q, n); //транспонированная матрица 
            double** R_1 = inverse_matrix(R, n); //обратная матрица

            double* x = (double*)malloc(n * sizeof(double*)); //вектор решений
            x = multiplyMatrixByColumn(multiplyMatrixByMatrix(R_1, QT, n), b_real, n);

            double norm_xx = NormDiff(x, x_real, n); //вычисляет норму абсолютной погрешности решения
            double* b = multiplyMatrixByColumn(A, x, n); //получает новый столбец свободных чисел
            double norm_AXB = NormDiff(b, b_real, n); //вычисляет норму невязки

            //печатает решение
            for (int s = 0; s < n; s++) {
                printf("%.20lf\n", x[s]);
            }

            /* заполняет таблицу Excel - ';' для перехода к следующей ячейке */
            fprintf(f1, "%.20lf;%.20lf;%lf\n", norm_xx, norm_AXB, cond);

            //чистит память
            freeMatrix(QT, n);
            freeMatrix(R_1, n);
            free(x);
            free(b);
        }
        else printf("Определитель матрицы равен 0"); //проверка условий применимости метода

        //чистит память
        freeMatrix(A, n);
        freeMatrix(Q, n);
        freeMatrix(R, n);
        free(b_real);
    }
    fclose(fp); //закрывает файл
    fclose(f1); //закрывает файл

    return 0;
}
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

//считает евклидову норму разности двух векторов
double NormDiff(double* x1, double* x2, int n) {
    double norm = 0.0;
    for (int i = 0; i < n; i++) { //находит сумму квадратов
        norm += (x1[i] - x2[i]) * (x1[i] - x2[i]);
    }
    norm = sqrt(norm); //извлекает корень
    return norm;
}

double conjugate_gradient(double** A, double* b, double* x, int n, double eps, int* iter) {
    LARGE_INTEGER freq, start, end;
    QueryPerformanceFrequency(&freq);
    QueryPerformanceCounter(&start); //запоминает время начала выполнения метода

    double* r = malloc(n * sizeof(double));
    double* p = malloc(n * sizeof(double));
    double* Ap = malloc(n * sizeof(double));

    //Инициализация начального приближения
    for (int i = 0; i < n; i++) {
        x[i] = 5.0;
    }

    for (int i = 0; i < n; i++) {
        r[i] = b[i];
        for (int j = 0; j < n; j++) {
            r[i] -= A[i][j] * x[j]; //вектор невязки: r = b - A * x
        }
        p[i] = r[i]; //Инициализация начального направления
    }

    int Iter = 0;
    double r_new = 1000.0; //квадрат нормы нового вектора невязки
    //Итерационный процесс
    while (sqrt(r_new) > eps) {
        Iter++;
        double r_old = 0.0; //квадрат нормы старого вектора невязки
        for (int i = 0; i < n; i++) {
            r_old += r[i] * r[i]; //знаменатель коэффициента бетта и числитель коэффициента альфа
        }

        //Умножение матрицы A на вектор p: Ap = A * p
        for (int i = 0; i < n; i++) {
            Ap[i] = 0.0;
            for (int j = 0; j < n; j++) {
                Ap[i] += A[i][j] * p[j];
            }
        }

        double alpha = 0.0;
        double pAp = 0.0;
        for (int i = 0; i < n; i++) {
            pAp += p[i] * Ap[i]; //знаменатель коэффициента альфа
        }
        alpha = r_old / pAp; //коэффициент для получения следующего приближения x

        for (int i = 0; i < n; i++) {
            x[i] += alpha * p[i]; //получает следующее приближение: x = x + alfa * p
            r[i] -= alpha * Ap[i]; //получает новый вектор невязки: r = r - alfa * Ap
        }

        r_new = 0.0; //квадрат нормы нового вектора невязки
        for (int i = 0; i < n; i++) {
            r_new += r[i] * r[i]; //числитель коэффициента бетта
        }

        double beta = r_new / r_old; //коэффициент для получения сопряженного направления p
        for (int i = 0; i < n; i++) {
            p[i] = r[i] + beta * p[i]; //получает сопряженное направление: p = r + beta * p
        }
    }
    *iter = Iter;

    free(r);
    free(p);
    free(Ap);

    QueryPerformanceCounter(&end); //запоминает время конца выполнения метода
    //находит время в секундах на исполнение метода
    double t = (double)(end.QuadPart - start.QuadPart) / (double)freq.QuadPart;
    return t;
}

int main() {
    setlocale(LC_CTYPE, "Russian");
    int n = 10; //размерность СЛАУ
    int a = 1;
    int c = 100; //границы диапазона значений элементов матрицы
    double eps = 0.0; //точность сходимости
    int Iter = 0; //количество итераций

    /* ИССЛЕДОВАНЕ ЗАВИСИМОСТИ ВРЕМЕНИ ВЫПОЛНЕНИЯ МЕТОДА ОТ РАЗМЕРНОСТИ МАТРИЦЫ */

    //FILE* f1; //файл для записи
    //fopen_s(&f1, "time4.csv", "w"); //открывает файл для записи времени
    //if (!f1) { //проверка успешного открытия файла
    //    printf("Ошибка при открытии файла\n");
    //    exit(-1); //выходит из программы с кодом -1
    //}

    //for (int d = 0; d < 8; d++) {
    //    double** A = createMatrix(n, n + 1);
    //    for (int i = 0; i < n; i++) { //заполняет матрицу случайными числами
    //        A[i][i] = (double)rand() / RAND_MAX * (c - a) + a;
    //        for (int j = i + 1; j < n; j++) { //получает симметричную положительно-определенную матрицу
    //            A[i][j] = (double)rand() / RAND_MAX * (c - a) + a;
    //            A[j][i] = A[i][j];
    //        }
    //        A[i][n] = (double)rand() / RAND_MAX * (c - a) + a;
    //    }
    //    double* b = malloc(n * sizeof(double)); //столбец свободных членов исходной матрицы
    //    for (int z = 0; z < n; z++) { //копирует все элементы столбца свободных членов
    //        b[z] = A[z][n];
    //    }

    //    double* x = malloc(n * sizeof(double)); //столбец решений
    //    eps = 1e-10;

    //    //Вызов метода сопряженных градиентов
    //    double t = conjugate_gradient(A, b, x, n, eps, &Iter);
    //    /* заполняет таблицу Excel - ';' для перехода к следующей ячейке */
    //    fprintf(f1, "%lf;%d\n", t, n);
    //    printf("%lf\n", t);

    //    //Освобождает память
    //    freeMatrix(A, n);
    //    free(b);
    //    free(x);
    //    n *= 2;
    //}

    //fclose(f1);

    /* ИССЛЕДОВАНИЕ ЗАВИСИМОСТИ АБСОЛЮТНОЙ ПОГРЕШНОСТИ И ЧИСЛА ИТЕРАЦИЙ ОТ ТОЧНОСТИ */
    
    FILE* fp; //файл с матрицей
    fopen_s(&fp, "file.txt", "r"); //открывает файл из MatLab
    if (!fp) { //проверка успешного открытия файла
        printf("Ошибка при открытии файла\n");
        exit(-1); //выходит из программы с кодом -1
    }

    FILE* f2; //файл для записи
    fopen_s(&f2, "norm4.csv", "w"); //открывает файл для записи норм
    if (!f2) { //проверка успешного открытия файла
        printf("Ошибка при открытии файла\n");
        exit(-1); //выходит из программы с кодом -1
    }

    n = 10;
    double x_real[10] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 }; //столбец настоящих решений
    double** A = createMatrix(n, n + 1);
    A = readMatrix(fp, A, n, n + 1); //инициализирует матрицу

    double* b = malloc(n * sizeof(double)); //столбец свободных членов исходной матрицы
    for (int z = 0; z < n; z++) { //копирует все элементы столбца свободных членов
        b[z] = A[z][n];
    }

    for (int k = 1; k <= 10; k++) {
        double* x = malloc(n * sizeof(double)); //столбец решений
        eps = pow(10, -k);
        
        //Вызов метода сопряженных градиентов
        double t = conjugate_gradient(A, b, x, n, eps, &Iter);

        double norm_xx = NormDiff(x, x_real, n); //вычисляет норму абсолютной погрешности решения
        
        /* заполняет таблицу Excel - ';' для перехода к следующей ячейке */
        fprintf(f2, "%d;%.20lf;%.*lf\n", Iter, norm_xx, k, eps);
        printf("%d\n", Iter);

        //Освобождает память
        free(x);
    }

    freeMatrix(A, n);
    free(b);
    fclose(fp);
    fclose(f2);

    return 0;
}
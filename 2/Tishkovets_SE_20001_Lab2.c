#define _USE_MATH_DEFINES
#pragma warning(disable:4996)
#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <locale.h>
#include <math.h>
#include <stdlib.h>
#include <windows.h>
#include <time.h>
#define DET_ERROR -2

//динамически выделяет память для матрицы размера n x m
//n - длина столбца (кол-во строк); m - длина строки (кол-во столбцов)
double** createMatrix(int n, int m) {
    double** matrix = (double**)malloc(n * sizeof(double*));
    if (matrix == NULL) { //проверка на ошибку с выделением памяти
        printf("Ошибка веделения памяти в функции createMatrix");
        exit(2); //выходит из программы с кодом 2
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

//реализует метод Гаусса с поиском ведущего элемента по строке
double gaussMethod(int n, int m, double** a, double* x, int* number) {
    LARGE_INTEGER freq, start, end;
    QueryPerformanceFrequency(&freq);
    QueryPerformanceCounter(&start); //запоминает время начала выполнения метода
    
    int i, j, k;
    int s; //переменная для запоминания номера столбца
    
    for (i = 0; i < n - 1; i++) { //идёт по строкам
        double max = 0.0;
        s = i;
        //поиск ведущего элемента по строке
        for (int z = i; z < n; z++) { //сначала пробегает по всей строке
            if (fabs(a[i][z]) > fabs(max)) {
                max = a[i][z]; //находит максимальный элемент
                s = z; //запоминает номер столбца
            }
        }
        if (s != i) { // если ведущий элемент не на диагонали
            for (j = 0; j < n; j++) { //меняет столбцы местами
                double temp = a[j][i]; //временная переменная
                a[j][i] = a[j][s];
                a[j][s] = temp;
            }
            double temp1 = number[i]; //запоминает изменение порядка следования переменных
            number[i] = number[s];
            number[s] = temp1;
        }
        //прямой ход
        for (k = i + 1; k < n; k++) {
            if (a[i][i] == 0) { //если диагональный элемент равен 0
                return DET_ERROR;
            }
            double  term = a[k][i] / a[i][i]; //отношение элемента под диагональю к диагональному
            for (j = 0; j < m; j++) { //получает нули в столбце под диагональным элементом
                a[k][j] = a[k][j] - term * a[i][j];
            }
        }
    }
    //обратный ход
    for (i = n - 1; i >= 0; i--) {
        x[i] = a[i][m - 1]; //нижний правый элемент
        for (j = i + 1; j < m - 1; j++) { //идет снизу вверх
            x[i] = x[i] - a[i][j] * x[j]; //получает остальные решения, используя полученный ранее значения
        }
        //чтобы не делить на 0
        if (a[i][i] != 0) x[i] = x[i] / a[i][i]; //коэффициент перед иксом должен быть 1
    }
    QueryPerformanceCounter(&end); //запоминает время конца выполнения метода
    //находит время в секундах на исполнение метода
    double t = (double)(end.QuadPart - start.QuadPart) / (double)freq.QuadPart;
    return t;
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

//считает евклидову норму разности двух векторов
double Norm(double* x1, double* x2, int n) {
    double norm = 0.0;
    for (int i = 0; i < n; i++) { //находит сумму квадратов
        norm += (x1[i] - x2[i]) * (x1[i] - x2[i]);
    }
    norm = sqrt(norm); //извлекает корень
    return norm;
}

int main(void) {
    setlocale(LC_CTYPE, "Russian");
    int a = -100, b = 100; //границы диапазона значений элементов матрицы
    int n = 10; //длина столбца (кол-во строк)
    int m = n + 1; //длина строки (кол-во столбцов)
    int* number = (int*)malloc(n * sizeof(int*)); //массив для запоминания порядка переменных
    for (int i = 0; i < n; i++) {
        number[i] = i;
    }

    /* ИССЛЕДОВАНЕ ЗАВИСИМОСТИ ВРЕМЕНИ ИСПОЛНЕНИЯ МЕТОДА ОТ РАЗМЕРНОСТИ МАТРИЦЫ */

    //FILE* f; //файл для записи
    //fopen_s(&f, "time.csv", "w"); //открывает файл для записи времени
    //if (!f) { //проверка успешного открытия файла
    //    printf("Ошибка при открытии файла\n");
    //    exit(-1); //выходит из программы с кодом -1
    //}
    //
    //for (int z = 1; z <= 10; z++) {
    //    //srand(time(NULL)); //чтобы каждый раз генерировались другие числа
    //    double** matrix = createMatrix(n, m);
    //    for (int i = 0; i < n; i++) { //заполняет матрицу случайными числами
    //        for (int j = 0; j < m; j++) {
    //            matrix[i][j] = (double)rand() / RAND_MAX * (b - a) + a;
    //        }
    //    }
    //    printf("До преобразований:\n");
    //    print_matrix(matrix, n, m);
    //    double* x = (double*)malloc(n * sizeof(double*)); //вектор решений
    //    double t = gaussMethod(n, m, matrix, x, number);
    //    if (t != DET_ERROR) { //обрабатывает ошибку на случай вырожденности матрицы
    //        printf("После преобразований:\n");
    //        print_matrix(matrix, n, m);
    //        printf("%.15lf\n", t);
    //        /* заполняет таблицу Excel - ';' для перехода к следующей ячейке */
    //        //fprintf(f, "%lf;%d\n", t, n);
    //        for (int p = 0; p < n; p++) { //восстанавливает порядок следования неизвестных
    //            int c = number[p];
    //            /*в случае изменения порядка следования переменных находит порядковый номер l
    //            элемента, в который записан нужный номер p, и меняет местами элементы с номерами l и p*/
    //            if (c != p) {
    //                int l = p;
    //                while (number[l] != p) {
    //                    l++;
    //                } //меняет элементы местами в массиве для запоминания порядка переменных
    //                int tmp = number[l];
    //                number[l] = c;
    //                number[p] = tmp;
    //                //меняет элементы местами в массиве полученных решений
    //                double tmp1 = x[l];
    //                x[l] = x[p];
    //                x[p] = tmp1;
    //            }
    //        }
    //    }
    //    else printf("Определитель матрицы равен 0"); //проверка условий применимости метода

    //    n *= 2;
    //    m = n + 1;
    //}

    /* ИССЛЕДОВАНИЕ ЗАВИСИМОСТИ АБСОЛЮТНОЙ ПОГРЕШНОСТИ И НЕВЯЗКИ ОТ ЧИСЛА ОБУСЛОВЛЕННОСТИ */

    FILE* fp; //файл с матрицами
    fopen_s(&fp, "file.txt", "r"); //открывает файл из MatLab
    if (!fp) { //проверка успешного открытия файла
        printf("Ошибка при открытии файла\n");
        exit(-1); //выходит из программы с кодом -1
    }

    FILE* f1; //файл для записи
    fopen_s(&f1, "norm.csv", "w"); //открывает файл для записи норм
    if (!f1) { //проверка успешного открытия файла
        printf("Ошибка при открытии файла\n");
        exit(-1); //выходит из программы с кодом -1
    }

    double x_real[10] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 }; //столбец настоящих решений
    for (int k = 0; k < 10; k++) {
        double** A = createMatrix(n, m);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) { //считывает элементы матрицы в двумерный массив
                fscanf(fp, "%lf", &A[i][j]);
            }
        }

        double* b_real = createMatrix(n, 1); //столбец свободных членов исходной матрицы
        double** old = createMatrix(n, m); //исходная матрица
        //запоминает исходную матрицу, до преобразований
        for (int z = 0; z < n; z++) { //копирует все элементы в старую матрицу
            for (int l = 0; l < m; l++) {
                old[z][l] = A[z][l];
            }
            b_real[z] = A[z][n];
        }

        double cond = pow(10, k + 1); //число обусловленностей
        printf("До преобразований:\n");
        print_matrix(A, n, m);
        double* x = (double*)malloc(n * sizeof(double*)); //вектор решений
        double t = gaussMethod(n, m, A, x, number);
        if (t != DET_ERROR) { //обрабатывает ошибку на случай вырожденности матрицы
            printf("После преобразований:\n");
            print_matrix(A, n, m);

            for (int p = 0; p < n; p++) { //восстанавливает порядок следования неизвестных
                int c = number[p];
                /*в случае изменения порядка следования переменных находит порядковый номер l
                элемента, в который записан нужный номер p, и меняет местами элементы с номерами l и p*/
                if (c != p) {
                    int l = p;
                    while (number[l] != p) {
                        l++;
                    }
                    //меняет элементы местами в массиве для запоминания порядка переменных
                    int tmp = number[l];
                    number[l] = c;
                    number[p] = tmp;
                    //меняет элементы местами в массиве полученных решений
                    double tmp1 = x[l];
                    x[l] = x[p];
                    x[p] = tmp1;
                }
            }

            double norm_xx = Norm(x, x_real, n); //вычисляет норму абсолютной погрешности решения
            double* b = multiplyMatrixByColumn(old, x, n); //получает новый столбец свободных чисел
            double norm_AXB = Norm(b, b_real, n); //вычисляет норму невязки

            /* заполняет таблицу Excel - ';' для перехода к следующей ячейке */
            fprintf(f1, "%.20lf;%.20lf;%lf\n", norm_xx, norm_AXB, cond);
        }
        else printf("Определитель матрицы равен 0"); //проверка условий применимости метода

        freeMatrix(A, n); //очищает память
        freeMatrix(old, n); //очищает память
    }

    //fclose(f); //закрывает файл
    fclose(fp); //закрывает файл
    fclose(f1); //закрывает файл

    return 0;
}
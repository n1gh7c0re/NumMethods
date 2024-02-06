#define _USE_MATH_DEFINES
#pragma warning(disable:4996)
#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <locale.h>
#include <math.h>
#include <stdlib.h>
#include <windows.h>
#define MEMORY_ERROR 2 //��� ������ ��������� ������

//����������� �������� ������ ��� ������� ������� n x m
//n - ����� ������� (���-�� �����); m - ����� ������ (���-�� ��������)
double** createMatrix(int n, int m) {
    double** matrix = (double**)malloc(n * sizeof(double*));
    if (matrix == NULL) { //�������� �� ������ � ���������� ������
        return MEMORY_ERROR;
    }
    for (int i = 0; i < n; i++) {
        matrix[i] = (double*)malloc(m * sizeof(double));
    }
    return matrix;
}

//������� ������, ���������� ��� ������� 
void freeMatrix(double** matrix, int n) {
    for (int i = 0; i < n; i++) {
        free(matrix[i]);
    }
    free(matrix);
}

//��������� �� ����� ������� � ��������� ������
double** readMatrix(FILE* file, double** A, int n, int m) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            fscanf(file, "%lf", &A[i][j]);
        }
    }
    return A;
}

//�������� �������
void print_matrix(double** matrix, int n, int m) {
    printf("\n�������:\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            if (j == n) printf("| "); //��� ��������� ������� ��������� �����
            printf("%15.8lf ", matrix[i][j]); //������ ���� 5, 1 ����� ����� �������, ������������ �� ������� ����
        }
        printf("\n");
    }
    printf("\n");
}

//������� ��������� ����� �������� ���� ��������
double NormDiff(double* x1, double* x2, int n) {
    double norm = 0.0;
    for (int i = 0; i < n; i++) { //������� ����� ���������
        norm += (x1[i] - x2[i]) * (x1[i] - x2[i]);
    }
    norm = sqrt(norm); //��������� ������
    return norm;
}

double conjugate_gradient(double** A, double* b, double* x, int n, double eps, int* iter) {
    LARGE_INTEGER freq, start, end;
    QueryPerformanceFrequency(&freq);
    QueryPerformanceCounter(&start); //���������� ����� ������ ���������� ������

    double* r = malloc(n * sizeof(double));
    double* p = malloc(n * sizeof(double));
    double* Ap = malloc(n * sizeof(double));

    //������������� ���������� �����������
    for (int i = 0; i < n; i++) {
        x[i] = 5.0;
    }

    for (int i = 0; i < n; i++) {
        r[i] = b[i];
        for (int j = 0; j < n; j++) {
            r[i] -= A[i][j] * x[j]; //������ �������: r = b - A * x
        }
        p[i] = r[i]; //������������� ���������� �����������
    }

    int Iter = 0;
    double r_new = 1000.0; //������� ����� ������ ������� �������
    //������������ �������
    while (sqrt(r_new) > eps) {
        Iter++;
        double r_old = 0.0; //������� ����� ������� ������� �������
        for (int i = 0; i < n; i++) {
            r_old += r[i] * r[i]; //����������� ������������ ����� � ��������� ������������ �����
        }

        //��������� ������� A �� ������ p: Ap = A * p
        for (int i = 0; i < n; i++) {
            Ap[i] = 0.0;
            for (int j = 0; j < n; j++) {
                Ap[i] += A[i][j] * p[j];
            }
        }

        double alpha = 0.0;
        double pAp = 0.0;
        for (int i = 0; i < n; i++) {
            pAp += p[i] * Ap[i]; //����������� ������������ �����
        }
        alpha = r_old / pAp; //����������� ��� ��������� ���������� ����������� x

        for (int i = 0; i < n; i++) {
            x[i] += alpha * p[i]; //�������� ��������� �����������: x = x + alfa * p
            r[i] -= alpha * Ap[i]; //�������� ����� ������ �������: r = r - alfa * Ap
        }

        r_new = 0.0; //������� ����� ������ ������� �������
        for (int i = 0; i < n; i++) {
            r_new += r[i] * r[i]; //��������� ������������ �����
        }

        double beta = r_new / r_old; //����������� ��� ��������� ������������ ����������� p
        for (int i = 0; i < n; i++) {
            p[i] = r[i] + beta * p[i]; //�������� ����������� �����������: p = r + beta * p
        }
    }
    *iter = Iter;

    free(r);
    free(p);
    free(Ap);

    QueryPerformanceCounter(&end); //���������� ����� ����� ���������� ������
    //������� ����� � �������� �� ���������� ������
    double t = (double)(end.QuadPart - start.QuadPart) / (double)freq.QuadPart;
    return t;
}

int main() {
    setlocale(LC_CTYPE, "Russian");
    int n = 10; //����������� ����
    int a = 1;
    int c = 100; //������� ��������� �������� ��������� �������
    double eps = 0.0; //�������� ����������
    int Iter = 0; //���������� ��������

    /* ����������� ����������� ������� ���������� ������ �� ����������� ������� */

    //FILE* f1; //���� ��� ������
    //fopen_s(&f1, "time4.csv", "w"); //��������� ���� ��� ������ �������
    //if (!f1) { //�������� ��������� �������� �����
    //    printf("������ ��� �������� �����\n");
    //    exit(-1); //������� �� ��������� � ����� -1
    //}

    //for (int d = 0; d < 8; d++) {
    //    double** A = createMatrix(n, n + 1);
    //    for (int i = 0; i < n; i++) { //��������� ������� ���������� �������
    //        A[i][i] = (double)rand() / RAND_MAX * (c - a) + a;
    //        for (int j = i + 1; j < n; j++) { //�������� ������������ ������������-������������ �������
    //            A[i][j] = (double)rand() / RAND_MAX * (c - a) + a;
    //            A[j][i] = A[i][j];
    //        }
    //        A[i][n] = (double)rand() / RAND_MAX * (c - a) + a;
    //    }
    //    double* b = malloc(n * sizeof(double)); //������� ��������� ������ �������� �������
    //    for (int z = 0; z < n; z++) { //�������� ��� �������� ������� ��������� ������
    //        b[z] = A[z][n];
    //    }

    //    double* x = malloc(n * sizeof(double)); //������� �������
    //    eps = 1e-10;

    //    //����� ������ ����������� ����������
    //    double t = conjugate_gradient(A, b, x, n, eps, &Iter);
    //    /* ��������� ������� Excel - ';' ��� �������� � ��������� ������ */
    //    fprintf(f1, "%lf;%d\n", t, n);
    //    printf("%lf\n", t);

    //    //����������� ������
    //    freeMatrix(A, n);
    //    free(b);
    //    free(x);
    //    n *= 2;
    //}

    //fclose(f1);

    /* ������������ ����������� ���������� ����������� � ����� �������� �� �������� */
    
    FILE* fp; //���� � ��������
    fopen_s(&fp, "file.txt", "r"); //��������� ���� �� MatLab
    if (!fp) { //�������� ��������� �������� �����
        printf("������ ��� �������� �����\n");
        exit(-1); //������� �� ��������� � ����� -1
    }

    FILE* f2; //���� ��� ������
    fopen_s(&f2, "norm4.csv", "w"); //��������� ���� ��� ������ ����
    if (!f2) { //�������� ��������� �������� �����
        printf("������ ��� �������� �����\n");
        exit(-1); //������� �� ��������� � ����� -1
    }

    n = 10;
    double x_real[10] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 }; //������� ��������� �������
    double** A = createMatrix(n, n + 1);
    A = readMatrix(fp, A, n, n + 1); //�������������� �������

    double* b = malloc(n * sizeof(double)); //������� ��������� ������ �������� �������
    for (int z = 0; z < n; z++) { //�������� ��� �������� ������� ��������� ������
        b[z] = A[z][n];
    }

    for (int k = 1; k <= 10; k++) {
        double* x = malloc(n * sizeof(double)); //������� �������
        eps = pow(10, -k);
        
        //����� ������ ����������� ����������
        double t = conjugate_gradient(A, b, x, n, eps, &Iter);

        double norm_xx = NormDiff(x, x_real, n); //��������� ����� ���������� ����������� �������
        
        /* ��������� ������� Excel - ';' ��� �������� � ��������� ������ */
        fprintf(f2, "%d;%.20lf;%.*lf\n", Iter, norm_xx, k, eps);
        printf("%d\n", Iter);

        //����������� ������
        free(x);
    }

    freeMatrix(A, n);
    free(b);
    fclose(fp);
    fclose(f2);

    return 0;
}
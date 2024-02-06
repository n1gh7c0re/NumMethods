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

//����������� �������� ������ ��� ������� ������� n x m
//n - ����� ������� (���-�� �����); m - ����� ������ (���-�� ��������)
double** createMatrix(int n, int m) {
    double** matrix = (double**)malloc(n * sizeof(double*));
    if (matrix == NULL) { //�������� �� ������ � ���������� ������
        printf("������ ��������� ������ � ������� createMatrix");
        exit(2); //������� �� ��������� � ����� 2
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

//�������� �������
void print_matrix(double** matrix, int n, int m) {
    printf("\n�������:\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            if (j == n) printf("| "); //��� ��������� ������� ��������� �����
            printf("%15.10lf ", matrix[i][j]); //������ ���� 5, 1 ����� ����� �������, ������������ �� ������� ����
        }
        printf("\n");
    }
    printf("\n");
}

//��������� ����� ������ � ������� �������� �������� �� ������
double gaussMethod(int n, int m, double** a, double* x, int* number) {
    LARGE_INTEGER freq, start, end;
    QueryPerformanceFrequency(&freq);
    QueryPerformanceCounter(&start); //���������� ����� ������ ���������� ������
    
    int i, j, k;
    int s; //���������� ��� ����������� ������ �������
    
    for (i = 0; i < n - 1; i++) { //��� �� �������
        double max = 0.0;
        s = i;
        //����� �������� �������� �� ������
        for (int z = i; z < n; z++) { //������� ��������� �� ���� ������
            if (fabs(a[i][z]) > fabs(max)) {
                max = a[i][z]; //������� ������������ �������
                s = z; //���������� ����� �������
            }
        }
        if (s != i) { // ���� ������� ������� �� �� ���������
            for (j = 0; j < n; j++) { //������ ������� �������
                double temp = a[j][i]; //��������� ����������
                a[j][i] = a[j][s];
                a[j][s] = temp;
            }
            double temp1 = number[i]; //���������� ��������� ������� ���������� ����������
            number[i] = number[s];
            number[s] = temp1;
        }
        //������ ���
        for (k = i + 1; k < n; k++) {
            if (a[i][i] == 0) { //���� ������������ ������� ����� 0
                return DET_ERROR;
            }
            double  term = a[k][i] / a[i][i]; //��������� �������� ��� ���������� � �������������
            for (j = 0; j < m; j++) { //�������� ���� � ������� ��� ������������ ���������
                a[k][j] = a[k][j] - term * a[i][j];
            }
        }
    }
    //�������� ���
    for (i = n - 1; i >= 0; i--) {
        x[i] = a[i][m - 1]; //������ ������ �������
        for (j = i + 1; j < m - 1; j++) { //���� ����� �����
            x[i] = x[i] - a[i][j] * x[j]; //�������� ��������� �������, ��������� ���������� ����� ��������
        }
        //����� �� ������ �� 0
        if (a[i][i] != 0) x[i] = x[i] / a[i][i]; //����������� ����� ����� ������ ���� 1
    }
    QueryPerformanceCounter(&end); //���������� ����� ����� ���������� ������
    //������� ����� � �������� �� ���������� ������
    double t = (double)(end.QuadPart - start.QuadPart) / (double)freq.QuadPart;
    return t;
}

//�������� ������� �� ������ �������
double* multiplyMatrixByColumn(double** matrix, double* column, int n) {
    double* res = createMatrix(n, 1);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            res[i] += matrix[i][j] * column[j];
        }
    }
    return res;
}

//������� ��������� ����� �������� ���� ��������
double Norm(double* x1, double* x2, int n) {
    double norm = 0.0;
    for (int i = 0; i < n; i++) { //������� ����� ���������
        norm += (x1[i] - x2[i]) * (x1[i] - x2[i]);
    }
    norm = sqrt(norm); //��������� ������
    return norm;
}

int main(void) {
    setlocale(LC_CTYPE, "Russian");
    int a = -100, b = 100; //������� ��������� �������� ��������� �������
    int n = 10; //����� ������� (���-�� �����)
    int m = n + 1; //����� ������ (���-�� ��������)
    int* number = (int*)malloc(n * sizeof(int*)); //������ ��� ����������� ������� ����������
    for (int i = 0; i < n; i++) {
        number[i] = i;
    }

    /* ����������� ����������� ������� ���������� ������ �� ����������� ������� */

    //FILE* f; //���� ��� ������
    //fopen_s(&f, "time.csv", "w"); //��������� ���� ��� ������ �������
    //if (!f) { //�������� ��������� �������� �����
    //    printf("������ ��� �������� �����\n");
    //    exit(-1); //������� �� ��������� � ����� -1
    //}
    //
    //for (int z = 1; z <= 10; z++) {
    //    //srand(time(NULL)); //����� ������ ��� �������������� ������ �����
    //    double** matrix = createMatrix(n, m);
    //    for (int i = 0; i < n; i++) { //��������� ������� ���������� �������
    //        for (int j = 0; j < m; j++) {
    //            matrix[i][j] = (double)rand() / RAND_MAX * (b - a) + a;
    //        }
    //    }
    //    printf("�� ��������������:\n");
    //    print_matrix(matrix, n, m);
    //    double* x = (double*)malloc(n * sizeof(double*)); //������ �������
    //    double t = gaussMethod(n, m, matrix, x, number);
    //    if (t != DET_ERROR) { //������������ ������ �� ������ ������������� �������
    //        printf("����� ��������������:\n");
    //        print_matrix(matrix, n, m);
    //        printf("%.15lf\n", t);
    //        /* ��������� ������� Excel - ';' ��� �������� � ��������� ������ */
    //        //fprintf(f, "%lf;%d\n", t, n);
    //        for (int p = 0; p < n; p++) { //��������������� ������� ���������� �����������
    //            int c = number[p];
    //            /*� ������ ��������� ������� ���������� ���������� ������� ���������� ����� l
    //            ��������, � ������� ������� ������ ����� p, � ������ ������� �������� � �������� l � p*/
    //            if (c != p) {
    //                int l = p;
    //                while (number[l] != p) {
    //                    l++;
    //                } //������ �������� ������� � ������� ��� ����������� ������� ����������
    //                int tmp = number[l];
    //                number[l] = c;
    //                number[p] = tmp;
    //                //������ �������� ������� � ������� ���������� �������
    //                double tmp1 = x[l];
    //                x[l] = x[p];
    //                x[p] = tmp1;
    //            }
    //        }
    //    }
    //    else printf("������������ ������� ����� 0"); //�������� ������� ������������ ������

    //    n *= 2;
    //    m = n + 1;
    //}

    /* ������������ ����������� ���������� ����������� � ������� �� ����� ��������������� */

    FILE* fp; //���� � ���������
    fopen_s(&fp, "file.txt", "r"); //��������� ���� �� MatLab
    if (!fp) { //�������� ��������� �������� �����
        printf("������ ��� �������� �����\n");
        exit(-1); //������� �� ��������� � ����� -1
    }

    FILE* f1; //���� ��� ������
    fopen_s(&f1, "norm.csv", "w"); //��������� ���� ��� ������ ����
    if (!f1) { //�������� ��������� �������� �����
        printf("������ ��� �������� �����\n");
        exit(-1); //������� �� ��������� � ����� -1
    }

    double x_real[10] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 }; //������� ��������� �������
    for (int k = 0; k < 10; k++) {
        double** A = createMatrix(n, m);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) { //��������� �������� ������� � ��������� ������
                fscanf(fp, "%lf", &A[i][j]);
            }
        }

        double* b_real = createMatrix(n, 1); //������� ��������� ������ �������� �������
        double** old = createMatrix(n, m); //�������� �������
        //���������� �������� �������, �� ��������������
        for (int z = 0; z < n; z++) { //�������� ��� �������� � ������ �������
            for (int l = 0; l < m; l++) {
                old[z][l] = A[z][l];
            }
            b_real[z] = A[z][n];
        }

        double cond = pow(10, k + 1); //����� ����������������
        printf("�� ��������������:\n");
        print_matrix(A, n, m);
        double* x = (double*)malloc(n * sizeof(double*)); //������ �������
        double t = gaussMethod(n, m, A, x, number);
        if (t != DET_ERROR) { //������������ ������ �� ������ ������������� �������
            printf("����� ��������������:\n");
            print_matrix(A, n, m);

            for (int p = 0; p < n; p++) { //��������������� ������� ���������� �����������
                int c = number[p];
                /*� ������ ��������� ������� ���������� ���������� ������� ���������� ����� l
                ��������, � ������� ������� ������ ����� p, � ������ ������� �������� � �������� l � p*/
                if (c != p) {
                    int l = p;
                    while (number[l] != p) {
                        l++;
                    }
                    //������ �������� ������� � ������� ��� ����������� ������� ����������
                    int tmp = number[l];
                    number[l] = c;
                    number[p] = tmp;
                    //������ �������� ������� � ������� ���������� �������
                    double tmp1 = x[l];
                    x[l] = x[p];
                    x[p] = tmp1;
                }
            }

            double norm_xx = Norm(x, x_real, n); //��������� ����� ���������� ����������� �������
            double* b = multiplyMatrixByColumn(old, x, n); //�������� ����� ������� ��������� �����
            double norm_AXB = Norm(b, b_real, n); //��������� ����� �������

            /* ��������� ������� Excel - ';' ��� �������� � ��������� ������ */
            fprintf(f1, "%.20lf;%.20lf;%lf\n", norm_xx, norm_AXB, cond);
        }
        else printf("������������ ������� ����� 0"); //�������� ������� ������������ ������

        freeMatrix(A, n); //������� ������
        freeMatrix(old, n); //������� ������
    }

    //fclose(f); //��������� ����
    fclose(fp); //��������� ����
    fclose(f1); //��������� ����

    return 0;
}
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

//����������� ��� ���������� ������� ����� �����������
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

//��������� ����� ������� ��� ��������� ����������� �����
double leverrier(double* p, double** A, int n) {
    LARGE_INTEGER freq, start, end;
    QueryPerformanceFrequency(&freq);
    QueryPerformanceCounter(&start); //���������� ����� ������ ���������� ������

    double* S = malloc(n * sizeof(double)); //������ ������ �������
    for (int i = 0; i < n; i++) {
        S[i] = 0.0; //����������
        p[i] = 0.0; //����������
    }

    double** B = createMatrix(n, n); //��������������� ������� ��� ���������� �������� ������� �
    for (int d = 0; d < n; d++) { //�������� �������� ������� � � ������� �
        for (int l = 0; l < n; l++) {
            B[d][l] = A[d][l];
        }
    }

    for (int k = 0; k < n; k++) { //������� ���� ������� � ������� k = 1, 2, ..., n
        for (int j = 0; j < n; j++) {
            S[k] += B[j][j];
        }
        B = multiplyMatrixByMatrix(B, A, n);
    }

    //�������� ������������ �������� p(k) = 1/k * [ S(k) - p(1) * S(k-1) - ... - p(k-1) * S(1) ]
    for (int i = 0; i < n; i++) {
        p[i] = S[i] / (i + 1);
        for (int j = 0; j < i; j++) {
            p[i] -= p[j] * S[i - j - 1] / (i + 1);
        }
    }

    QueryPerformanceCounter(&end); //���������� ����� ����� ���������� ������
    //������� ����� � �������� �� ���������� ������
    double t = (double)(end.QuadPart - start.QuadPart) / (double)freq.QuadPart;
    return t;
}

int main(void) {
    setlocale(LC_CTYPE, "Russian");
    int n = 10; //����������� ����
    int a = 1;
    int c = 10; //������� ��������� �������� ��������� �������

    /* ����������� ����������� ������� ���������� ������ �� ����������� ������� */

    //FILE* f2; //���� ��� ������
    //fopen_s(&f2, "time5.csv", "w"); //��������� ���� ��� ������ �������
    //if (!f2) { //�������� ��������� �������� �����
    //    printf("������ ��� �������� �����\n");
    //    exit(-1); //������� �� ��������� � ����� -1
    //}

    //for (int d = 0; d < 6; d++) {
    //    double** A = createMatrix(n, n);
    //    for (int i = 0; i < n; i++) { //��������� ������� ���������� �������
    //        for (int j = 0; j < n; j++) {
    //            A[i][j] = (double)rand() / RAND_MAX * (c - a) + a;
    //        }
    //    }
    //    double* p = malloc(n * sizeof(double)); //������ ������������� ������������� ��������

    //    double t = leverrier(p, A, n); //�������� ����� �������

    //    /* ��������� ������� Excel - ';' ��� �������� � ��������� ������ */
    //    fprintf(f2, "%lf;%d\n", t, n);
    //    printf("%lf\n", t);

    //    freeMatrix(A, n);
    //    free(p);
    //    n *= 2;
    //}
    //fclose(f2);
    
    /* ������������ ������������ ������ */

    //FILE* f1; //���� ��� ������
    //fopen_s(&f1, "polynomial1.txt", "w"); //��������� ���� ��� ������ ������������� ������������ �������
    //if (!f1) { //�������� ��������� �������� �����
    //    printf("������ ��� �������� �����\n");
    //    exit(-1); //������� �� ��������� � ����� -1
    //}

    //FILE* fp; //���� � ���������
    //fopen_s(&fp, "file1.txt", "r"); //��������� ���� �� MatLab
    //if (!fp) { //�������� ��������� �������� �����
    //    printf("������ ��� �������� �����\n");
    //    exit(-1); //������� �� ��������� � ����� -1
    //}

    //for (int d = 0; d < 3; d++) {
    //    double** A = createMatrix(n, n);
    //    A = readMatrix(fp, A, n, n);
    //    print_matrix(A, n, n);

    //    double* p = malloc(n * sizeof(double)); //������ ������������� ������������� ��������

    //    double t = leverrier(p, A, n); //�������� ����� �������

    //    double* p_true = malloc((n + 1) * sizeof(double)); //������ ������������� ������������� ��������
    //    p_true[0] = 1; //�������� ������� ������ ������������ ����� 1


    //    for (int i = 0; i < n; i++) {
    //        p_true[i + 1] = -p[i]; //������ ����
    //    }

    //    for (int i = 0; i < n + 1; i++) { //���������� � ��������� ���� ������������
    //        fprintf(f1, "%.15lf ", p_true[i]);
    //    }
    //    fprintf(f1, "\n");

    //    freeMatrix(A, n);
    //    free(p);
    //    free(p_true);
    //}
    //fclose(f1);
    //fclose(fp);

    /* ������������ �� ������������ ������ */

    FILE* f3; //���� ��� ������
    fopen_s(&f3, "polynomial2.txt", "w"); //��������� ���� ��� ������ ������������� �������������� �������
    if (!f3) { //�������� ��������� �������� �����
        printf("������ ��� �������� �����\n");
        exit(-1); //������� �� ��������� � ����� -1
    }

    FILE* fp2; //���� � ���������
    fopen_s(&fp2, "file2.txt", "r"); //��������� ���� �� MatLab
    if (!fp2) { //�������� ��������� �������� �����
        printf("������ ��� �������� �����\n");
        exit(-1); //������� �� ��������� � ����� -1
    }

    for (int d = 0; d < 3; d++) {
        double** A = createMatrix(n, n);
        A = readMatrix(fp2, A, n, n);
        print_matrix(A, n, n);

        double* p = malloc(n * sizeof(double)); //������ ������������� ������������� ��������

        double t = leverrier(p, A, n); //�������� ����� �������

        double* p_true = malloc((n + 1) * sizeof(double)); //������ ������������� ������������� ��������
        p_true[0] = 1; //�������� ������� ������ ������������ ����� 1


        for (int i = 0; i < n; i++) {
            p_true[i + 1] = -p[i]; //������ ����
        }

        for (int i = 0; i < n + 1; i++) { //���������� � ��������� ���� ������������
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
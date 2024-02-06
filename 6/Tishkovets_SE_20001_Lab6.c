#define _USE_MATH_DEFINES
#pragma warning(disable:4996)
#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <locale.h>
#include <math.h>
#include <stdlib.h>
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

//��������� �� ����� ������� � ��������� ������
double** readMatrix(FILE* file, double** A, int n, int m) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            fscanf(file, "%lf", &A[i][j]);
        }
    }
    return A;
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

//�������� ���� ������ ������� (�����������������) �� ������
double multiplyColumnByColumn(double* x1, double* x2, int n) {
    double res = 0.0;
    for (int i = 0; i < n; i++) {
        res += x1[i] * x2[i];
    }
    return res;
}

//������� ��������� ����� �������
double Norm(double* x, int n) {
    double norm = 0.0;
    for (int i = 0; i < n; i++) { //������� ����� ���������
        norm += x[i] * x[i];
    }
    norm = sqrt(norm); //��������� ������
    return norm;
}

//��������� ��������� ����� � ����������� �������
double power_method_shift(double** A, double* l_real, double eps, int* iter, int n) {
    double* x_prev = (double*)malloc(n * sizeof(double*)); //��������� �����������
    double* x = (double*)malloc(n * sizeof(double*)); //��������� �����������
    for (int i = 0; i < n; i++) {
        x_prev[i] = 1.0; //������������� ���������� �����������
    }

    double shift = (l_real[1] + l_real[n - 1]) / 2; //����������� �����

    double** B = createMatrix(n, n); //����� �������, ��� �� �������
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i != j) B[i][j] = A[i][j];
            else B[i][j] = A[i][j] - shift; // B = A - shift * E
        }
    }

    double lambda = 0; //������
    double lambda_prev = 10; //������ �� ���������� ��������
    int Iter = 0; //����� ��������

    while (fabs(lambda - lambda_prev) > eps) {
        x = multiplyMatrixByColumn(B, x_prev, n); //��������� �����������

        double temp1 = multiplyColumnByColumn(x, x_prev, n); //���������
       // double temp2 = multiplyColumnByColumn(x_prev, x_prev, n); //�����������

        lambda_prev = lambda;
        lambda = temp1; //��������� ������ lambda_k+1 = (x_k+1, x_k) / (x_k, x_k)
            // ��� ��� ������� �����������, temp2 = 1

        double norm = Norm(x, n); //������� ��������� ����� �������
        for (int i = 0; i < n; i++) {
            x[i] /= norm; //��������� ������
        }

        for (int i = 0; i < n; i++) { //��������� �����������
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

    /* ������������ ����������� ���������� ����������� � ����� �������� �� �������� */

    FILE* fp; //���� � 3 ���������
    fopen_s(&fp, "file1.txt", "r"); //��������� ���� �� MatLab
    if (!fp) { //�������� ��������� �������� �����
        printf("������ ��� �������� �����\n");
        exit(-1); //������� �� ��������� � ����� -1
    }

    FILE* f1; //���� ��� ������
    fopen_s(&f1, "norm6.csv", "w"); //��������� ���� ��� ������ ����
    if (!f1) { //�������� ��������� �������� �����
        printf("������ ��� �������� �����\n");
        exit(-1); //������� �� ��������� � ����� -1
    }

    double lambda_real[10] = { 100, 10, 9, 8, 7, 6, 5, 4, 3, 2 }; //������� ��������� ��

    for (int s = 0; s < 3; s++) {
        double** A = createMatrix(n, n);
        A = readMatrix(fp, A, n, n);
        print_matrix(A, n, n);

        for (int d = 1; d <= 10; d++) {
            eps = pow(10, -d);

            //����� ���������� ������ �� �������
            double lambda = power_method_shift(A, lambda_real, eps, &Iter, n);

            double error = fabs(lambda - lambda_real[0]);

            printf("������������ ��: %.15lf\n", lambda);

            /* ��������� ������� Excel - ';' ��� �������� � ��������� ������ */
            fprintf(f1, "%.10lf;%d;%.15lf\n", eps, Iter, error);
        }
        fprintf(f1, "\n");
    }

    fclose(f1);
    fclose(fp);

    return 0;
}
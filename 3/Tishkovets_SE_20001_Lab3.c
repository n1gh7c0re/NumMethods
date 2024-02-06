#define _USE_MATH_DEFINES
#pragma warning(disable:4996)
#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <locale.h>
#include <math.h>
#include <stdlib.h>
#include <windows.h>
#include <time.h>
#define MEMORY_ERROR 2 //��� ������ ��������� ������
#define DET_ERROR -2   //��� ��������� ������� ������������ ������

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

//������� ��������� ����� �������
double Norm(double* x, int n) {
    double norm = 0.0;
    for (int i = 0; i < n; i++) { //������� ����� ���������
        norm += x[i] * x[i];
    }
    norm = sqrt(norm); //��������� ������
    return norm;
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

//������������� ������� � ����������� n
double** transpose_matrix(double** A, int n) {
    double** B = createMatrix(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            B[j][i] = A[i][j];
        }
    }
    return B;
}

/* �������� �������� ������� � ������� ������ ������ - ��������:
   1 2 3 | 1 0 0     1 0 0 | 2 8 9
   4 5 6 | 0 1 0 --> 0 1 0 | 3 6 5
   7 8 9 | 0 0 1     0 0 1 | 1 7 4 */
double** inverse_matrix(double** A, int n) {
    double r, a;
    double** matrix = createMatrix(n, 2 * n);
    for (int i = 0; i < n; i++) { //�������� �������� �������� �������
        for (int j = 0; j < n; j++) {
            matrix[i][j] = A[i][j];
        }
    }

    for (int i = 0; i < n; i++) {
        //�������� ��������� ������� ������ �� ��������
        for (int j = n; j < 2 * n; j++) {
            if (i == (j - n)) matrix[i][j] = 1;
            else matrix[i][j] = 0;
        }
    }

    //������ ��� ������ ������
    for (int i = 0; i < n; i++) { //���� �� ��������
        for (int j = 0; j < n; j++) { //�� �������
            if (i != j) { //�������� ���� ��� � ��� ����������, �.�. j ���������� � 0
                r = matrix[j][i] / matrix[i][i];
                for (int k = 0; k < 2 * n; k++) {
                    matrix[j][k] -= r * matrix[i][k];
                }
            }
        }
    }

    //�������� ��������� ������� �����
    for (int i = 0; i < n; i++) {
        a = matrix[i][i];
        for (int j = 0; j < 2 * n; j++) {
            matrix[i][j] /= a;
        }
    }

    double** A_1 = createMatrix(n, n); //�������� �������
    for (int i = 0; i < n; i++) { //�������� �������� ���������� ������� ����������� n
        for (int j = 0; j < n; j++) {
            A_1[i][j] = matrix[i][n + j];
        }
    }
    freeMatrix(matrix, n);

    return A_1;
}

//��������� ���������������� ����� �����-������
double gramSchmidt(double** A, double** Q, double** R, int n) {
    LARGE_INTEGER freq, start, end;
    QueryPerformanceFrequency(&freq);
    QueryPerformanceCounter(&start); //���������� ����� ������ ���������� ������

    for (int j = 0; j < n; j++) {
        for (int s = 0; s < n; s++) { //�������� j-� ������� ������� A � j-� ������� ������� Q
            Q[s][j] = A[s][j];
        }

        double* Q_j = (double*)malloc(n * sizeof(double*));
        for (int k = 0; k < n; k++) { //�������� j-� ������� ������� Q
            Q_j[k] = Q[k][j];
        }

        for (int i = 0; i < j; i++) {
            double* Q_i = (double*)malloc(n * sizeof(double*));
            for (int z = 0; z < n; z++) { //�������� i-� ������� ������� Q
                Q_i[z] = Q[z][i];
            }

            R[i][j] = multiplyColumnByColumn(Q_i, Q_j, n);
            for (int t = 0; t < n; t++) {
                Q_j[t] = Q_j[t] - R[i][j] * Q_i[t];
            }
            free(Q_i);
        }

        R[j][j] = Norm(Q_j, n);
        if (R[j][j] == 0) return DET_ERROR; //������� ������� � �� �������� ������� ������������

        for (int l = 0; l < n; l++) {
            Q_j[l] = Q_j[l] / R[j][j];
        }

        for (int h = 0; h < n; h++) { //��������� j-� ��� ����� ������� ������� Q
            Q[h][j] = Q_j[h];
        }
        free(Q_j);
    }

    QueryPerformanceCounter(&end); //���������� ����� ����� ���������� ������
    //������� ����� � �������� �� ���������� ������
    double t = (double)(end.QuadPart - start.QuadPart) / (double)freq.QuadPart;
    return t;
}

int main(void) {
    setlocale(LC_CTYPE, "Russian");
    int n = 10;
    int m = n + 1;
    int a = -100, b = 100; //������� ��������� �������� ��������� �������

    /* ����������� ����������� ������� ���������� ������ �� ����������� ������� */

    //FILE* f2; //���� ��� ������
    //fopen_s(&f2, "time3.csv", "w"); //��������� ���� ��� ������ ����
    //if (!f2) { //�������� ��������� �������� �����
    //    printf("������ ��� �������� �����\n");
    //    exit(-1); //������� �� ��������� � ����� -1
    //}

    //for (int d = 0; d < 8; d++) {
    //    double** matrix = createMatrix(n, m);
    //    for (int i = 0; i < n; i++) { //��������� ������� ���������� �������
    //        for (int j = 0; j < m; j++) {
    //            matrix[i][j] = (double)rand() / RAND_MAX * (b - a) + a;
    //        }
    //    }

    //    double* b_real = createMatrix(n, 1); //������� ��������� ������ �������� �������
    //    for (int z = 0; z < n; z++) { //�������� ��� ��������
    //        b_real[z] = matrix[z][n];
    //    }

    //    double** Q = createMatrix(n, n); //������������� �������
    //    double** R = createMatrix(n, n); //������� �������� �� ������������������ ������ � ���������
    //    for (int i = 0; i < n; i++) { //��������� ������� R � Q ������
    //        for (int j = 0; j < n; j++) {
    //            R[i][j] = 0.0;
    //            Q[i][j] = 0.0;
    //        }
    //    }

    //    double t = gramSchmidt(matrix, Q, R, n);
    //    if (t != DET_ERROR) {
    //        double** QT = transpose_matrix(Q, n); //����������������� ������� 
    //        double** R_1 = inverse_matrix(R, n); //�������� �������

    //        double* x = (double*)malloc(n * sizeof(double*)); //������ �������
    //        x = multiplyMatrixByColumn(multiplyMatrixByMatrix(R_1, QT, n), b_real, n);

    //        /* ��������� ������� Excel - ';' ��� �������� � ��������� ������ */
    //        fprintf(f2, "%lf;%d\n", t, n);
    //    }
    //    else printf("������������ ������� ����� 0"); //�������� ������� ������������ ������

    //    n *= 2;
    //    m = n + 1;
    //}
    //fclose(f2); //��������� ����

    /* ������������ ����������� ���������� ����������� � ������� �� ����� ��������������� */

    FILE* fp; //���� � ���������
    fopen_s(&fp, "file.txt", "r"); //��������� ���� �� MatLab
    if (!fp) { //�������� ��������� �������� �����
        printf("������ ��� �������� �����\n");
        exit(-1); //������� �� ��������� � ����� -1
    }

    FILE* f1; //���� ��� ������
    fopen_s(&f1, "norm3.csv", "w"); //��������� ���� ��� ������ ����
    if (!f1) { //�������� ��������� �������� �����
        printf("������ ��� �������� �����\n");
        exit(-1); //������� �� ��������� � ����� -1
    }

    double x_real[10] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 }; //������� ��������� �������

    for (int d = 0; d < 8; d++) {
        double** A = createMatrix(n, m);
        A = readMatrix(fp, A, n, m);

        double* b_real = createMatrix(n, 1); //������� ��������� ������ �������� �������
        for (int z = 0; z < n; z++) { //�������� ��� ��������
            b_real[z] = A[z][n];
        }

        double cond = pow(10, d); //����� ���������������
        printf("��������:\n");
        print_matrix(A, n, m);

        double** Q = createMatrix(n, n); //������������� �������
        double** R = createMatrix(n, n); //������� ����������� �������
        for (int i = 0; i < n; i++) { //��������� ������� R � Q ������
            for (int j = 0; j < n; j++) {
                R[i][j] = 0.0;
                Q[i][j] = 0.0;
            }
        }

        double t = gramSchmidt(A, Q, R, n);
        if (t != DET_ERROR) {
            double** QT = transpose_matrix(Q, n); //����������������� ������� 
            double** R_1 = inverse_matrix(R, n); //�������� �������

            double* x = (double*)malloc(n * sizeof(double*)); //������ �������
            x = multiplyMatrixByColumn(multiplyMatrixByMatrix(R_1, QT, n), b_real, n);

            double norm_xx = NormDiff(x, x_real, n); //��������� ����� ���������� ����������� �������
            double* b = multiplyMatrixByColumn(A, x, n); //�������� ����� ������� ��������� �����
            double norm_AXB = NormDiff(b, b_real, n); //��������� ����� �������

            //�������� �������
            for (int s = 0; s < n; s++) {
                printf("%.20lf\n", x[s]);
            }

            /* ��������� ������� Excel - ';' ��� �������� � ��������� ������ */
            fprintf(f1, "%.20lf;%.20lf;%lf\n", norm_xx, norm_AXB, cond);

            //������ ������
            freeMatrix(QT, n);
            freeMatrix(R_1, n);
            free(x);
            free(b);
        }
        else printf("������������ ������� ����� 0"); //�������� ������� ������������ ������

        //������ ������
        freeMatrix(A, n);
        freeMatrix(Q, n);
        freeMatrix(R, n);
        free(b_real);
    }
    fclose(fp); //��������� ����
    fclose(f1); //��������� ����

    return 0;
}
#define _USE_MATH_DEFINES
#pragma warning(disable:4996)
#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <locale.h>
#include <math.h>

//���������� �������� �������, n - ��� ������������� ���� ��� ����� ���������
//1 - ��������������, 2 - ���������������
double f(double x, int n) {
	if (n == 1) {
		return pow(x, 4) - 18 * x - 10; //�������������� ���������
	}
	else {
		return pow(5, x) - 6 * x - 7; //��������������� ��������
	}
}

//���������� �������� 1-�� ����������� �������, n - ��� ������������� ���� ��� ����� ���������
//1 - ��������������, 2 - ���������������
double Diff(double x, int n) {
	if (n == 1) {
		return 4 * pow(x, 3) - 18; //����������� ��������������� ���������
	}
	else {
		return pow(5, x) * log(5) - 6; //����������� ���������������� ��������
	}
}

//����� ����������� �������
double BisectionMethod(double a, double b, double epsilon, int n, int* count) {
	double fa = 0.0, fc = 0.0;
	double c = 0.0;
	int cnt = 0; //�������
	while (fabs(b - a) > 2 * epsilon) { //������� ��������� �����
		c = (a + b) / 2; //�������� �������
		fa = f(a, n);
		fc = f(c, n);
		if (fc == 0) { //������� ������ �� �����
			return c;
		}
		if (fa * fc < 0) { //������� ������ ����������� ����������
			b = c;
		}
		else a = c;
		cnt++; //������� ���������� ��������
	}
	*count = cnt;
	return (a + b) / 2;
}

//����� ������� (����� �����������)
double NewtonsMethod(double x0, double epsilon, int n, int* count) {
	int cnt = 0; //�������
	double f1 = f(x0, n);
	double f2 = Diff(x0, n);
	while (fabs(f1 / f2) > epsilon) { //������� ��������� �����
		x0 = x0 - f1 / f2;
		f1 = f(x0, n);
		f2 = Diff(x0, n);
		cnt++; //������� ���������� ��������
	}
	*count = cnt;
	return x0;
}

int main(void) {
	setlocale(LC_CTYPE, "Russian");

	FILE* f1; //���� ��� ������
	fopen_s(&f1, "norm1.csv", "w"); //��������� ���� ��� ������ ����
	if (!f1) { //�������� ��������� �������� �����
		printf("������ ��� �������� �����\n");
		exit(-1); //������� �� ��������� � ����� -1
	}

	double e1 = 0.0, e2 = pow(10, -15); //�����������
	int count1 = 0, count2 = 0, count3 = 0, count4 = 0; //�������� ��� ���������� ��������
	//����� �������������� �������
	double x1_real = 2.78457554401786;
	double x2_real = -0.550455030617542;
	//����� ��������������� �������
	double x3_real = 1.786009501751085;
	double x4_real = -1.140060605624304;

	printf("\t����� ����������� �������\n");

	printf("\n�������������� ��������� x^4 - 18 * x - 10 = 0\n");
	for (int i = 1; i <= 10; i++) {
		double x1 = 0.0;
		double x2 = 0.0;
		e1 = pow(10, -i);
		x1 = BisectionMethod(0.641, 3.621, e1, 1, &count1);
		x2 = BisectionMethod(-1.778, -0.357, e1, 1, &count2);
		printf("%.10lf   %.10lf\n", x1, x2);
		/* ��������� ������� Excel - ';' ��� �������� � ��������� ������ */
		fprintf(f1, "%.10lf;%d;%.15lf;;%d;%.15lf\n", e1, count1, fabs(x1 - x1_real), count2, fabs(x2 - x2_real));
	}
	fprintf(f1, "\n");

	printf("\n��������������� ��������� 5^x - 6 * x = 7\n");
	for (int i = 1; i <= 10; i++) {
		double x3 = 0.0;
		double x4 = 0.0;
		e1 = pow(10, -i);
		x3 = BisectionMethod(1.5, 2.5, e1, 2, &count1);
		x4 = BisectionMethod(-2, -0.5, e1, 2, &count2);
		printf("%.10lf   %.10lf\n", x3, x4);
		/* ��������� ������� Excel - ';' ��� �������� � ��������� ������ */
		fprintf(f1, "%.10lf;%d;%.15lf;;%d;%.15lf\n", e1, count1, fabs(x3 - x3_real), count2, fabs(x4 - x4_real));
	}
	fprintf(f1, "\n");

	printf("\n\t\t����� �������\n");

	printf("\n�������������� ��������� x^4 - 18 * x - 10 = 0\n");
	
	for (int i = 1; i <= 10; i++) {
		double x5 = 0.0;
		double x6 = 0.0;
		e1 = pow(10, -i);
		x5 = NewtonsMethod(3, e1, 1, &count1);
		x6 = NewtonsMethod(-1, e1, 1, &count2);
		printf("%.20lf   %.20lf\n", x5, x6);
		/* ��������� ������� Excel - ';' ��� �������� � ��������� ������ */
		fprintf(f1, "%.10lf;%d;%.15lf;;%d;%.15lf\n", e1, count1, fabs(x5 - x1_real), count2, fabs(x6 - x2_real));
	}
	fprintf(f1, "\n");

	printf("\n��������������� ��������� 5^x - 6 * x = 7\n");
	
	for (int i = 1; i <= 10; i++) {
		double x7 = 0.0;
		double x8 = 0.0;
		e1 = pow(10, -i);
		x7 = NewtonsMethod(2, e1, 2, &count1);
		x8 = NewtonsMethod(-2, e1, 2, &count2);
		printf("%.20lf   %.20lf\n", x7, x8);
		/* ��������� ������� Excel - ';' ��� �������� � ��������� ������ */
		fprintf(f1, "%.10lf;%d;%.15lf;;%d;%.15lf\n", e1, count1, fabs(x7 - x3_real), count2, fabs(x8 - x4_real));
	}

	fclose(f1);
	return 0;
}
#define _USE_MATH_DEFINES
#pragma warning(disable:4996)
#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <locale.h>
#include <math.h>
#include <stdlib.h>
#define MEMORY_ERROR 2 //код ошибки выделения памяти

typedef struct {
	double** matrix;
	int height;
	int weight;
}matrix_t;

typedef struct {
	matrix_t* matrix;
	double betta;
	matrix_t* amega;//здесь амега на амега транспонир (сочет св-во)
}reflection_matrix_t;

double norma(matrix_t* matrix) {
	double norm = 0;
	for (int i = 0; i < matrix->height; i++) {
		for (int j = 0; j < matrix->weight; j++) {
			norm += matrix->matrix[i][j] * matrix->matrix[i][j];
		}
	}
	return pow(norm, 0.5);
}

void creatematrix(matrix_t* matrix, int i, int j) {
	matrix->height = i;
	matrix->weight = j;
	matrix->matrix = malloc(sizeof(double*) * i);
	for (int c = 0; c < i; c++) {
		double* ptr = malloc(sizeof(double) * j);
		matrix->matrix[c] = ptr;
	}
}

//считывает из файла матрицу в двумерный массив
matrix_t** read(FILE* file, matrix_t* A, int n, int m) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			fscanf(file, "%lf", &A->matrix[i][j]);
		}
	}
	return A;
}

matrix_t* naturalle_matrix(int i, int j) {
	matrix_t* matrix = malloc(sizeof(matrix_t));
	if (i != j) {
		char massage[57] = "matrix is not squary so you cant create naturalle_matrix";
	}
	creatematrix(matrix, i, j);
	for (int i0 = 0; i0 < i; i0++) {
		for (int j0 = 0; j0 < j; j0++) {
			if (i0 == j0) matrix->matrix[i0][j0] = 1;
			else matrix->matrix[i0][j0] = 0;
		}
	}
	return matrix;
}

matrix_t* transposematrix(matrix_t* matrix) {
	matrix_t* result = malloc(sizeof(matrix_t));
	creatematrix(result, matrix->weight, matrix->height);
	for (int i = 0; i < matrix->height; i++) {
		for (int j = 0; j < matrix->weight; j++) {
			result->matrix[j][i] = matrix->matrix[i][j];
		}
	}
	return result;
}

matrix_t* multiplydigit(matrix_t* matrix, double digit) {
	matrix_t* result = malloc(sizeof(matrix_t));
	creatematrix(result, matrix->height, matrix->weight);
	for (int i = 0; i < matrix->height; i++) {
		for (int j = 0; j < matrix->weight; j++) {
			result->matrix[i][j] = matrix->matrix[i][j] * digit;
		}
	}
	return result;
}

matrix_t* sum_matrix(matrix_t* first, matrix_t* second, int sign) {//0 = +; 1 = -
	matrix_t* result = malloc(sizeof(matrix_t));
	creatematrix(result, first->height, first->weight);
	for (int i = 0; i < first->height; i++) {
		for (int j = 0; j < first->weight; j++) {
			result->matrix[i][j] = first->matrix[i][j] + pow(-1, sign) * second->matrix[i][j];
		}
	}
	return result;
}

matrix_t* multiplymatrix(matrix_t* first, matrix_t* second) {
	matrix_t* result = malloc(sizeof(matrix_t));
	creatematrix(result, first->height, second->weight);
	if (first->weight != second->height) return result;
	for (int string = 0; string < first->height; string++) {
		for (int column = 0; column < second->weight; column++) {
			double sum = 0;
			for (int i = 0; i < first->weight; i++) {
				sum += first->matrix[string][i] * second->matrix[i][column];
			}
			result->matrix[string][column] = sum;
		}
	}
	return result;
}

matrix_t* getting_column(matrix_t* matrix, int count) {
	matrix_t* new_matrix = malloc(sizeof(matrix_t));
	creatematrix(new_matrix, matrix->height - count, 1);
	for (int i = count; i < matrix->height; i++) {
		new_matrix->matrix[i - count][0] = matrix->matrix[i][count - 1];
	}
	return new_matrix;
}

reflection_matrix_t getting_reflection(matrix_t* matrix, int count) {
	reflection_matrix_t reflection;
	matrix = getting_column(matrix, count);
	count = 0;
	double alfa;
	if (matrix->matrix[count][count] != 0) {
		alfa = -matrix->matrix[count][count] / matrix->matrix[0][0] * norma(matrix);
	}
	else alfa = norma(matrix);
	matrix->matrix[count][count] -= alfa;
	reflection.amega = multiplymatrix(matrix, transposematrix(matrix));
	reflection.betta = 2 / (pow(norma(matrix), 2));
	reflection.matrix = sum_matrix(naturalle_matrix(matrix->height - count, matrix->height - count), multiplydigit(reflection.amega, reflection.betta), 1);
	return reflection;
}

reflection_matrix_t expand_reflection(reflection_matrix_t reflection, int size) {
	reflection_matrix_t expand_reflect;
	int difference = size - reflection.amega->weight;
	expand_reflect.matrix = malloc(sizeof(matrix_t));
	creatematrix(expand_reflect.matrix, size, size);
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			if (i < size - reflection.matrix->height) {
				if (i == j) expand_reflect.matrix->matrix[i][j] = 1;
				else expand_reflect.matrix->matrix[i][j] = 0;
			}
			else if (j < size - reflection.matrix->weight) expand_reflect.matrix->matrix[i][j] = 0;
			else {
				expand_reflect.matrix->matrix[i][j] = reflection.matrix->matrix[i - difference][j - difference];
			}
		}
	}
	expand_reflect.amega = sum_matrix(naturalle_matrix(size, size), expand_reflect.matrix, 1);
	expand_reflect.betta = 1;
	return expand_reflect;
}

matrix_t* get_matrix_H(matrix_t* matrix) {
	for (int count = 1; count < matrix->height - 1; count++) {
		reflection_matrix_t reflection = expand_reflection(getting_reflection(matrix, count), matrix->height);
		matrix = multiplymatrix(multiplymatrix(reflection.matrix, matrix), reflection.matrix);
	}
	return matrix;
}

double count_coeff(matrix_t* matrix, int start, int n) {
	double res = 1;
	for (int i = start; i < n; i++) {
		res *= matrix->matrix[i][i - 1];
	}
	return res;
}

double* determinant(matrix_t* matrix, int n) {
	double* coeff = malloc(sizeof(double) * (n + 1));
	if (n == 0) {
		coeff[0] = 0;
		return coeff;
	}
	for (int i = 0; i < n; i++) {
		coeff[i] = pow(-1, n + 1 + i) * matrix->matrix[i][n - 1] * count_coeff(matrix, i + 1, n);
		double* tmp_coeff = determinant(matrix, i);
		for (int k = 0; k < i; k++) {
			double tmp = tmp_coeff[k] * coeff[i];//1 * 3 = 4 степень
			coeff[k] += tmp;
			if (k == i - 1) coeff[i] *= tmp_coeff[i];
		}
		if (i == n - 1) {//надо учесть что умножается скобка с минус лямда
			for (int k = 1; k < n; k++) {
				coeff[k] -= tmp_coeff[k - 1];
			}
		}
	}
	coeff[n] = pow(-1, n);//тк -лямда
	return coeff;
}

/* А ТЕПЕРЬ НОРМАЛЬНЫЙ КОД */

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

//реализует модифицированный метод Грама-Шмидта
void gramSchmidt(double** A, double** Q, double** R, int n) {
	for (int i = 0; i < n; i++) { //заполняет матрицы R и Q нулями
		for (int j = 0; j < n; j++) {
			R[i][j] = 0.0;
			Q[i][j] = 0.0;
		}
	}

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
		if (R[j][j] == 0) return -1; //столбцы матрицы А не являются линейно независимыми

		for (int l = 0; l < n; l++) {
			Q_j[l] = Q_j[l] / R[j][j];
		}

		for (int h = 0; h < n; h++) { //заполняем j-й уже новый столбец матрицы Q
			Q[h][j] = Q_j[h];
		}
		free(Q_j);
	}
}

//реализует метод QR со сдвигом для нахождения собственных чисел
double** QR_shift(double** A, int n, double eps, int* iter, int flag) {
	double** B = createMatrix(n, n);
	for (int i = 0; i < n; i++) { //копирует элементы матрицы
		for (int j = 0; j < n; j++) {
			B[i][j] = A[i][j];
		}
	}

	double* x = (double*)malloc((n - 1) * sizeof(double*));
	for (int i = 0; i < n - 1; i++) {
		x[i] = B[i + 1][i]; //поддиагональные элементы
	}
	int Iter = 0;

	double shift = 0;
	while (Norm(x, n - 1) > eps) { //пока поддиагональные элементы не близки к нулю
		if (flag == 0) shift = B[n - 1][n - 1]; //сдвиг
		else if (flag == 1) shift = 10;
		else if (flag == 2) shift = -10;

		double** C = createMatrix(n, n); //временная матрица
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (i == j) C[i][j] = B[i][j] + shift; // C = B + shift * E
				else C[i][j] = B[i][j];
			}
		}

		double** Q = createMatrix(n, n); //ортогональная матрица
		double** R = createMatrix(n, n); //верхняя треугольная матрица
		//получает разложение матрицы C на Q и R
		gramSchmidt(C, Q, R, n);

		//получает матрицу для следующей итерации A = R * Q - shift * E
		B = multiplyMatrixByMatrix(R, Q, n);
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (i == j) B[i][j] -= shift; // B = B - shift * E
			}
		}

		//обновляет поддиагональные элементы
		for (int i = 0; i < n - 1; i++) {
			x[i] = B[i + 1][i];
		}

		Iter++;

		//освобождает память
		freeMatrix(Q, n);
		freeMatrix(R, n);
		freeMatrix(C, n);
	}
	printf("\nСдвиг %lf", shift);

	*iter = Iter;

	free(x);
	return B;
}

int main(void) {
	setlocale(LC_CTYPE, "Russian");
	int n = 10; //размерность
	double eps = 10e-10; //точность
	int Iter = 0;

	FILE* fp; //файл с матрицами
	fopen_s(&fp, "file.txt", "r"); //открывает файл из MatLab
	if (!fp) { //проверка успешного открытия файла
		printf("Ошибка при открытии файла\n");
		exit(-1); //выходит из программы с кодом -1
	}

	FILE* f1; //файл для записи
	fopen_s(&f1, "norm7.csv", "w"); //открывает файл для записи норм
	if (!f1) { //проверка успешного открытия файла
		printf("Ошибка при открытии файла\n");
		exit(-1); //выходит из программы с кодом -1
	}

	//приведение исходной матрицы к матрице Хессенберга
	matrix_t* matrix = malloc(sizeof(matrix_t));
	matrix_t* matrix_h = malloc(sizeof(matrix_t));
	creatematrix(matrix, n, n);
	creatematrix(matrix_h, n, n);
	matrix = read(fp, matrix, n, n);
	matrix_h = get_matrix_H(matrix); //получает матрицу Хессенберга

	//копирует элементы старой матрицы в новую (без тупых структур - через двойной указатель)
	double** A = createMatrix(n, n);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			A[i][j] = matrix_h->matrix[i][j];
		}
	}

	double lambda_real[10] = { 10, 9, 8, 7, 6, 5, 4, 3, 2, 1 }; //столбец настоящих сч
	
	for (int flag = 0; flag < 4; flag++) {
		eps = 1;
		for (int d = 1; d <= 10; d++) {
			eps /= 10;

			double** B = createMatrix(n, n);
			//вызывает QR метод со сдвигом
			B = QR_shift(A, n, eps, &Iter, flag); //на диагонали стоят сч

			printf("Число итераций %d\n", Iter);

			double* lambda = (double*)malloc(n * sizeof(double*));
			if (flag == 2) {
				for (int i = 0; i < n; i++) {
					lambda[n - 1 - i] = B[i][i];
				}
			}
			else {
				for (int i = 0; i < n; i++) {
					lambda[i] = B[i][i];
				}
			}

			double norm = NormDiff(lambda, lambda_real, n);

			/* заполняет таблицу Excel - ';' для перехода к следующей ячейке */
			fprintf(f1, "%.10lf;%d;%.15lf\n", eps, Iter, norm);

			free(lambda);
			freeMatrix(B, n);
		}
		fprintf(f1, "\n");
	}
	
	

	//очищает память
	freeMatrix(A, n);
	free(matrix);
	free(matrix_h);
	fclose(fp);

	return 0;
}
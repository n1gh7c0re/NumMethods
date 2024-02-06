FileName3=fopen('file.txt','w'); %открывает файл для записи данных
format longEng %изменяет формат выходного отображения на longEng (15 цифр после запятой)
X = [1; 2; 3; 4; 5; 6; 7; 8; 9; 10]; %столбец решений
[Q, r]=qr(rand(10)); %генерирует случайную матрицу 10 на 10, разложенную на Q*R
QT = transpose(Q); %транспонирует матрицу Q
D = zeros(10); %создает матрицу 10 на 10, заполненную нулями
det_curr = 1; %определитель матрицы
for i=1:9 %считает определитель диагональной матрицы
   D(i, i) = rand(1); %генерирует случайное число
   det_curr = det_curr * D(i, i);
end
D(10,10) = 10 / det_curr; %получает определитель, равный 10
A = QT * D * Q; %получает симметричную матрицу 
b = A * X; %получает столбец свободных чисел
FullMatrix = horzcat(A,b); %получает расширенную матрицу
for j = 1:10
    for m = 1:11 %записывает в файл элементы матрицы
        fprintf(FileName3,'%.15f ', FullMatrix(j,m));
    end
end

%график зависимости времени выполнения метода t от размерности матрицы N
x1 = xlsread('time4.csv', 'B1:B8'); %считывает из таблицы значения для оси Х
y1 = xlsread('time4.csv', 'A1:A8'); %считывает из таблицы значения для оси Y

%крайние точки аппроксимирующей прямой для графика времени
x12 = x1(1);
x13 = x1(8);
y12 = y1(1)/1.3;
y13 = y1(8)/1.3;

figure
loglog(x1, y1, '-b','LineWidth', 2); %строит график
grid on; %строит сетку
hold on; %сохраняет график
loglog([x12 x13], [y12 y13], '--r','LineWidth', 2); %строит аппроксимирующую прямую
hold on; %сохраняет график
title('График зависимости времени выполнения метода t от размерности матрицы N') %заголовок
xlabel('Размерность матрицы N') %подпись оси х
ylabel('Время выполнения метода t, сек') %подпись оси у
%легенда
legend('Зависимость времени выполнения метода от размерности матрицы', ...
    'Аппроксимирующая прямая');

%график зависимости нормы погрешности от точности
x2 = xlsread('norm4.csv', 'C1:C10'); %считывает из таблицы значения для оси Х
y3 = xlsread('norm4.csv', 'B1:B10'); %считывает из таблицы значения для оси Y

figure
loglog(x2, y3, '-r','LineWidth', 2); %строит график для абсолютной погрешности
grid on; %строит сетку
hold on; %сохраняет график
title('График зависимости нормы абсолютной погрешности от точности решения'); %заголовок
xlabel('Точность решения'); %подпись оси x
ylabel('Норма абсолютной погрешности'); %подпись оси у

%график зависимости числа итераций от точности решения
y2 = xlsread('norm4.csv', 'A1:A10'); %считывает из таблицы значения для оси Y
figure
loglog(x2, y2, '-r','LineWidth', 2); %строит график для абсолютной погрешности
grid on; %строит сетку
hold on; %сохраняет график
title('График зависимости числа итераций от точности решения'); %заголовок
xlabel('Точность решения'); %подпись оси x
ylabel('Число итераций'); %подпись оси у
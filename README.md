#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
using namespace std;

// Точное решение: y(x) = (x(2+x)- 2(1+x)*ln(1+x)) /(x*(-1 + ln(4)))
double exactSolution(double x) {
    if (x == 0) return 0;
    return (x *(2+x)-2*(1 + x)*log(1 + x))/(x *(-1 + log(4)));
}

int main() {
    setlocale(LC_ALL, "Russian");

    int N=100;
    double ax=0;
    double bx =1;
    double h =(bx-ax) / N;
    double* x =new double[N + 3];
    double* alpha=new double[N + 3];
    double*beta =new double[N + 3];
    double*y =new double[N + 2];

    cout << "h=" << h << endl;
    cout << "N=" << N << endl;


    // Сетка
    for (int i = 1; i <=N + 1; i++) {
        x[i] = (i -1) * h;
    }

    // Начальные условия для прогоночных коэффициентов
    alpha[2]= 0;
    beta[2] =0;  //y(0) = 0

    // Прямая прогонка
    for (int i =2; i <= N; i++) {
        double A=x[i] *x[i] *(x[i] + 1);
        double C=2 * x[i] *x[i] * (x[i] +1) + 2 * h * h;
        double B=x[i]*x[i]*(x[i] + 1);
        double F =0;

        alpha[i + 1]= B/(C - A * alpha[i]);
        beta[i + 1]=(A *beta[i]+F) / (C -A * alpha[i]);
    }

    // Граничное условие: y'(1) = 1 => y(N+1) = y_N + h
    y[N + 1] = (beta[N + 1]+h) /(1-alpha[N + 1]);

    // Обратная прогонка
    for (int i = N; i >= 1; i--) {
        y[i] =alpha[i + 1] * y[i +1]+beta[i + 1];
    }
    ofstream file("C:\\Users\\rustamova.43011\\Desktop\\solution.txt");
    if (!file) {
        cout << "Не удалось создать файл" << endl;
        return 1;
    }

    file << "x\t\t y(числ решение)\t\t Точное решение\t\tПогрешность" << endl;

    // Точное решение вывод результатов 
    for (int i = 1; i <= N + 1; i++) {
        double t_resh=exactSolution(x[i]); //точное решение
        double pogr =fabs(y[i] - t_resh);
        file << x[i] << "\t\t" << y[i] << "\t\t\t" << t_resh << "\t\t\t" << pogr << endl;
        cout << "x=" << x[i] << " num=" << y[i] << " exact=" << t_resh << endl;
    }
    file.close();
    delete[] x;
    delete[] alpha;
    delete[] beta;
    delete[] y;

    return 0;
}


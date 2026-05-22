#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <stdexcept>

// ====== КЛАСС MATRIX ======
class Matrix {
private:
    std::vector<std::vector<double>> data;
    int rows, cols;

    // Создание единичной матрицы
    static Matrix identity(int n) {
        Matrix result(n, n);
        for (int i = 0; i < n; ++i)
            result.data[i][i] = 1.0;
        return result;
    }

public:
    // Конструкторы
    Matrix() : rows(0), cols(0) {}
    explicit Matrix(int r, int c) : rows(r), cols(c) {
        data.assign(rows, std::vector<double>(cols, 0.0));
    }
    explicit Matrix(int r, int c, double init) : rows(r), cols(c) {
        data.assign(rows, std::vector<double>(cols, init));
    }
    Matrix(const Matrix& other) = default;

    // Оператор присваивания
    Matrix& operator=(const Matrix& other) {
        if (this != &other) {
            rows = other.rows;
            cols = other.cols;
            data = other.data;
        }
        return *this;
    }

    // Доступ к элементам с проверкой границ
    double& operator()(int i, int j) {
        if (i < 0 || i >= rows || j < 0 || j >= cols)
            throw std::out_of_range("Индекс выходит за пределы матрицы");
        return data[i][j];
    }
    const double& operator()(int i, int j) const {
        if (i < 0 || i >= rows || j < 0 || j >= cols)
            throw std::out_of_range("Индекс выходит за пределы матрицы");
        return data[i][j];
    }

    int getRows() const { return rows; }
    int getCols() const { return cols; }

    // Арифметические операции
    Matrix operator+(const Matrix& other) const {
        if (rows != other.rows || cols != other.cols)
            throw std::invalid_argument("Сложение: размеры матриц не совпадают");
        Matrix result(rows, cols);
        for (int i = 0; i < rows; ++i)
            for (int j = 0; j < cols; ++j)
                result.data[i][j] = data[i][j] + other.data[i][j];
        return result;
    }

    Matrix operator-(const Matrix& other) const {
        if (rows != other.rows || cols != other.cols)
            throw std::invalid_argument("Вычитание: размеры матриц не совпадают");
        Matrix result(rows, cols);
        for (int i = 0; i < rows; ++i)
            for (int j = 0; j < cols; ++j)
                result.data[i][j] = data[i][j] - other.data[i][j];
        return result;
    }

    Matrix operator*(const Matrix& other) const {
        if (cols != other.rows)
            throw std::invalid_argument("Умножение: число столбцов A не равно числу строк B");
        Matrix result(rows, other.cols);
        for (int i = 0; i < rows; ++i)
            for (int j = 0; j < other.cols; ++j)
                for (int k = 0; k < cols; ++k)
                    result.data[i][j] += data[i][k] * other.data[k][j];
        return result;
    }

    Matrix operator*(double scalar) const {
        Matrix result(rows, cols);
        for (int i = 0; i < rows; ++i)
            for (int j = 0; j < cols; ++j)
                result.data[i][j] = data[i][j] * scalar;
        return result;
    }

    Matrix operator/(double scalar) const {
        if (scalar == 0) throw std::invalid_argument("Деление на ноль");
        Matrix result(rows, cols);
        for (int i = 0; i < rows; ++i)
            for (int j = 0; j < cols; ++j)
                result.data[i][j] = data[i][j] / scalar;
        return result;
    }

    // Комбинированные операции
    Matrix& operator+=(const Matrix& other) { *this = *this + other; return *this; }
    Matrix& operator-=(const Matrix& other) { *this = *this - other; return *this; }
    Matrix& operator*=(const Matrix& other) { *this = *this * other; return *this; }
    Matrix& operator*=(double scalar) { *this = *this * scalar; return *this; }
    Matrix& operator/=(double scalar) { *this = *this / scalar; return *this; }

    // Сравнение (прямое равенство)
    bool operator==(const Matrix& other) const {
        if (rows != other.rows || cols != other.cols) return false;
        for (int i = 0; i < rows; ++i)
            for (int j = 0; j < cols; ++j)
                if (data[i][j] != other.data[i][j])
                    return false;
        return true;
    }
    bool operator!=(const Matrix& other) const { return !(*this == other); }

    // Транспонирование
    Matrix transpose() const {
        Matrix result(cols, rows);
        for (int i = 0; i < rows; ++i)
            for (int j = 0; j < cols; ++j)
                result.data[j][i] = data[i][j];
        return result;
    }

    // Определитель (сравнение с нулём напрямую)
    double determinant() const {
        if (rows != cols) //Проверка квадратности. Определитель существует только у квадратных матриц.
            throw std::invalid_argument("Определитель только для квадратных матриц");
        Matrix temp = *this; //копия, чтобы не испортить оригинальную матрицу.
        double det = 1.0; //нач.знач
        for (int i = 0; i < rows; ++i) {
            int pivot = i; //Ищется строка с максимальным элементом.Чтобы избежать деления на маленькие числа и повысить устойчивость.
            for (int j = i + 1; j < rows; ++j)
                if (std::fabs(temp.data[j][i]) > std::fabs(temp.data[pivot][i]))
                    pivot = j;
            if (temp.data[pivot][i] == 0.0) //Проверка вырожденности. Если ведущий элемент 0 то определитель 0.
                return 0.0;
            if (pivot != i) {
                std::swap(temp.data[i], temp.data[pivot]); //перестановка строк
                det = -det; //перестановка строк меняет знак определителя.
            }
            det *= temp.data[i][i]; //умножение диагонали. После метода Гаусса determinant равен произведению диагональных элементов.
            for (int j = i + 1; j < rows; ++j) {
                double factor = temp.data[j][i] / temp.data[i][i]; //обнуление под диагональю(глав шаг), factor-на сколько нужно вычесть текущую строку.
                for (int k = i + 1; k < cols; ++k)
                    temp.data[j][k] -= factor * temp.data[i][k];
            }
        }
        return det;
    }

    // Обратная матрица
    Matrix inverse() const { //Создаётся расширенная матрица:
        if (rows != cols)
            throw std::invalid_argument("Обратная матрица только для квадратных матриц");
        double det = determinant();
        if (det == 0)
            throw std::invalid_argument("Матрица вырожденная, обратной не существует");

        int n = rows;
        Matrix aug(n, 2 * n);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j)
                aug.data[i][j] = data[i][j];
            aug.data[i][n + i] = 1.0;
        }

        for (int i = 0; i < n; ++i) {
            double pivot = aug.data[i][i];
            if (pivot == 0.0) {
                int swapRow = -1;
                for (int j = i + 1; j < n; ++j) {
                    if (aug.data[j][i] != 0.0) {
                        swapRow = j;
                        break;
                    }
                }
                if (swapRow == -1) throw std::invalid_argument("Обратная матрица не существует");
                std::swap(aug.data[i], aug.data[swapRow]);
                pivot = aug.data[i][i];
            }
            for (int j = 0; j < 2 * n; ++j)
                aug.data[i][j] /= pivot;
            for (int j = 0; j < n; ++j) {
                if (j == i) continue;
                double factor = aug.data[j][i];
                for (int k = 0; k < 2 * n; ++k)
                    aug.data[j][k] -= factor * aug.data[i][k];
            }
        }

        Matrix result(n, n);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                result.data[i][j] = aug.data[i][n + j];
        return result;
    }

    // Возведение в целую степень
    Matrix power(int exponent) const {
        if (rows != cols)
            throw std::invalid_argument("Степень только для квадратных матриц");
        if (exponent == 0) return identity(rows);
        if (exponent < 0) return inverse().power(-exponent);
        Matrix result = identity(rows);
        Matrix base = *this;
        while (exponent > 0) {
            if (exponent & 1) result = result * base;
            base = base * base;
            exponent >>= 1;
        }
        return result;
    }

    // Норма (максимум суммы модулей по столбцам)
    double norm() const {
        double maxSum = 0.0;
        for (int j = 0; j < cols; ++j) {
            double colSum = 0.0;
            for (int i = 0; i < rows; ++i)
                colSum += std::fabs(data[i][j]);
            if (colSum > maxSum) maxSum = colSum;
        }
        return maxSum;
    }

    // Проверки типов (прямое сравнение)
    bool isSquare() const { return rows == cols; }

    bool isZero() const {
        for (int i = 0; i < rows; ++i)
            for (int j = 0; j < cols; ++j)
                if (data[i][j] != 0.0) return false;
        return true;
    }

    bool isIdentity() const {
        if (!isSquare()) return false;
        for (int i = 0; i < rows; ++i)
            for (int j = 0; j < cols; ++j) {
                if (i == j && data[i][j] != 1.0) return false;
                if (i != j && data[i][j] != 0.0) return false;
            }
        return true;
    }

    bool isDiagonal() const {
        if (!isSquare()) return false;
        for (int i = 0; i < rows; ++i)
            for (int j = 0; j < cols; ++j)
                if (i != j && data[i][j] != 0.0) return false;
        return true;
    }

    bool isSymmetric() const {
        if (!isSquare()) return false;
        for (int i = 0; i < rows; ++i)
            for (int j = i + 1; j < cols; ++j)
                if (data[i][j] != data[j][i]) return false;
        return true;
    }

    bool isUpperTriangular() const {
        if (!isSquare()) return false;
        for (int i = 1; i < rows; ++i)
            for (int j = 0; j < i; ++j)
                if (data[i][j] != 0.0) return false;
        return true;
    }

    bool isLowerTriangular() const {
        if (!isSquare()) return false;
        for (int i = 0; i < rows; ++i)
            for (int j = i + 1; j < cols; ++j)
                if (data[i][j] != 0.0) return false;
        return true;
    }

    // Ввод / вывод
    void input() {
        std::cout << "  Введите элементы матрицы " << rows << "x" << cols << ":\n";
        for (int i = 0; i < rows; ++i)
            for (int j = 0; j < cols; ++j)
                std::cin >> data[i][j];
    }

    void print() const {
        for (int i = 0; i < rows; ++i) {
            std::cout << "  ";
            for (int j = 0; j < cols; ++j)
                std::cout << std::setw(10) << std::setprecision(4) << data[i][j] << " ";
            std::cout << std::endl;
        }
    }

    friend std::ostream& operator<<(std::ostream& os, const Matrix& m) {
        m.print();
        return os;
    }
    friend std::istream& operator>>(std::istream& is, Matrix& m) {
        for (int i = 0; i < m.rows; ++i)
            for (int j = 0; j < m.cols; ++j)
                is >> m.data[i][j];
        return is;
    }
};

// ====== ФУНКЦИИ ОФОРМЛЕНИЯ ======
void printHeader(const std::string& title) {
    std::cout << "\n  === " << title << " ===\n";
}

void printMenu() {
    std::cout << "\n  === ГЛАВНОЕ МЕНЮ ===\n";
    std::cout << "   1. Создать матрицу A\n";
    std::cout << "   2. Создать матрицу B\n";
    std::cout << "   3. A + B\n";
    std::cout << "   4. A - B\n";
    std::cout << "   5. A * B\n";
    std::cout << "   6. A * скаляр\n";
    std::cout << "   7. A / скаляр\n";
    std::cout << "   8. A += B\n";
    std::cout << "   9. A -= B\n";
    std::cout << "  10. A *= B\n";
    std::cout << "  11. A *= скаляр\n";
    std::cout << "  12. A /= скаляр\n";
    std::cout << "  13. A == B\n";
    std::cout << "  14. A != B\n";
    std::cout << "  15. Транспонирование A\n";
    std::cout << "  16. Определитель A\n";
    std::cout << "  17. Обратная матрица A\n";
    std::cout << "  18. A^N (степень)\n";
    std::cout << "  19. Норма матрицы A\n";
    std::cout << "  20. Проверка типов A\n";
    std::cout << "  21. Вывести A\n";
    std::cout << "  22. Вывести B\n";
    std::cout << "  23. Ввести A заново\n";
    std::cout << "  24. Ввести B заново\n";
    std::cout << "   0. Выход\n";
    std::cout << "  Ваш выбор: ";
}

void printTypeCheck(const Matrix& m, const std::string& name) {
    std::cout << "\n  === Проверка типов матрицы " << name << " ===\n";
    std::cout << "  Квадратная:          " << (m.isSquare() ? "Да" : "Нет") << "\n";
    if (m.isSquare()) {
        std::cout << "  Диагональная:        " << (m.isDiagonal() ? "Да" : "Нет") << "\n";
        std::cout << "  Единичная:           " << (m.isIdentity() ? "Да" : "Нет") << "\n";
        std::cout << "  Симметрическая:      " << (m.isSymmetric() ? "Да" : "Нет") << "\n";
        std::cout << "  Верхнетреугольная:   " << (m.isUpperTriangular() ? "Да" : "Нет") << "\n";
        std::cout << "  Нижнетреугольная:    " << (m.isLowerTriangular() ? "Да" : "Нет") << "\n";
    }
    std::cout << "  Нулевая:             " << (m.isZero() ? "Да" : "Нет") << "\n";
}

// ====== ГЛАВНАЯ ФУНКЦИЯ ======
int main() {
    setlocale(LC_ALL, "Russian");

    Matrix A, B;
    bool hasA = false, hasB = false;
    int choice, exponent;
    double scalar;

     std::cout << "\n  === ЛАБОРАТОРНАЯ РАБОТА №6 ===\n";
    std::cout << "  === Класс \"Вещественная матрица\" ===\n";

    do {
        printMenu();
        std::cin >> choice;

        switch (choice) {
        case 1: {
            int r, c;
            std::cout << "  Введите размеры матрицы A (строки столбцы): ";
            std::cin >> r >> c;
            A = Matrix(r, c);
            std::cout << "  Ввод матрицы A:\n";
            A.input();
            hasA = true;
            printHeader("Матрица A");
            std::cout << A;
        } break;

        case 2: {
            int r, c;
            std::cout << "  Введите размеры матрицы B (строки столбцы): ";
            std::cin >> r >> c;
            B = Matrix(r, c);
            std::cout << "  Ввод матрицы B:\n";
            B.input();
            hasB = true;
            printHeader("Матрица B");
            std::cout << B;
        } break;

        case 3:
            if (hasA && hasB) {
                printHeader("Сложение матриц A + B");
                try { std::cout << (A + B) << std::endl; }
                catch (const std::exception& e) { std::cout << "  Ошибка: " << e.what() << std::endl; }
            }
            else std::cout << "  Сначала создайте обе матрицы!\n";
            break;

        case 4:
            if (hasA && hasB) {
                printHeader("Вычитание матриц A - B");
                try { std::cout << (A - B) << std::endl; }
                catch (const std::exception& e) { std::cout << "  Ошибка: " << e.what() << std::endl; }
            }
            else std::cout << "  Сначала создайте обе матрицы!\n";
            break;

        case 5:
            if (hasA && hasB) {
                printHeader("Умножение матриц A * B");
                try { std::cout << (A * B) << std::endl; }
                catch (const std::exception& e) { std::cout << "  Ошибка: " << e.what() << std::endl; }
            }
            else std::cout << "  Сначала создайте обе матрицы!\n";
            break;

        case 6:
            if (hasA) {
                printHeader("Умножение А на скаляр");
                std::cout << "  Введите скаляр: ";
                std::cin >> scalar;
                std::cout << (A * scalar) << std::endl;
            }
            else std::cout << "  Сначала создайте матрицу A!\n";
            break;

        case 7:
            if (hasA) {
                printHeader("Деление A на скаляр");
                std::cout << "  Введите скаляр: ";
                std::cin >> scalar;
                try { std::cout << (A / scalar) << std::endl; }
                catch (const std::exception& e) { std::cout << "  Ошибка: " << e.what() << std::endl; }
            }
            else std::cout << "  Сначала создайте матрицу A!\n";
            break;

        case 8:
            if (hasA && hasB) {
                printHeader("Операция A += B");
                try { A += B; std::cout << "  Результат:\n" << A; }
                catch (const std::exception& e) { std::cout << "  Ошибка: " << e.what() << std::endl; }
            }
            else std::cout << "  Сначала создайте обе матрицы!\n";
            break;

        case 9:
            if (hasA && hasB) {
                printHeader("Операция A -= B");
                try { A -= B; std::cout << "  Результат:\n" << A; }
                catch (const std::exception& e) { std::cout << "  Ошибка: " << e.what() << std::endl; }
            }
            else std::cout << "  Сначала создайте обе матрицы!\n";
            break;

        case 10:
            if (hasA && hasB) {
                printHeader("Операция A *= B");
                try { A *= B; std::cout << "  Результат:\n" << A; }
                catch (const std::exception& e) { std::cout << "  Ошибка: " << e.what() << std::endl; }
            }
            else std::cout << "  Сначала создайте обе матрицы!\n";
            break;

        case 11:
            if (hasA) {
                printHeader("Операция A *= скаляр");
                std::cout << "  Введите скаляр: ";
                std::cin >> scalar;
                A *= scalar;
                std::cout << "  Результат:\n" << A;
            }
            else std::cout << "  Сначала создайте матрицу A!\n";
            break;

        case 12:
            if (hasA) {
                printHeader("Операция A /= скаляр");
                std::cout << "  Введите скаляр: ";
                std::cin >> scalar;
                try { A /= scalar; std::cout << "  Результат:\n" << A; }
                catch (const std::exception& e) { std::cout << "  Ошибка: " << e.what() << std::endl; }
            }
            else std::cout << "  Сначала создайте матрицу A!\n";
            break;

        case 13:
            if (hasA && hasB) {
                printHeader("Сравнение A == B");
                std::cout << "  Результат: " << (A == B ? "true (матрицы равны)" : "false (матрицы не равны)") << std::endl;
            }
            else std::cout << "  Сначала создайте обе матрицы!\n";
            break;

        case 14:
            if (hasA && hasB) {
                printHeader("Сравнение A != B");
                std::cout << "  Результат: " << (A != B ? "true (матрицы не равны)" : "false (матрицы равны)") << std::endl;
            }
            else std::cout << "  Сначала создайте обе матрицы!\n";
            break;

        case 15:
            if (hasA) {
                printHeader("Транспонирование матрицы A");
                std::cout << A.transpose() << std::endl;
            }
            else std::cout << "  Сначала создайте матрицу A!\n";
            break;

        case 16:
            if (hasA) {
                printHeader("Определитель матрицы A");
                try { std::cout << "  det(A) = " << std::setprecision(6) << A.determinant() << std::endl; }
                catch (const std::exception& e) { std::cout << "  Ошибка: " << e.what() << std::endl; }
            }
            else std::cout << "  Сначала создайте матрицу A!\n";
            break;

        case 17:
            if (hasA) {
                printHeader("Обратная матрица A⁻¹");
                try { std::cout << A.inverse() << std::endl; }
                catch (const std::exception& e) { std::cout << "  Ошибка: " << e.what() << std::endl; }
            }
            else std::cout << "  Сначала создайте матрицу A!\n";
            break;

        case 18:
            if (hasA) {
                printHeader("Возведение А в степень");
                std::cout << "  Введите степень: ";
                std::cin >> exponent;
                try { std::cout << (A.power(exponent)) << std::endl; }
                catch (const std::exception& e) { std::cout << "  Ошибка: " << e.what() << std::endl; }
            }
            else std::cout << "  Сначала создайте матрицу A!\n";
            break;

        case 19:
            if (hasA) {
                printHeader("Норма матрицы A");
                std::cout << "  ||A|| = " << std::setprecision(6) << A.norm() << std::endl;
            }
            else std::cout << "  Сначала создайте матрицу A!\n";
            break;

        case 20:
            if (hasA) printTypeCheck(A, "A");
            else std::cout << "  Сначала создайте матрицу A!\n";
            break;

        case 21:
            if (hasA) {
                printHeader("Матрица A");
                std::cout << A;
            }
            else std::cout << "  Матрица A не создана!\n";
            break;

        case 22:
            if (hasB) {
                printHeader("Матрица B");
                std::cout << B;
            }
            else std::cout << "  Матрица B не создана!\n";
            break;

        case 23:
            if (hasA) {
                int r, c;
                std::cout << "  Введите размеры матрицы A (строки столбцы): ";
                std::cin >> r >> c;
                A = Matrix(r, c);
                std::cout << "  Ввод матрицы A:\n";
                A.input();
                std::cout << "  Матрица A обновлена:\n" << A;
            }
            else std::cout << "  Сначала создайте матрицу A!\n";
            break;

        case 24:
            if (hasB) {
                int r, c;
                std::cout << "  Введите размеры матрицы B (строки столбцы): ";
                std::cin >> r >> c;
                B = Matrix(r, c);
                std::cout << "  Ввод матрицы B:\n";
                B.input();
                std::cout << "  Матрица B обновлена:\n" << B;
            }
            else std::cout << "  Сначала создайте матрицу B!\n";
            break;

        case 0:
            std::cout << "\n  === Работа завершена. До свидания! ===\n";
            break;

        default:
            std::cout << "  Неверный выбор! Попробуйте снова.\n";
        }
    } while (choice != 0);

    return 0;
}

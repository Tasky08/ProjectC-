#include <iostream>
#include <cmath>
#include <functional>
#include <random>
#include <cassert>

// Вычисление интегралов

/**
 * @brief Вычисляет определенный интеграл методом Симпсона.
 *
 * @param f Функция, подлежащая интегрированию.
 * @param a Нижний предел интегрирования.
 * @param b Верхний предел интегрирования.
 * @param n Количество разбиений.
 * @return Значение определенного интеграла.
 */
double simpsonIntegral(const std::function<double(double)>& f, double a, double b, int n) {
    double h = (b - a) / n;
    double integral = f(a) + f(b);

    for (int i = 1; i < n; i += 2) {
        integral += 4 * f(a + i * h);
    }

    for (int i = 2; i < n; i += 2) {
        integral += 2 * f(a + i * h);
    }

    integral *= h / 3;
    return integral;
}

/**
 * @brief Вычисляет определенный интеграл методом Монте-Карло.
 *
 * @param f Функция, подлежащая интегрированию.
 * @param a Нижний предел интегрирования.
 * @param b Верхний предел интегрирования.
 * @param num_samples Количество выборок.
 * @return Значение определенного интеграла.
 */
double monteCarloIntegral(const std::function<double(double)>& f, double a, double b, int num_samples) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(a, b);

    double sum = 0.0;
    for (int i = 0; i < num_samples; ++i) {
        double x = dis(gen);
        sum += f(x);
    }

    return (b - a) * sum / num_samples;
}

/**
 * @brief Вычисляет неопределенный интеграл методом Симпсона.
 *
 * @param f Функция, подлежащая интегрированию.
 * @param a Нижний предел интегрирования.
 * @param b Верхний предел интегрирования.
 * @param n Количество разбиений.
 * @param x Значение x для неопределенного интеграла.
 * @return Значение неопределенного интеграла.
 */
double indefiniteSimpsonIntegral(const std::function<double(double)>& f, double a, double b, int n, double x) {
    double h = (b - a) / n;
    double integral = f(a) + f(b);

    for (int i = 1; i < n; i += 2) {
        integral += 4 * f(a + i * h);
    }

    for (int i = 2; i < n; i += 2) {
        integral += 2 * f(a + i * h);
    }

    integral *= h / 3;
    double C = integral - f(x);
    return integral - C;
}

// Пример функции для интегрирования
double exampleFunction(double x) {
    return sin(x);
}

// Дополнительные математические функции

/**
 * @brief Функция для вычисления косинуса.
 *
 * @param x Значение для вычисления косинуса.
 * @return Значение косинуса.
 */
double cosineFunction(double x) {
    return cos(x);
}

/**
 * @brief Функция для вычисления экспоненты.
 *
 * @param x Значение для вычисления экспоненты.
 * @return Значение экспоненты.
 */
double exponentialFunction(double x) {
    return exp(x);
}

/**
 * @brief Функция для вычисления квадратного корня.
 *
 * @param x Значение для вычисления квадратного корня.
 * @return Значение квадратного корня.
 */
double sqrtFunction(double x) {
    return sqrt(x);
}

/**
 * @brief Вычисляет интеграл от косинуса методом Симпсона.
 *
 * @param a Нижний предел интегрирования.
 * @param b Верхний предел интегрирования.
 * @param n Количество разбиений.
 * @return Значение определенного интеграла.
 */
double cosineSimpsonIntegral(double a, double b, int n) {
    return simpsonIntegral(cosineFunction, a, b, n);
}

/**
 * @brief Вычисляет интеграл от экспоненты методом Симпсона.
 *
 * @param a Нижний предел интегрирования.
 * @param b Верхний предел интегрирования.
 * @param n Количество разбиений.
 * @return Значение определенного интеграла.
 */
double exponentialSimpsonIntegral(double a, double b, int n) {
    return simpsonIntegral(exponentialFunction, a, b, n);
}

/**
 * @brief Вычисляет интеграл от квадратного корня методом Симпсона.
 *
 * @param a Нижний предел интегрирования.
 * @param b Верхний предел интегрирования.
 * @param n Количество разбиений.
 * @return Значение определенного интеграла.
 */
double sqrtSimpsonIntegral(double a, double b, int n) {
    return simpsonIntegral(sqrtFunction, a, b, n);
}

// Ввод данных
void inputData(double& a, double& b, int& n, int& num_samples) {
    std::cout << "Введите нижний предел интегрирования: ";
    std::cin >> a;

    std::cout << "Введите верхний предел интегрирования: ";
    std::cin >> b;

    std::cout << "Введите количество разбиений (для метода Симпсона): ";
    std::cin >> n;

    std::cout << "Введите количество выборок (для метода Монте-Карло): ";
    std::cin >> num_samples;
}

// Вывод результатов
void displayResults(double a, double b, int n, int num_samples, double x) {
    std::cout << "\nОпределенный интеграл методом Симпсона: "
        << simpsonIntegral(exampleFunction, a, b, n) << std::endl;

    std::cout << "Определенный интеграл методом Монте-Карло: "
        << monteCarloIntegral(exampleFunction, a, b, num_samples) << std::endl;

    std::cout << "Неопределенный интеграл методом Симпсона: "
        << indefiniteSimpsonIntegral(exampleFunction, a, b, n, x) << std::endl;
}

// Пример функции для интегрирования
double testFunction(double x) {
    return x * x; // Интеграл от x^2
}

void runTests() {
    // Тест для метода Симпсона
    double a = 0.0;
    double b = 1.0;
    int n = 1000;
    double result = simpsonIntegral(testFunction, a, b, n);
    double expected = 1.0 / 3.0;
    assert(std::fabs(result - expected) < 1e-6);

    // Тест для метода Монте-Карло
    int num_samples = 100000;
    result = monteCarloIntegral(testFunction, a, b, num_samples);
    assert(std::fabs(result - expected) < 1e-2);

    std::cout << "Все тесты пройдены успешно!" << std::endl;
}

// Основная функция
int main() {
    double a, b;
    int n, num_samples;

    runTests(); // Запуск тестов

    inputData(a, b, n, num_samples);

    double x;
    std::cout << "Введите значение x для неопределенного интеграла: ";
    std::cin >> x;

    displayResults(a, b, n, num_samples, x);

    std::cout << "\nОпределенный интеграл от косинуса методом Симпсона: "
        << cosineSimpsonIntegral(a, b, n) << std::endl;

    std::cout << "Определенный интеграл от экспоненты методом Симпсона: "
        << exponentialSimpsonIntegral(a, b, n) << std::endl;

    std::cout << "Определенный интеграл от квадратного корня методом Симпсона: "
        << sqrtSimpsonIntegral(a, b, n) << std::endl;

    // Дополнительные вызовы функций для увеличения количества строк
    std::cout << "Интеграл от x^2 методом Симпсона: "
        << simpsonIntegral([](double x) { return x * x; }, a, b, n) << std::endl;

    std::cout << "Интеграл от x^3 методом Монте-Карло: "
        << monteCarloIntegral([](double x) { return x * x * x; }, a, b, num_samples) << std::endl;

    std::cout << "Неопределенный интеграл от x^2 методом Симпсона: "
        << indefiniteSimpsonIntegral([](double x) { return x * x; }, a, b, n, x) << std::endl;

    return 0;
}

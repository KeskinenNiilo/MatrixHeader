#ifndef MATRIXNIILOK_H
#define MATRIXNIILOK_H

#include <vector>
#include <iostream>
#include <stdexcept>
#include <random>
#include <cmath>
#include <iomanip>

class Matrix {
public:
    std::vector<std::vector<double>> matrix;
    Matrix(std::initializer_list<std::vector<double>> init) : matrix(init) {
        if (!matrix.empty()) {
            size_t width = matrix[0].size();
            for (const auto& row : matrix) {
                if (row.size() != width) throw std::runtime_error("Inconsistent row sizes");
            }
        }
    }
    Matrix() = default;

    // -------------------------------| NEW |-----------------------------------------
    static Matrix random(int rows = 4, int cols = 4, double low = 0, double high = 1.0f) {
        Matrix result;
        static std::mt19937 gen(std::random_device{}());
        std::uniform_real_distribution<double> dist(low, high);
        for (int i = 0; i < rows; ++i) {
            std::vector<double> temp(cols);
            for (int j = 0; j < cols; ++j) temp[j] = dist(gen);
            result.push(temp);
        }
        return result;
    }

    static Matrix full(double value =0, int rows =4, int cols =4) {
        Matrix result;
        for (int i = 0; i < rows; ++i) {
            std::vector<double> temp(cols, value);
            result.push(temp);
        }
        return result;
    }

    // -------------------------------| UTIL |----------------------------------------
    void print(int precision = 4) const {
        for (const auto& row : matrix) {
            for (double val : row)
                std::cout << std::fixed << std::setprecision(precision) << val << "\t";
            std::cout << "\n";
        }
    }

    Matrix apply(double (*func)(double)) const {
        Matrix result = *this;
        for (auto& row : result.matrix) {
            for (auto& val : row) {
                val = func(val);
            }
        }
        return result;
    }

    Matrix transpose() const {
        if (matrix.empty()) return Matrix();
        int rows = matrix.size();
        int cols = matrix[0].size();
        Matrix result;
        result.matrix.resize(cols, std::vector<double>(rows));
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                result.matrix[j][i] = matrix[i][j];
            }
        }
        return result;
    }

    double sum() const {
        if (matrix.empty()) return 0;
        double sum = 0;
        for (auto& row : matrix) {
            for (double val : row) {
                sum += val;
            }
        }
        return sum;
    }

    double mean() const {
        if (matrix.empty() || matrix[0].empty()) throw std::runtime_error("Matrix empty");
        int s = sum();
        return s / (matrix.size() * matrix[0].size());
    }

    std::vector<double> flatten() const {
        std::vector<double> result;
        for (auto& row : matrix) {
            result.insert(result.end(), row.begin(), row.end());
        }
        return result;
    }

    std::pair<int, int> shape() const {
        if (matrix.empty()) return {0,0};
        return {matrix.size(), matrix[0].size()};
    }

    Matrix reshape(int rows, int cols) const {
        std::vector<double> flat = flatten();
        if (rows * cols != flat.size()) throw std::runtime_error("Size mismatch");
        Matrix result;
        result.matrix.resize(rows, std::vector<double>(cols));
        int index = 0;
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                result.matrix[i][j] = flat[index++];
            }
        }
        return result;
    }

    Matrix getRow(int index) const {
        if (index < 0 || index >= matrix.size()) throw std::runtime_error("Row index out of range");
        Matrix result;
        result.push(matrix[index]);
        return result;
    }

    Matrix getCol(int index) const {
        if (matrix.empty() || index < 0 || index >= matrix[0].size()) throw std::runtime_error("Column index out of range");
        Matrix result;
        for (const auto& row : matrix) {
            result.push({row[index]});
        }
        return result;
    }

    double max() const {
        if (matrix.empty()) throw std::runtime_error("Matrix empty");
        double mx = matrix[0][0];
        for (const auto& row : matrix) {
            for (double val : row) {
                if (val > mx) mx = val;
            }
        }
        return mx;
    }

    double min() const {
        if (matrix.empty()) throw std::runtime_error("Matrix empty");
        double mn = matrix[0][0];
        for (const auto& row : matrix) {
            for (double val : row) {
                if (val < mn) mn = val;
            }
        }
        return mn;
    }

    // ------------------------------| INSERT |---------------------------------------
    void insert(const std::vector<double> value, int index) {
        if(!matrix.empty() && value.size() != matrix[0].size()) throw std::runtime_error("Row size mismatch");
        if(index < 0 || index > matrix.size()) throw std::runtime_error("Index out of matrix range");
        matrix.insert(matrix.begin() + index, value);
    }

    void push(const std::vector<double>& r) {
        if (matrix.empty() || r.size() == matrix[0].size()) matrix.push_back(r);
        else throw std::runtime_error("Row size mismatch");
    }

    // --------------------------------| * |------------------------------------------
    Matrix dot(const Matrix& other) const {
        if(matrix.empty() || other.matrix.empty()) throw std::runtime_error("One of the matrices is empty");
        int aRows = matrix.size();
        int aCols = matrix[0].size();
        int bRows = other.matrix.size();
        int bCols = other.matrix[0].size();
        if (aCols != bRows) throw std::runtime_error("Matrix size mismatch");
        Matrix result;
        result.matrix.resize(aRows, std::vector<double>(bCols, 0));
        for (int i = 0; i < aRows; ++i) {
            for (int j = 0; j < bCols; ++j) {
                for (int k = 0; k < aCols; ++k) {
                    result.matrix[i][j] += matrix[i][k] * other.matrix[k][j];
                }
            }
        }
        return result;
    }

    Matrix operator*(const Matrix& other) const {
        return this->dot(other);
    }
    
    Matrix multiply(const Matrix& other) const {
        if (matrix.empty() || other.matrix.empty()) throw std::runtime_error("One of the matrices is empty");
        if (matrix.size() != other.matrix.size() || matrix[0].size() != other.matrix[0].size()) throw std::runtime_error("Matrix size mismatch");
        Matrix result = *this;
        int rows = matrix.size();
        int cols = matrix[0].size();
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                result.matrix[i][j] *= other.matrix[i][j];
            }
        }
        return result;
    }

    Matrix& operator*=(const Matrix& other) {
        if (matrix.empty() || other.matrix.empty()) throw std::runtime_error("One of the matrices is empty");
        if (matrix.size() != other.matrix.size() || matrix[0].size() != other.matrix[0].size()) throw std::runtime_error("Matrix size mismatch");
        for (int i = 0; i < matrix.size(); ++i) {
            for (int j = 0; j < matrix[0].size(); ++j) {
                matrix[i][j] *= other.matrix[i][j];
            }
        }
        return *this;
    }

    Matrix operator*(double scalar) const {
        Matrix result = *this;
        for (auto& row : result.matrix) {
            for (auto& val : row) {
                val *= scalar;
            }
        }
        return result;
    }
    
    // --------------------------------| / |------------------------------------------
    Matrix operator/(double scalar) const {
        if (scalar == 0) throw std::runtime_error("Division by zero");
        Matrix result = *this;
        for (auto& row : result.matrix) {
            for (auto& val : row) {
                val /= scalar;
            }
        }
        return result;
    }

    // --------------------------------| + |------------------------------------------
    void add(const Matrix& other) {
        if(matrix.empty() || other.matrix.empty()) throw std::runtime_error("One of the matrices is empty");
        if (matrix.size() != other.matrix.size() || matrix[0].size() != other.matrix[0].size()) throw std::runtime_error("Matrix size mismatch");
        int rows = matrix.size();
        int cols = matrix[0].size();
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                matrix[i][j] += other.matrix[i][j];
            }
        }
    }

    Matrix& operator+=(const Matrix& other) {
        add(other);
        return *this;
    }

    Matrix operator+(const Matrix& other) const {
        Matrix result = *this;
        result.add(other);
        return result;
    }

    Matrix operator+(double scalar) const {
        Matrix result = *this;
        for (auto& row : result.matrix) {
            for (auto& val : row) {
                val += scalar;
            }
        }
        return result;
    }

    // --------------------------------| - |------------------------------------------
    void subtract(const Matrix& other) {
        if(matrix.empty() || other.matrix.empty()) throw std::runtime_error("One of the matrices is empty");
        if (matrix.size() != other.matrix.size() || matrix[0].size() != other.matrix[0].size()) throw std::runtime_error("Matrix size mismatch");
        int rows = matrix.size();
        int cols = matrix[0].size();
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                matrix[i][j] -= other.matrix[i][j];
            }
        }
    }

    Matrix& operator-=(const Matrix& other) {
        subtract(other);
        return *this;
    }

    Matrix operator-(const Matrix& other) const {
        Matrix result = *this;
        result.subtract(other);
        return result;
    }

    Matrix operator-(double scalar) const {
        Matrix result = *this;
        for (auto& row : result.matrix) {
            for (auto& val : row) {
                val -= scalar;
            }
        }
        return result;
    }
};

namespace MatrixActivations {
    double sigmoid(double x) {
        return 1.0 / (1.0 + std::exp(-x));
    }
    double sigmoidDerivative(double x) {
        return sigmoid(x) * (1 - sigmoid(x));
    }
    double relu(double x) {
        return x > 0 ? x : 0;
    }
    double reluDerivative(double x) {
        return x > 0 ? 1 : 0;
    }
    double tanh(double x) {
        return std::tanh(x);
    }
    double tanhDerivative(double x) {
        return 1 - std::tanh(x) * std::tanh(x);
    }
}

#endif

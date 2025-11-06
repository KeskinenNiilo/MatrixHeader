#ifndef MATRIXNIILOK_H
#define MATRIXNIILOK_H

#include <vector>
#include <iostream>
#include <stdexcept>
#include <random>

class Matrix {
public:
    std::vector<std::vector<float>> matrix;
    Matrix(std::initializer_list<std::vector<float>> init) : matrix(init) {
        if (!matrix.empty()) {
            size_t width = matrix[0].size();
            for (const auto& row : matrix) {
                if (row.size() != width) throw std::runtime_error("Inconsistent row sizes");
            }
        }
    }
    Matrix() = default;

    void insert(std::vector<float> value, int index) {
        if(!matrix.empty() && value.size() != matrix[0].size()) throw std::runtime_error("Row size mismatch");
        if(index < 0 || index > static_cast<int>(matrix.size())) throw std::runtime_error("Index out of matrix scope");
        matrix.insert(matrix.begin() + index, value);
    }

    void print() const {
        for (const auto& row : matrix) {
            for (float val : row) std::cout << val << " ";
            std::cout << "\n";
        }
    }

    void push(const std::vector<float>& r) {
        if (matrix.empty() || r.size() == matrix[0].size()) matrix.push_back(r);
        else throw std::runtime_error("Row size mismatch");
    }

    Matrix dot(const Matrix& other) const {
        if(matrix.empty() || other.matrix.empty()) throw std::runtime_error("One of the matrices is empty");
        int aRows = matrix.size();
        int aCols = matrix[0].size();
        int bRows = other.matrix.size();
        int bCols = other.matrix[0].size();
        if (aCols != bRows) throw std::runtime_error("Matrix size mismatch");
        Matrix result;
        result.matrix.resize(aRows, std::vector<float>(bCols, 0));
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
    Matrix operator-(const Matrix& other) const {
        Matrix result = *this;
        result.subtract(other);
        return result;
    }

    static Matrix random(int rows = 4, int cols = 4, float low = 0, float high = 1.0f) {
        Matrix result;
        static std::mt19937 gen(std::random_device{}());
        std::uniform_real_distribution<float> dist(low, high);
        for (int i = 0; i < rows; ++i) {
            std::vector<float> temp(cols);
            for (int j = 0; j < cols; ++j) temp[j] = dist(gen);
            result.push(temp);
        }
        return result;
    }


    static Matrix full(float value =0, int rows =4, int cols =4) {
        Matrix result;
        for (int i = 0; i < rows; ++i) {
            std::vector<float> temp(cols, value);
            result.push(temp);
        }
        return result;
    }

    Matrix apply(float (*func)(float)) const {
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
        result.matrix.resize(cols, std::vector<float>(rows));
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                result.matrix[j][i] = matrix[i][j];
            }
        }
        return result;
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

    Matrix operator*(float scalar) const {
        Matrix result = *this;
        for (auto& row : result.matrix) {
            for (auto& val : row) {
                val *= scalar;
            }
        }
        return result;
    }

    Matrix operator/(float scalar) const {
        if (scalar == 0) throw std::runtime_error("Division by zero");
        Matrix result = *this;
        for (auto& row : result.matrix) {
            for (auto& val : row) {
                val /= scalar;
            }
        }
        return result;
    }

    Matrix operator+(float scalar) const {
        Matrix result = *this;
        for (auto& row : result.matrix) {
            for (auto& val : row) {
                val += scalar;
            }
        }
        return result;
    }

    Matrix operator-(float scalar) const {
        Matrix result = *this;
        for (auto& row : result.matrix) {
            for (auto& val : row) {
                val -= scalar;
            }
        }
        return result;
    }

    float sum() const {
        if (matrix.empty()) return 0;
        float sum = 0;
        for (auto& row : matrix) {
            for (float val : row) {
                sum += val;
            }
        }
        return sum;
    }

    float mean() const {
        if (matrix.empty() || matrix[0].empty()) throw std::runtime_error("Matrix empty");
        return sum() / (matrix.size() * matrix[0].size());
    }

    std::vector<float> flatten() const {
        std::vector<float> result;
        for (auto& row : matrix) {
            result.insert(result.end(), row.begin(), row.end());
        }
        return result;
    }

};

#endif

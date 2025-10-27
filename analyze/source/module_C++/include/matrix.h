#include <vector>
#include <stdexcept>

namespace itp{
template <typename T>
class Matrix {
public:
    using Type = std::vector<std::vector<T>>;

private:
    Type data;

public:

    Matrix() = default;

    // 构造指定维度并初始化为 value
    Matrix(size_t dim1, size_t dim2, const T& value = T{})
        : data(dim1, std::vector<T>(dim2, value)) {}

    // 从已有二维矩阵复制结构，用新值初始化
    Matrix(const Type& ori_matrix, const T& value = T{}) {
        if (!ori_matrix.empty()) {
            size_t dim1 = ori_matrix.size();
            size_t dim2 = ori_matrix[0].size();
            data = Type(dim1, std::vector<T>(dim2, value));
        }
        // 如果 ori_matrix 为空，则 data 保持为空
    }

    // 提供对内部数据的访问
    const Type& get() const { return data; }
    Type& get() { return data; }

    // 可选：重载括号操作符实现访问（推荐用于性能敏感场景）
    T& operator()(size_t i, size_t j) {
        if (i >= data.size() || j >= (data[i].size())) {
            throw std::out_of_range("Matrix index out of range");
        }
        return data[i][j];
    }

    const T& operator()(size_t i, size_t j) const {
        if (i >= data.size() || j >= (data[i].size())) {
            throw std::out_of_range("Matrix index out of range");
        }
        return data[i][j];
    }

    // 重载下标操作符实现行访问
    std::vector<T>& operator[](size_t i) {
        if (i >= data.size()) {
            throw std::out_of_range("Matrix index out of range");
        }
        return data[i];
    }

    const std::vector<T>& operator[](size_t i) const {
        if (i >= data.size()) {
            throw std::out_of_range("Matrix index out of range");
        }
        return data[i];
    }

    // 正确的 + 运算符重载（成员函数）
    Matrix<T> operator+(const Matrix<T>& other) const {
        if (this->rows() != other.rows() || this->cols() != other.cols()) {
            throw std::invalid_argument("Matrix addition dimension mismatch");
        }
        
        Matrix<T> result(this->rows(), this->cols());
        for (size_t i = 0; i < this->rows(); i++) {
            for (size_t j = 0; j < this->cols(); j++) {
                result(i, j) = (*this)(i, j) + other(i, j);
            }
        }
        return result;
    }   

    Matrix<T> operator+=(const Matrix<T>& other) {
        if (this->rows() != other.rows() || this->cols() != other.cols()) {
            throw std::invalid_argument("Matrix addition dimension mismatch");
        }
        
        for (size_t i = 0; i < this->rows(); i++) {
            for (size_t j = 0; j < this->cols(); j++) {
                (*this)(i, j) += other(i, j);
            }
        }
        return (*this);
    }   

    // 获取尺寸
    size_t rows() const { return data.size(); }
    size_t cols() const { return data.empty() ? 0 : data[0].size(); }
    size_t size() const { return rows(); }
};




template <typename T>
class Matrix_3d {
public:
    using Type = std::vector<Matrix<T>>;

private:
    Type data;

public:
    // 构造指定维度并初始化为 value
    Matrix_3d(size_t dim1, size_t dim2, size_t dim3, const T& value = T{})
        : data(dim1, Matrix(dim2, dim3, value)) {}

    // 从已有三维矩阵复制结构，用新值初始化
    Matrix_3d(const Type& ori_matrix, const T& value = T{}) {
        if (!ori_matrix.empty() && !ori_matrix[0].empty() && !ori_matrix[0][0].empty()) {
            size_t dim1 = ori_matrix.size();
            size_t dim2 = ori_matrix[0].size();
            size_t dim3 = ori_matrix[0][0].size();
            data = Type(dim1, Matrix(dim2, dim3, value));
        }
        // 处理空矩阵情况
    }

    // 访问内部数据
    const Type& get() const { return data; }
    Type& get() { return data; }

    // 重载括号操作符，实现三维访问
    T& operator()(size_t i, size_t j, size_t k) {
        if (i >= data.size() || j >= data[i].size() || k >= data[i][j].size()) {
            throw std::out_of_range("Matrix3D index out of range");
        }
        return data[i][j][k];
    }

    const T& operator()(size_t i, size_t j, size_t k) const {
        if (i >= data.size() || j >= data[i].size() || k >= data[i][j].size()) {
            throw std::out_of_range("Matrix3D index out of range");
        }
        return data[i][j][k];
    }

    Matrix<T>& operator[](size_t i) {
        if (i >= data.size()) {
            throw std::out_of_range("Matrix3D index out of range");
        }
        return data[i];
    }

    const Matrix<T>& operator[](size_t i) const {
        if (i >= data.size()) {
            throw std::out_of_range("Matrix3D index out of range");
        }
        return data[i];
    }

    std::vector<T>& operator()(size_t i, size_t j) {
        if (i >= data.size() || j >= data[i].size()) {
            throw std::out_of_range("Matrix3D index out of range");
        }
        return data[i][j];
    }

    const std::vector<T>& operator()(size_t i, size_t j) const {
        if (i >= data.size() || j >= data[i].size()) {
            throw std::out_of_range("Matrix3D index out of range");
        }
        return data[i][j];
    }

    // 获取维度
    size_t dim1() const { return data.size(); }
    size_t dim2() const { return data.empty() ? 0 : data[0].size(); }
    size_t dim3() const { 
        return (data.empty() || data[0].empty()) ? 0 : data[0][0].size(); 
    }

    size_t size() const { return dim1(); }
};
}//namespace itp

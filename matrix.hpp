#pragma once

#include<iostream>
#include<type_traits>
#include<vector>

namespace m_calc {
    template<typename T>
    class Matrix {
        public:
            static_assert(std::is_arithmetic<T>::value, "type T must be arithmetic!");
            using Mat = typename std::vector<std::vector<T>>;

        private:
            Mat data;

        public:
            /*member funcions*/
            bool operator==(const Matrix &m) const;
            inline bool operator!=(const Matrix &m) const;
            Matrix &operator=(const Matrix &m);
            Matrix &operator+=(const Matrix &m);
            Matrix &operator-=(const Matrix &m);
            Matrix &operator*=(const T c);
            Matrix &operator*=(const Matrix &m);

            inline const Mat get_data() const;
            inline const T get_element(int i, int j) const;
            inline const size_t row_num() const;
            inline const size_t column_num() const;
            inline const bool is_square() const;
            const Matrix transpose() const;
            const Matrix echelon() const;
            const Matrix inverse() const;
            const int rank() const;
            const T determinant() const;
            void print() const;
            static Matrix identity(size_t size);
            static Matrix elementary(size_t size, int i, T c);
            static Matrix elementary(size_t size, int i, int j, T c);
            static Matrix elementary_switch(size_t size, int i, int j);

            /*constrctors*/
            Matrix(const Mat v);
            Matrix(const size_t size);
            Matrix(const size_t row, const size_t col);
    };
    template<typename T> const Matrix<T> operator+(const Matrix<T> &m1, const Matrix<T> &m2);
    template <typename T> const Matrix<T> operator-(const Matrix<T> &m1, const Matrix<T> &m2);
    template <typename T> const Matrix<T> operator*(const Matrix<T> &m, const T c);
    template <typename T> const Matrix<T> operator*(const T c, const Matrix<T> &m);
    template <typename T> const Matrix<T> operator*(const Matrix<T> &m1, const Matrix<T> &m2);
}

/*constrctors*/
template <typename T> m_calc::Matrix<T>::Matrix(Mat v) : data(v) {}
template <typename T> m_calc::Matrix<T>::Matrix(size_t size) : data(size, std::vector<T>(size, 0)) {}
template <typename T> m_calc::Matrix<T>::Matrix(size_t row, size_t col) : data(row, std::vector<T>(col, 0)) {}

/*definition of member functions*/
template <typename T> bool m_calc::Matrix<T>::operator==(const Matrix &m) const {
    if (row_num() != m.row_num() || column_num() != m.column_num()) {
        return false;
    }
    for (int i = 0; i < row_num(); i++) {
        for (int j = 0; j < column_num(); j++) {
            if (data.at(i).at(j) != m.data.at(i).at(j)) {
                return false;
            }
        }
    }
    return true;
}
template <typename T> inline bool m_calc::Matrix<T>::operator!=(const Matrix &m) const {
    return !(*this == m);
}
template <typename T> m_calc::Matrix<T> &m_calc::Matrix<T>::operator=(const Matrix &m) {
    data = m.data;
    return *this;
}
template <typename T> m_calc::Matrix<T> &m_calc::Matrix<T>::operator+=(const Matrix &m) {
    if (row_num() != m.row_num() && column_num() != m.column_num()) {
        std::cerr << "The pair of the matrices are an error!";
        std::abort();
    }
    for (int i = 0; i < row_num(); i++) {
        for (int j = 0; j < column_num(); j++) {
            data.at(i).at(j) -= m.data.at(i).at(j);
        }
    }
    return *this;
}
template <typename T> m_calc::Matrix<T> &m_calc::Matrix<T>::operator-=(const Matrix &m) {
    if (row_num() != m.row_num() && column_num() != m.column_num()) {
        std::cerr << "The pair of the matrices are an error !";
        std::abort();
    }
    for (int i = 0; i < row_num(); i++) {
        for (int j = 0; j < column_num(); j++) {
            data.at(i).at(j) -= m.data.at(i).at(j);
        }
    }
    return *this;
}
template <typename T> m_calc::Matrix<T> &m_calc::Matrix<T>::operator*=(const T c) {
    for (int i = 0; i < row_num(); i++) {
        for (int j = 0; j < column_num(); j++) {
            data.at(i).at(j) *= c;
        }
    }
    return *this;
    }
template <typename T> m_calc::Matrix<T> &m_calc::Matrix<T>::operator*=(const Matrix &m) {
    if (column_num() != m.row_num()) {
        std::cerr << "The pair of the matrices are an error!";
        std::abort();
    }
    Mat data_result(row_num(), std::vector<T>(m.column_num(), 0));
    for (int i = 0; i < row_num(); i++) {
        for (int j = 0; j < m.column_num(); j++) {
            for (int k = 0; k < column_num(); k++) {
                data_result.at(i).at(j) += data.at(i).at(k) * m.data.at(k).at(j);
            }
        }
    }
    data = data_result;
    return *this;
}

template <typename T> inline const m_calc::Matrix<T>::Mat m_calc::Matrix<T>::get_data() const {
    return data;
}
template <typename T> inline const T m_calc::Matrix<T>::get_element(int i, int j) const {
    return data.at(i).at(j);
}
template <typename T> inline const size_t m_calc::Matrix<T>::row_num() const {
    return data.size();
}
template <typename T> inline const size_t m_calc::Matrix<T>::column_num() const {
    return data.at(0).size();
}
template <typename T> inline const bool m_calc::Matrix<T>::is_square() const {
    return row_num() == column_num();
}
template <typename T> const m_calc::Matrix<T> m_calc::Matrix<T>::transpose() const {
m_calc::Matrix<T> m_result(column_num(), row_num());
for (int i = 0; i < row_num(); i++) {
    for (int j = 0; j < column_num(); j++) {
        m_result.data.at(j).at(i) = data.at(i).at(j);
    }
}
return m_result;
}
template <typename T> const m_calc::Matrix<T> m_calc::Matrix<T>::echelon() const {
    Matrix<T> m_result = *this;
    bool skip = false;
    for (int i = 0; i < row_num(); i++) {
        if (m_result.data.at(i).at(i) == 0) {
            for (int k = i; k < row_num(); k++) {
                if (m_result.data.at(k).at(i) != 0) {
                    m_result = elementary_switch(row_num(), k, i) * m_result;
                    break;
                } else if (k == row_num() - 1) {
                    skip = true;
                }
            }
        }
        if (!skip) {
            T s = 1 / m_result.data.at(i).at(i);
            m_result = elementary(row_num(), i, s) * m_result;
            for (int k = 0; k < row_num(); k++) {
                if (k != i) {
                    T c = (-1) * (m_result.data.at(k).at(i) /
                                    m_result.data.at(i).at(i));
                    m_result = elementary(row_num(), k, i, c) * m_result;
                }
            }
        }
    }
    return m_result;
}
template <typename T> const m_calc::Matrix<T> m_calc::Matrix<T>::inverse() const {
    Matrix<T> m_result = identity(row_num());
    Matrix<T> m_echelon = *this;
    if (!is_square()) {
        std::cerr << "This is not a square matrix!" << std::endl;
        std::abort();
    }
    bool skip = false;
    for (int i = 0; i < row_num(); i++) {
        if (m_echelon.data.at(i).at(i) == 0) {
            for (int k = i; k < row_num(); k++) {
                if (m_echelon.data.at(k).at(i) != 0) {
                    m_result = Matrix<T>::elementary_switch(row_num(), k, i) * m_result;
                    m_echelon = elementary_switch(row_num(), k, i) * m_echelon;
                    break;
                } else if (k == row_num() - 1) {
                    skip = true;
                }
            }
        }
        if (!skip) {
            T s = 1 / m_echelon.data.at(i).at(i);
            m_result = elementary(row_num(), i, s) * m_result;
            m_echelon = elementary(row_num(), i, s) * m_echelon;
            for (int k = 0; k < row_num(); k++) {
                if (k != i) {
                    T c = (-1) * (m_echelon.data.at(k).at(i) / m_echelon.data.at(i).at(i));
                    m_result = elementary(row_num(), k, i, c) * m_result;
                    m_echelon = elementary(row_num(), k, i, c) * m_echelon;
                }
            }
        }
    }
    return m_result;
}
template <typename T> const int m_calc::Matrix<T>::rank() const {
int rank = 0;
Mat e = echelon().data;
for (int i = 0; i < row_num(); i++) {
    for (int j = i; j < column_num(); j++) {
        if (e.at(i).at(j) != 0) {
            rank++;
            break;
        }
    }
}
return rank;
}
template <typename T> const T m_calc::Matrix<T>::determinant() const {
m_calc::Matrix<T> m_tri = *this;
T det;
if (is_square()) {
    std::cerr << "This is not a square matrix!" << std::endl;
}
for (int i = 0; i < row_num(); i++) {
    for (int j = i + 1; j < column_num(); j++) {
        T c = m_tri.data.at(j).at(i) / m_tri.data.at(i).at(i);
        m_tri = elementary(row_num(), j, i, (-1 * c)) * m_tri;
    }
}
det = m_tri.data.at(0).at(0);
for (int i = 1; i < row_num(); i++) {
    det *= m_tri.data.at(i).at(i);
}
return det;
}
template <typename T> void m_calc::Matrix<T>::print() const {
for (int i = 0; i < row_num(); i++) {
    std::cout << data.at(i).at(0);
    for (int j = 1; j < column_num(); j++) {
        std::cout << " " << data.at(i).at(j);
    }
    std::cout << std::endl;
}
return;
}
template <typename T> m_calc::Matrix<T> m_calc::Matrix<T>::identity(size_t size) {
Mat data(size, std::vector<T>(size, 0));
for (int i = 0; i < size; i++) {
    data.at(i).at(i) = 1;
}
m_calc::Matrix<T> m(data);
return m;
}
template <typename T> m_calc::Matrix<T> m_calc::Matrix<T>::elementary(size_t size, int i, T c) {
m_calc::Matrix<T> m = identity(size);
m.data.at(i).at(i) = c;
return m;
}
template <typename T> m_calc::Matrix<T> m_calc::Matrix<T>::elementary(size_t size, int i, int j, T c) {
m_calc::Matrix<T> m = identity(size);
m.data.at(i).at(j) = c;
return m;
}
template <typename T> m_calc::Matrix<T> m_calc::Matrix<T>::elementary_switch(size_t size, int i, int j) {
m_calc::Matrix<T> p = identity(size);
p.data.at(i).at(i) = 0;
p.data.at(j).at(j) = 0;
p.data.at(i).at(j) = 1;
p.data.at(j).at(i) = 1;
return p;
}

/*definition of functions*/
template <typename T> const m_calc::Matrix<T> m_calc::operator+(const Matrix<T> &m1, const Matrix<T> &m2) {
    static_assert(std::is_arithmetic<T>::value, "type T must be arithmetic!");
    m_calc::Matrix<T> m_result(m1);
    return m_result += m2;
}
template <typename T> const m_calc::Matrix<T> m_calc::operator-(const Matrix<T> &m1, const Matrix<T> &m2) {
    static_assert(std::is_arithmetic<T>::value, "type T must be arithmetic!");
    m_calc::Matrix<T> m_result(m1);
    return m_result -= m2;
}
template <typename T> const m_calc::Matrix<T> m_calc::operator*(const Matrix<T> &m, const T c) {
    static_assert(std::is_arithmetic<T>::value, "type T must be arithmetic!");
    m_calc::Matrix<T> m_result(m);
    return m_result *= c;
}
template <typename T> const m_calc::Matrix<T> m_calc::operator*(const T c, const Matrix<T> &m) {
    static_assert(std::is_arithmetic<T>::value, "type T must be arithmetic!");
    m_calc::Matrix<T> m_result(m);
    return m_result *= c;
}
template <typename T> const m_calc::Matrix<T> m_calc::operator*(const Matrix<T> &m1, const Matrix<T> &m2) {
    static_assert(std::is_arithmetic<T>::value, "type T must be arithmetic!");
    m_calc::Matrix<T> m_result(m1);
    return m_result *= m2;
}
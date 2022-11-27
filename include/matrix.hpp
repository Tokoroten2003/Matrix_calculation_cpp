#pragma once

#include<iostream>
#include<type_traits>
#include<vector>

namespace Mcalc {
    template <typename T> using Mat = typename std::vector<std::vector<T>>;
    template <typename T> class Matrix {
        public:
            static_assert(std::is_arithmetic<T>::value, "type T must be arithmetic!");

            /*constrctors*/
            Matrix(const Mat<T> v);
            Matrix(const size_t size);
            Matrix(const size_t row, const size_t col);

        private:
            Mat<T> data;

        public:
            bool operator==(const Matrix &m) const;
            inline bool operator!=(const Matrix &m) const;
            Matrix &operator=(const Matrix &m);
            Matrix &operator+=(const Matrix &m);
            Matrix &operator-=(const Matrix &m);
            Matrix &operator*=(const T c);
            Matrix &operator*=(const Matrix &m);

            /*member funcions*/
            inline const Mat<T> get_data() const;                       //return the matrix
            inline const T get_component(int i, int j) const;           //return the component of the matrix
            inline const size_t row_num() const;                        //return the number of rows of the matrix
            inline const size_t column_num() const;                     //return the mumber of columns of the matrix
            inline const bool is_square() const;                        //return if the matrix is square
            const Matrix transpose() const;                             //return the transposed matrix
            const Matrix echelon() const;                               //return the echelon matrix
            const Matrix inverse() const;                               //return the inverse matrix
            const int rank() const;                                     //return the rank of the matrix
            const T determinant() const;                                //return the determinant of the matrix
            void print() const;                                         //print out the matrix
            static Matrix identity(size_t size);                        //return identity matrix
            static Matrix elementary(size_t size, int i, T c);          //return elementary matrix
            static Matrix elementary(size_t size, int i, int j, T c);   //return elementary matrix
            static Matrix elementary_swap(size_t size, int i, int j);   //return elementary matrix
    };

    template<typename T> const Matrix<T> operator+(const Matrix<T> &m1, const Matrix<T> &m2);
    template <typename T> const Matrix<T> operator-(const Matrix<T> &m1, const Matrix<T> &m2);
    template <typename T> const Matrix<T> operator*(const Matrix<T> &m, const T c);
    template <typename T> const Matrix<T> operator*(const T c, const Matrix<T> &m);
    template <typename T> const Matrix<T> operator*(const Matrix<T> &m1, const Matrix<T> &m2);
}

/*constrctors*/
template <typename T> Mcalc::Matrix<T>::Matrix(Mat<T> v) : data(v) {}
template <typename T> Mcalc::Matrix<T>::Matrix(size_t size) : data(size, std::vector<T>(size, 0)) {}
template <typename T> Mcalc::Matrix<T>::Matrix(size_t row, size_t col) : data(row, std::vector<T>(col, 0)) {}

/*definition of member functions*/
template <typename T> bool Mcalc::Matrix<T>::operator==(const Matrix &m) const {
    if (row_num() != m.row_num() || column_num() != m.column_num()) {       //if the type of matrices are different, return false
        return false;
    }
    for (int i = 0; i < row_num(); i++) {                                   //if a component is different to the other one, return false
        for (int j = 0; j < column_num(); j++) {
            if (data.at(i).at(j) != m.data.at(i).at(j)) {
                return false;
            }
        }
    }
    return true;
}
template <typename T> inline bool Mcalc::Matrix<T>::operator!=(const Matrix &m) const {
    return !(*this == m);
}
template <typename T> Mcalc::Matrix<T> &Mcalc::Matrix<T>::operator=(const Matrix &m) {
    data = m.data;
    return *this;
}
template <typename T> Mcalc::Matrix<T> &Mcalc::Matrix<T>::operator+=(const Matrix &m) {
    if (row_num() != m.row_num() && column_num() != m.column_num()) {       //if the type of matrices are different, stop running
        std::cerr << "The pair of the matrices are an error!";
        std::abort();
    }
    for (int i = 0; i < row_num(); i++) {                                   //add all every component to each component corresponding
        for (int j = 0; j < column_num(); j++) {
            data.at(i).at(j) += m.data.at(i).at(j);
        }
    }
    return *this;
}
template <typename T> Mcalc::Matrix<T> &Mcalc::Matrix<T>::operator-=(const Matrix &m) {
    if (row_num() != m.row_num() && column_num() != m.column_num()) {       //if the type of matrices are different, stop running
        std::cerr << "The pair of the matrices are an error !";
        std::abort();
    }
    for (int i = 0; i < row_num(); i++) {                                   //substract every component from each component corresponding
        for (int j = 0; j < column_num(); j++) {
            data.at(i).at(j) -= m.data.at(i).at(j);
        }
    }
    return *this;
}
template <typename T> Mcalc::Matrix<T> &Mcalc::Matrix<T>::operator*=(const T c) {
    for (int i = 0; i < row_num(); i++) {           //multiply every component by c
        for (int j = 0; j < column_num(); j++) {
            data.at(i).at(j) *= c;
        }
    }
    return *this;
    }
template <typename T> Mcalc::Matrix<T> &Mcalc::Matrix<T>::operator*=(const Matrix &m) {
    if (column_num() != m.row_num()) {                          //if the type of matrices is not appropreate, stop running
        std::cerr << "The pair of the matrices are an error!";
        std::abort();
    }
    Mat<T> data_result(row_num(), std::vector<T>(m.column_num(), 0));
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

template <typename T> inline const Mcalc::Mat<T> Mcalc::Matrix<T>::get_data() const {
    return data;
}
template <typename T> inline const T Mcalc::Matrix<T>::get_component(int i, int j) const {
    return data.at(i).at(j);
}
template <typename T> inline const size_t Mcalc::Matrix<T>::row_num() const {
    return data.size();
}
template <typename T> inline const size_t Mcalc::Matrix<T>::column_num() const {
    return data.at(0).size();
}
template <typename T> inline const bool Mcalc::Matrix<T>::is_square() const {
    return row_num() == column_num();
}
template <typename T> const Mcalc::Matrix<T> Mcalc::Matrix<T>::transpose() const {
Mcalc::Matrix<T> m_result(column_num(), row_num());        //make a matrix whose column size and row size is swaped
for (int i = 0; i < row_num(); i++) {
    for (int j = 0; j < column_num(); j++) {                //swap the component of column and row
        m_result.data.at(j).at(i) = data.at(i).at(j);
    }
}
return m_result;
}
template <typename T> const Mcalc::Matrix<T> Mcalc::Matrix<T>::echelon() const {
    Matrix<T> m_result = *this;
    bool skip = false;
    for (int i = 0; i < row_num(); i++) {
        if (m_result.data.at(i).at(i) == 0) {           //if the diagonal component of the row is "0", swap the order with the one with a diagonal component which is not "0"
            for (int k = i; k < row_num(); k++) {
                if (m_result.data.at(k).at(i) != 0) {
                    m_result = elementary_swap(row_num(), k, i) * m_result;
                    break;
                } else if (k == row_num() - 1) {        //if there is no appropreate rows, skip the next step
                    skip = true;
                }
            }
        }
        if (!skip) {                                    //if there is no need to skip, do elementary row operation
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
template <typename T> const Mcalc::Matrix<T> Mcalc::Matrix<T>::inverse() const {
    Matrix<T> m_result = identity(row_num());
    Matrix<T> m_echelon = *this;
    if (!is_square()) {     //if the matrix is not square, stop running
        std::cerr << "This is not a square matrix!" << std::endl;
        std::abort();
    }
    bool skip = false;
    for (int i = 0; i < row_num(); i++) {       //do same operation to echelon() to the identity matrix
        if (m_echelon.data.at(i).at(i) == 0) {
            for (int k = i; k < row_num(); k++) {
                if (m_echelon.data.at(k).at(i) != 0) {
                    m_result = Matrix<T>::elementary_swap(row_num(), k, i) * m_result;
                    m_echelon = elementary_swap(row_num(), k, i) * m_echelon;
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
template <typename T> const int Mcalc::Matrix<T>::rank() const {
int rank = 0;
Mat<T> e = echelon().data;
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
template <typename T> const T Mcalc::Matrix<T>::determinant() const {
if (!is_square()) {                                              //if the matrix is not square, stop running
    std::cerr << "This is not a square matrix!" << std::endl;
}
Mcalc::Matrix<T> m_tri = *this;
T det;
for (int i = 0; i < row_num(); i++) {                           //make an upper triangular matrix
    if (m_tri.data.at(i).at(i) == 0) {
        for (int k = i + 1; k < row_num(); k++) {
            if (m_tri.data.at(k).at(i) != 0) {
                m_tri = elementary_swap(row_num(), i, k) * m_tri;
                break;
            }
            if (k == row_num() - 1) {return 0;}
        }
    }
    for (int k = i + 1; k < row_num(); k++) {
        T c = m_tri.data.at(k).at(i) / m_tri.data.at(i).at(i);
        m_tri = elementary(row_num(), k, i, (-1 * c)) * m_tri;
    }
}
det = m_tri.data.at(0).at(0);
for (int i = 1; i < row_num(); i++) {                           //multiple all the diagonal components of the upper triangular matirx
    det *= m_tri.data.at(i).at(i);
}
return det;
}
template <typename T> void Mcalc::Matrix<T>::print() const {
for (int i = 0; i < row_num(); i++) {
    std::cout << data.at(i).at(0);
    for (int j = 1; j < column_num(); j++) {
        std::cout << " " << data.at(i).at(j);
    }
    std::cout << std::endl;
}
return;
}
template <typename T> Mcalc::Matrix<T> Mcalc::Matrix<T>::identity(size_t size) {
Mat<T> data(size, std::vector<T>(size, 0));    //all components are "0"
for (int i = 0; i < size; i++) {            //change the diagonal components to "1"
    data.at(i).at(i) = 1;
}
Mcalc::Matrix<T> m(data);
return m;
}
template <typename T> Mcalc::Matrix<T> Mcalc::Matrix<T>::elementary(size_t size, int i, T c) {
Mcalc::Matrix<T> m = identity(size);
m.data.at(i).at(i) = c;
return m;
}
template <typename T> Mcalc::Matrix<T> Mcalc::Matrix<T>::elementary(size_t size, int i, int j, T c) {
Mcalc::Matrix<T> m = identity(size);
m.data.at(i).at(j) = c;
return m;
}
template <typename T> Mcalc::Matrix<T> Mcalc::Matrix<T>::elementary_swap(size_t size, int i, int j) {
Mcalc::Matrix<T> p = identity(size);
p.data.at(i).at(i) = 0;
p.data.at(j).at(j) = 0;
p.data.at(i).at(j) = 1;
p.data.at(j).at(i) = 1;
return p;
}

/*definition of functions*/
template <typename T> const Mcalc::Matrix<T> Mcalc::operator+(const Matrix<T> &m1, const Matrix<T> &m2) {
    static_assert(std::is_arithmetic<T>::value, "type T must be arithmetic!");
    Mcalc::Matrix<T> m_result(m1);
    return m_result += m2;
}
template <typename T> const Mcalc::Matrix<T> Mcalc::operator-(const Matrix<T> &m1, const Matrix<T> &m2) {
    static_assert(std::is_arithmetic<T>::value, "type T must be arithmetic!");
    Mcalc::Matrix<T> m_result(m1);
    return m_result -= m2;
}
template <typename T> const Mcalc::Matrix<T> Mcalc::operator*(const Matrix<T> &m, const T c) {
    static_assert(std::is_arithmetic<T>::value, "type T must be arithmetic!");
    Mcalc::Matrix<T> m_result(m);
    return m_result *= c;
}
template <typename T> const Mcalc::Matrix<T> Mcalc::operator*(const T c, const Matrix<T> &m) {
    static_assert(std::is_arithmetic<T>::value, "type T must be arithmetic!");
    Mcalc::Matrix<T> m_result(m);
    return m_result *= c;
}
template <typename T> const Mcalc::Matrix<T> Mcalc::operator*(const Matrix<T> &m1, const Matrix<T> &m2) {
    static_assert(std::is_arithmetic<T>::value, "type T must be arithmetic!");
    Mcalc::Matrix<T> m_result(m1);
    return m_result *= m2;
}
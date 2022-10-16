#pragma once

#include<iostream>
#include<vector>

namespace m_calc {
    template<typename T>
    class Matrix {
        public:
            using Mat = typename std::vector<std::vector<T>>;

        private:
            Mat data;
            int row_size = data.size();
            int column_size = data.at(0).size();
            bool is_square = (row_size == column_size);

        public:
            /*member funcions*/
            bool operator==(const Matrix &m) const;
            inline bool operator!=(const Matrix &m) const;
            Matrix &operator=(const Matrix &m);
            Matrix &operator+=(const Matrix &m);
            Matrix &operator-=(const Matrix &m);
            Matrix &operator*=(const T c);
            Matrix &operator*=(const Matrix &m);

            const Matrix transpose() const;
            const Matrix echelon() const;
            const Matrix inverse() const;
            const int rank() const;
            const T determinant() const;
            void print() const;
            static Matrix identity(int size);
            static Matrix elementary(int size, int i, T c);
            static Matrix elementary(int size, int i, int j, T c);
            static Matrix elementary_switch(int size, int i, int j);

            /*constrctors*/
            Matrix(int size) : data(size, std::vector<T>(size, 0)) {}
            Matrix(int row_size, int column_size) : data(row_size, std::vector<T>(column_size, 0)) {}
            Matrix(Mat v) : data(v) {}
    };

    /*definition of member functions*/
    template<typename T>
    bool Matrix<T>::operator==(const Matrix &m) const {
            bool equality = true;
            if(row_size == m.row_size && column_size == m.column_size) {
                for(int i=0; i<row_size; i++) {
                    for(int j=0; j<column_size; j++) {
                        if(data.at(i).at(j) != m.data.at(i).at(j)) {
                            equality = false;
                            break;
                        }
                    }
                }
            }
            else {
                equality = false;
            }
            return equality;
    }
    template<typename T>
    inline bool Matrix<T>::operator!=(const Matrix &m) const {
        return !(*this == m);
    }
    template<typename T>
    Matrix<T> &Matrix<T>::operator=(const Matrix &m) {
        data = m.data;
        return *this;
    }
    template<typename T>
    Matrix<T> &Matrix<T>::operator+=(const Matrix &m) {
        if(row_size == m.row_size && column_size == m.column_size) {
            for(int i=0; i<row_size; i++) {
                for(int j=0; j<column_size; j++) {
                    data.at(i).at(j) += m.data.at(i).at(j);
                }
            }
        }
        else {
            std::cout << "The pair of the matrices are an error!" << std::endl;
        }
        return *this;
    }
    template<typename T>
    Matrix<T> &Matrix<T>::operator-=(const Matrix &m) {
        if(row_size == m.row_size && column_size == m.column_size) {
            for(int i=0; i<row_size; i++) {
                for(int j=0; j<column_size; j++) {
                    data.at(i).at(j) -= m.data.at(i).at(j);
                }
            }
        }
        else {
            std::cout << "The pair of the matrices are an error!" << std::endl;
        }
        return *this;
    }
    template<typename T>
    Matrix<T> &Matrix<T>::operator*=(const T c) {
        for(int i=0; i<row_size; i++) {
            for(int j=0; j<column_size; j++) {
                data.at(i).at(j) *= c;
            }
        }
        return *this;
    }
    template<typename T>
    Matrix<T> &Matrix<T>::operator*=(const Matrix &m) {
        Mat data_result(row_size, std::vector<T>(column_size, 0));
        if(column_size == m.row_size) {
            for (int i=0; i<row_size; i++) {
                for(int j=0; j<m.column_size; j++) {
                    for(int k=0; k<column_size; k++) {
                        data_result.at(i).at(j) += data.at(i).at(k) * m.data.at(k).at(j);
                    }
                }
            }
        }
        else {
            std::cout << "The pair of the matrices are an error!" << std::endl;
        }
        data = data_result;
        return *this;
    }

    template<typename T>
    const Matrix<T> Matrix<T>::transpose() const {
        Matrix<T> m_result(column_size, row_size);
        for(int i=0; i<row_size; i++) {
            for(int j=0; j<column_size; j++) {
                m_result.data.at(j).at(i) = data.at(i).at(j);
            }
        }
        return m_result;
    }
    template<typename T>
    const Matrix<T> Matrix<T>::echelon() const {
        Matrix<T> m_result = *this;
        bool skip = false;
        for(int i=0; i<row_size; i++) {
            if(m_result.data.at(i).at(i) == 0) {
                for(int k=i; k<row_size; k++) {
                    if(m_result.data.at(k).at(i) != 0) {
                        m_result = elementary_switch(row_size, k, i) * m_result;
                        break;
                    }
                    else if(k == row_size - 1) {
                        skip = true;
                    }
                }
            }
            if(!skip) {
                T s = 1 / m_result.data.at(i).at(i);
                m_result = elementary(row_size, i, s) * m_result;
                for(int k=0; k<row_size; k++) {
                    if(k != i) {
                        T c = (-1) * (m_result.data.at(k).at(i) / m_result.data.at(i).at(i));
                        m_result = elementary(row_size, k, i, c) * m_result;
                    }
                    else {}
                }
            }
            else {}
        }
        return m_result;
    }
    template<typename T>
    const Matrix<T> Matrix<T>::inverse() const {
        Matrix<T> m_result = identity(row_size);
        Matrix<T> m_echelon = *this;
        if(is_square) {
            bool skip = false;
            for(int i=0; i<row_size; i++) {
                if(m_echelon.data.at(i).at(i) == 0) {
                    for(int k=i; k<row_size; k++) {
                        if(m_echelon.data.at(k).at(i) != 0) {
                            m_result = elementary_switch(row_size, k, i) * m_result;
                            m_echelon = elementary_switch(row_size, k, i) * m_echelon;
                            break;
                        }
                        else if(k == row_size - 1) {
                            skip = true;
                        }
                    }
                }
                if(!skip) {
                    T s = 1 / m_echelon.data.at(i).at(i);
                    m_result = elementary(row_size, i, s) * m_result;
                    m_echelon = elementary(row_size, i, s) * m_echelon;
                    for(int k=0; k<row_size; k++) {
                        if(k != i) {
                            T c = (-1) * (m_echelon.data.at(k).at(i) / m_echelon.data.at(i).at(i));
                            m_result = elementary(row_size, k, i, c) * m_result;
                            m_echelon = elementary(row_size, k, i, c) * m_echelon;
                        }
                        else {}
                    }
                }
                else {}
            }
        }
        else {
            std::cout << "This is not a square matrix!" << std::endl;
        }
        return m_result;
    }
    template<typename T>
    const int Matrix<T>::rank() const {
        int rank = 0;
        Mat e = echelon().data;
        for(int i=0; i<row_size; i++) {
            for(int j=i; j<column_size; j++) {
                if(e.at(i).at(j) != 0) {
                    rank++;
                    break;
                }
            }
        }
        return rank;
    }
    template<typename T>
    const T Matrix<T>::determinant() const {
        Matrix<T> m_tri = *this;
        T det;
        if(is_square) {
            for(int i=0; i<row_size; i++) {
                for(int j=i+1; j<column_size; j++) {
                    T c = m_tri.data.at(j).at(i) / m_tri.data.at(i).at(i);
                    m_tri = elementary(row_size, j, i, (-1 * c)) * m_tri;
                }
            }
            det = m_tri.data.at(0).at(0);
            for(int i=1; i<row_size; i++) {
                det *= m_tri.data.at(i).at(i);
            }
        }
        else {
            std::cout << "This is not a square matrix!" << std::endl;
        }
        return det;
    }
    template<typename T>
    void Matrix<T>::print() const {
        for(int i=0; i<row_size; i++) {
            std::cout << data.at(i).at(0);
            for(int j=1; j<column_size; j++) {
                std::cout << " " << data.at(i).at(j);
            }
            std::cout << std::endl;
        }
        return;
    }
    template<typename T>
    Matrix<T> Matrix<T>::identity(int size) {
        Mat data(size, std::vector<T>(size, 0));
        for(int i=0; i<size; i++) {
            data.at(i).at(i) = 1;
        }
        Matrix<T> m(data);
        return m;
    }
    template<typename T>
    Matrix<T> Matrix<T>::elementary(int size, int i, T c) {
        Matrix<T> m = identity(size);
        m.data.at(i).at(i) = c;
        return m;
    }
    template<typename T>
    Matrix<T> Matrix<T>::elementary(int size, int i, int j, T c) {
        Matrix<T> m = identity(size);
        m.data.at(i).at(j) = c;
        return m;
    }
    template<typename T>
    Matrix<T> Matrix<T>::elementary_switch(int size, int i, int j) {
        Matrix<T> p = identity(size);
        p.data.at(i).at(i) = 0;
        p.data.at(j).at(j) = 0;
        p.data.at(i).at(j) = 1;
        p.data.at(j).at(i) = 1;
        return p;
    }

    /*daclaration of functions*/
    template<typename T>
    const Matrix<T> operator+(const Matrix<T> &m1, const Matrix<T> &m2) {
        Matrix<T> m_result = m1;
        m_result += m2;
        return m_result;
    }
    template<typename T>
    const Matrix<T> operator-(const Matrix<T> &m1, const Matrix<T> &m2) {
        Matrix<T> m_result = m1;
        m_result -= m2;
        return m_result;
    }
    template<typename T>
    const Matrix<T> operator*(const Matrix<T> &m, const T c) {
        Matrix m_result = m;
        m_result *= c;
        return m_result;
    }
    template<typename T>
    const Matrix<T> operator*(const T c, const Matrix<T> &m) {
        return m * c;
    }
    template<typename T>
    const Matrix<T> operator*(const Matrix<T> &m1, const Matrix<T> &m2) {
        Matrix<T> m_result = m1;
        m_result *= m2;
        return m_result;
    }

}

#pragma once

#include<iostream>
#include<vector>

namespace l_alg {
    template<typename T>
    class Matrix {
        public:
            using Mat = typename std::vector<std::vector<T>>;
            Mat data;
            int rsize = data.size();
            int csize = data.at(0).size();
            bool is_square = (rsize == csize);

            /*member funcions*/
            bool operator==(const Matrix &m) const;
            inline bool operator!=(const Matrix &m) const;
            Matrix &operator=(const Matrix &m);
            Matrix &operator+=(const Matrix &m);
            Matrix &operator-=(const Matrix &m);
            Matrix &operator*=(const T c);
            Matrix &operator*=(const Matrix &m);

            const Matrix transpose() const;
            const T determinant() const;
            void print() const;
            static Matrix makeUnit(int n);

            /*constrctors*/
            Matrix(int size) : data(size, std::vector<T>(size, 0)) {}
            Matrix(int rsize, int csize) : data(rsize, std::vector<T>(csize, 0)) {}
            Matrix(Mat v) : data(v) {}
    };

    /*definition of member functions*/
    template<typename T>
    bool Matrix<T>::operator==(const Matrix &m) const {
            bool equality = true;
            if(rsize == m.rsize && csize == m.csize) {
                for(int i=0; i<rsize; i++) {
                    for(int j=0; j<csize; j++) {
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
        if(rsize == m.rsize && csize == m.csize) {
            for(int i=0; i<rsize; i++) {
                for(int j=0; j<csize; j++) {
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
        if(rsize == m.rsize && csize == m.csize) {
            for(int i=0; i<rsize; i++) {
                for(int j=0; j<csize; j++) {
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
        for(int i=0; i<rsize; i++) {
            for(int j=0; j<csize; j++) {
                data.at(i).at(j) *= c;
            }
        }
        return *this;
    }
    template<typename T>
    Matrix<T> &Matrix<T>::operator*=(const Matrix &m) {
        Mat data_result(rsize, std::vector<T>(csize, 0));
        if(csize == m.rsize) {
            for (int i=0; i<rsize; i++) {
                for(int j=0; j<m.csize; j++) {
                    for(int k=0; k<csize; k++) {
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
        Matrix<T> m_result(csize, rsize);
        for(int i=0; i<rsize; i++) {
            for(int j=0; j<csize; j++) {
                m_result.data.at(j).at(i) = data.at(i).at(j);
            }
        }
        return m_result;
    }
    template<typename T>
    const T Matrix<T>::determinant() const {
        Matrix<T> m_tri = *this;
        T det;
        if(is_square) {
            for(int i=0; i<rsize; i++) {
                for(int j=i+1; j<csize; j++) {
                    T c = m_tri.data.at(j).at(i) / m_tri.data.at(i).at(i);
                    for(int k=0; k<csize; k++) {
                        m_tri.data.at(j).at(k) -= c * m_tri.data.at(i).at(k);
                    }
                }
            }
            det = m_tri.data.at(0).at(0);
            for(int i=1; i<rsize; i++) {
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
        for(int i=0; i<rsize; i++) {
            std::cout << data.at(i).at(0);
            for(int j=1; j<csize; j++) {
                std::cout << " " << data.at(i).at(j);
            }
            std::cout << std::endl;
        }
        return;
    }
    template<typename T>
    Matrix<T> Matrix<T>::makeUnit(int size) {
        Mat data(size, std::vector<T>(size, 0));
        for(int i=0; i<size; i++) {
            data.at(i).at(i) = 1;
        }
        Matrix<T> m(data);
        return m;
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

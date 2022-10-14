#pragma once
#include<iostream>
#include<vector>

namespace l_alg {
    template<typename T>
    struct Matrix {
        using Mat = typename std::vector<std::vector<T>>;
        Mat data;
        int rsize = data.size();
        int csize = data.at(0).size();
        bool is_square = (rsize == csize);

        bool operator==(const Matrix &m) const;
        inline bool operator!=(const Matrix &m) const;
        Matrix &operator=(const Matrix &m);
        inline Matrix &operator+=(const Matrix &m);
        inline Matrix &operator-=(const Matrix &m);
        inline Matrix &operator*=(const T c);
        inline Matrix &operator*=(const Matrix &m);

        const Matrix transpose(const Matrix &m) const;
        const T determinant(const Matrix &m) const;

        void print() {
            for(int i=0; i<rsize; i++) {
                std::cout << data.at(i).at(0);
                for(int j=1; j<csize; j++) {
                    std::cout << " " << data.at(i).at(j);
                }
                std::cout << std::endl;
            }
            return;
        }
        Matrix(int m, int n) : data(m, std::vector<T>(n, 0)) {}
        Matrix(Mat v) : data(v) {}
    };
    template<typename T> 
    bool Matrix<T>::operator==(const Matrix &m) const {
            bool equality = true;
            for(int i=0; i<rsize; i++) {
                for(int j=0; j<csize; j++) {
                    if(data.at(i).at(j) != m.data.at(i).at(j)) {
                        equality != equality;                       //all entries equals to each other
                        break;
                    }
                }
                if(!equality) {
                    break;
                }
            }
            return equality;
    }
    template<typename T> 
    inline bool Matrix<T>::operator!=(const Matrix &m) const {
        return !(*this == m);
    }
    template<typename T> 
    Matrix<T> &Matrix<T>::operator=(const Matrix &m) {
        if(rsize == m.rsize && csize == m.csize) {
            for(int i=0; i<rsize; i++) {
                for(int j=0; j<csize; j++) {
                    data.at(i).at(j) = m.data.at(i).at(j);
                }
            }
        }
        return *this;
    }
    template<typename T> 
    const Matrix<T> operator+(const Matrix<T> &m1, const Matrix<T> &m2) {
        Matrix m3(m1.rsize, m1.csize);
        if(m1.rsize == m2.rsize && m1.csize == m2.csize) {
            for(int i=0; i<m1.rsize; i++) {
                for(int j=0; j<m1.csize; j++) {
                    m3.data.at(i).at(j) = m1.data.at(i).at(j) + m2.data.at(i).at(j);
                }
            }
        }
        else {
            m3 = m1;
        }
        return m3;
    }
    template<typename T> 
    const Matrix<T> operator-(const Matrix<T> &m1, const Matrix<T> &m2) {
        Matrix m3(m1.rsize, m1.csize);
        if(m1.rsize == m2.rsize && m1.csize == m2.csize) {
            for(int i=0; i<m1.rsize; i++) {
                for(int j=0; j<m1.csize; j++) {
                    m3.data.at(i).at(j) = m1.data.at(i).at(j) - m2.data.at(i).at(j);
                }
            }
        }
        else {
            m3 = m1;
        }
        return m3;
    }
    template<typename T> 
    const Matrix<T> operator*(const Matrix<T> &m, const T c) {
        Matrix m_result(m.rsize, m.csize);
        for(int i=0; i<m.rsize; i++) {
            for(int j=0; j<m.csize; j++) {
                m_result.data.at(i).at(j) = m.data.at(i).at(j) * c;
            }
        }
        return m_result;
    }
    template<typename T>
    const Matrix<T> operator*(const T c, const Matrix<T> &m) {
        return m * c;
    }
    template<typename T> 
    const Matrix<T> operator*(const Matrix<T> &m1, const Matrix<T> &m2) {
        Matrix m3(m1.rsize, m2.csize);
        if(m1.csize == m2.rsize) {
            for (int i=0; i<m1.rsize; i++) {
                for(int j=0; j<m2.csize; j++) {
                    for(int k=0; k<m1.csize; k++) {
                        m3.data.at(i).at(j) = m1.data.at(i).at(k) + m2.data.at(k).at(j);
                    }
                }
            }
        }
        else {
            m3 = m1;
        }
        return m3;
    }
    template<typename T> 
    inline Matrix<T> &Matrix<T>::operator+=(const Matrix &m) {
        *this = *this + m;
        return *this;
    }
    template<typename T> 
    inline Matrix<T> &Matrix<T>::operator-=(const Matrix &m) {
        *this = *this - m;
        return *this;
    }
    template<typename T>
    inline Matrix<T> &Matrix<T>::operator*=(const T c) {
        *this = *this * c;
        return *this;
    }
    template<typename T> 
    inline Matrix<T> &Matrix<T>::operator*=(const Matrix &m) {
        *this = *this * m;
        return *this;
    }
    template<typename T> 
    const Matrix<T> Matrix<T>::transpose(const Matrix &m) const {
        Matrix m_result(m.csize, m.rsize);
        for(int i=0; i<m.rsize; i++) {
            for(int j=0; j<m.csize; j++) {
                m_result.data.at(j).at(i) = m.data.at(i).at(j);
            }
        }
        return m_result;
    }
    template<typename T>
    const T Matrix<T>::determinant(const Matrix &m) const {
        Matrix m_tri(m.rsize, m.csize);
        T det;
        if(m.is_square) {
            for(int i=0; i<m.rsize; i++) {
                for(int j=0; j<i; j++) {
                    T c = m.data.at(j).at(i) / m.data.at(i).at(i);
                    for(int k=0; k<m.csize; k++) {
                        m_tri.data.at(j).at(k) = m.data.at(j).at(k) - (c * m.data.at(i).at(k));
                    }
                }
            }
            det = m_tri.at(0).at(0);
            for(int i=1; i<m.rsize; i++) {
                det *= m_tri.data.at(i).at(i);
            }
        }
        return det;
    }
}
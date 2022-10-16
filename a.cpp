#include<iostream>
#include<vector>
#include"matrix.hpp"
using namespace m_calc;

int main() {
    Matrix<float> a({{1, 2}, {3, 4}});
    Matrix<float> e = Matrix<float>::identity(2);
    Matrix<float> b = a + e;
    float c = 5;
    std::cout << "a = " << std::endl;
    a.print();
    std::cout << "c = " << std::endl << c << std::endl;
    std::cout << std::endl;
    std::cout << "a + e =" << std::endl;
    b.print();
    std::cout << "a - e =" << std::endl;
    (a - e).print();
    std::cout << "a * e =" << std::endl;
    (a * e).print();
    std::cout << "a * c =" << std::endl;
    (a * c).print();
    std::cout << "transposed a =" << std::endl;
    a.transpose().print();
    std::cout << "determinant a =" << std::endl;
    std::cout << a.determinant() << std::endl;
    std::cout << "echelon form of a =" << std::endl;
    a.echelon().print();
    std::cout << "rank of a =" << std::endl << a.rank() << std::endl;
}
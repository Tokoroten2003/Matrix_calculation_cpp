#include<iostream>
#include<vector>
#include"linear_algebra.hpp"
using namespace l_alg;

int main() {
    Matrix<int> a({{1, 2}, {3, 4}});
    Matrix<int> e = Matrix<int>::makeUnit(2);
    Matrix<int> b = a + e;
    int c = 5;
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
    (a.transpose()).print();
    std::cout << "determinant a =" << std::endl;
    std::cout << a.determinant() << std::endl;
}
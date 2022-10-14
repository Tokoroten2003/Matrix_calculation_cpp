#include<iostream>
#include<vector>
#include"linear_algebra.hpp"
using namespace l_alg;

int main() {
    Matrix<int> a({{0, 0}, {0, 0}});
    Matrix<int> e = Matrix<int>::makeUnit(2);
    (a + e).print();
}
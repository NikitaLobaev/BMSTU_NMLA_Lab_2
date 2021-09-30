#include <Eigen/Dense>
#include <iostream>
#include "Gauss.cpp"
//#include "Matrix.cpp"

int main() {
    const Lobaev::Math::Matrix<double> matrix({
        {0, 1, 3},
        {2, 4, 6}
    });

    Lobaev::Math::Vector<double> solution(1); //stub
    try {
        solution = Lobaev::Math::Gauss::gauss(matrix);
    } catch (const char *exception) {
        std::cerr << exception << std::endl;
        return 1;
    }

    Eigen::Matrix<int, 2, 3> eigen_matrix;
    eigen_matrix(0, 0) = 0;
    eigen_matrix(0, 1) = 1;
    eigen_matrix(0, 2) = 3;
    eigen_matrix(1, 0) = 2;
    eigen_matrix(1, 1) = 4;
    eigen_matrix(1, 2) = 6;

    return 0;
}


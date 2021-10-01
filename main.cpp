#include <Eigen/Dense>
#include "Gauss.cpp"

int main() {
    const Lobaev::Math::Matrix<double> matrix_a({
        {0, 1},
        {2, 4}
    });
    const Lobaev::Math::Vector<double> vector_f({
        3,
        6
    });

    const Lobaev::Math::Vector<double> solution = Lobaev::Math::Gauss::gauss(matrix_a, vector_f);

    Eigen::Matrix<int, 2, 2> eigen_matrix_a;
    eigen_matrix_a(0, 0) = 0;
    eigen_matrix_a(0, 1) = 1;
    eigen_matrix_a(1, 0) = 2;
    eigen_matrix_a(1, 1) = 4;

    Eigen::Matrix<int, 2, 1> eigen_vector_f;
    eigen_vector_f(0, 0) = 3;
    eigen_vector_f(1, 0) = 6;

    const Eigen::Matrix<int, 2, 1> eigen_vector_x = eigen_matrix_a * eigen_vector_f;

    return 0;
}


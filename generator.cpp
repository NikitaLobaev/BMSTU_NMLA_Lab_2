#include "Matrix.h"

template <class T>
Lobaev::Math::Matrix<T> gen_random_matrix(const size_t dimension, const T from, const T to) {
    using namespace Lobaev::Math;

    Matrix<T> matrix(dimension, dimension);
    for (size_t row_index = 0; row_index < dimension; row_index++) {
        for (size_t column_index = 0; column_index < dimension; column_index++) {
            matrix(row_index, column_index) = (T) (std::rand() % (int) to) - from;
        }
    }
    return matrix;
}

template <class T>
void make_matrix_diagonally_dominant(Lobaev::Math::Matrix<T> &matrix, const T from, const T to) {
    if (!matrix.is_square()) {
        throw "make_matrix_diagonally_dominant: only square matrix may become diagonally dominant.";
    }

    for (size_t row_index = 0; row_index < matrix.rows_count(); row_index++) {
        T sum = (T) 0;

        for (size_t column_index = 0; column_index < matrix.columns_count(); column_index++) {
            if (row_index == column_index) {
                continue;
            }

            sum += matrix(row_index, column_index);
        }

        matrix(row_index, row_index) = sum + (T) (std::rand() % (int) to) - from + (T) 1;
    }
}

template <class T>
Lobaev::Math::Vector<T> gen_random_vector(const size_t dimension, const T from, const T to) {
    using namespace Lobaev::Math;

    Vector<T> vector(dimension);
    for (size_t index = 0; index < dimension; index++) {
        vector(index) = (T) (std::rand() % (int) to) - from;
    }
    return vector;
}

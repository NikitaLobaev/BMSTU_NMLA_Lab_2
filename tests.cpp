#include <cmath>
#include "gtest/gtest.h"
#include "Gauss.cpp"

using namespace Lobaev::Math;
using namespace Lobaev::Math::Gauss;

template <class T>
Vector<T> round(const Vector<T>&);

TEST(gauss_basic, test_1_ok) {
    const Matrix<double> matrix_a({
        {2, 3},
        {5, 7}
    });
    const Vector<double> vector_f({
        11,
        13
    });

    const Vector<double> vector_x = round(gauss_basic(matrix_a, vector_f));

    const Vector<double> expected_vector_x({
        -38,
        29
    });

    ASSERT_EQ(vector_x, expected_vector_x);
}

TEST(gauss_basic, test_2_ok) {
    const Matrix<double> matrix_a({
        {2, 1, 1},
        {1, -1, 0},
        {3, -1, 2}
    });
    const Vector<double> vector_f({
        2,
        -2,
        2
    });

    const Vector<double> vector_x = round(gauss_basic(matrix_a, vector_f));

    const Vector<double> expected_vector_x({
        -1,
        1,
        3
    });

    ASSERT_EQ(vector_x, expected_vector_x);
}

TEST(gauss_basic, test_matrix_a_dimension_fail) {
    const Matrix<double> matrix_a({
        {2, 1, 1},
        {1, -1, 0}
    });
    const Vector<double> vector_f({
        2,
        -2,
        0
    });

    ASSERT_ANY_THROW(gauss_basic(matrix_a, vector_f));
}

TEST(gauss_basic, test_vector_f_dimension_fail) {
    const Matrix<double> matrix_a({
        {2, 1},
        {1, -1}
    });
    const Vector<double> vector_f({
        2,
        -2,
        0
    });

    ASSERT_ANY_THROW(gauss_basic(matrix_a, vector_f));
}

TEST(gauss_basic, test_divide_by_0_fail) {
    const Matrix<double> matrix_a({
        {0, 1},
        {2, 4}
    });
    const Vector<double> vector_f({
        3,
        6
    });

    ASSERT_ANY_THROW(gauss_basic(matrix_a, vector_f));
}

TEST(gauss_by_row, test_1_ok) {
    const Matrix<double> matrix_a({
        {2, 3},
        {5, 7}
    });
    const Vector<double> vector_f({
        11,
        13
    });

    const Vector<double> vector_x = round(gauss_by_row(matrix_a, vector_f));

    const Vector<double> expected_vector_x({
        -38,
        29
    });

    ASSERT_EQ(vector_x, expected_vector_x);
}

TEST(gauss_by_row, test_2_ok) {
    const Matrix<double> matrix_a({
        {2, 1, 1},
        {1, -1, 0},
        {3, -1, 2}
    });
    const Vector<double> vector_f({
        2,
        -2,
        2
    });

    const Vector<double> vector_x = round(gauss_by_row(matrix_a, vector_f));

    const Vector<double> expected_vector_x({
        -1,
        1,
        3
    });

    ASSERT_EQ(vector_x, expected_vector_x);
}

TEST(gauss_by_row, test_3_ok) {
    const Matrix<double> matrix_a({
        {0, 1},
        {2, 4}
    });
    const Vector<double> vector_f({
        3,
        6
    });

    const Vector<double> vector_x = round(gauss_by_row(matrix_a, vector_f));

    const Vector<double> expected_vector_x({
        -3,
        3
    });

    ASSERT_EQ(vector_x, expected_vector_x);
}

TEST(gauss_by_row, test_matrix_a_dimension_fail) {
    const Matrix<double> matrix_a({
        {2, 1, 1},
        {1, -1, 0}
    });
    const Vector<double> vector_f({
        2,
        -2,
        0
    });

    ASSERT_ANY_THROW(gauss_by_row(matrix_a, vector_f));
}

TEST(gauss_by_row, test_vector_f_dimension_fail) {
    const Matrix<double> matrix_a({
        {2, 1},
        {1, -1}
    });
    const Vector<double> vector_f({
        2,
        -2,
        0
    });

    ASSERT_ANY_THROW(gauss_by_row(matrix_a, vector_f));
}

TEST(gauss_by_columns, test_1_ok) {
    const Matrix<double> matrix_a({
        {2, 3},
        {5, 7}
    });
    const Vector<double> vector_f({
        11,
        13
    });

    const Vector<double> vector_x = round(gauss_by_columns(matrix_a, vector_f));

    const Vector<double> expected_vector_x({
        -38,
        29
    });

    ASSERT_EQ(vector_x, expected_vector_x);
}

TEST(gauss_by_columns, test_2_ok) {
    const Matrix<double> matrix_a({
        {2, 1, 1},
        {1, -1, 0},
        {3, -1, 2}
    });
    const Vector<double> vector_f({
        2,
        -2,
        2
    });

    const Vector<double> vector_x = round(gauss_by_columns(matrix_a, vector_f));

    const Vector<double> expected_vector_x({
        -1,
        1,
        3
    });

    ASSERT_EQ(vector_x, expected_vector_x);
}

TEST(gauss_by_columns, test_3_ok) {
    const Matrix<double> matrix_a({
        {0, 1},
        {2, 4}
    });
    const Vector<double> vector_f({
        3,
        6
    });

    const Vector<double> vector_x = round(gauss_by_columns(matrix_a, vector_f));

    const Vector<double> expected_vector_x({
        -3,
        3
    });

    ASSERT_EQ(vector_x, expected_vector_x);
}

TEST(gauss_by_columns, test_matrix_a_dimension_fail) {
    const Matrix<double> matrix_a({
        {2, 1, 1},
        {1, -1, 0}
    });
    const Vector<double> vector_f({
        2,
        -2,
        0
    });

    ASSERT_ANY_THROW(gauss_by_columns(matrix_a, vector_f));
}

TEST(gauss_by_columns, test_vector_f_dimension_fail) {
    const Matrix<double> matrix_a({
        {2, 1},
        {1, -1}
    });
    const Vector<double> vector_f({
        2,
        -2,
        0
    });

    ASSERT_ANY_THROW(gauss_by_row(matrix_a, vector_f));
}

TEST(gauss, test_1_ok) {
    const Matrix<double> matrix_a({
        {2, 3},
        {5, 7}
    });
    const Vector<double> vector_f({
        11,
        13
    });

    const Vector<double> vector_x = round(gauss(matrix_a, vector_f));

    const Vector<double> expected_vector_x({
        -38,
        29
    });

    ASSERT_EQ(vector_x, expected_vector_x);
}

TEST(gauss, test_2_ok) {
    const Matrix<double> matrix_a({
        {2, 1, 1},
        {1, -1, 0},
        {3, -1, 2}
    });
    const Vector<double> vector_f({
        2,
        -2,
        2
    });

    const Vector<double> vector_x = round(gauss(matrix_a, vector_f));

    const Vector<double> expected_vector_x({
        -1,
        1,
        3
    });

    ASSERT_EQ(vector_x, expected_vector_x);
}

TEST(gauss, test_3_ok) {
    const Matrix<double> matrix_a({
        {0, 1},
        {2, 4}
    });
    const Vector<double> vector_f({
        3,
        6
    });

    const Vector<double> vector_x = round(gauss(matrix_a, vector_f));

    const Vector<double> expected_vector_x({
        -3,
        3
    });

    ASSERT_EQ(vector_x, expected_vector_x);
}

TEST(gauss, test_matrix_a_dimension_fail) {
    const Matrix<double> matrix_a({
        {2, 1, 1},
        {1, -1, 0}
    });
    const Vector<double> vector_f({
        2,
        -2,
        0
    });

    ASSERT_ANY_THROW(gauss(matrix_a, vector_f));
}

TEST(gauss, test_vector_f_dimension_fail) {
    const Matrix<double> matrix_a({
        {2, 1},
        {1, -1}
    });
    const Vector<double> vector_f({
        2,
        -2,
        0
    });

    ASSERT_ANY_THROW(gauss(matrix_a, vector_f));
}

template <class T>
Vector<T> round(const Vector<T> &vector) {
    Vector<T> result(vector);
    for (size_t index = 0; index < result.size(); index++) {
        result(index) = std::round(result(index));
    }
    return result;
}

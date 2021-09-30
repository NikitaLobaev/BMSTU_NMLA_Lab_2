#include "Matrix.cpp"

namespace Lobaev::Math::Gauss {

    template <class T>
    Vector<T> gauss_basic(Matrix<T> matrix_a, Vector<T> vector_f) {
        if (!matrix_a.is_square()) {
            throw "gauss_basic: only square matrices are supported.";
        }
        if (vector_f.size() != matrix_a.columns_count()) {
            throw "gauss_basic: cannot solve SLAU, matrix_a columns count isn't equal to vector_f size.";
        }

        for (size_t diagonal_index = 0; diagonal_index + 1 < matrix_a.rows_count(); diagonal_index++) {
            if (matrix_a(diagonal_index, diagonal_index) == 0) {
                throw "gauss_basic: cannot solve SLAU, divide by 0.";
            }

            for (size_t row_index = diagonal_index + 1; row_index < matrix_a.rows_count(); row_index++) {
                for (size_t column_index = diagonal_index + 1; column_index < matrix_a.columns_count(); column_index++) {
                    matrix_a(row_index, column_index) -= matrix_a(diagonal_index, column_index) *
                                                         matrix_a(row_index, diagonal_index) /
                                                         matrix_a(diagonal_index, diagonal_index);
                }
                vector_f(row_index) -= vector_f(diagonal_index) * matrix_a(row_index, diagonal_index) /
                                       matrix_a(diagonal_index, diagonal_index);
            }
        }

        Vector<T> vector_x = vector_f;
        for (size_t row_index = vector_x.size(); row_index-- > 0;) {
            for (size_t column_index = row_index + 1; column_index < matrix_a.columns_count(); column_index++) {
                vector_x(row_index) -= matrix_a(row_index, column_index) * vector_x(column_index);
            }
            vector_x(row_index) /= matrix_a(row_index, row_index);
        }

        return vector_x;
    }

    template <class T>
    Vector<T> gauss_by_row(Matrix<T> matrix_a, Vector<T> vector_f) {
        if (!matrix_a.is_square()) {
            throw "gauss_by_row: only square matrices are supported.";
        }
        if (vector_f.size() != matrix_a.columns_count()) {
            throw "gauss_by_row: cannot solve SLAU, matrix_a columns count isn't equal to vector_f size.";
        }

        std::vector<size_t> rows_permutation(matrix_a.columns_count());
        for (size_t index = 0; index < rows_permutation.size(); index++) {
            rows_permutation[index] = index;
        }

        for (size_t diagonal_index = 0; diagonal_index + 1 < matrix_a.rows_count(); diagonal_index++) {
            size_t best_column_index = diagonal_index;
            for (size_t column_index = diagonal_index + 1; column_index < matrix_a.columns_count(); column_index++) {
                if (matrix_a(diagonal_index, column_index) != 0 &&
                        (matrix_a(diagonal_index, best_column_index) < matrix_a(diagonal_index, column_index)) ||
                        matrix_a(diagonal_index, best_column_index) == 0) {
                    best_column_index = column_index;
                }
            }

            matrix_a.swap_columns(diagonal_index, best_column_index);
            std::swap(rows_permutation[diagonal_index], rows_permutation[best_column_index]);

            if (matrix_a(diagonal_index, diagonal_index) == 0) {
                throw "gauss_by_row: infinitely many solutions.";
            }

            for (size_t row_index = diagonal_index + 1; row_index < matrix_a.rows_count(); row_index++) {
                for (size_t column_index = diagonal_index + 1; column_index < matrix_a.columns_count(); column_index++) {
                    matrix_a(row_index, column_index) -= matrix_a(diagonal_index, column_index) *
                                                         matrix_a(row_index, diagonal_index) /
                                                         matrix_a(diagonal_index, diagonal_index);
                }
                vector_f(row_index) -= vector_f(diagonal_index) * matrix_a(row_index, diagonal_index) /
                                       matrix_a(diagonal_index, diagonal_index);
            }
        }

        Vector<T> vector_x(vector_f.size());
        for (size_t row_index = vector_x.size(); row_index-- > 0;) {
            vector_x(rows_permutation[row_index]) = vector_f(row_index);

            for (size_t column_index = row_index + 1; column_index < matrix_a.columns_count(); column_index++) {
                vector_x(rows_permutation[row_index]) -= matrix_a(row_index, column_index) *
                                                         vector_x(rows_permutation[column_index]);
            }

            vector_x(rows_permutation[row_index]) /= matrix_a(row_index, row_index);
        }

        return vector_x;
    }

    template <class T>
    Vector<T> gauss_by_columns(Matrix<T> matrix_a, Vector<T> vector_f) {
        if (!matrix_a.is_square()) {
            throw "gauss_by_columns: only square matrices are supported.";
        }
        if (vector_f.size() != matrix_a.columns_count()) {
            throw "gauss_by_columns: cannot solve SLAU, matrix_a columns count isn't equal to vector_f size.";
        }

        for (size_t diagonal_index = 0; diagonal_index + 1 < matrix_a.rows_count(); diagonal_index++) {
            size_t best_row_index = diagonal_index;
            for (size_t row_index = diagonal_index + 1; row_index < matrix_a.rows_count(); row_index++) {
                if (matrix_a(row_index, diagonal_index) != 0 &&
                        (matrix_a(best_row_index, diagonal_index) < matrix_a(row_index, diagonal_index)) ||
                        matrix_a(best_row_index, diagonal_index) == 0) {
                    best_row_index = row_index;
                }
            }

            matrix_a.swap_rows(diagonal_index, best_row_index);
            vector_f.swap(diagonal_index, best_row_index);

            if (matrix_a(diagonal_index, diagonal_index) == 0) {
                throw "gauss_by_columns: infinitely many solutions.";
            }

            for (size_t row_index = diagonal_index + 1; row_index < matrix_a.rows_count(); row_index++) {
                for (size_t column_index = diagonal_index + 1; column_index < matrix_a.columns_count(); column_index++) {
                    matrix_a(row_index, column_index) -= matrix_a(diagonal_index, column_index) *
                                                         matrix_a(row_index, diagonal_index) /
                                                         matrix_a(diagonal_index, diagonal_index);
                }
                vector_f(row_index) -= vector_f(diagonal_index) * matrix_a(row_index, diagonal_index) /
                                       matrix_a(diagonal_index, diagonal_index);
            }
        }

        Vector<T> vector_x = vector_f;
        for (size_t row_index = vector_x.size(); row_index-- > 0;) {
            for (size_t column_index = row_index + 1; column_index < matrix_a.columns_count(); column_index++) {
                vector_x(row_index) -= matrix_a(row_index, column_index) * vector_x(column_index);
            }
            vector_x(row_index) /= matrix_a(row_index, row_index);
        }

        return vector_x;
    }

    template <class T>
    Vector<T> gauss(Matrix<T> matrix_a, Vector<T> vector_f) {
        if (!matrix_a.is_square()) {
            throw "gauss: only square matrices are supported.";
        }
        if (vector_f.size() != matrix_a.columns_count()) {
            throw "gauss: cannot solve SLAU, matrix_a columns count isn't equal to vector_f size.";
        }

        std::vector<size_t> rows_permutation(matrix_a.columns_count());
        for (size_t index = 0; index < rows_permutation.size(); index++) {
            rows_permutation[index] = index;
        }

        for (size_t diagonal_index = 0; diagonal_index + 1 < matrix_a.rows_count(); diagonal_index++) {
            size_t best_row_index = diagonal_index;
            size_t best_column_index = diagonal_index;

            for (size_t row_index = diagonal_index + 1; row_index < matrix_a.rows_count(); row_index++) {
                if (matrix_a(row_index, diagonal_index) != 0 &&
                        (matrix_a(best_row_index, diagonal_index) < matrix_a(row_index, diagonal_index)) ||
                        matrix_a(best_row_index, diagonal_index) == 0) {
                    best_row_index = row_index;
                }
            }

            for (size_t column_index = diagonal_index + 1; column_index < matrix_a.columns_count(); column_index++) {
                if (matrix_a(diagonal_index, column_index) != 0 &&
                        (matrix_a(best_row_index, best_column_index) < matrix_a(diagonal_index, column_index)) ||
                        matrix_a(best_row_index, best_column_index) == 0) {
                    best_row_index = diagonal_index;
                    best_column_index = column_index;
                }
            }

            matrix_a.swap_rows(diagonal_index, best_row_index);
            vector_f.swap(diagonal_index, best_row_index);

            matrix_a.swap_columns(diagonal_index, best_column_index);
            std::swap(rows_permutation[diagonal_index], rows_permutation[best_column_index]);

            if (matrix_a(diagonal_index, diagonal_index) == 0) {
                throw "gauss: infinitely many solutions.";
            }

            for (size_t row_index = diagonal_index + 1; row_index < matrix_a.rows_count(); row_index++) {
                for (size_t column_index = diagonal_index + 1; column_index < matrix_a.columns_count(); column_index++) {
                    matrix_a(row_index, column_index) -= matrix_a(diagonal_index, column_index) *
                                                         matrix_a(row_index, diagonal_index) /
                                                         matrix_a(diagonal_index, diagonal_index);
                }
                vector_f(row_index) -= vector_f(diagonal_index) * matrix_a(row_index, diagonal_index) /
                                       matrix_a(diagonal_index, diagonal_index);
            }
        }

        Vector<T> vector_x(vector_f.size());
        for (size_t row_index = vector_x.size(); row_index-- > 0;) {
            vector_x(rows_permutation[row_index]) = vector_f(row_index);

            for (size_t column_index = row_index + 1; column_index < matrix_a.columns_count(); column_index++) {
                vector_x(rows_permutation[row_index]) -= matrix_a(row_index, column_index) *
                                                         vector_x(rows_permutation[column_index]);
            }

            vector_x(rows_permutation[row_index]) /= matrix_a(row_index, row_index);
        }

        return vector_x;
    }

}

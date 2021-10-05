#include <Eigen/Dense>
#include <ctime>
#include <iostream>
#include <fstream>
#include "Gauss.cpp"
#include "generator.cpp"

template <class T, class T2>
T2 delta(const Lobaev::Math::Vector<T>&, const Lobaev::Math::Vector<T>&);

enum GnuplotTask9Graph {
    GAUSS_BASIC,
    GAUSS_BY_ROW,
    GAUSS_BY_COLUMNS,
    GAUSS_BASIC_WITH_DD,
    GAUSS_BY_ROW_WITH_DD,
    GAUSS_BY_COLUMNS_WITH_DD
};

const size_t gnuplot_task_9_graphs_size = 6;
const GnuplotTask9Graph gnuplot_task_9_graphs[gnuplot_task_9_graphs_size]{
    GAUSS_BASIC,
    GAUSS_BY_ROW,
    GAUSS_BY_COLUMNS,
    GAUSS_BASIC_WITH_DD,
    GAUSS_BY_ROW_WITH_DD,
    GAUSS_BY_COLUMNS_WITH_DD
};

const char *usage = "Usage: lab2 <gnuplot task 9 command file> <gnuplot task 9 output file> <selected graph index>";

int main(int argc, char **argv) {
    if (argc != 4) {
        std::cerr << usage << std::endl;
        return 1;
    }

    const std::string gnuplot_task_9_command_filename = argv[1];

    const std::string gnuplot_task_9_output_filename = argv[2];

    const size_t selected_gnuplot_task_9_graph_index = std::stoul(argv[3]);
    if (selected_gnuplot_task_9_graph_index >= gnuplot_task_9_graphs_size) {
        std::cerr << usage << std::endl;
        return 1;
    }

    const GnuplotTask9Graph selected_gnuplot_task_9_graph = gnuplot_task_9_graphs[selected_gnuplot_task_9_graph_index];

    std::srand(std::time(nullptr)); //для генерации случайных матриц с диагональным преобладанием

    { //проверка примера
        const Lobaev::Math::Matrix<int> matrix_a({
            {0, 1},
            {2, 4}
        });
        const Lobaev::Math::Vector<int> vector_f({
            3,
            6
        });

        const Lobaev::Math::Vector<int> solution = Lobaev::Math::Gauss::gauss(matrix_a, vector_f);
    }

    { //проверка на сторонней библиотеке (Eigen3)
        Eigen::Matrix<int, 2, 2> eigen_matrix_a;
        eigen_matrix_a(0, 0) = 0;
        eigen_matrix_a(0, 1) = 1;
        eigen_matrix_a(1, 0) = 2;
        eigen_matrix_a(1, 1) = 4;

        Eigen::Matrix<int, 2, 1> eigen_vector_f;
        eigen_vector_f(0, 0) = 3;
        eigen_vector_f(1, 0) = 6;

        const Eigen::Matrix<int, 2, 1> eigen_vector_x = eigen_matrix_a * eigen_vector_f;
    }

    { //вычисление погрешности вычислений методов Гаусса простого, "по строкам", и "по столбцам"
        std::cout << "--- Task 7 ---" << std::endl;

        const Lobaev::Math::Matrix<int> matrix_a({
            {3, 2, 1, 1},
            {1, -1, 4, -1},
            {-2, -2, -3, 1},
            {1, 5, -1, 2}
        });
        const Lobaev::Math::Vector<int> vector_f({
            -2,
            -1,
            9,
            4
        });

        //эталонный правильный ответ: {-3, -1, 2, 7}

        const Lobaev::Math::Vector<int> solution = Lobaev::Math::Gauss::gauss(matrix_a, vector_f);

        {
            const Lobaev::Math::Vector<int> solution_basic = Lobaev::Math::Gauss::gauss_basic(matrix_a, vector_f);

            std::cout << delta<int, double>(solution_basic, solution) << '%' << std::endl;
        }
        {
            const Lobaev::Math::Vector<int> solution_by_row = Lobaev::Math::Gauss::gauss_by_row(matrix_a, vector_f);

            std::cout << delta<int, double>(solution_by_row, solution) << '%' << std::endl;
        }
        {
            const Lobaev::Math::Vector<int> solution_by_columns = Lobaev::Math::Gauss::gauss_by_columns(matrix_a, vector_f);

            std::cout << delta<int, double>(solution_by_columns, solution) << '%' << std::endl;
        }

        std::cout << "--------------" << std::endl;
    }

    { //генерация случайных матриц с диагональным преобладанием
        std::cout << "--- Task 9 ---" << std::endl;

        Lobaev::Math::Vector<double> (*selected_gauss)(Lobaev::Math::Matrix<double>, Lobaev::Math::Vector<double>);
        switch (selected_gnuplot_task_9_graph) {
            case GAUSS_BASIC:
            case GAUSS_BASIC_WITH_DD:
                selected_gauss = Lobaev::Math::Gauss::gauss_basic;
                break;
            case GAUSS_BY_ROW:
            case GAUSS_BY_ROW_WITH_DD:
                selected_gauss = Lobaev::Math::Gauss::gauss_by_row;
                break;
            default: //GAUSS_BY_COLUMNS, GAUSS_BY_COLUMNS_WITH_DD
                selected_gauss = Lobaev::Math::Gauss::gauss_by_columns;
                break;
        }

        std::ofstream gnuplot_task_9_ofstream(gnuplot_task_9_output_filename);

        const double from = -1000;
        const double to = 1000;

        for (size_t dimension = 2; dimension <= 150; dimension++) {
            if (dimension % 20 == 0) {
                std::cout << "\tdimension = " << dimension << std::endl;
            }

            Lobaev::Math::Matrix<double> matrix_a = gen_random_matrix<double>(dimension, from, to);
            const Lobaev::Math::Vector<double> vector_f = gen_random_vector<double>(dimension, from, to);

            if (selected_gnuplot_task_9_graph == GAUSS_BASIC_WITH_DD ||
                    selected_gnuplot_task_9_graph == GAUSS_BY_ROW_WITH_DD ||
                    selected_gnuplot_task_9_graph == GAUSS_BY_COLUMNS_WITH_DD) {
                make_matrix_diagonally_dominant(matrix_a, from, to);
            }

            try {
                const Lobaev::Math::Vector<double> ref_solution = Lobaev::Math::Gauss::gauss(matrix_a, vector_f);
                const Lobaev::Math::Vector<double> solution = selected_gauss(matrix_a, vector_f);

                const auto delta_value = delta<double, double>(solution, ref_solution);

                gnuplot_task_9_ofstream << dimension << ' ' << delta_value << std::endl;
            } catch (const char *exception) {
            }
        }

        gnuplot_task_9_ofstream.close();

        gnuplot_task_9_ofstream = std::ofstream(gnuplot_task_9_command_filename);
        gnuplot_task_9_ofstream << "set xlabel \"Matrix dimension (n)\"" << std::endl << "set ylabel \"delta\"" <<
        std::endl << "plot \"" << gnuplot_task_9_output_filename << "\" with lines" << std::endl;
        gnuplot_task_9_ofstream.close();

        system(("gnuplot --persist \"" + gnuplot_task_9_command_filename + "\"").c_str());

        std::cout << "--------------" << std::endl;
    }

    return 0;
}

template <class T, class T2>
T2 delta(const Lobaev::Math::Vector<T> &result, const Lobaev::Math::Vector<T> &expected_result) {
    const Lobaev::Math::Vector<T> diff_vector = result - expected_result;

    const T2 diff_norm_euclidean = diff_vector.template norm_euclidean<T2>();

    const T2 result_norm_euclidean = expected_result.template norm_euclidean<T2>();

    return diff_norm_euclidean * 100.0 / result_norm_euclidean;
}

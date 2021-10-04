#include <Eigen/Dense>
#include <ctime>
#include <iostream>
#include <fstream>
#include "Gauss.cpp"
#include "generator.cpp"

template <class T, class T2>
T2 delta(const Lobaev::Math::Vector<T>&, const Lobaev::Math::Vector<T>&);

const size_t gnuplot_task_9_graph_1 = 1,
gnuplot_task_9_graph_2 = 2,
gnuplot_task_9_graph_3 = 3,
gnuplot_task_9_graph_4 = 4,
gnuplot_task_9_graph_5 = 5,
gnuplot_task_9_graph_6 = 6;

const char *usage = "Usage: lab2 <gnuplot task 9 command file> <gnuplot task 9 output file> <selected graph number>";

int main(int argc, char **argv) {
    if (argc != 4) {
        std::cerr << usage << std::endl;
        return 1;
    }

    const std::string gnuplot_task_9_command_filename = argv[1];

    const std::string gnuplot_task_9_output_filename = argv[2];

    const size_t selected_gnuplot_task_9_graph_number = std::stoul(argv[3]);
    if (selected_gnuplot_task_9_graph_number > gnuplot_task_9_graph_6) {
        std::cerr << usage << std::endl;
        return 1;
    }

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
    }

    { //генерация случайных матриц с диагональным преобладанием
        std::ofstream gnuplot_task_9_ofstream(gnuplot_task_9_output_filename);

        const double from = -1000;
        const double to = 1000;
        for (size_t dimension = 2; dimension <= 100; dimension++) {
            Lobaev::Math::Matrix<double> matrix_a = gen_random_matrix<double>(dimension, from, to);
            const Lobaev::Math::Vector<double> vector_f = gen_random_vector<double>(dimension, from, to);

            try {
                const Lobaev::Math::Vector<double> solution_basic = Lobaev::Math::Gauss::gauss_basic(matrix_a, vector_f);
                const Lobaev::Math::Vector<double> solution_by_row = Lobaev::Math::Gauss::gauss_by_row(matrix_a, vector_f);
                const Lobaev::Math::Vector<double> solution_by_columns = Lobaev::Math::Gauss::gauss_by_columns(matrix_a, vector_f);
                const Lobaev::Math::Vector<double> solution = Lobaev::Math::Gauss::gauss(matrix_a, vector_f);

                const double delta_basic = delta<double, double>(solution_basic, solution);
                const double delta_by_row = delta<double, double>(solution_by_row, solution);
                const double delta_by_columns = delta<double, double>(solution_by_columns, solution);

                switch (selected_gnuplot_task_9_graph_number) {
                    case gnuplot_task_9_graph_1:
                        gnuplot_task_9_ofstream << dimension << ' ' << delta_basic << std::endl;
                        break;
                    case gnuplot_task_9_graph_2:
                        gnuplot_task_9_ofstream << dimension << ' ' << delta_by_row << std::endl;
                        break;
                    case gnuplot_task_9_graph_3:
                        gnuplot_task_9_ofstream << dimension << ' ' << delta_by_columns << std::endl;
                        break;
                }
            } catch (const char *exception) {
            }

            make_matrix_diagonally_dominant(matrix_a, from, to);

            try {
                const Lobaev::Math::Vector<double> solution_basic = Lobaev::Math::Gauss::gauss_basic(matrix_a, vector_f);
                const Lobaev::Math::Vector<double> solution_by_row = Lobaev::Math::Gauss::gauss_by_row(matrix_a, vector_f);
                const Lobaev::Math::Vector<double> solution_by_columns = Lobaev::Math::Gauss::gauss_by_columns(matrix_a, vector_f);
                const Lobaev::Math::Vector<double> solution = Lobaev::Math::Gauss::gauss(matrix_a, vector_f);

                const double delta_basic = delta<double, double>(solution_basic, solution);
                const double delta_by_row = delta<double, double>(solution_by_row, solution);
                const double delta_by_columns = delta<double, double>(solution_by_columns, solution);

                switch (selected_gnuplot_task_9_graph_number) {
                    case gnuplot_task_9_graph_4:
                        gnuplot_task_9_ofstream << dimension << ' ' << delta_basic << std::endl;
                        break;
                    case gnuplot_task_9_graph_5:
                        gnuplot_task_9_ofstream << dimension << ' ' << delta_by_row << std::endl;
                        break;
                    case gnuplot_task_9_graph_6:
                        gnuplot_task_9_ofstream << dimension << ' ' << delta_by_columns << std::endl;
                        break;
                }
            } catch (const char *exception) {
            }
        }

        gnuplot_task_9_ofstream.close();

        gnuplot_task_9_ofstream = std::ofstream(gnuplot_task_9_command_filename);
        gnuplot_task_9_ofstream << "set xlabel \"Matrix dimension (n)\"" << std::endl << "set ylabel \"delta\"" <<
        std::endl << "plot \"" << gnuplot_task_9_output_filename << "\" with lines" << std::endl;
        gnuplot_task_9_ofstream.close();

        system(("gnuplot --persist \"" + gnuplot_task_9_command_filename + "\"").c_str());
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

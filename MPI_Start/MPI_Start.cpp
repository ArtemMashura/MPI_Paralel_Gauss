//#include <mpi.h>
//#include <stdio.h>
//#include <stdlib.h>
//#include <math.h>
//#include <iostream>
//#include <locale.h>
//using namespace std;
//#define N 4 // Размер системы уравнений (число неизвестных)
//
//void print_matrix(double matrix[N][N + 1]) {
//    for (int i = 0; i < N; i++) {
//        for (int j = 0; j < N + 1; j++) {
//            printf("%10.4f ", matrix[i][j]);
//        }
//        printf("\n");
//    }
//    printf("\n");
//}
//
//void gauss_elimination(double matrix[N][N + 1], int rank, int size) {
//    for (int k = 0; k < N; k++) {
//        if (rank == 0) {
//            // Нормализуем ведущую строку
//            for (int j = k + 1; j < N + 1; j++) {
//                matrix[k][j] /= matrix[k][k];
//            }
//            matrix[k][k] = 1.0;
//        }
//
//        // Передаем ведущую строку всем процессам
//        MPI_Bcast(matrix[k], N + 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//
//        // Обновляем строки, принадлежащие текущему процессу
//        for (int i = k + 1 + rank; i < N; i += size) {
//            double factor = matrix[i][k];
//            for (int j = k; j < N + 1; j++) {
//                matrix[i][j] -= factor * matrix[k][j];
//            }
//        }
//
//        // Синхронизируем процессы
//        for (int i = 0; i < N; i++) {
//            MPI_Bcast(matrix[i], N + 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//        }
//    }
//}
//
//void back_substitution(double matrix[N][N + 1], double result[N], int rank, int size) {
//    for (int i = N - 1; i >= 0; i--) {
//        if (rank == 0) {
//            result[i] = matrix[i][N];
//            for (int j = i + 1; j < N; j++) {
//                result[i] -= matrix[i][j] * result[j];
//            }
//        }
//        MPI_Bcast(&result[i], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//    }
//}
//
//int main(int argc, char** argv) {
//    setlocale(LC_ALL, "Russian");
//    MPI_Init(&argc, &argv);
//    int rank, size;
//    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//    MPI_Comm_size(MPI_COMM_WORLD, &size);
//
//    double matrix[N][N + 1] = {
//        {2, -1, 1, 3, 8},
//        {1, 3, 2, -4, 4},
//        {4, -1, -2, 2, 2},
//        {3, 1, 3, -1, 6}
//    };
//
//    if (rank == 0) {
//        printf("Ishodnaya matrica:\n");
//        print_matrix(matrix);
//    }
//
//    // Прямой ход метода Гаусса
//    gauss_elimination(matrix, rank, size);
//
//    if (rank == 0) {
//        printf("Matrica posle pryamogo:\n");
//        print_matrix(matrix);
//    }
//
//    // Обратный ход для нахождения неизвестных
//    double result[N] = { 0 };
//    back_substitution(matrix, result, rank, size);
//
//    if (rank == 0) {
//        printf("Solution:\n");
//        for (int i = 0; i < N; i++) {
//            printf("x%d = %f\n", i + 1, result[i]);
//        }
//    }
//
//    MPI_Finalize();
//    return 0;
//}


#include <iostream>
#include <mpi.h>
#include <chrono>
using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::milliseconds;
using namespace std;

int main(int argc, char** argv)
{
    const int size = 100; // ця змінна має дорівнювати число процесів * num_iter, інакше програма не працює

    int size_proc, rank;
    MPI_Status status;
    MPI_Request request;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    const int num_iter = 10; // на жаль, я не можу ставити цю змінну за користувача. Деталі: https://shorturl.at/vw6is
                             // коротко: я міг би ставити значення розмірів масивів на константу яка явно не визначена, якби проект був на с, а не с++
                             // а я не хочу створювати новий проект тому що не хочу заново підключати mpi.
    double matrix[size][size];
    double B[size];
    double E[size][size];

    double segmentM[num_iter][size];
    double segmentE[num_iter][size];

    double segmentMK[size];
    double segmentEK[size];

    double mini_segmentM[num_iter];
    double mini_segmentE[num_iter];

    // Заповнення матриць
    if (rank == 0)
    {
        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < size; j++)
            {
                matrix[i][j] = rand() % 10;
                cout << matrix[i][j] <<"x" << j+1 ;
                if (j + 1 < size) {
                    cout << " + ";
                }
                if (i == j) E[i][j] = 1.0;
                else E[i][j] = 0.0;
            }
            B[i] = rand() % 10;
            cout << " = " << B[i] << endl;
        }
    }

    // Прямий хід
    auto start = high_resolution_clock::now();

    double div, multi;
    for (int k = 0; k < size; k++)
    {
        if (rank == 0)
        {
            if (matrix[k][k] == 0.0)
            {
                bool changed = false;
                for (int i = k + 1; i < size; i++)
                {
                    if (matrix[i][k] != 0)
                    {
                        swap(matrix[k], matrix[i]);
                        swap(E[k], E[i]);
                        changed = true;
                        break;
                    }
                }
                if (!changed)
                {
                    cout << endl << "Error: матрица не может быть найдена" << endl;
                    return -1;
                }
            }
            div = matrix[k][k];
        }

        MPI_Scatter(matrix[k], num_iter, MPI_DOUBLE, mini_segmentM, num_iter, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Scatter(E[k], num_iter, MPI_DOUBLE, mini_segmentE, num_iter, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&div, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        for (int j = 0; j < num_iter; j++)
        {
            mini_segmentM[j] /= div;
            mini_segmentE[j] /= div;
        }

        MPI_Gather(mini_segmentM, num_iter, MPI_DOUBLE, matrix[k], num_iter, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Gather(mini_segmentE, num_iter, MPI_DOUBLE, E[k], num_iter, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        if (rank == 0)
        {
            for (int i = 0; i < size; i++)
            {
                segmentMK[i] = matrix[k][i];
                segmentEK[i] = E[k][i];
            }
        }

        MPI_Bcast(segmentMK, size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(segmentEK, size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        MPI_Scatter(matrix, num_iter * size, MPI_DOUBLE, segmentM, num_iter * size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Scatter(E, num_iter * size, MPI_DOUBLE, segmentE, num_iter * size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        for (int i = 0; i < num_iter; i++)
        {
            if ((rank * num_iter) + i <= k)
                continue;

            multi = segmentM[i][k];
            for (int j = 0; j < size; j++)
            {
                segmentM[i][j] -= multi * segmentMK[j];
                segmentE[i][j] -= multi * segmentEK[j];
            }
        }

        MPI_Gather(segmentM, size * num_iter, MPI_DOUBLE, matrix, size * num_iter, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Gather(segmentE, size * num_iter, MPI_DOUBLE, E, size * num_iter, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    //Зворотній хід
    for (int k = size - 1; k > 0; k--)
    {
        if (rank == 0)
        {
            for (int i = 0; i < size; i++)
            {
                segmentMK[i] = matrix[k][i];
                segmentEK[i] = E[k][i];
            }
        }

        MPI_Bcast(segmentMK, size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(segmentEK, size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Scatter(matrix, num_iter * size, MPI_DOUBLE, segmentM, num_iter * size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Scatter(E, num_iter * size, MPI_DOUBLE, segmentE, num_iter * size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        for (int i = num_iter - 1; i > -1; i--)
        {
            if ((rank * num_iter) + i >= k)
                continue;

            multi = segmentM[i][k];
            for (int j = 0; j < size; j++)
            {
                segmentM[i][j] -= multi * segmentMK[j];
                segmentE[i][j] -= multi * segmentEK[j];
            }
        }

        MPI_Gather(segmentM, size * num_iter, MPI_DOUBLE, matrix, size * num_iter, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Gather(segmentE, size * num_iter, MPI_DOUBLE, E, size * num_iter, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
    /*if (rank == 0)
    {
        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < size; j++)
                cout << E[i][j] << ' ';

            cout << '\n';
        }
    }*/

    double X[size];
    double segmentX[num_iter];

    MPI_Bcast(B, size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(E, num_iter * size, MPI_DOUBLE, segmentE, num_iter * size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    for (int i = 0; i < num_iter; i++)
    {
        segmentX[i] = 0;
        for (int j = 0; j < size; j++)
            segmentX[i] += segmentE[i][j] * B[j];
    }

    MPI_Gather(segmentX, num_iter, MPI_DOUBLE, X, num_iter, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    /*if (rank == 0)
    {
        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < size; j++)
                cout << matrix[i][j] << " ";
            cout << endl;
        }
    }*/

    if (rank == 0)
    {
        for (int i = 0; i < size; i++)
            cout << "\nx" << i + 1 << " = " << X[i];

        auto end = high_resolution_clock::now();
        duration<double, std::milli> ms_double = end - start;
        std::cout << endl << ms_double.count() << "ms\n";
    }
    MPI_Finalize();
    return 0;
}
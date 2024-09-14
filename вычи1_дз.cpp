#include <fstream>
#include <iostream>
using namespace std;

//часть а - трехмерный массив
void task_a() {
    const int size1 = 20;
    const int size2 = 30;
    const int size3 = 40;
    int last_iter = 100;
    double h1 = 1.0;
    double h2 = 1.0;
    double h3 = 1.0;
    double tau = 0.01;
    double D = 1.0;
    double k1 = (D * tau) / (h1 * h1);
    double k2 = (D * tau) / (h2 * h2);
    double k3 = (D * tau) / (h3 * h3);
    double curr_conc[size1][size2][size3];
    double next_conc[size1][size2][size3];
    //чтение файла
    ifstream initial("input.raw", ios_base::binary);
    initial.read((char*)curr_conc, sizeof(double) * size1 * size2 * size3);
    initial.close();
    //граничные условия
    for (size_t i = 0; i < size1; ++i) {
        for (size_t j = 0; j < size2; ++j) {
            next_conc[i][j][0] = 0.0;
            next_conc[i][j][size3 - 1] = 0.0;
        }
    }
    for (size_t j = 0; j < size2; ++j) {
        for (size_t l = 1; l < size3 - 1; ++l) {
            next_conc[0][j][l] = 0.0;
            next_conc[size1 - 1][j][l] = 0.0;
        }
    }
    for (size_t l = 1; l < size3 - 1; ++l) {
        for (size_t i = 1; i < size1 - 1; ++i) {
            next_conc[i][0][l] = 0.0;
            next_conc[i][size2 - 1][l] = 0.0;
        }
    }
    //численное решение
    for (int iter = 0; iter < last_iter; ++iter) {
        for (size_t i = 1; i < size1 - 1; ++i) {
            for (size_t j = 1; j < size2 - 1; ++j) {
                for (size_t l = 1; l < size3 - 1; ++l) {
                    next_conc[i][j][l] = (1 - 2 * k1 - 2 * k2 - 2 * k3) * curr_conc[i][j][l] 
                        + k1 * (curr_conc[i + 1][j][l] + curr_conc[i - 1][j][l])
                        + k2 * (curr_conc[i][j + 1][l] + curr_conc[i][j - 1][l])
                        + k3 * (curr_conc[i][j][l + 1] + curr_conc[i][j][l - 1]);
                }
            }
        }
        swap(next_conc, curr_conc);
    }
    //запись в файл
    ofstream fin("output.raw", ios_base::binary);
    fin.write((char*)curr_conc, sizeof(double) * size1 * size2 * size3);
    fin.close();

    return;
}
//часть б - одномерный массив
void task_b() {
    const int size1 = 20;
    const int size2 = 30;
    const int size3 = 40;
    int last_iter = 100;
    double h1 = 1.0;
    double h2 = 1.0;
    double h3 = 1.0;
    double tau = 0.01;
    double D = 1.0;
    double k1 = (D * tau) / (h1 * h1);
    double k2 = (D * tau) / (h2 * h2);
    double k3 = (D * tau) / (h3 * h3);
    double curr_conc[size1 * size2 * size3];
    double next_conc[size1 * size2 * size3];
    //чтение файла
    ifstream initial("input.raw", ios_base::binary);
    initial.read((char*)curr_conc, sizeof(double) * size1 * size2 * size3);
    initial.close();
    //граничные условия
    size_t ind;
    for (size_t i = 0; i < size1; ++i) {
        for (size_t j = 0; j < size2; ++j) {
            ind = (i * size2 + j) * size3;
            next_conc[ind] = 0.0;
            next_conc[ind + size3 - 1] = 0.0;
        }
    }
    for (size_t j = 0; j < size2; ++j) {
        for (size_t l = 1; l < size3 - 1; ++l) {
            ind = j * size3 + l;
            next_conc[ind] = 0.0;
            next_conc[(size1 - 1) * size2 * size3 + ind] = 0.0;
        }
    }
    for (size_t l = 1; l < size3 - 1; ++l) {
        for (size_t i = 1; i < size1 - 1; ++i) {
            ind = i * size2 * size3 + l;
            next_conc[ind] = 0.0;
            next_conc[ind + (size2 - 1) * size3] = 0.0;
        }
    }
    //численное решение
    for (int iter = 0; iter < last_iter; ++iter) {
        for (size_t i = 1; i < size1 - 1; ++i) {
            for (size_t j = 1; j < size2 - 1; ++j) {
                for (size_t l = 1; l < size3 - 1; ++l) {
                    ind = (i * size2 + j) * size3 + l;
                    next_conc[ind] = (1 - 2 * k1 - 2 * k2 - 2 * k3) * curr_conc[ind]
                        + k1 * (curr_conc[ind + size2 * size3] + curr_conc[ind - size2 * size3])
                        + k2 * (curr_conc[ind + size3] + curr_conc[ind - size3])
                        + k3 * (curr_conc[ind + 1] + curr_conc[ind - 1]);
                }
            }
        }
        swap(next_conc, curr_conc);
    }
    //запись в файл
    ofstream fin("output.raw", ios_base::binary);
    fin.write((char*)curr_conc, sizeof(double) * size1 * size2 * size3);
    fin.close();

    return;
}

int main() {

    //часть а
    task_a();
    
    //часть б
    //task_b();

    return 0;
} 
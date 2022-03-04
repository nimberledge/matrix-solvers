#include "Matrix.h"
#include "Matrix.cpp"
#include "CSRMatrix.h"
#include "CSRMatrix.cpp"
#include <iostream>
#include <memory>
#include <string>
#include <chrono>

using std::cout;
using std::endl;
using std::string;

// Time a run of jacobi based on input matrix in filename with size n
void runJacobiBenchmark(string filename, unsigned int n) {
    cout << "Running Jacobi benchmark for square matrix size: " << n << endl;
    Matrix<double> m(filename);
    double b[m.rows], res[m.rows];
    double tol = 1e-6;
    for (unsigned int i = 0; i < m.rows; i++) {
        b[i] = (double) (i+1) / (n*n);
        res[i] = b[i];
    }
    auto start = std::chrono::steady_clock::now();
    m.jacobiSolve(b, tol, res);
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    cout << "elapsed time: " << elapsed_seconds.count() << "s" << endl;
}

// Time a run of cholesky based on input matrix in filename with size n
void runCholeskyBenchmark(string filename, unsigned int n) {
    cout << "Running dense Cholesky benchmark for square matrix size: " << n << endl;
    Matrix<double> m(filename);
    double b[m.rows], res[m.rows];
    for (unsigned int i = 0; i < m.rows; i++) {
        b[i] = i+1;
    }
    auto start = std::chrono::steady_clock::now();
    m.choleskySolve(b, res);
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    cout << "elapsed time: " << elapsed_seconds.count() << "s" << endl;
}

// Time a run of sparse Jacobi based on input matrix in filename with size n
void runSparseJacobiBenchmark(string filename, unsigned int n) {
    cout << "Running sparse Jacobi benchmark for sparse matrix size: " << n << endl;
    CSRMatrix<double> m(filename);
    double b[m.rows], res[m.rows];
    for (unsigned int i = 0; i < m.rows; i++) {
        b[i] = i+1;
        res[i] = b[i];
    }
    auto start = std::chrono::steady_clock::now();
    double tol = 1e-6;
    m.jacobiSolve(b, tol, res);
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    cout << "elapsed time: " << elapsed_seconds.count() << "s" << endl;
}

// Time a run of the sparse CG Method 
void runSparseCGBenchmark(string filename, unsigned int n) {
    cout << "Running sparse Conjugate Gradient benchmark for sparse matrix size: " << n << endl;
    CSRMatrix<double> m(filename);
    double b[m.rows], res[m.rows];
    for (unsigned int i = 0; i < m.rows; i++) {
        b[i] = i+1;
        res[i] = b[i];
    }
    auto start = std::chrono::steady_clock::now();
    double tol = 1e-6;
    m.conjGradSolve(b, tol, res);
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    cout << "elapsed time: " << elapsed_seconds.count() << "s" << endl;
}

int main() {
    // Use matrices generated using the python script in tools/
    // Benchmark cholesky
    string cfNames[] = {"data/choltest10by10.txt",
        "data/choltest15by15.txt",
        "data/choltest20by20.txt",
        "data/choltest40by40.txt",
        "data/choltest60by60.txt",
        "data/choltest100by100.txt",
        "data/choltest200by200.txt",
        "data/choltest400by400.txt"};
    unsigned int Ns[] = {10, 15, 20, 40, 60, 100, 200, 400};
    for (unsigned int i = 0; i < 8; i++) {
        runCholeskyBenchmark(cfNames[i], Ns[i]);
    }

    string jfNames[] = {"data/jactest10by10.txt",
        "data/jactest20by20.txt",
        "data/jactest40by40.txt",
        "data/jactest100by100.txt",
        "data/jactest200by200.txt"};
    unsigned int jNs[] = {10, 20, 40, 100, 200};
    for (unsigned int i = 0; i < 5; i++) {
        runJacobiBenchmark(jfNames[i], jNs[i]);
    }
    string sfNames[] = {"data/sparsetest10by10.txt",
        "data/sparsetest20by20.txt",
        "data/sparsetest40by40.txt",
        "data/sparsetest100by100.txt",
        "data/sparsetest200by200.txt"};
    unsigned int sNs[] = {10, 20, 40, 100, 200};
    for (unsigned int i = 0; i < 5; i++) {
        runSparseJacobiBenchmark(sfNames[i], sNs[i]);
    }

    string cgfNames[] = {"data/sparsetest10by10.txt",
        "data/sparsetest20by20.txt",
        "data/sparsetest40by40.txt",
        "data/sparsetest100by100.txt",
        "data/sparsetest200by200.txt"};
    unsigned int cgNs[] = {10, 20, 40, 100, 200};
    for (unsigned int i = 0; i < 5; i++) {
        runSparseCGBenchmark(cgfNames[i], cgNs[i]);
    }

    return 0;
}
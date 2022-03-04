#include "Matrix.h"
#include "Matrix.cpp"
#include "CSRMatrix.h"
#include "CSRMatrix.cpp"
#include <iostream>
#include <string>

using std::cout;
using std::endl;
using std::string;

int main() {
    cout << "Here's a demonstration of the LADS solverTM(patent pending) working on 20x20 matrices" << endl;
    Matrix<double> m1("data/choltest20by20.txt");
    cout << "Reading a matrix from data/choltest20by20.txt" << endl;
    cout << "Wow, would you look at that? It's symmetric and positive definite" << endl;
    cout << "We solve the linear system Ax=b for x, where b[i] = i+1" << endl;
    double b[m1.rows];
    for (unsigned int i = 0; i < m1.rows; i++) {
        b[i] = i+1;
    }
    double res[m1.rows];
    m1.choleskySolve(b, res);
    cout << "Using a Cholesky solver, we get result vector x = " << endl;
    cout << "[";
    for (unsigned int i = 0; i < m1.rows; i++) {
        if (i != m1.rows - 1) {
            cout << res[i] << ", ";
        }
        else {
            cout << res[i];
        }
    }
    cout << "]" << endl;
    
    cout << "\nLet's see another one, then shall we?" << endl;
    cout << "Reading a matrix from data/jactest20by20.txt" << endl;
    Matrix<double> m2("data/jactest20by20.txt");
    cout << "Using a Jacobi solver now, same b: " << endl;
    double tol = 1e-7;
    for (unsigned int i = 0; i < m1.rows; i++) {
        b[i] = i+1;
        res[i] = b[i];
    }
    m2.jacobiSolve(b, tol, res);
    cout << "[";
    for (unsigned int i = 0; i < m1.rows; i++) {
        if (i != m2.rows - 1) {
            cout << res[i] << ", ";
        }
        else {
            cout << res[i];
        }
    }
    cout << "]" << endl;
    
    cout << "\nI know what you're thinking, what if my matrix is mostly zeros?" << endl;
    cout << "We got you covered, reading now from data/sparsetest20by20.txt" << endl;
    CSRMatrix<double> m3("data/sparsetest20by20.txt");
    cout << "It's only got " << m3.nnz << " non-zeros out of " << m3.rows * m3.cols << " entries!" << endl;
    for (unsigned int i = 0; i < m3.rows; i++) {
        b[i] = i+1;
        res[i] = b[i];
    }
    cout << "Solving Ax=b, b[i] = i+1:" << endl;
    m3.jacobiSolve(b, tol, res);
    cout << "[";
    for (unsigned int i = 0; i < m1.rows; i++) {
        if (i != m3.rows - 1) {
            cout << res[i] << ", ";
        }
        else {
            cout << res[i];
        }
    }
    cout << "]" << endl;

    cout << "\nBut wait, what's that I heard about a Krylov subspace?" << endl;
    cout << "Yeah no worries, we got a conjugate gradient solver too, it'll do the job" << endl;
    cout << "Solving the same system with the same nice SPD matrix: " << endl;
    m3.conjGradSolve(b, tol, res);
    cout << "[";
    for (unsigned int i = 0; i < m1.rows; i++) {
        if (i != m3.rows - 1) {
            cout << res[i] << ", ";
        }
        else {
            cout << res[i];
        }
    }
    cout << "]" << endl;
    return 0;
}
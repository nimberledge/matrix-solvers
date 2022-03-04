#include "Matrix.h"
#include "Matrix.cpp"
#include "CSRMatrix.h"
#include "CSRMatrix.cpp"
#include <iostream>
#include <memory>

using std::cout;
using std::endl;

// Test if we can initialize a CSR Matrix correctly
int testCSRMatrixInit() {
    unsigned int rows = 4;
    unsigned int cols = 4;
    Matrix<double> m(rows, cols, false);
    for (unsigned int i = 0; i < rows; i++) {
        for (unsigned int j = 0; j < cols; j++) {
            if (i == j) {
                m.values[i*cols + j] = i+1;
            } else {
                m.values[i*cols + j] = 0.0;
            }
        }
    }
    CSRMatrix<double> csrM(m);
    double tolerance = 1e-9;
    for (unsigned int i = 0; i < rows; i++) {
        for (unsigned int j = 0; j < cols; j++) {
            if (std::abs(csrM.at(i, j) - m.values[i*cols + j]) > tolerance) {
                return 1;
            }
        }
    }
    cout << "testCSRMatrixInit() " << endl;
    csrM.printMatrix();
    return 0;
}

// Test lower triangular solve
int testCSRLowerTriangularSolve() {
    cout << "testCSRLowerTriangularSolve(): ";
    double b[] = {1.0, 2.0, 3.0};
    unsigned int rowarr[] = {0, 1, 2, 4};
    unsigned int colarr[] = {0, 1, 0, 2};
    double valarr[] = {1.0, 2.0, 3.0, 4.0};
    CSRMatrix<double> csrM(3, 3, 4, rowarr, colarr, valarr);
    double expected[] = {1.0, 1.0, 0.0};
    double tolerance = 1e-8;
    double res[3];
    double res2[3];
    csrM.solveDenseTriangular(b, res, true);
    csrM.solveSparseTriangular(b, res2, true);
    for (unsigned int i = 0; i < 3; i++) {
        if (std::abs(expected[i] - res[i]) > tolerance) {
            return 1;
        }
        if (std::abs(expected[i] - res2[i]) > tolerance) {
            return 1;
        }
    }
    cout << "It's good!" << endl;
    return 0;
}

// Test upper triangular solve
int testCSRDenseUpperTriangularSolve() {
    cout << "testCSRUpperTriangularSolve(): ";
    double b[] = {1.0, 2.0, 3.0};
    unsigned int rowarr[] = {0, 2, 3, 4};
    unsigned int colarr[] = {0, 2, 1, 2};
    double valarr[] = {1.0, 3.0, 2.0, 4.0};
    CSRMatrix<double> csrM(3, 3, 4, rowarr, colarr, valarr);
    double expected[] = {-1.25, 1.0, 0.75};
    double tolerance = 1e-8;
    double res[3];
    csrM.solveDenseTriangular(b, res, false);
    for (unsigned int i = 0; i < 3; i++) {
        if (std::abs(expected[i] - res[i]) > tolerance) {
            return 1;
        }
    }
    cout << "It's good!" << endl;
    return 0;
}

// Test upper triangular solve
int testCSRSparseUpperTriangularSolve() {
    cout << "testCSRSparseUpperTriangularSolve(): ";
    double b[] = {1.0, 1.0, 0.0, 0.0};
    unsigned int rowarr[] = {0, 2, 3, 4, 5};
    unsigned int colarr[] = {0, 3, 1, 2, 3};
    double valarr[] = {1.0, 2.0, 3.0, 5.0, 6.0};
    CSRMatrix<double> csrM(4, 4, 5, rowarr, colarr, valarr);
    double expected[] = {1.0, 0.33333333333, 0.0, 0.0};
    double tolerance = 1e-8;
    double res[4];
    csrM.solveSparseTriangular(b, res, false);
    for (unsigned int i = 0; i < 4; i++) {
        if (std::abs(expected[i] - res[i]) > tolerance) {
            return 1;
        }
    }
    cout << "It's good!" << endl;
    return 0;
}

// Test upper triangular solve
int testCSRSparseUpperTriangularSolve2() {
    cout << "testCSRSparseUpperTriangularSolve2(): ";
    double b[] = {1.0, 0.0, 3.0};
    unsigned int rowarr[] = {0, 2, 3, 4};
    unsigned int colarr[] = {0, 2, 1, 2};
    double valarr[] = {1.0, 3.0, 2.0, 4.0};
    CSRMatrix<double> csrM(3, 3, 4, rowarr, colarr, valarr);
    double expected[] = {-1.25, 0.0, 0.75};
    double tolerance = 1e-8;
    double res[3];
    csrM.solveSparseTriangular(b, res, false);
    for (unsigned int i = 0; i < 3; i++) {
        if (std::abs(expected[i] - res[i]) > tolerance) {
            return 1;
        }
    }
    cout << "It's good!" << endl;
    return 0;
}

int testCSRMatrixWrite() {
    cout << "testCSRMatrixWrite(): ";
    unsigned int rowarr[] = {0, 2, 3, 4};
    unsigned int colarr[] = {0, 2, 1, 2};
    double valarr[] = {1.0, 3.0, 2.0, 4.0};
    CSRMatrix<double> csrM(3, 3, 4, rowarr, colarr, valarr);

    if (csrM.writeToFile("data/testCSRMatrix.txt")) {
        cout << "It's good!" << endl;
        return 0;
    }
    return 1;
}

int testCSRMatrixRead() {
    cout << "testCSRMatrixRead()" << endl;
    CSRMatrix<double> csrM("data/testCSRMatrix.txt");
    Matrix<double> m(3, 3, false);
    for (unsigned int i = 0; i < 3*3; i++) {
        m.values[i] = 0.0;
    }
    m.values[0] = 1.0;
    m.values[2] = 3.0;
    m.values[4] = 2.0;
    m.values[8] = 4.0;
    double tolerance = 1e-8;
    for (unsigned int i = 0; i < 3; i++) {
        for (unsigned int j = 0; j < 3; j++) {
            if (std::abs(m.values[i*3 + j] - csrM.at(i, j)) > tolerance) {
                return 1;
            }
        }
    }
    if (csrM.nnz != 4) {
        return 1;
    }
    cout << "Printing matrix that was read from: " << "data/testCSRMatrix.txt" << endl;
    csrM.printMatrix();
    cout << "It's good!" << endl;
    return 0;
}

int testJacobiSolve() {
    cout << "testJacobiSolve(): ";
    double b[] = {1.0, 2.0, 3.0};
    double expected[] = {-5.0, 1.0, 2.0};
    double res[] = {0, 0, 0};
    unsigned int rowarr[] = {0, 2, 3, 5};
    unsigned int colarr[] = {0, 2, 1, 0, 2};
    double valarr[] = {1.0, 3.0, 2.0, 1.0, 4.0};
    CSRMatrix<double> csrM(3, 3, 5, rowarr, colarr, valarr);
    double tol = 1e-6;
    double testTol = 1e-5;
    csrM.jacobiSolve(b, tol, res);
    for (unsigned int i = 0; i < 3; i++) {
        if (std::abs(res[i] - expected[i]) > testTol) {
            return 1;
        }
    }
    cout << "It's good!" << endl;
    return 0;
}

int testCholeskySolve() {
    cout << "testCholeskySolve(): " << endl;;
    double b[] = {1.0, 2.0, 3.0};
    unsigned int rowarr[] = {0, 1, 3, 5};
    unsigned int colarr[] = {0, 1, 2, 1, 2};
    double valarr[] = {3.0, 3.0, -1.0, -1.0, 3.0};
    CSRMatrix<double> csrM(3, 3, 5, rowarr, colarr, valarr);
    double expected[] = {1.75, 2.0, 2.25};
    double res[3];
    csrM.choleskySolve(b, res);
    cout << "it fuckin sucks " << endl;
    return 0;
}

int testConjGradSolve() {
    cout << "testConjGradSolve(): ";
    double b[] = {1.0, 2.0, 3.0};
    double expected[] = {0.75, (double) (2) / 3, 1.25};
    double res[] = {0, 0, 0};
    unsigned int rowarr[] = {0, 2, 3, 5};
    unsigned int colarr[] = {0, 2, 1, 0, 2};
    double valarr[] = {3.0, -1.0, 3.0, -1.0, 3.0};
    CSRMatrix<double> csrM(3, 3, 5, rowarr, colarr, valarr);
    double tol = 1e-6;
    double testTol = 1e-4;
    csrM.conjGradSolve(b, tol, res);
    for (unsigned int i = 0; i < 3; i++) {
        // cout << "res[" << i << "]: " << res[i] << endl;
        if (std::abs(res[i] - expected[i]) > testTol) {
            return 1;
        }
    }
    cout << "It's good!" << endl;
    return 0;
}

int testConjGradSolve2() {
    cout << "testConjGradSolve2(): ";
    double b[10];
    for (unsigned int i = 0; i < 10; i++) {
        b[i] = i+1;
    }
    CSRMatrix<double> csrM("data/sparsetest10by10.txt");
    double res[10];
    double tol = 1e-4;
    double testTol = 1e-3;
    csrM.conjGradSolve(b, tol, res);
    csrM.matVecMult(res, b);
    for (unsigned int i = 0; i < 10; i++) {
        if (std::abs(b[i] - i - 1) > testTol) {
            return 1;
        }
    }
    cout << "It's good!" << endl;
    return 0;
}

int main() {
    int err = testCSRMatrixInit();
    err = err || testCSRLowerTriangularSolve();
    err = err || testCSRDenseUpperTriangularSolve();
    err = err || testCSRSparseUpperTriangularSolve();
    err = err || testCSRSparseUpperTriangularSolve2();
    err = err || testJacobiSolve();
    err = err || testCSRMatrixWrite();
    err = err || testCSRMatrixRead();
    err = err || testConjGradSolve();
    err = err || testConjGradSolve2();
    if (err == 0) {
        cout << "All tests passed" << endl;
    } else {
        cout << "Errorrr" << endl;
    }
    return err;
}
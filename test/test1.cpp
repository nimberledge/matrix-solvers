#include "Matrix.h"
#include "Matrix.cpp"
#include <iostream>
#include <memory>

using std::cout;
using std::endl;

// Test if we can print when matrix is not pre-alloc'd
int testPrint() {
    unsigned int rows = 4;
    unsigned int cols = 4;
    Matrix<double> *m = new Matrix<double>(rows, cols, false);
    for (unsigned int i = 0; i < rows * cols; i++) {
        m->values[i] = i+1;
    }
    cout << "testPrint()" << endl;
    m->printMatrix();
    m->printValues();
    delete m;
    return 0;
}

// Test if we can print when matrix is pre-alloc'd
int testPrintPrealloc() {
    unsigned int rows = 4;
    unsigned int cols = 4;
    double *vals = new double[rows * cols];
    Matrix<double> *m = new Matrix<double>(rows, cols, vals);
    for (unsigned int i = 0; i < rows * cols; i++) {
        m->values[i] = i+1;
    }
    cout << "testPrintPrealloc()" << endl;
    m->printMatrix();
    m->printValues();
    delete m;
    return 0;
}

// Test if we can print when matrix is statically instantiated
int testPrintStack() {
    unsigned int rows = 4;
    unsigned int cols = 4;
    Matrix<double> m(rows, cols, false);
    for (unsigned int i = 0; i < rows * cols; i++) {
        m.values[i] = i+1;
    }
    cout << "testPrintStack()" << endl;
    m.printMatrix();
    m.printValues();
    return 0;
}
 
 // Test solving Ax=b using our Jacobi solver
int testJacobiSolver() {
    unsigned int rows = 3;
    unsigned int cols = 3;
    double b[rows];
    double res[rows];

    Matrix<double> m(rows, cols, false);
    m.values[0] = 1.0;
    m.values[1] = 2.0;
    m.values[2] = 3.0;
    m.values[3] = 0;
    m.values[4] = 5.0;
    m.values[5] = 6.0;
    m.values[6] = 0;
    m.values[7] = 0;
    m.values[8] = 7.0;

    for (int i = 0; i < rows; i++) {
        b[i] = i+1;
        res[i] = 0;
    }
    m.jacobiSolve(b, 1e-6, res);
    double expected[] = {-0.05714286, -0.114286, 0.4285714};
    double testTolerance = 1e-6;
    cout << "testJacobi()" << endl;
    for (int i = 0; i < rows; i++) {
        if (std::abs(res[i] - expected[i]) > testTolerance) {
            return 1;
        }
    }
    return 0;
}  

// Test a simple 3x3 matVecMult
int testMatVecMult() {
    double b[] = {1.0, 2.0, 3.0};
    double expected[] = {-8.0, 4.0, 11.4};
    Matrix<double> m(3, 3, false);
    m.values[0] = 1.0;
    m.values[1] = 0.0;
    m.values[2] = -3.0;
    m.values[3] = 2.0;
    m.values[4] = 1.0;
    m.values[5] = 0.0;
    m.values[6] = 3.2;
    m.values[7] = 4.1;
    m.values[8] = 0.0;
    double res[3];
    m.matVecMult(b, res);
    double tolerance = 1e-8;
    cout << "testMatVecMult()" << endl;
    for (unsigned int i = 0; i < 3; i++) {
        if (std::abs(expected[i] - res[i]) > tolerance) {
            return 1;
        }
    }
    return 0;
}

// Test lower triangular solve
int testLowerTriangularSolve() {
    double b[] = {1.0, 2.0, 3.0};
    Matrix<double> m(3, 3, false);
    m.values[0] = 1.0;
    m.values[1] = 0.0;
    m.values[2] = 0.0;
    m.values[3] = 2.0;
    m.values[4] = 3.0;
    m.values[5] = 0.0;
    m.values[6] = 4.0;
    m.values[7] = 5.0;
    m.values[8] = 6.0;
    double expected[] = {1.0, 0.0, -0.1666666666666};
    double tolerance = 1e-8;
    double res[3];
    m.triangularSolve(b, res, true);
    cout << "testLowerTriangularSolve()" << endl;
    for (unsigned int i = 0; i < 3; i++) {
        if (std::abs(expected[i] - res[i]) > tolerance) {
            return 1;
        }
    }
    return 0;
}

// Test upper triangular solve
int testUpperTriangularSolve() {
    double b[] = {1.0, 2.0, 3.0};
    Matrix<double> m(3, 3, false);
    m.values[0] = 1.0;
    m.values[1] = 2.0;
    m.values[2] = 4.0;
    m.values[3] = 0.0;
    m.values[4] = 3.0;
    m.values[5] = 5.0;
    m.values[6] = 0.0;
    m.values[7] = 0.0;
    m.values[8] = 6.0;
    double expected[] = {-0.66666666, -0.166666666, 0.5};
    double tolerance = 1e-8;
    double res[3];
    m.triangularSolve(b, res, false);
    cout << "testUpperTriangularSolve()" << endl;
    for (unsigned int i = 0; i < 3; i++) {
        if (std::abs(expected[i] - res[i]) > tolerance) {
            return 1;
        }
    }
    return 0;
}

// Test the cholesky solver on a SPD matrix
int testCholeskySolve() {
    double b[] = {1.0, 2.0, 3.0};
    Matrix<double> m(3, 3, false);
    m.values[0] = 3.0;
    m.values[1] = -1.0;
    m.values[2] = -1.0;
    m.values[3] = -1.0;
    m.values[4] = 3.0;
    m.values[5] = -1.0;
    m.values[6] = -1.0;
    m.values[7] = -1.0;
    m.values[8] = 3.0;
    double expected[] = {1.75, 2.0, 2.25};
    double tolerance = 1e-8;
    double res[3];
    cout << "testCholeskySolve()" << endl;
    m.choleskySolve(b, res);
    for (unsigned int i = 0; i < 3; i++) {
        if (std::abs(expected[i] - res[i]) > tolerance) {
            return 1;
        }
    }
    return 0;
}

int testMatrixRead() {
    std::string filename = "data/testMatrix.txt";
    cout << "testMatrixRead()" << endl;
    Matrix<double> m(filename);
    double tolerance = 1e-10;
    for (unsigned int i = 0; i < 16; i++) { 
        if (std::abs(m.values[i] - i - 1) > tolerance) {
            return 1;
        }
    }
    if (m.rows != 4 || m.cols != 4) {
        return 1;
    }
    return 0;
}

int testMatrixWrite() {
    cout << "testMatrixWrite(): ";
    Matrix<double> m(4, 4, false);
    for (unsigned int i = 0; i < 16; i++) {
        m.values[i] = i+1;
    }

    if (m.writeToFile("data/testMatrix.txt")) {
        cout << "It's good!" << endl;
        return 0;
    }
    return 1;
}

int main() {
    int err = testPrint();
    err = err || testPrintPrealloc();
    err = err || testPrintStack();
    err = err || testJacobiSolver();
    err = err || testMatVecMult();
    err = err || testLowerTriangularSolve();
    err = err || testUpperTriangularSolve();
    err = err || testCholeskySolve();
    err = err || testMatrixWrite();
    err = err || testMatrixRead();
    if (err == 0) {
        cout << "All tests passed" << endl;
    } else {
        cout << "Errorrr" << endl;
    }
    return err;
}
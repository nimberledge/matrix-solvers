#pragma once
#include <iostream>
#include "Matrix.h"
#include "Matrix.cpp"
#include "CSRMatrix.h"
#include "CSRMatrix.cpp"

using std::string;


class SolverInterface {
    public:
        SolverInterface();
        ~SolverInterface();

        int start();

    private:
        void printMessage(string message);
        int getDenseOrSparse();
        void getMatrixSize(int *rowptr, int *colptr);
        // Returns 0 for Jacobi, 1 for Cholesky, -1 otherwise
        int getMethodDense();
        int getMethodSparse();
        void getVectorb(int rows, double *b);
        void getElement(int i, int j, Matrix<double> *m);
        string getFileName();
        Matrix<double> *getMatrixDense();
        CSRMatrix<double> *getMatrixSparse();
};
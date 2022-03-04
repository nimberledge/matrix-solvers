#pragma once
#include <iostream>
#include <vector>
#include "Matrix.h"
#include <string>
#include <sstream>
#include <memory>

using std::vector;
using std::string;

// Reference - https://epubs.siam.org/doi/pdf/10.1137/1.9780898718881.fm
template<class T>
class CSRMatrix : public Matrix<T> {
    public:
        // Number of non-zero entries
        unsigned int nnz;

        // Row and column indices as specified by any standard CSR format
        std::unique_ptr<unsigned int[]> rowIdxs;
        std::unique_ptr<unsigned int[]> colIdxs;
        /*
        --------Constructors / Destructor-------------
        */
        // Default one
        CSRMatrix(unsigned int rows, unsigned int cols, unsigned int nnz);
        // If you want to pass in the specified format
        // PLEASE NOTE THAT THIS DOES NOT USE THE POINTERS PROVIDED
        // WE USE SMART POINTERS, SO WE DO A DEEP COPY OF THE ARRAYS PROVIDED
        // PLEASE MANAGE YOUR OWN MEMORY, THANKS
        CSRMatrix(unsigned int rows, unsigned int cols, unsigned int nnz,
                  unsigned int *rowVals, unsigned int *colVals,
                  T *values);
        // Construct from a dense matrix
        CSRMatrix(const Matrix<T> &m);
        // Read from a file
        CSRMatrix(string filename);
        // Destructor
        ~CSRMatrix();
        
        /*
        -------Accessor-------------------
        DON'T USE THIS IT'S HELLA SLOW
        */
        const T at(unsigned int row, unsigned col);

        /*
        -------Printing / File I/O---------------------
        */
        // This guy hella slow mate
        virtual void printMatrix();
        virtual void printValues();
        bool writeToFile(string filename);

        /*
        -------Linear Algebra Functionality-----
        */
        // Solve a triangular system on a sparse matrix but with dense 
        // RHS vector b
        // set lower=false for upper triangular matrix.
        // Assumes the matrix is non-singular
        void solveDenseTriangular(T *b, T *res, bool lower);
        // Solve a sparse upper triangular system with sparse RHS
        // This should be very efficient
        void solveSparseTriangular(T *b, T *res, bool lower);
        // For a strictly diagonal dominant Sparse Matrix
        void jacobiSolve(T *b, T tolerance, T *res);
        // For an SPD Sparse Matrix
        void conjGradSolve(T *b, T tol, T *res);
        void matVecMult(T *b, T *res);
        // Needed for inheritance, but due to the above 2, let's do
        // it slightly differently. We check b's sparsity, and then call
        // the appropriate method
        void triangularSolve(T *b, T *res, bool lower);
        // Also needed for inheritance purposes
        void choleskySolve(T *b, T *res);
    
    private:
        // Numbers below this (absolute value) tolerance are treated as 0
        // Can change this to something more useful but for practical purposes
        // 1e-9 is not bad
        double kZeroTolerance = 1e-9;
        T dotProduct(T *a, T *b);
};
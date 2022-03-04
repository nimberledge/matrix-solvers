#include "CSRMatrix.h"
#include <fstream>
#include <cmath>
#include <iostream>
#include <memory>

using std::cerr;
using std::cout;
using std::endl;

template<class T>
CSRMatrix<T>::CSRMatrix(unsigned int rows, unsigned int cols, unsigned int nnz)
:Matrix<T>(rows, cols, true), nnz(nnz) {
    rowIdxs = std::unique_ptr<unsigned int[]> (new unsigned int[rows+1]);
    colIdxs = std::unique_ptr<unsigned int[]> (new unsigned int[nnz]);
    this->values = std::shared_ptr<T[]> (new T[nnz]);
}

template<class T>
CSRMatrix<T>::CSRMatrix(unsigned int rows, unsigned int cols, unsigned int nnz,
                  unsigned int *rowVals, unsigned int *colVals,
                  T *values) : Matrix<T>(rows, cols, true), nnz(nnz) { 
    // The following lines required a bit of thought, because we have to now create
    // unique pointers from the values passed
    // got memory errors when this was not deep-copied
    this->values = std::shared_ptr<T[]> (new T[nnz]);
    rowIdxs = std::unique_ptr<unsigned int[]> (new unsigned int[this->rows+1]);
    colIdxs = std::unique_ptr<unsigned int[]> (new unsigned int[nnz]);
    for (unsigned int i = 0; i < nnz; i++) {
        this->values[i] = values[i];
        this->colIdxs[i] = colVals[i];
        if (i <= rows) {
            rowIdxs[i] = rowVals[i];
        }
    }
}

template<class T>
CSRMatrix<T>::CSRMatrix(const Matrix<T> &m): Matrix<T>(m.rows, m.cols, true) {
    nnz = 0;
    for (unsigned int i = 0; i < m.rows * m.cols; i++) {
        if (std::abs(m.values[i]) > kZeroTolerance) {
            nnz++;
        }
    }
    rowIdxs = std::unique_ptr<unsigned int[]> (new unsigned int[this->rows+1]);
    colIdxs = std::unique_ptr<unsigned int[]> (new unsigned int[nnz]);
    this->values = std::shared_ptr<T[]> (new T[nnz]);
    nnz = 0;
    for (unsigned int i = 0; i < this->rows; i++) {
        // index of where row i starts
        rowIdxs[i] = nnz;
        for (unsigned int j = 0; j < this->cols; j++) {
            // Check if it is zero, if not you gots to do it
            if (std::abs(m.values[i*this->cols + j]) > kZeroTolerance) {
                colIdxs[nnz] = j;
                this->values[nnz] = m.values[i*this->cols + j];
                nnz++;
            }
        }
    }
    rowIdxs[this->rows] = nnz;
}

template<class T>
CSRMatrix<T>::CSRMatrix(string filename): Matrix<T>(0, 0, true) {
    string line;
    std::ifstream inFile(filename.c_str());
    if (!inFile.good()) {
        cerr << "Bad input file, Initializing empty CSR Matrix" << endl;
        cerr << "filename: " << filename << endl;
        this->rows = 0;
        this->cols = 0;
        nnz = 0;
        rowIdxs = nullptr;
        colIdxs = nullptr;
        this->values = nullptr;
        return;
    }
    // First line is rows, cols, nnz
    std::getline(inFile, line);
    std::istringstream strm1(line);
    strm1 >> this->rows;
    strm1 >> this->cols;
    strm1 >> nnz;
    rowIdxs = std::unique_ptr<unsigned int[]> (new unsigned int[this->rows+1]);
    colIdxs = std::unique_ptr<unsigned int[]> (new unsigned int[nnz]);
    this->values = std::shared_ptr<T[]> (new T[nnz]);
    unsigned int idx = 0;
    unsigned int row_idx = 0;
    rowIdxs[0] = 0;
    while (std::getline(inFile, line)) {
        // We need a new stream
        std::istringstream strm(line);
        unsigned int row;
        string temp;
        // row, then col, then value
        strm >> row;
        strm >> colIdxs[idx];
        strm >> temp;
        // Check for empty rows
        while (row_idx < row) {
            rowIdxs[++row_idx] = idx;
        }
        
        // Check if it's a floating point type because atof
        if (std::is_same<T, float>::value || std::is_same<T, double>::value) {
            this->values[idx] = std::atof(temp.c_str());
        } else {
            this->values[idx] = std::atoi(temp.c_str());
        }
        ++idx;
    }
    rowIdxs[this->rows] = nnz;
}

template<class T>
CSRMatrix<T>::~CSRMatrix() {  }

template<class T>
const T CSRMatrix<T>::at(unsigned int row, unsigned int col) {
    if (row >= this->rows || col >= this->cols) {
        cerr << "Invalid row / col returning 0.0" << endl;
        return 0.0;
    }
    // If there are no non-zeros, easy optimization
    if (nnz == 0)  return 0.0;
    // Deal with last row edge case
    unsigned int searchIdx = rowIdxs[row];
    while (searchIdx < rowIdxs[row+1]) {
        if (colIdxs[searchIdx] == col) {
            return this->values[searchIdx];
        }
        searchIdx++;
    }
    return 0.0;
}

template<class T>
void CSRMatrix<T>::printMatrix() {
    for (unsigned int i = 0; i < this->rows; i++) {
        for (unsigned int j = 0; j < this->cols; j++) {
            cout << this->at(i, j) << "\t\t";
        }
        cout << endl;
    }
}

template<class T>
void CSRMatrix<T>::printValues() {
    cout << "RowIdxs: ";
    for (unsigned int i = 0; i <= this->rows; i++) {
        cout << rowIdxs[i] << ", ";
    }
    cout << "\ncolIdxs: ";
    for (unsigned int i = 0; i < nnz; i++) {
        cout << colIdxs[i] << ", ";
    }
    cout << "\nvalues: ";
    for (unsigned int i = 0; i < nnz; i++) {
        cout << this->values[i] << ", ";
    }
}

template<class T>
void CSRMatrix<T>::solveDenseTriangular(T *b, T *res, bool lower) {
    if (lower) {
        for (unsigned int row = 0; row < this->rows; row++) {
            unsigned int colIdx = rowIdxs[row];
            T tmp = b[row];
            while (colIdx < rowIdxs[row+1] && colIdxs[colIdx] < row) {
                tmp -= this->values[colIdx] * res[colIdxs[colIdx]];
                colIdx++;
            }
            // If we have a non-singular triangular matrix, then 
            // the last colVec[colIdx] value in the row will be 
            // colVec[colIdx] = row, the diagonal entry
            res[row] = tmp / this->values[colIdx];
        }
        return;
    }
    // Upper triangular matrix
    for (int row = this->rows-1; row >= 0; row--) {
        // If non-singular upper triangular, we are assured that
        // rowVec[row] corresponds to the diagonal entry
        // So instead, we now update colIdx by 1
        unsigned int colIdx = rowIdxs[row] + 1;
        T tmp = b[row];
        while (colIdx < rowIdxs[row+1]) {
            tmp -= this->values[colIdx] * res[colIdxs[colIdx]];
            colIdx++;
        }
        // In theory, dataValues[rowVec[row]] should be M[row, row]
        // But requires a non-singular triangular matrix
        res[row] = tmp / this->values[rowIdxs[row]];
    }
}

template<class T>
void CSRMatrix<T>::solveSparseTriangular(T *b, T *res, bool lower) {
    if (lower) {
        // There is not an efficient way to do lower triangular solves 
        // using a CSR Matrix. You would need a CSC Matrix for that
        // The reason being, the efficiency gained here is through the 
        // checking for non-zeros in x step. If we were to do this for 
        // a lower triangular matrix, we would have an expensive check
        // at which cost, we might as well just do a dense solve.
        solveDenseTriangular(b, res, lower);
        return;
    }
    // Only do the clever stuff for an Upper Triangular Matrix
    // visitStack = non-zero entries in x, in a logical order
    int visitStack[this->rows];
    int stackUpdateIndex = this->rows;
    bool marked[this->rows];
    // Set marked to false, result to 0
    std::fill(res, res+this->rows, 0);
    std::fill(marked, marked+this->rows, 0);
    // Compute indices i of x that are non-zero
    // Rule 1: if b[i] is non-zero, then x[i] is also
    // Rule 2: if there is j such that x[j] is non-zero and 
    // A[i, j] is non-zero, then x[i] is non-zero
    for (int i = this->rows-1; i >= 0; i--) {
        // if b[i] is non-zero, x[i] is non-zero so we gotta visit this
        if (!marked[i] && std::abs(b[i]) > kZeroTolerance) {
            visitStack[--stackUpdateIndex] = i;
            marked[i] = true;
        } else {
            // Otherwise, if there's a column j in this row
            // such that A_ij is non-zero and x_j is non-zero
            // then x_i is non-zero
            unsigned int colIdx = rowIdxs[i];
            while (colIdx < rowIdxs[i+1]) {
                if (marked[colIdxs[colIdx]]) {
                    marked[i] = true;
                    visitStack[--stackUpdateIndex] = i;
                }
                colIdx++;
            }
        }
    }
    // Now that we have the non-zeros in x, for each of those indices
    // Calculate the back-substitution
    // The i-loop loops through our visitStack, which holds non-zero indices
    // of x. The colIdx loop loops over non-zeros in our sparse matrix
    // Overall, this gives us a cheap back substitution
    for (int i = this->rows-1; i >= stackUpdateIndex; i--) {
        unsigned int row = visitStack[i];
        unsigned int colIdx = rowIdxs[row] + 1;
        T tmp = b[row];
        while (colIdx < rowIdxs[row+1]) {
            tmp = tmp - this->values[colIdx] * res[colIdxs[colIdx]];
            colIdx++;
        }
        res[row] = tmp / this->values[rowIdxs[row]];
    }
}

template<class T>
void CSRMatrix<T>::jacobiSolve(T *b, T tolerance, T *res) {
    int maxIter = 1000;
    int iter = 0;
    T x_new[this->rows];
    double difference = 0.0;
    while (iter < maxIter) {
        difference = 0.0;
        for (unsigned int i = 0; i < this->rows; i++) {
            unsigned int diagIdx = rowIdxs[i];
            double sigma = 0.0;
            for (unsigned int j = rowIdxs[i]; j < rowIdxs[i+1]; j++) {
                if (colIdxs[j] != i) { 
                    // Non-diagonal element, do computation
                    sigma += this->values[j] * res[colIdxs[j]];
                } else {
                    // Diagonal element, store the index
                    diagIdx = j;
                }
            }
            x_new[i] = (b[i] - sigma) / (this->values[diagIdx]);
            difference += pow(x_new[i] - res[i], 2);
        }
        // Check for stopping condition using L^2 norm
        for (unsigned int i = 0; i < this->rows; i++) {
            res[i] = x_new[i];
        }
        if (sqrt(difference) < tolerance) {
            return;
        }
        iter++;
    }
    cerr << "jacobiSolve() did not terminate early, solution could be inaccurate" << endl;
}

template<class T>
bool CSRMatrix<T>::writeToFile(string filename) {
    std::ofstream outFile;
    outFile.open(filename, std::ios_base::out);
    if (!outFile.good()) {
        return false;
    }
    // Defining the file format
    outFile << this->rows << "\t" << this->cols << "\t" << nnz << "\n";
    for (unsigned int i = 0; i < this->rows; i++) {
        unsigned int colIdx = rowIdxs[i];
        while (colIdx < rowIdxs[i+1]) {
            outFile << i << "\t" << colIdxs[colIdx] << "\t" << this->values[colIdx] << "\n";
            colIdx++;
        }
    }
    outFile.close();
    return true;
}

template <class T>
T CSRMatrix<T>::dotProduct(T *a, T *b) {
    T res = 0;
    for (unsigned int i = 0; i < this->rows; i++) {
        res += a[i] * b[i];
    }
    return res;
}

template <class T>
void CSRMatrix<T>::matVecMult(T *b, T *res) {
    for (unsigned int i = 0; i < this->rows; i++) {
        res[i] = 0.0;
        for (unsigned int colIdx = rowIdxs[i]; colIdx < rowIdxs[i+1]; colIdx++) {
            unsigned int j = colIdxs[colIdx];
            res[i] += this->values[colIdx] * b[j];
        }
    }
}

template <class T>
void CSRMatrix<T>::conjGradSolve(T *b, T tol, T *res) {
    // Sparse Krylov Solver
    // Reference: http://adl.stanford.edu/aa222/Lecture_Notes_files/CG_lecture.pdf
    T r[this->rows], p[this->rows], Ap[this->rows], Ax[this->rows];
    this->matVecMult(res, Ax);
    for (unsigned int i = 0; i < this->rows; i++) {
        p[i] = r[i] = b[i] - Ax[i];
    }
    T oldResidual = this->dotProduct(r, r);
    for (unsigned int i = 0; i < this->rows; i++) {
        this->matVecMult(p, Ap);
        T alpha = oldResidual / this->dotProduct(p, Ap);
        // Could easily be vectorized, hope the compiler gods know this
        for (unsigned int j = 0; j < this->rows; j++) {
            res[j] += alpha * p[j];
            r[j] -= alpha * Ap[j];
        }
        // Check stopping condition
        T newResidual = this->dotProduct(r, r);
        if (sqrt(newResidual) < tol) {
            return;
        }
        // Update p
        T resRatio = newResidual / oldResidual;
        for (unsigned int j = 0; j < this->rows; j++) {
            p[j] = r[j] + resRatio * p[j];
        }
        oldResidual = newResidual;
    }
    cerr << "maxIter reached? " << endl;
}

template <class T>
void CSRMatrix<T>::triangularSolve(T *b, T *res, bool lower) {
    if (lower) {
        this->solveDenseTriangular(b, res, lower);
        return;
    }
    unsigned int nnzsB = 0;
    for (unsigned int i = 0; i < this->rows; i++) {
        if (std::abs(b[i]) < kZeroTolerance) {
            nnzsB++;
        }
    }
    double RHSSparsityLimit = 0.5;
    if (((double) nnzsB / this->rows) < RHSSparsityLimit) {
        this->solveSparseTriangular(b, res, lower);
        return;
    }
    this->solveDenseTriangular(b, res, lower);
}

template <class T>
void CSRMatrix<T>::choleskySolve(T *b, T *res) {
    cerr << "CSRMatrix choleskySolve() not supported, sorry" << endl;
}
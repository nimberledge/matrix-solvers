#include "Matrix.h"
#include <cmath>

using std::cerr;
using std::endl;
using std::cout;

// Use initializer list for performance (sorta)
template<class T>
Matrix<T>::Matrix(unsigned int rows, unsigned int cols, bool preallocated)
: rows(rows), cols(cols), preallocated(preallocated) {
    if (!preallocated) {
        values = (std::shared_ptr<T[]>) new T[rows * cols];
    }
}

template<class T>
Matrix<T>::Matrix(unsigned int rows, unsigned int cols, T *values)
: rows(rows), cols(cols), values(std::shared_ptr<T[]> (values)) {
    preallocated = true;
}

template<class T>
Matrix<T>::Matrix(std::string filename) {
    preallocated = false;
    std::ifstream inFile(filename.c_str());
    if (!inFile.good()) {
        cerr << "Bad input file, empty matrix incoming" << endl;
        rows = 0;
        cols = 0;
        values = nullptr;
        return;
    }
    // Now read the file
    inFile >> rows >> cols;
    values = std::shared_ptr<T[]> (new T[rows * cols]);
    for(int i = 0; i < rows*cols; i++) {
        std::string temp;
        inFile >> temp;
        // Check if it's a floating point type because atof
        if (std::is_same<T, float>::value || std::is_same<T, double>::value) {
            values[i] = std::atof(temp.c_str());
        } else {
            values[i] = std::atoi(temp.c_str());
        }
    }
}

template<class T>
void Matrix<T>::printMatrix() {
    for (unsigned i = 0; i < rows; i++) {
        for (unsigned j = 0; j < cols; j++) {
            // print a tab so everything is nicely aligned
            cout << values[i*cols + j] << "\t\t";
        }
        cout << endl;
    }
}

template<class T>
void Matrix<T>::printValues() {
    // print the values and a space after
    for (unsigned i = 0; i < rows * cols; i++) {
        cout << values[i] << " ";
    }
    cout << endl;
}

template<class T>
Matrix<T>::~Matrix() { }

template<class T>
void Matrix<T>::jacobiSolve(T* b, T tolerance, T* res) {
    double difference = 0;
    double sig;
    T new_solution[cols];
    int maxIter = 1000;
    for(int iterations = 0; iterations < maxIter; iterations++) {   
        difference = 0.0;
        for(int i = 0; i < rows; i++) {
            sig = 0;
            for(int j = 0; j < cols; j++) {
                if( j != i) {
                    sig += values[i*cols + j] * res[j];
                }
            }
            new_solution[i] = (b[i] - sig)/(values[i*cols + i]);
            difference += pow(res[i] - new_solution[i], 2);
        }
        // Copy new_sol -> solution
        for (int i = 0; i < rows; i++) {
            res[i] = new_solution[i];
        }
        // Check stopping condition - L2 norm
        if(sqrt(difference) < tolerance) {
            return;
        }
    }
    cerr << "maxIter reached, solution may not be good" << endl;
}

template<class T>
void Matrix<T>::matVecMult(T *b, T *res) {
    for (unsigned int i = 0; i < rows; i++) {
        res[i] = 0.0;
        for (unsigned int j = 0; j < cols; j++) {
            res[i] += this->values[i*cols + j] * b[j];
        }
    }
}

template<class T>
void Matrix<T>::triangularSolve(T *b, T *res, bool lower) {
    if (lower) {
        for (int i = 0; i < rows; i++) {
            T tmp = b[i];
            for (int j = 0; j < i; j++) {
                tmp -= this->values[i*cols + j] * res[j];
            }
            res[i] = tmp / this->values[i*cols + i];
        }
        return;
    }
    // Upper triangular
    for (int i = rows-1; i >= 0; i--) {
        T tmp = b[i];
        for (int j = i+1; j < cols; j++) {
            tmp -= this->values[i*cols + j] * res[j];
        }
        res[i] = tmp / this->values[i*cols + i];
    }
}

template<class T>
void Matrix<T>::choleskySolve(T *b, T *res) {
    Matrix L(rows, cols, false);
    this->_copyValues(&L);
    for (unsigned int k = 0; k < cols; k++) {
        L.values[k*cols + k] = sqrt(L.values[k*cols + k]);
        T Lkk = L.values[k*cols + k];
        // Scale current column
        for (unsigned int i = k+1; i < cols; i++) {
            L.values[i*cols + k] = L.values[i*cols + k] / Lkk;
        }
        // Column operation subtract a multiple of column k
        for (unsigned int j = k+1; j < cols; j++) {
            for (unsigned int i = j; i < cols; i++) {
                L.values[i*cols + j] = L.values[i*cols + j] - L.values[i*cols + k] * L.values[j*cols + k];
            }
        }
    }
    
    // Compute L.T
    Matrix LT(rows, cols, false);
    for (unsigned int i = 0; i < rows; i++) {
        for (unsigned int j = 0; j < cols; j++) {
            LT.values[j*cols + i] = L.values[i*cols + j];
        }
    }
    // Solve Ly = b
    // then L.T x = y
    T y[rows];
    L.triangularSolve(b, y, true);
    LT.triangularSolve(y, res, false);
}

template<class T>
void Matrix<T>::_copyValues(Matrix<T> *other) {
    for (unsigned int i = 0; i < rows * cols; i++) {
        other->values[i] = values[i];
    }
}

template <class T>
bool Matrix<T>::writeToFile(string filename) {
    std::ofstream file (filename);
    if (!file.good()) {
	return false;
    }

    file << rows << "\t" << cols << endl;
    
    for (int i = 0;i < rows * cols; i++)
    {
        file << values[i] << "\t";
    }    
    
    file << endl;

    return true;
}

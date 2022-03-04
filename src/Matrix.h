#pragma once
#include <iostream>
#include <memory>
#include <sstream>
#include <fstream>
#include <string>

using std::string;

template<class T>
class Matrix {
    public:
        /*
        ---------Class Members-----------
        */
        // Take unsigned int for rows and columns
        // so that they can't be negative
        unsigned int rows;
        unsigned int cols;
        /*
        Matrix stores data in ROW-MAJOR FORM.
        ROW MAJOR
        */
        std::shared_ptr<T[]> values;


        /*
        ---------Constructors / Destructors-----------
        */

        // Constructor without pre-allocated memory
        Matrix(unsigned int rows, unsigned int cols, bool preallocated);
        // Constructor with pre-allocated memory
        Matrix(unsigned int rows, unsigned int cols, T* values);
        // Destructor
        ~Matrix();
        // Constructor to read file
        Matrix(string filename);

        /*
        ---------Printing--------------
        */
        virtual void printMatrix();
        virtual void printValues();
        bool writeToFile(string filename);

        /*
        --------Linear Algebra Functionality------
        */
        // Jacobi solver that solves Ax=b iteratively to a given tolerance
        virtual void jacobiSolve(T* b, T tolerance, T* res);
        // Assumes that b and res are the correct shapes, segfaults are
        // your problem
        void matVecMult(T *b, T *res);
        // Use lower=false for upper triangular solve
        void triangularSolve(T *b, T *res, bool lower);
        // Assumes Matrix is SPD (and therefore square of course)
        void choleskySolve(T *b, T *res);
    
    private:
        // Determines whether or not to dynamically allocate memory
        // This is opposite to the course code given to us. The way I read
        // it is preallocated means that the memory was allocated 
        // pre-instantiation, and thus not allocated by us. 
        bool preallocated;

        // Do a deep copy of this into other
        void _copyValues(Matrix *other);
};
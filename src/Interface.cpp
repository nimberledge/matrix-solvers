#include "Interface.h"
#include <iostream>
#include <fstream>

using std::cout;
using std::endl;
using std::cerr;
using std::cin;



SolverInterface::SolverInterface() { }

SolverInterface::~SolverInterface() { }

int SolverInterface::start() {
    string message = "Welcome to the best linear solver in the world; THE LADS' SOLVER.";
    this->printMessage(message);
    int matrixType = this->getDenseOrSparse();
    // Dense solvers
    if (matrixType == 0) {
        Matrix<double> *m = this->getMatrixDense();
        cout << "Printing the matrix: " << endl;
        m->printMatrix();
        int solverChoice = this->getMethodDense();
        double b[m->rows], res[m->rows];
        this->getVectorb(m->rows, b);
        // Jacobi Solver
        double tol = 1e-6;
        if (solverChoice == 0) {
            m->jacobiSolve(b, tol, res);
        } else if (solverChoice == 1) {
            m->choleskySolve(b, res);
        }
        cout << "Solution vector: " << endl;
        cout << "[";
        for (unsigned int i = 0; i < m->rows; i++) {
            if (i != m->rows - 1) {
                cout << res[i] << ", ";
            }
            else {
                cout << res[i];
            }
        }
        cout << "]" << endl;
        delete m;
    }
    // Sparse solvers
    else if (matrixType == 1) {
        CSRMatrix<double> *m = this->getMatrixSparse();
        cout << "Printing the matrix: " << endl;
        m->printMatrix();
        int solverChoice = this->getMethodSparse();
        double b[m->rows], res[m->rows];
        std::fill(res, res+m->rows, 0);
        this->getVectorb(m->rows, b);
        // Jacobi Solver
        double tol = 1e-6;
        if (solverChoice == 0) {
            m->jacobiSolve(b, tol, res);
        } else if (solverChoice == 1) {
            m->conjGradSolve(b, tol, res);
        }

        cout << "Solution vector: " << endl;
        cout << "[";
        for (unsigned int i = 0; i < m->rows; i++) {
            if (i != m->rows - 1) {
                cout << res[i] << ", ";
            }
            else {
                cout << res[i];
            }
        }
        cout << "]" << endl;
        delete m;
    }

    return 0;
}

int SolverInterface::getDenseOrSparse() {
    int option;
    cout << "Will your input matrix be in Dense or Sparse form? Press 0 for Dense and 1 for Sparse: ";
    cin >> option;

    if (option == 0) {
        cout << "You have chosen to input a Dense matrix!" << endl;
        return 0;
    }

    if (option == 1) {
        cout << "You have chosen to input a Sparse matrix!" << endl;
        return 1;
    }

    cout << "That's not an option I gave you...." << endl;
    return -1;
}


Matrix<double> *SolverInterface::getMatrixDense() {
    // Ask user for file or input on CLI
    int choice;
    cout << "Press 0 to input matrix values manually or press 1 to read in a text file: ";
    cin >> choice;
    do {
        if (choice == 0) {
            int rows, cols;
            this->getMatrixSize(&rows, &cols);
            // create matrix m
            Matrix<double> *m = new Matrix<double>(rows, cols, false);
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < cols; j++) {
                    this->getElement(i, j, m);
                }
            }
            return m;
        }

        if (choice == 1) {
            string filename = this->getFileName();
            Matrix<double> *m = new Matrix<double>(filename);
            cout << "Exiting getMatrixDense()" << endl;
            return m;
        }
    } while (choice != 0 && choice != 1);
    // Should never happen
    return nullptr;
}


CSRMatrix<double> *SolverInterface::getMatrixSparse() {
    // Ask user for filename
    string filename = this->getFileName();
    CSRMatrix<double> *m = new CSRMatrix<double>(filename);
    return m;
}

void SolverInterface::getMatrixSize(int *rowptr, int *colptr) {
    int rows;
    int cols;
    cout << "\nPlease enter the dimensions of your matrix below"<< endl;
    cout << "\nNumber of rows: ";
    cin >> rows;
    while (rows < 1) {
        cout << "\n Sorry! the number of rows cannot be less than 1! Please enter the number of rows; ";
        cin >> rows;
    }
    
    cout << "\nNumber of columns: ";
    cin >> cols;
    while(cols < 1) {
        cout << "\n Sorry! the number of columns cannot be less than 1! Please enter the number of columns; ";
        cin >> cols;
    } 

    *rowptr = rows;
    *colptr = cols;
}

void SolverInterface::getElement(int i, int j, Matrix<double> *m) {
    double element;
    cout << "\nProvide Element (" << i << ", " << j << "): ";
    cin >> element;
    m->values[i * m->cols + j] = element;
}

string SolverInterface::getFileName() {
    string filename;
    cout << "\nProvide a filename (eg. somefile.txt): ";
    cin >> filename;
    return filename;
}

int SolverInterface::getMethodDense() {
    int method;
    cout << "Press 0 to use the Jacobi Solver, Press 1 to use the Cholesky Solver: ";
    cin >> method;
    if (method == 0) {
        cout << "You have chosen the Greatest Jacobi Solver in the world!" << endl;
        return 0;
    }
    if(method == 1) {
        cout << "Good Choice! You have chosen the Greatest Cholesky Solver in the world!" << endl;
        return 1;
    }
    cout << "That's not an option I gave you...." << endl;
    return -1;
}

int SolverInterface::getMethodSparse() {
    int method;
    cout << "Press 0 to use the Jacobi Sparse Solver, Press 1 to use the Sparse Krylov (Conjugate Gradient) Solver: ";
    cin >> method;
    if (method == 0) {
        cout << "You have chosen the Greatest Jacobi Sparse Solver in the world!" << endl;
        return 0;
    }
    if(method == 1) {
        cout << "Good Choice! You have chosen the Greatest Krylov Solver in the world!" << endl;
        return 1;
    }
    cout << "That's not an option I gave you...." << endl;
    return -1;
}

void SolverInterface::getVectorb(int rows, double *b) {
    for (unsigned int i = 0; i < rows; i++) {
        cout << "Enter b[" << i << "]:" << endl;
        string temp;
        cin >> temp;
        b[i] = std::atof(temp.c_str());
    }
    cout << "b as entered: " << endl;
    cout << "[";
    for (unsigned int i = 0; i < rows; i++) {
        if (i != rows - 1) {
            cout << b[i] << ", ";
        }
        else {
            cout << b[i];
        }
    }
    cout << "]" << endl;
}

void SolverInterface::printMessage(string message) {
    cout << message << endl;
}

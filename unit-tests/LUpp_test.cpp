#include <iostream>
#include <math.h>
#include <ctime>
#include "../Matrix.h"
#include "../Matrix.cpp"

using namespace std;

// test LU decomposition with partial pivoting by ensuring norm is sufficiently small

int main()
{
    // generate random dense diagonally dominant real positive definite matrix of given size

    int rows = rand() % 100 + 10;
    int cols = rows;

    // initialise matrix A
    unique_ptr<Matrix<double> > A(new Matrix<double>(rows, cols, true));
    
    A->GenerateRandomSPDMatrix(rows, cols, *A);

    // create vector and fill with randomly generated numbers
    unique_ptr<double []> b(new double[cols]);
    for (int i=0; i<cols; i++)
    {
        b[i] = rand() % 100 + 1;
    }

    // initialise solution x
    unique_ptr<double []> x(new double[cols]);

    // solve linear system to find x
    A->LU_pp_solver(b, x);

    // calculate norm
    double norm = A->getNorm(x, b);

    // give warning if norm is greater than 1e-10
    if (norm > 1e-10)
    {
        cerr << "Norm is greater than 1e-10:" << norm << endl;
    }

}

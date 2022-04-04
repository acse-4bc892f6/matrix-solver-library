#include <iostream>
#include <math.h>
#include <ctime>
#include <algorithm>
#include "../Matrix.h"
#include "../Matrix.cpp"
#include "../CSRMatrix.h"
#include "../CSRMatrix.cpp"

using namespace std;

// TEST JACOBI SOLVER WORKS FOR DENSE MATRICES AND CSR MATRICES
// AND THE L2 NORM OF RESIDUAL IS SUFFICIENTLY SMALL FOR BOTH

int main()
{
    srand(time(0));

    int rows = rand() % 10 + 10;
    int cols = rows;

    // number of non zero rows below diagonal
    int bandwidth = rand() % 10 + 0;

    // initialise dense matrix
    unique_ptr<Matrix<double> > dense_A(new Matrix<double>(rows, cols, true));
    // fill it with randomly generated numbers based on criteria above
    dense_A->GenerateRandomSPDMatrix(rows, cols, *dense_A);

    // create vector and fill with randomly generated numbers
    unique_ptr<double []> b(new double[cols]);
    for (int i=0; i<cols; i++)
    {
        b[i] = rand() % 100 + 1;
    }

    // initialise solution x to A_dense and fill it with zeroes
    unique_ptr<double []> x(new double[cols]);

    for (int i=0; i<cols; i++)
    {
        x[i] = 0.0;
    }

    // use Jacobi to find x
    dense_A->Jacobi(b,x,10);
    // calculate L2 norm of residual
    double norm = dense_A->getNorm(x, b);

    unique_ptr<CSRMatrix<double> > A(new CSRMatrix<double>(rows, cols, 0, false));
    A->GenerateRandomSPDMatrix(rows, cols, bandwidth, *A);

    // initialise solution x2 and fill it with zeroes
    unique_ptr<double []> x2(new double[cols]);
    for (int i=0; i<cols; i++)
    {
        x2[i] = 0.0;
    }

    // use Jacobi to find x2
    A->Jacobi(b,x2,500);
    // calculate L2 norm of residual
    double norm2 = A->getNorm(x2, b);

    // give out warnings if L2 norm of either residual is greater than 0.001
    
    if (norm > 1e-3)
    {
        cerr << "Norm of residual from dense matrix is greater than 1e-3: " << norm << endl;
        cerr << "rows, cols" << endl;
        cerr << rows << ", " << cols << endl;
    }

    if (norm2 > 1e-3)
    {
        cerr << "Norm of residual from sparse matrix is greater than 1e-3: " << norm2 << endl;
        cerr << "rows, cols, bandwidth" << endl;
        cerr << rows << ", " << cols << ", " << bandwidth  << endl;
    }

}
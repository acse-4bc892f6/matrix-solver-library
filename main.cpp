#include <iostream>
#include <math.h>
#include <ctime>
#include <algorithm>
#include "Matrix.h"
#include "Matrix.cpp"
#include "CSRMatrix.h"
#include "CSRMatrix.cpp"
#include <fstream>

using namespace std;

// MAIN FILE THAT RUNS DIFFERENT FUNCTIONS AND 
// SAVES OUTPUT SOLVER AND NORM TO FILE


int main(){

    srand(time(0));

    // generate a random number
    int m = rand() % 10 + 10;
    int rows = m;
    int cols = m;

    // number of non zero rows below diagonal
    int bandwidth = rand() % 10 + 0;

    // generate a sparse matrix
    unique_ptr<CSRMatrix<double> > sparse_A(new CSRMatrix<double>(rows, cols, 0, false));
    sparse_A->GenerateRandomSPDMatrix(rows, cols, bandwidth,  *sparse_A);



    unique_ptr<Matrix<double> > A(new Matrix<double>(rows, cols, true));

    // convert sparse matrix to dense
    // to perform the same functions on the dense solvers
    sparse_A->CSRtoDense(*A);


    // create vector and fill with randomly generated numbers
    unique_ptr<double []> b(new double[cols]);
    for (int i=0; i<cols; i++)
    {
    b[i] = rand() % 10 + 1;
    }

    // initialise solution x
    // dense solutions
    unique_ptr<double []> LUpp_x(new double[cols]);
    unique_ptr<double []> jacobi_x(new double[cols]);
    unique_ptr<double []> gauss_seidel_x(new double[cols]);

    // sparse solution
    unique_ptr<double []> gauss_seidel_x_CSR(new double[cols]);
    unique_ptr<double []> jacobi_x_CSR(new double[cols]);

    // solve linear system to find x and time each solver

    // LU Solver
    A->LU_pp_solver(b, LUpp_x);
    // Jacobi Solver
    A->Jacobi(b, jacobi_x, 100);
       // Gauss-Seidel solver
    A->GaussSeidel(b, gauss_seidel_x, 100);

    // Jacobi sparse solver
    sparse_A->Jacobi(b, jacobi_x_CSR, 100);
    // Gauss seidel solver
    sparse_A->GaussSeidel(b, gauss_seidel_x_CSR, 100);


    // check that the b and A * x match
    // compute the norms
    double LUpp_norm = A->getNorm(LUpp_x, b);

    double jacobi_norm = A->getNorm(jacobi_x, b);

    double gauss_seidel_norm = A->getNorm(gauss_seidel_x, b);

    double jacobi_norm_CSR = sparse_A->getNorm(jacobi_x_CSR, b);

    double gauss_seidel_norm_CSR = sparse_A->getNorm(gauss_seidel_x_CSR, b);


    // print norms to file 
    // print solutions to file
    ofstream myfile;
    myfile.open("Matrix-solutions-test.txt");
    // Matrix Size
    myfile << "Matrix Size:" << rows << "x" << rows <<endl;

    // Dense solutions
    myfile << "Dense solution LUpp norms:" << endl;
    myfile << LUpp_norm << " ";
    myfile << endl;

    myfile << "LUpp Dense solution:" << endl;
    for (int i=0; i<cols; i++)
        {
            myfile << LUpp_x[i] << " ";
        }
    myfile << endl;


    myfile << "Dense solution Jacobi norms:" << endl;
    myfile << jacobi_norm << " ";
    myfile << endl;

    
    myfile << "Jacobi Dense solution:" << endl;
    for (int i=0; i<cols; i++)
        {
            myfile << jacobi_x[i] << " ";
        }
    myfile << endl;


    myfile << "Dense solution Gauss-Seidel norms:" << endl;
    myfile << gauss_seidel_norm << " ";
    myfile << endl;

    myfile << "Gauss-Seidel Dense solution:" << endl;
    for (int i=0; i<cols; i++)
    {
        myfile << gauss_seidel_x[i] << " ";
    }
    myfile << endl;
    
    // Sparse solutions
    myfile << "Jacobi Sparse solution:" << endl;
    myfile << jacobi_norm_CSR << " ";
    myfile << endl;

    myfile << "Jacobi norms Sparse solution:" << endl;
    for (int i=0; i<cols; i++)
    {
        myfile << jacobi_x_CSR[i] << " ";
    }
    myfile << endl;


    myfile << "Gauss-Seidel Sparse solution norm:" << endl;
    myfile << gauss_seidel_norm_CSR << " ";
    myfile << endl;

    myfile << "Gauss-Seidel Sparse solution:" << endl;
    for (int i=0; i<cols; i++)
    {
        myfile << gauss_seidel_x_CSR[i] << " ";
    }
    myfile << endl;
    myfile.close();


}

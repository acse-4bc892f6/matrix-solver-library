#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "Matrix.h"
#include "Matrix.cpp"
#include "CSRMatrix.h"
#include "CSRMatrix.cpp"
#include <stdlib.h> 

using namespace std;

int main()
{
    int rows, cols;
    
    // allow use to choose between dense matrix and csr sparse matrix
    int matrix_type = -1;

    // in case use types in anything other than 1 and 2
    while (matrix_type != 1 && matrix_type != 2)
    {
        cout << "Linear Solvers" << endl;
        cout << " Choose which matrix you would like to solve for. Enter 1 or 2 for:" << endl;
        cout << "1. Dense matrix" << endl;
        cout << "2. Sparse matrix" << endl;
    
        cin >> matrix_type;
    }

    // allow user to choose between importing matrix from text file or randomly generate matrix
    int import_matrix = -1;

    if (matrix_type == 1)
    {
        // define text file format
        cout << "If you choose to import matrix values from a text file, this program assumes" << endl;
        cout << "the text file is in the formats shown below" << endl;
        cout << endl;
        cout << "1) For dense matrix format input a series of single space separated values" << endl;
        cout << "in row-major ordering, for example for the 3 x 3 matrix:" << endl;
        cout << "1 2 3" << endl; 
        cout << "4 5 6" << endl; 
        cout << "7 8 9" << endl; 
        cout << "The input text file should read:" << endl;
        cout << "1 2 3 4 5 6 7 8 9" << endl;
        cout << "If both the vector and the matrix is imported, they must be provided in separate files" << endl;
        cout << endl;
        // in case use types in anything other than 1 or 2
        while (import_matrix != 1 && import_matrix != 2)
        {
            cout << "Do you want to import matrix values from text file? Enter 1 or 2 for:" << endl;
            cout << "1. Yes" << endl;
            cout << "2. No" << endl;
            cout << "If no, matrix values will be randomly generated" << endl;
            cin >> import_matrix;
        }
    }

    else
    {
        // csr sparse matrix has different text file format required
        cout << "If you choose to import matrix values from a text file, this program assumes" << endl;
        cout << "the text file is in the formats shown below" << endl;
        cout << endl;
        cout << "2) For sparse matrix, it is assumed that the matrix is stored in CSR format with non-zero values" << endl;
        cout << "followed by row position array and then column index array. For example, for the  4 x 4 matrix:" << endl;
        cout << "0 0 0 0" << endl;
        cout << "5 8 0 0" << endl;
        cout << "0 0 3 0" << endl;
        cout << "0 6 0 0" << endl;
        cout << "The input text file should read:" << endl;
        cout << "5 8 3 6" << endl;
        cout << "0 0 2 3 4" << endl;
        cout << "0 1 2 1" << endl;
        cout << "If both the vector and the matrix is imported, they must be provided in separate files" << endl;

        while (import_matrix != 1 && import_matrix != 2)
        {
            cout << "Do you want to import matrix values from text file? Enter 1 or 2 for" << endl;
            cout << "1. Yes" << endl;
            cout << "2. No" << endl;
            cout << "If no, matrix values will be randomly generated" << endl;
            cin >> import_matrix;
        }   
    }

    // ask use for matrix size i.e. rows/columns
    int size = -1;
    // in case use type in 0 or negative numbers
    while (size < 1)
    {
        cout << "Enter the size (n) of the square Matrix (n x n): ";
        cin >> size;
    }
        
    rows = size;
    cols = size;

    // allow use to choose which solver to use
    int solver = -1;

    if (matrix_type == 1)
    {
        // three solvers implemented for dense matrix
        // in case user type in anything other than 1, 2 or 3
        while (solver != 1 && solver != 2 && solver != 3)
        {
            cout << "What solver would you like to solve your problem with? Enter 1, 2 or 3 for:" << endl;
            cout << "1. Gauss-Seidel" << endl;
            cout << "2. Jacobi" << endl;
            cout << "3. LU Decomposition" << endl;
            cin >> solver;
        }
    }

    else
    {
        // two solvers implemented for csr sparse matrix
        // in case user type in anything other than 1 or 2
        while (solver != 1 && solver != 2)
        {
            cout << "What solver would you like to solve your problem with? Enter 1 or 2 for" << endl;
            cout << "1. Gauss-Seidel" << endl;
            cout << "2. Jacobi" << endl;
            cin >> solver;
        }
    }

    // allow use to input the number of iterations if iterative solver is chosen
    int iter;
    if (solver == 1 or solver == 2)
    {
        cout << "Number of iterations: ";
        cin >> iter;
        cout << endl;
    }

    string solver_type;
    double norm_resid;

    // Initialize array to store solver results
    unique_ptr<double[]> x(new double[cols]);
    for (int i = 0; i<cols; i++)
    {
        x[i] = 0.0;
    }

    // dense matrix
    if (matrix_type == 1)
    {       
        // declare matrix as unique pointer
        unique_ptr<Matrix<double> > A;
        // define vector as unique pointer and allocate column number of doubles
        unique_ptr<double[]> b(new double[cols]);
            
        // user chose to randomly generate matrix values
        if (import_matrix == 2)
        {
            A.reset(new Matrix<double>(rows, cols, false));
            A->GenerateRandomSPDMatrix(rows, cols, *A);
        }
        
        // use chose to import matrix values from text file
        else
        {
            // ask for text file name
            string matrix_fname;
            cout << "Enter file name that stores matrix data: ";
            cin >> matrix_fname;
            fstream f_matrix_in;
            f_matrix_in.open(matrix_fname, fstream::in);
            if (f_matrix_in.fail()) return -1;

            // import from text file
            double *value_ptr = new double[rows * cols];
            for (int i=0; i<rows*cols; i++)
            {
                f_matrix_in >> value_ptr[i];
            }
            f_matrix_in.close();
            cout << "Matrix read in from " << matrix_fname << endl;

            A.reset(new Matrix<double>(rows, cols, value_ptr));
            delete[] value_ptr;
        }

        // allow user to choose between importing vector and randomly generating vector
        int import_vector = -1;

        // in case user input anything other than 1 or 2
        while (import_vector != 1 && import_vector != 2)
        {
            cout << "Do you want to import vector values from text file? Enter 1 or 2 for:" << endl;
            cout << "1. Yes" << endl;
            cout << "2. No" << endl;
            cout << "If no, vector values will be randomly generated" << endl;
            cin >> import_vector;
        }

        // user chose to randomly generate vector
        if (import_vector == 2)
        {
            for (int i=0; i<cols; i++)
            {
                b[i] = rand() % 100 + 1;
            }
        }

        // user chose to import vector from text file
        else
        {
            // ask user for text file name
            string vector_fname;
            cout << "Enter file name that stores vector data: ";
            cin >> vector_fname;

            // import vector
            fstream f_vector_in;
            f_vector_in.open(vector_fname, fstream::in);
            if (f_vector_in.fail()) return -1;

            for (int i=0; i<cols; i++)
            {
                f_vector_in >> b[i];
            }
            f_vector_in.close();
            cout << "Vector read in from " << vector_fname << endl;
        }

        // user chose Gauss Seidel
        if (solver == 1)
        {
            A->GaussSeidel(b, x, iter);
            solver_type = "Gauss-Seidel";
            norm_resid = A->getNorm(x,b);
        }

        // user chose Jacobi
        else if (solver == 2)
        {
            A->Jacobi(b, x, iter);
            solver_type = "Jacobi";
            norm_resid = A->getNorm(x,b);
        }

        // user chose LU decomposition with partial pivoting
        else
        {
            A->LU_pp_solver(b, x);
            solver_type = "LU decomposition";
            norm_resid = A->getNorm(x,b);
        }

        cout << "Solver: " << solver_type << endl;
        cout << "L2 norm of the residual: " << norm_resid << endl;

        // allow user the option to save matrix to text file if matrix is randomly generated
        int store_matrix = -1;
        if (import_matrix == 2)
        {
            // in case the input is not 1 or 2
            while (store_matrix != 1 && store_matrix!= 2)
            {
                cout << "Would you like to store the matrix to text file? Enter 1 or 2 for:" << endl;
                cout << "1. Yes" << endl;
                cout << "2. No" << endl;
                cin >> store_matrix;
            }
        }

        // user wants to store randomly generated matrix
        if (store_matrix == 1)
        {
            // ask for file name
            string store_matrix_fname;
            cout << "Enter the file name where you would like to store the matrix: ";
            cin >> store_matrix_fname;
            cout << endl;

            // write to file
            // append in case user wants to add to an existing file
            fstream matrix_out;
            matrix_out.open(store_matrix_fname, fstream::out | fstream::app);

            for (int i=0; i<A->rows * A->cols; i++)
            {
                matrix_out << A->values[i] << "\t";
            }
            matrix_out << endl;
            
            matrix_out.close();
            cout << "Matrix stored in " << store_matrix_fname << endl;
        }
        
        
        // allow user the option to save vector to text file if matrix is randomly generated
        int store_vector = -1;
        if (import_vector == 2)
        {
            // in case the input is not 1 or 2
            while (store_vector != 1 && store_vector != 2)
            {
                cout << "Would you like to store the vector to text file? Enter 1 or 2 for:" << endl;
                cout << "1. Yes" << endl;
                cout << "2. No" << endl;
                cin >> store_vector;
            }
        }

        // user wants to store randomly generated vector
        if (store_vector == 1)
        {
            // ask for file name
            string store_vector_fname;
            cout << "Enter the file name where you would like to store the vector: ";
            cin >> store_vector_fname;
            cout << endl;

            // write to file
            // append in case user wants to add to an existing file
            fstream vector_out;
            vector_out.open(store_vector_fname, fstream::out | fstream::app);

            for (int i=0; i<A->cols; i++)
            {
                vector_out << b[i] << "\t";
            }
            vector_out << endl;

            vector_out.close();
            cout << "Vector stored in " << store_vector_fname << endl;
        }

        // allow user the option to store the solution
        int store_sol = -1;
        // in case input is anything other than 1 or 2
        while (store_sol != 1 && store_sol != 2)
        {
            cout << "Would you like to store the solution to text file? Enter 1 or 2 for:" << endl;
            cout << "1. Yes" << endl;
            cout << "2. No" << endl;
            cin >> store_sol;
        }

        // user wish to store solution
        if (store_sol == 1)
        {
            // ask for file name
            string store_sol_fname;
            cout << "Enter the file name where you would like to store the solution: ";
            cin >> store_sol_fname;
            cout << endl;
            
            // write to file
            // append in case user wants to add to an existing file
            fstream sol_out;
            sol_out.open(store_sol_fname, fstream::out | fstream::app);

            for (int i=0; i<A->cols; i++)
            {
                sol_out << x[i] << "\t";
            }
            sol_out << endl;

            sol_out.close();
            cout << "Solution stored in " << store_sol_fname << endl;
        }
        if (store_sol == 2)
        {
         cout << "Printing solution:";
         for (int i=0; i<A->cols; i++)
            {
                cout << x[i] << " ";
            }
            cout << endl;
        }
    }

    // csr sparse matrix
    else
    {
        // declare csr sparse matrix as unique pointer
        unique_ptr<CSRMatrix<double> > A;
        // define vector b as unique pointer and allocate column number of doubles
        unique_ptr<double[]> b(new double[cols]);

        // user wish to randomly generate matrix values
        if (import_matrix == 2)
        {
            A.reset(new CSRMatrix<double>(rows, cols, 1, true));
            A->GenerateRandomSPDMatrix(rows, cols, 1, *A);
        }

        // user wish to import matrix values
        else
        {
            // ask for file name
            string csr_fname;
            cout << "Enter file name that stores matrix data: ";
            cin >> csr_fname;

            // import data from file
            ifstream f_csr_in;
            f_csr_in.open(csr_fname);
            if (f_csr_in.fail()) return -1;

            // Initialize vector to store CSR values in string format
            vector<string> result;

            vector<double> v_values;
            vector<int> v_row_pos;
            vector<int> v_col_index;

            // Read the lines from the file as a string
            string input;
            while (getline (f_csr_in, input))
                result.push_back(input);

            f_csr_in.close();
            cout << "Matrix data read in from " << csr_fname << endl;

            // first row contains the non-zero values
            string str_values = result[0];
            stringstream ss1(str_values);
            string values_temp;

            // pass in as string and convert to double
            while (!ss1.eof())
            {
                ss1 >> values_temp;
                v_values.push_back(stod(values_temp));
            }

            // second row contains row position
            string str_row_pos = result[1];
            stringstream ss2(str_row_pos);
            string row_pos_temp;

            // pass in as string and convert to integer
            while (!ss2.eof())
            {
                ss2 >> row_pos_temp;
                v_row_pos.push_back(stoi(row_pos_temp));
            }

            // third row contains column index
            string str_col_index = result[2];
            stringstream ss3(str_col_index);
            string col_index_temp;

            // pass in as string and convert to integer
            while (!ss3.eof())
            {
                ss3 >> col_index_temp;
                v_col_index.push_back(stoi(col_index_temp));
            }

            // define matrix parameters from the data
            int nnzs = v_values.size();
            double *value_ptr = &v_values[0];
            int *row_pos = &v_row_pos[0];
            int *col_index = &v_col_index[0];

            // define the csr sparse matrix
            A.reset(new CSRMatrix<double>(rows, cols, nnzs, value_ptr, row_pos, col_index));
        }

        // allow user the option to import vector from text file
        int import_vector = -1;
        // in case the input is anything other than 1 or 2
        while (import_vector != 1 && import_vector != 2)
        {
            cout << "Do you want to import vector values from text file? Enter 1 or 2 for:" << endl;
            cout << "1. Yes" << endl;
            cout << "2. No" << endl;
            cout << "If no, vector values will be randomly generated" << endl;
            cin >> import_vector;
        }

        // user wish to randomly generate vector
        if (import_vector == 2)
        {
            for (int i=0; i<cols; i++)
            {
                b[i] = rand() % 100 + 1;
            }
        }

        // user wish to import vector from text file
        else
        {
            // ask for file name
            string csr_vector_fname;
            cout << "Enter file name that stores vector data: ";
            cin >> csr_vector_fname;

            // import data
            ifstream f_csr_vector_in;
            f_csr_vector_in.open(csr_vector_fname);
            if (f_csr_vector_in.fail()) return -1;

            for (int i=0; i<cols; i++)
            {
                f_csr_vector_in >> b[i];
            }
            f_csr_vector_in.close();

            cout << "Vector data read in from " << csr_vector_fname << endl;
        }

        // user chose Gauss-Seidel 
        if (solver == 1)
        {
            A->GaussSeidel(b, x, iter);
            solver_type = "Gauss-Seidel";
            norm_resid = A->getNorm(x,b);
        }

        // user chose Jacobi
        else
        {
            A->Jacobi(b, x, iter);
            solver_type = "Jacobi";
            norm_resid = A->getNorm(x,b);
        }

        cout << "Solver: " << solver_type << endl;
        cout << "L2 norm of the residual: " << norm_resid << endl;

        // allow user the option to save the matrix if matrix is randomly generated
        int store_matrix = -1;
        if (import_matrix==2)
        {
            // in case input is anything other than 1 or 2
            while (store_matrix != 1 && store_matrix!= 2)
            {
                cout << "Would you like to store the matrix to text file? Enter 1 or 2 for:" << endl;
                cout << "1. Yes" << endl;
                cout << "2. No" << endl;
                cin >> store_matrix;
            }
        }

        // user wish to save the randomly generated matrix
        if (store_matrix == 1)
        {
            // ask for file name
            string store_matrix_fname;
            cout << "Enter the file name where you would like to store the matrix: ";
            cin >> store_matrix_fname;
            cout << endl;

            // write to file
            // append in case user wants to add to an existing file
            fstream matrix_out;
            matrix_out.open(store_matrix_fname, fstream::out | fstream::app);

            for (int i=0; i<A->nnzs; i++)
            {
                matrix_out << A->values[i] << "\t";
            }
            matrix_out << endl;

            for (int i=0; i<A->rows+1; i++)
            {
                matrix_out << A->row_position[i] << "\t";
            }
            matrix_out << endl;

            for (int i=0; i<A->nnzs; i++)
            {
                matrix_out << A->col_index[i] << "\t";
            }
            matrix_out << endl;
            
            matrix_out.close();
            cout << "Matrix stored in " << store_matrix_fname << endl;
        }
        
        // allow user the option to save the randomly generated vector
        int store_vector = -1;
        if (import_vector == 2)
        {
            // in case input is anything other than 1 or 2
            while (store_vector != 1 && store_vector != 2)
            {
                cout << "Would you like to store the vector to text file? Enter 1 or 2 for:" << endl;
                cout << "1. Yes" << endl;
                cout << "2. No" << endl;
                cin >> store_vector;
            }
        }

        // user wish to save the randomly generated vector
        if (store_vector == 1)
        {
            // ask for file name
            string store_vector_fname;
            cout << "Enter the file name where you would like to store the vector: ";
            cin >> store_vector_fname;
            cout << endl;
            
            // save to file
            // append in case user wants to add to an existing file
            fstream vector_out;
            vector_out.open(store_vector_fname, fstream::out | fstream::app);

            for (int i=0; i<A->cols; i++)
            {
                vector_out << b[i] << "\t";
            }

            vector_out.close();
            cout << "Vector stored in " << store_vector_fname << endl;
        }

        // allow user the option to save the solution
        int store_sol = -1;
        // in case input is anything other than 1 and 2
        while (store_sol != 1 && store_sol != 2)
        {
            cout << "Would you like to store the solution to text file? Enter 1 or 2 for:" << endl;
            cout << "1. Yes" << endl;
            cout << "2. No" << endl;
            cin >> store_sol;
        }

        // user wish to save solution
        if (store_sol == 1)
        {
            // ask for file name
            string store_sol_fname;
            cout << "Enter the file name where you would like to store the solution: ";
            cin >> store_sol_fname;
            cout << endl;
            
            // save to file
            // append in case user wants to add to an existing file
            fstream sol_out;
            sol_out.open(store_sol_fname, fstream::out | fstream::app);

            for (int i=0; i<A->cols; i++)
            {
                sol_out << x[i] << "\t";
            }

            sol_out.close();
            cout << "Solution stored in " << store_sol_fname << endl;
        }

        // print out solution if user doesn't want to save solution
        if (store_sol == 2)
        {
            cout << "Printing solution:";
            for (int i=0; i<A->cols; i++)
            {
                cout << x[i] << " ";
            }
            cout << endl;
        }
    }
}
    
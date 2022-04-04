#pragma once
#include <vector>
#include "Matrix.h"

template <class T>
class CSRMatrix: public Matrix<T>
{
public:

   // constructor where we want to preallocate ourselves
   CSRMatrix(int rows, int cols, int nnzs, bool preallocate);
   // constructor where we already have allocated memory outside
   CSRMatrix(int rows, int cols, int nnzs, T *values_ptr, int *row_position, int *col_index);
   // destructor
   ~CSRMatrix();

   // Print out the values in our matrix
	virtual void printMatrix();

   // Perform some operations with our matrix
   void matMatMult(CSRMatrix<T>& mat_left, CSRMatrix<T>& result);
   // Perform some operations with our matrix
   void matVecMult(std::unique_ptr<T []> &vec, std::unique_ptr<T []> &output); 
   // Convert CSR matrix to dense format
   void CSRtoDense(Matrix<T>& output);
   // Transpose CSR matrix
   void transpose(CSRMatrix<T>& output);
   // Jacobi solver for CSR matrix
   void Jacobi(std::unique_ptr<T []> &b, std::unique_ptr<T []> &x, int iter);
   // Gauss-Seidel solver for CSR matrix
   void GaussSeidel(std::unique_ptr<T []> &b, std::unique_ptr<T []> &x, int iter);

   // get norm of the residual (A*x - b)
   double getNorm(std::unique_ptr<T []> &x, std::unique_ptr<T []> &b);

   // Generate a SPD (Symmetric Positive Definite) banded sparse matrix stored in CSR format
   void GenerateRandomSPDMatrix(int rows, int cols, int bandwidth, CSRMatrix &output);

   // overload the equal operator to assign one CSRMatrix to another CSRMatrix
   void operator = (const CSRMatrix<T> &mat)
   {
      this->rows = mat.rows;
      this->cols = mat.cols;
      this->nnzs = mat.nnzs;
      this->preallocated = mat.preallocated;
      this->size_of_values = mat.size_of_values;

      this->values.reset(new T[mat.nnzs]);
      for (int i=0; i<mat.nnzs; i++)
      {
         this->values[i] = mat.values[i];
      }

      this->row_position = new int[mat.rows + 1];
      for (int i=0; i<mat.rows+1; i++)
      {
         this->row_position[i] = mat.row_position[i];
      }

      this->col_index = new int[mat.nnzs];
      for (int i=0; i<mat.nnzs; i++)
      {
         this->col_index[i] = mat.col_index[i];
      }
   }

   // Explicitly using the C++11 nullptr here
   int *row_position = nullptr;
   int *col_index = nullptr;

   // How many non-zero entries we have in the matrix
   int nnzs=-1;

// Private variables - there is no need for other classes 
// to know about these variables
private:
   
};
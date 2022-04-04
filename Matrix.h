#pragma once
#include <memory>
#include <vector>

template <class T>
class Matrix
{
public:

   // constructor where we want to preallocate ourselves
   Matrix(int rows, int cols, bool preallocate);
   // constructor where we already have allocated memory outside
   Matrix(int rows, int cols, T *values_ptr);
   // destructor
   virtual ~Matrix();

   // Print out the values in our matrix
   void printValues();
	virtual void printMatrix();

   // Perform some operations with our matrix
   virtual void matMatMult(Matrix& mat_left, Matrix& output);

   // multiply a vector with our matrix
   virtual void matVecMult(std::unique_ptr<T []> &vec, std::unique_ptr<T []> &output);

   //transpose marix
   virtual void Transpose(Matrix &A);

   // Perform Gauss-Seidel to find result vector
   virtual void GaussSeidel(std::unique_ptr<T []> &b, std::unique_ptr<T []> &x, int iter, std::vector<T> *norm_resid=nullptr);

   // Perform Jacobi to find result vector
   virtual void Jacobi(std::unique_ptr<T []> &b, std::unique_ptr<T []> &x, int iter, std::vector<T> *norm_resid=nullptr);

   // solving linear system through LU decomposition with partial pivoting

   // swap jth and kth row in matrix for partial pivoting
   virtual void swap_rows(std::unique_ptr<T []> &A, int j, int k);
   // break down matrix to L and U with partial pivoting
   virtual void LU_decomposition_pp();
   // print lower triangular matrix
   virtual void printLowerTriangular();
   // print upper triangular matrix
   virtual void printUpperTriangular();
   // print permutation matrix
   virtual void printPermutation();

   // combine forward substitution and backward substitution to solve for x in the linear system
   // vector b and solution x in linear system is passed by reference
   virtual void LU_pp_solver(std::unique_ptr<T []> &b, std::unique_ptr<T []> &x);

   // swap_rows and LU_decomposition_pp could have been moved to private
   // but as a public method they can be accessed outside the class 
   // in the case where the user wants to break down matrix into P, L and U only

   // get norm of the residual (A*x - b)
   virtual double getNorm(std::unique_ptr<T []> &x, std::unique_ptr<T []> &b);

   // Generate a SPD (Symmetric Positive Definite) Matrix
   virtual void GenerateRandomSPDMatrix(int rows, int cols, Matrix& output);

   // overload + operator for matrix addition
   Matrix operator + (Matrix &mat)
   {
      Matrix res(mat.rows, mat.cols, true);
      for (int i=0; i<mat.size_of_values; i++)
      {
         res.values[i] = values[i] + mat.values[i];
      }
      return res;
   }

   // overload - operator for matrix subtraction
   Matrix operator - (Matrix &mat)
   {
      Matrix res(mat.rows, mat.cols, true);
      for (int i=0; i<mat.size_of_values; i++)
      {
         res.values[i] = values[i] - mat.values[i];
      }
      return res;
   }

   // Explicitly using the C++11 nullptr here
   // T *values = nullptr;   
   std::unique_ptr<T[]> values;
   int rows = -1;
   int cols = -1;

   // variables for LU decomposition
   std::unique_ptr<T[]> U; // pointer to store upper triangular matrix values
   std::unique_ptr<T[]> L; // pointer to store lower triangular matrix values
   std::unique_ptr<T[]> P_; // pointer to store permutation matrix values

   // variables for transpose
   std::unique_ptr<T[]> T_;

// Private variables - there is no need for other classes 
// to know about these variables
protected:
   int size_of_values = -1;
   bool preallocated = false;

private:
   bool decomposed = false; // true if matrix has been decomposed to L, U and P_
};
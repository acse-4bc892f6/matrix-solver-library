#include <iostream>
#include <cmath>
#include <vector>
#include "Matrix.h"

// Constructor - using an initialisation list here
template <class T>
Matrix<T>::Matrix(int rows, int cols, bool preallocate): rows(rows), cols(cols), size_of_values(rows * cols), preallocated(preallocate)
{
   // If we want to handle memory ourselves
   if (this->preallocated)
   {
      this->values.reset(new T[size_of_values]);
   }
}

// Constructor - now just setting the value of our double pointer
template <class T>
Matrix<T>::Matrix(int rows, int cols, T *values_ptr): rows(rows), cols(cols), size_of_values(rows * cols)
{
   this->values.reset(new T[this->size_of_values]);
   // copy values from values_ptr to values
   // don' use reset because memory is preallocated outside
   // we do not own the memory
   for (int i=0; i<this->size_of_values; i++)
   {
      values[i] = values_ptr[i];
   }
}

// destructor
template <class T>
Matrix<T>::~Matrix()
{}

// Just print out the values in our values array
template <class T>
void Matrix<T>::printValues() 
{ 
   std::cout << "Printing values" << std::endl;
	for (int i = 0; i< this->size_of_values; i++)
   {
      std::cout << this->values[i] << " ";
   }
   std::cout << std::endl;
}

// Explicitly print out the values in values array as if they are a matrix
template <class T>
void Matrix<T>::printMatrix() 
{ 
   std::cout << "Printing matrix" << std::endl;
   for (int j = 0; j< this->rows; j++)
   {  
      std::cout << std::endl;
      for (int i = 0; i< this->cols; i++)
      {
         // We have explicitly used a row-major ordering here
         std::cout << this->values[i + j * this->cols] << " ";
         // std::cout << this->identity_matrix_uptr[i + j * this->cols] << " ";
      }
   }
   std::cout << std::endl;
}

// Do matrix matrix multiplication
// output = mat_left * this
// m * k = m * n * n * k
template <class T>
void Matrix<T>::matMatMult(Matrix& mat_left, Matrix& output)
{

   // Check our dimensions match
   if (this->cols != output.cols)
   {
      std::cerr << "Input dimensions for matrices don't match" << std::endl;
      return;
   }

   // Check if our output matrix has had space allocated to it
   if (output.values != nullptr) 
   {
      // Check our dimensions match
      if (this->rows != mat_left.cols || mat_left.rows != output.rows)
      {
         std::cerr << "Input dimensions for matrices don't match" << std::endl;
         return;
      }      
   }
   // The output hasn't been preallocated, so we are going to do that
   else
   {
      T *tmp = new T[this->cols * mat_left.rows];
      output.values.reset(tmp);
      output.preallocated = true;
   }

   // Set values to zero before hand
   for (int i = 0; i < output.size_of_values; i++)
   {
      output.values[i] = 0;
   }

   // Now we can do our matrix-matrix multiplication
   // CHANGE THIS FOR LOOP ORDERING AROUND
   // AND CHECK THE TIME SPENT
   // Does the ordering matter for performance. Why??
   for(int i = 0; i < mat_left.rows; i++)
   {
      for(int j = 0; j < this->cols; j++)
      {
         for(int k = 0; k < mat_left.cols; k++)
         {            
               output.values[i * this->cols + j] += mat_left.values[i * mat_left.cols + k] * this->values[k * this->cols + j];
         }
      }
   }
}

// computes matrix vector multiplication
// output = this * vec
// n * 1 = n * k * k * 1
// pass vector and output by pointer
// vals is pointer for matrix values
template <class T>
void Matrix<T>::matVecMult(std::unique_ptr<T []> &vec, std::unique_ptr<T []> &output)
{
   if (vec.get() == nullptr)
   {
      std::cerr << "Vector hasn't been created" << std::endl;
      return;
   }
   if (output.get() == nullptr)
   {
      output.reset(new T[this->cols]);
   }

   for (int i=0; i<this->cols; i++)
   {
      output[i] = 0;
   }

   // carry out matrix vector multiplication

   // ROW-MAJOR ORDERING
   for (int i=0; i<this->rows; i++)
   {
      for (int j=0; j<this->cols; j++)
      {
         output[i] += this->values[i * cols + j] * vec[j];
      }
   }
}

template <class T>
void Matrix<T>::Transpose(Matrix &A)
{
   
    // Create a temporary matrix with reversed dimensions
   this->T_.reset(new T[this->size_of_values]);
    // Set the entries of B to be the same as those in A
      for (int i=0; i < this->size_of_values; i++)
      {
         for (int j=0; j < this->size_of_values; j++) 
         {
            this->T_[j + i * this->cols]=this->values[i + j * this->rows];
            }
         }
}

// Perform Gauss-Seidel method
// output x from Ax = b
template <class T>
void Matrix<T>::GaussSeidel(std::unique_ptr<T []> &b, std::unique_ptr<T []> &x, int iter, std::vector<T> *norm_resid)
{
   // check vector b has memory allocated
   if (b.get() == nullptr)
   {
      std::cerr << "b has no memory allocated" << std::endl;
      return;
   }
   // if x has no memory allocated, allocate here
   if (x.get() == nullptr)
   {
      x.reset(new T[this->cols]);
   }

   // iterate until residual becomes less than 1e-10 or until number of iterations is reached
   for (int k = 0; k < iter; k++)
   {
      for (int i = 0; i < this->rows; i++)
      {
         double dot_product = 0.0;
         for (int j = 0; j < this->cols; j++)
         {
            if (j != i) 
            {
               dot_product += this->values[i * this->cols + j] * x[j];
            }
         }
         x[i] = 1/this->values[i * this->cols + i] * (b[i] - dot_product);
      }

      double norm = this->getNorm(x, b);

      // Check if pointer to vector argument for storing the norm residual is a nullptr
      if (norm_resid) norm_resid->push_back(norm);

      if (norm < 1.0e-10) break;
   }
}

// perform Ax = b
// Jacobi method
template <class T>
void Matrix<T>::Jacobi(std::unique_ptr<T []> &b, std::unique_ptr<T []> &x, int iter, std::vector<T> *norm_resid)
{
   // check vector b has memory allocated
   if (b.get() == nullptr)
   {
      std::cerr << "b has no memory allocated" << std::endl;
      return;
   }
   // if x has no memory allocated, allocate here
   if (x.get() == nullptr)
   {
      x.reset(new T[this->cols]);
   }

   // iterate until residual becomes less than 1e-10 or until number of iterations is reached
   for (int k = 0; k < iter; k++)
   {
      std::unique_ptr<T []> x_new(new T[this->cols]);

      for (int i = 0; i < this->rows; i++)
      {
         double dot_product = 0.0;
         for (int j = 0; j < this->cols; j++)
         {
            if (j != i) 
            {
               dot_product += this->values[i * this->cols + j] * x[j];
            }
         }
         x_new[i] = 1/this->values[i * this->cols + i] * (b[i] - dot_product);
      }

      T *rawPtr = x_new.release();
      x.reset(rawPtr);

      double norm = this->getNorm(x, b);

      // Check if pointer to vector argument for storing the norm residual is a nullptr
      if (norm_resid) norm_resid->push_back(norm);

      if (norm < 1.0e-10) break;
      
   } 
}

template <class T>
void Matrix<T>::swap_rows(std::unique_ptr<T []> &A, int j, int k)
{
    // ensure matrix is not nullptr
    if (A.get() == nullptr)
    {
        std::cerr << "Matrix has no memory allocated to it" << std::endl;
        return;
    }

    // allocate memory to store rows as unique pointers
    std::unique_ptr<T[]> B(new T[this->cols]);
    std::unique_ptr<T[]> C(new T[this->cols]);

    // copy values in those rows to the memory allocated
    for (int i=0; i<this->cols; i++)
    {
        // copy jth row to B
        B[i] = A[j*this->rows + i];
        // copy kth row to C
        C[i] = A[k*this->rows + i];
    }

    // carry out the swap
    for (int i=0; i<this->cols; i++)
    {
        // kth row becomes jth row
        A[k*this->rows + i] = B[i];
        // jth row becomes kth row
        A[j*this->rows + i] = C[i];
    }
}

template <class T>
void Matrix<T>::LU_decomposition_pp()
{
    // allocate memory to store values for upper triangular matrix U
    // and lower triangular matrix L
    this->U.reset(new T[this->size_of_values]);
    this->L.reset(new T[this->size_of_values]);

    // initial values of U are matrix values
    // initial values of L are zero
    for (int i=0; i<this->size_of_values; i++)
    {
        this->U[i] = this->values[i];
        this->L[i] = 0;
    }

    // initial permutation matrix is the identity matrix
    this->P_.reset(new T[this->size_of_values]);
    for (int i=0; i<this->rows; i++)
    {
        for (int j=0; j<this->cols; j++)
        {
           if (i == j) P_[i*this->cols + j] = 1;
           else P_[i*this->cols + j] = 0;
        }
    }

    for (int k=0; k<(this->cols-1); k++)
    {
        // allocate memory to store values of upper triangular matrix 
        // in kth column, from kth row onwards to last row
        std::unique_ptr<T[]> temp(new T[this->rows - k]);

        for (int i=0; i<(this->rows - k); i++)
        {
            temp[i] = this->U[k + (i+k) * this->cols];
        }

        int j = 0; // index of the largest absolute value in temp

        // go over temp values
        for (int i=0; i<(this->rows - k); i++)
        {
            // replace j with i if the absolute value of the
            // ith element is larger than that of the jth element
            if (abs(temp[i]) > abs(temp[j]))
            {
                j = i;
            }
        }

        // j is counted from kth row onwards, so need to add k to j
        j += k;

        // swap jth and kth rows for upper triangular matrix, permutation matrix
        // and lower triangular matrix
        swap_rows(this->U, j, k);
        swap_rows(this->P_, j, k);
        swap_rows(this->L, j, k);

        // value of diagonal entry of upper triangular matrix
        T diag = this->U[k*this->rows+k];
        // go over entries from (k+1)th row onwards to the end row
        for (int i=k+1; i<this->rows; i++)
        {
            // divide entries from (k+1)th row onwards i.e. row i by the diagonal
            // these are the entries for lower triangular matrix
            T s = this->U[i * this->cols + k] / diag;

            // update entries to upper triangular matrix from column k onwards to last column at each row i
            for (int n=k; n<this->cols; n++)
            {
                this->U[i * this->cols + n] = this->U[i * this->cols + n] - s * this->U[k * this->cols + n];
            }

            // index remains the same as in upper triangular matrix
            this->L[i * this->cols + k] = s;
        }
    }

    // COULD USE DAXPY FOR MATRIX ADDITION???
    // try to optimise process by only performing operation on diagonal elements
   for (int i=0; i<this->rows; i++)
   {
      for (int j=0; j<this->cols; j++)
      {
         if (i == j) L[i*this->cols + j] += 1;
      }
   }

   // matrix is decomposed to L, U and P_
   this->decomposed = true;
}

template <class T>
void Matrix<T>::printLowerTriangular()
{
    // only able to print if matrix is decomposed
    if (!this->decomposed)
    {
       std::cerr << "Matrix is not deceomposed" << std::endl;
       return;
    }
    
    std::cout << "Printing lower triangular matrix:" << std::endl;
    
    // ROW-MAJOR ORDERING
    for (int j=0; j<this->rows; j++)
    {
        for (int i=0; i<this->cols; i++)
        {
            std::cout << this->L[i + j * this->cols] << " ";
        }
        std::cout << std::endl;
    }
}

template <class T>
void Matrix<T>::printUpperTriangular()
{
    // only able to print if matrix is decomposed
    if (!this->decomposed)
    {
       std::cerr << "Matrix is not deceomposed" << std::endl;
       return;
    }
    
    std::cout << "Printing upper triangular matrix:" << std::endl;
    
    // ROW-MAJOR ORDERING
    for (int j=0; j<this->rows; j++)
    {
        for (int i=0; i<this->cols; i++)
        {
            std::cout << this->U[i + j * this->cols] << " ";
        }
        std::cout << std::endl;
    }
}

template <class T>
void Matrix<T>::printPermutation()
{
    // only able to print if matrix is decomposed
    if (!this->decomposed)
    {
       std::cerr << "Matrix is not deceomposed" << std::endl;
       return;
    }
    
    std::cout << "Printing permutation matrix:" << std::endl;

    // ROW-MAJOR ORDERING
    for (int j=0; j<this->rows; j++)
    {
        for (int i=0; i<this->cols; i++)
        {
            std::cout << this->P_[i + j * this->cols] << " ";
        }
        std::cout << std::endl;
    }

}

template <class T>
void Matrix<T>::LU_pp_solver(std::unique_ptr<T []> &b, std::unique_ptr<T []> &x)
{
    // check vector b has memory allocated
   if (b.get() == nullptr)
   {
      std::cerr << "b has no memory allocated" << std::endl;
      return;
   }
   // if x has no memory allocated, allocate here
   if (x.get() == nullptr)
   {
      x.reset(new T[this->cols]);
   }
    
    this->LU_decomposition_pp();

   // linear system A * x = P_ * b
   // initialise P_ * b
    std::unique_ptr<T []> Pb(new T[this->cols]);

    for (int i=0; i<this->cols; i++)
    {
        Pb[i] = 0;
    }
    
    // compute P_ * b
    for (int i=0; i<this->rows; i++)
    {
        for (int j=0; j<this->cols; j++)
        {
            Pb[i] += this->P_[i * cols + j] * b[j];
        }
    }

    // A * x = b
    // L * (U * x) = P_ * b
    // compute for y in L * y = P_ * b where y = U * x

    // set initial y value to zero
    std::unique_ptr<T []> y(new T[this->cols]);
    for (int i=0; i<this->cols; i++)
    {
        y[i] = 0;
    }
   // forward substitution to solve for y where y = U * x
    for (int k=0; k<this->cols; k++)
    {
        T diag = this->L[k*this->rows+k];

        T s = 0;
        for (int j=0; j<k; j++)
        {
            s = s + this->L[k * this->cols + j] * y[j];
        }
        y[k] = (Pb[k] - s) / diag;
    }
    
    // compute for x in U * x = y
    // set initial x value to zero
    for (int i=0; i<this->cols; i++)
    {
        x[i] = 0;
    }
   // backward substitution
    for (int k=(this->cols-1); k>=0; k--)
    {
        T diag = this->U[k*this->cols+k];
        
        T s = 0;
        for (int j=0; j<this->cols; j++)
        {
            s = s + this->U[k * this->cols + j] * x[j];
        }
        x[k] = (y[k] - s) /  diag;
    }
}

template <class T>
double Matrix<T>::getNorm(std::unique_ptr<T []> &x, std::unique_ptr<T []> &b)
{
   // output = A * x, compute using matVecMult
   std::unique_ptr<T []> output(new T[this->cols]);

   this->matVecMult(x, output);

   // define residual and fill it with zeros initially
   std::unique_ptr<T []> residual(new T[this->cols]);

   for (int i=0; i<this->cols; i++)
   {
      residual[i] = 0.0;
   }

   // residual = A * x - b
   for (int i=0; i<this->cols; i++)
   {
      residual[i] = output[i] - b[i];
   }

   // compute L2 norm of residual
   double norm = 0.0;
   for (int i=0; i<this->rows; i++)
   {
      norm += pow(residual[i], 2);
   }
   norm = sqrt(norm);

   return norm;
}


template <class T>
void Matrix<T>::GenerateRandomSPDMatrix(int rows, int cols, Matrix& output)
{
   // use time as seed to generate different values each time this function is run
   srand(time(0));

   // Cholesky decomposition: real positive definite matrix = L * L.T
   // define L and L.T as unique pointers
   std::unique_ptr<Matrix<double> > L_temp(new Matrix<double>(this->rows, this->cols, true)); // L
   std::unique_ptr<Matrix<double> > LT_temp(new Matrix<double>(this->rows, this->cols, true)); // L.T

   // set initial values in L and L.T to zero
   for (int i=0; i<rows*cols; i++)
   {
      L_temp->values[i] = 0.0;
      LT_temp->values[i] = 0.0;
   }

   int count = 0;
   // loop over diagonals and fill it with randomly generaged positive number
   for (int i=0; i<rows*cols; i+=rows+1)
   {
      // generate positive large number to ensure a diagonally dominant matrix
      L_temp->values[i] = rand() % 100 + 100;
      // loop from beginning of row up to diagonal and 
      // fill them with smaller randomly generated positive numbers
      // to obtain diagonally dominant lower triangular matrix
      for (int j=(count*cols); j<i; j++)
      {
            L_temp->values[j] = rand() % 10 + 1;
      }

      count++;
   }

   // find L.T
   for (int i=0; i<cols; i++)
   {
      for (int j=0; j<rows; j++)
      {
            LT_temp->values[i*cols + j] = L_temp->values[i + j*cols];
      }
   }
   
   // A = L * L.T, computer from matMatMult
   LT_temp->matMatMult(*L_temp, output);

   // NOTE: MATRIX GENERATED IS LIKELY TO HAVE SPECTRAL RADIUS LARGER THAN ONE
   // SO WILL LEAD TO DIVERGENCE OF SOLUTION IF JACOBI SOLVER IS USED

   // attempted to generate matrix through eigendecomposition A = P * Lambda * P.T
   // where P is a matrix with columns as normalised eigenvectors
   // and Lambda is a matrix with diagonal filled with eigenvalues between 0 and 1 exclusive
   // SPD matrices are generated but solution still diverged when Jacobi solver is used
}


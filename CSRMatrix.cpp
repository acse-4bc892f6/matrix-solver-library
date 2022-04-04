#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>
#include "CSRMatrix.h"

using namespace std;

// constructor for csrmatrix
// we are calling the base contructor for matrix here with a false for the preallocate
// that way we don't accidently allocate rows*cols space for the values array
template <class T>
CSRMatrix<T>::CSRMatrix(int rows, int cols, int nnzs, bool preallocate): Matrix<T>(rows, cols, false), nnzs(nnzs)
{
   // If we don't pass false in the initialisation list base constructor, it would allocate values to be of size
   // rows * cols in our base matrix class
   // So then we need to set it to the real value we had passed in
   this->preallocated = preallocate;
   // this is set to be rows*cols in our parent class Matrix
   // which is not correct for our CSRMatrix
   this->size_of_values = nnzs;

   // If we want to handle memory ourselves
   if (this->preallocated)
   {
      this->values.reset(new T[this->nnzs]);

      // Must remember to delete this in the destructor
      this->row_position = new int[this->rows+1];
      this->col_index = new int[this->nnzs];
   }
}

// Constructor - now just setting the value of our T pointer
template <class T>
CSRMatrix<T>::CSRMatrix(int rows, int cols, int nnzs, T *values_ptr, int *row_position, int *col_index): Matrix<T>(rows, cols, values_ptr), nnzs(nnzs)
{
   // this is set to be rows*cols in our parent class Matrix
   // which is not correct for our CSRMatrix
   this->size_of_values = nnzs;  

   // copy values from the outside pointer to the local pointer instead of directly pointing to it
   // since we don't own the memory
   this->row_position = new int[rows+1];
   for (int i=0; i<rows+1; i++)
   {
      this->row_position[i] = row_position[i];
   }

   this->col_index = new int[nnzs];
   for (int i=0; i<nnzs; i++)
   {
      this->col_index[i] = col_index[i];
   }
}

// destructor
template <class T>
CSRMatrix<T>::~CSRMatrix()
{
   // Delete the values array
   if (this->preallocated){
      delete[] this->row_position;
      delete[] this->col_index;
   }
   // The super destructor is called after we finish here
   // This will delete this->values if preallocated is true
}

// Explicitly print out the values in values array as if they are a matrix
template <class T>
void CSRMatrix<T>::printMatrix() 
{ 
   std::cout << "Printing matrix" << std::endl;
   std::cout << "Values: ";
   for (int j = 0; j< this->nnzs; j++)
   {  
      std::cout << this->values[j] << " ";      
   }
   std::cout << std::endl;
   std::cout << "row_position: ";
   for (int j = 0; j< this->rows+1; j++)
   {  
      std::cout << this->row_position[j] << " ";      
   }
   std::cout << std::endl;   
   std::cout << "col_index: ";
   for (int j = 0; j< this->nnzs; j++)
   {  
      std::cout << this->col_index[j] << " ";      
   }
   std::cout << std::endl;   
}

// Do a matrix-vector product
// output = this * input
template<class T>
void CSRMatrix<T>::matVecMult(std::unique_ptr<T []> &vec, std::unique_ptr<T []> &output)
{
   if (vec == nullptr || output == nullptr)
   {
      std::cerr << "Input or output haven't been created" << std::endl;
      return;
   }

   // Set the output to zero
   for (int i = 0; i < this->rows; i++)
   {
      output[i] = 0.0;
   }

   // Loop over each row
   for (int i = 0; i < this->rows; i++)
   {
      // Loop over all the entries in this col
      for (int val_index = this->row_position[i]; val_index < this->row_position[i+1]; val_index++)
      {
         // This is an example of indirect addressing
         // Can make it harder for the compiler to vectorise!
         output[i] += this->values[val_index] * vec[this->col_index[val_index]];

      }
   }
}

// matmatmult
// preallocated in result must be false to prevent memory leak when overloaded equal operator is used
// to copy data from CSRMatrix object in this function, which will be deleted when it goes out of scope,
// to the result passed by reference 
template<class T>
void CSRMatrix<T>::matMatMult(CSRMatrix<T>& mat_left, CSRMatrix<T>& result)
{
   int output_val_index = 0;
   vector<T> output_val;
   vector<int> output_row_position = {0};
   vector<int> output_col_index;
   // Loop over each row in mat_left
   for (int i = 0; i < mat_left.rows; i++)
   {
      // output.row_position[i] = output_val_index;
      // Loop over each col in mat_right
      for (int j = 0; j < this->cols; j++)
      {
         // Initialize product for current matrix element
         T product = 0.0;
         // Loop through all non-zero elements in current row of mat_left
         for (int left_val_index = mat_left.row_position[i]; left_val_index < mat_left.row_position[i+1]; left_val_index++)
         {
            // Loop through column index of mat_right
            for (int k = 0; k < this->nnzs; k++)
            {
               // Check for current column in mat_right we are evaluating
               if (this->col_index[k] == j)
               {
                  // For each column index of the column of interest find the row it is located in mat_right
                  int* p = std::upper_bound( this->row_position, this->row_position + this->rows, k);

                  int row_pos_greater_than = p - this->row_position;
                  if (row_pos_greater_than - 1 == mat_left.col_index[left_val_index])
                  {
                     
                     product += mat_left.values[left_val_index] * this->values[k];
                  }
               }
            }
         }
         if (product > 0) 
         {
            output_val.push_back(product);
            output_col_index.push_back(j);
            output_val_index++;
         }
      }
      output_row_position.push_back(output_val_index);
   }

   int output_nnzs = output_val_index;
   // Convert vectors to arrays so they can be passed into CSRMatrix constructor
   T *val_arr = &output_val[0];
   int *col_index_arr = &output_col_index[0];
   int *row_position_arr = &output_row_position[0];

   // Initialize output CSRMatrix after determining nnzs within this function
   // CSRMatrix<T> output(mat_left.rows, this->cols, output_nnzs, true); 
   // Fill CSRMatrix values with values from vectors

   CSRMatrix<T> output(mat_left.rows, this->cols, output_nnzs, true); 
   output.nnzs = output_nnzs;
   for (int i = 0; i < output_nnzs; i++)
   {
      output.values[i] = output_val[i];
      output.col_index[i] = output_col_index[i];
   }
   for (int i = 0; i < output.rows + 1; i++)
   {
      output.row_position[i] = output_row_position[i];
   }

   result = output; // use overloaded equal operator to copy data from output to result
   return;
}
// function to convert CSR matrix to Dense format
template<class T>
void CSRMatrix<T>::CSRtoDense(Matrix<T>& output)
{
   // Set values in output array to zero by default
   for (int i = 0; i < output.rows * output.cols; i++)
   {
      output.values[i] = 0.0;
   }
   for (int i = 0; i < this->rows; i++)
   {
      for (int val_position = this->row_position[i]; val_position < this->row_position[i+1]; val_position++)
      {
         output.values[i * output.cols + this->col_index[val_position]] = this->values[val_position];
      }
   }
   return;
}

// Transpose CSR matrix to get CSR output
template<class T>
void CSRMatrix<T>::transpose(CSRMatrix<T>& output)
{
   // Set values in output array to zero by default
   for (int i = 0; i < output.rows * output.cols; i++)
   {
      output.values[i] = 0.0;
   }

   int output_value_index = 0;
   int output_row_index = 0;
   for (int i = 0; i < this->cols; i++)
   {
      output.row_position[output_row_index] = output_value_index;
      for (int j = 0; j < this->nnzs; j++)
      {
         int input_col = this->col_index[j];
         if (i == input_col)
         {
            T input_val = this->values[j];
            
            int* p = std::upper_bound(this->row_position, this->row_position + this->rows, j);

            int row_pos_greater_than = p - this->row_position;
            int input_row = row_pos_greater_than - 1;
            
            output.values[output_value_index] = input_val;
            output.col_index[output_value_index] = input_row;
            output_value_index++;
         }
      }
      output_row_index++;
   }
   output.row_position[output_row_index] = output_value_index;
   return;
}

// perform Ax = b
// Jacobi method
template <class T>
void CSRMatrix<T>::Jacobi(std::unique_ptr<T []> &b, std::unique_ptr<T []> &x, int iter)
{
   if (b.get() == nullptr)
   {
      std::cerr << "b has no memory allocated" << std::endl;
      return;
   }
   if (x.get() == nullptr)
   {
      x.reset(new T[this->cols]);
   }

   for (int k = 0; k < iter; k++)
   {
      std::unique_ptr<T []> x_new(new T[this->cols]);

      for (int i = 0; i < this->rows; i++)
      {
         T A_ii = 0.0;
         T dot_product = 0.0;
         for (int val_index = this->row_position[i]; val_index < this->row_position[i+1]; val_index++)
         {
            if (this->col_index[val_index] != i) 
            {
               dot_product += this->values[val_index] * x[this->col_index[val_index]];
            }
            else
            {
               A_ii = this->values[val_index];
            }
         }
         x_new[i] = 1/A_ii * (b[i] - dot_product);
      }

      T *rawPtr = x_new.release();
      x.reset(rawPtr);

      double norm = this->getNorm(x, b);

      // Check if pointer to vector argument for storing the norm residual is a nullptr
      // if (norm_resid) norm_resid->push_back(norm);

      if (norm < 1.0e-10) break;
      
   } 
}
// Perform Gauss-Seidel method
// output x from Ax = b
template <class T>
void CSRMatrix<T>::GaussSeidel(std::unique_ptr<T []> &b, std::unique_ptr<T []> &x, int iter)
{
   if (b.get() == nullptr)
   {
      std::cerr << "b has no memory allocated" << std::endl;
      return;
   }
   if (x.get() == nullptr)
   {
      x.reset(new T[this->cols]);
   }

   for (int k = 0; k < iter; k++)
   {
      for (int i = 0; i < this->rows; i++)
      {
         T A_ii = 0.0;
         T dot_product = 0.0;
         for (int val_index = this->row_position[i]; val_index < this->row_position[i+1]; val_index++)
         {
            if (this->col_index[val_index] != i) 
            {
               dot_product += this->values[val_index] * x[this->col_index[val_index]];
            }
            else
            {
               A_ii = this->values[val_index];
            }
         }
         x[i] = 1/A_ii * (b[i] - dot_product);
      }

      double norm = this->getNorm(x, b);

      // Check if pointer to vector argument for storing the norm residual is a nullptr
      // if (norm_resid) norm_resid->push_back(norm);

      if (norm < 1.0e-10) break;
   }
}

template <class T>
double CSRMatrix<T>::getNorm(std::unique_ptr<T []> &x, std::unique_ptr<T []> &b)
{
   // output = A * x
   std::unique_ptr<T []> output(new T[this->cols]);
   for (int i=0; i<this->cols; i++)
   {
      output[i] = 0;
   }

   this->matVecMult(x, output);

   std::unique_ptr<T []> residual(new T[this->cols]);

   for (int i=0; i<this->cols; i++)
   {
      residual[i] = 0.0;
   }

   for (int i=0; i<this->cols; i++)
   {
      residual[i] = output[i] - b[i];
   }

   double norm = 0.0;
   for (int i=0; i<this->rows; i++)
   {
      norm += pow(residual[i], 2);
   }
   norm = sqrt(norm);

   return norm;
}

template <class T>
void CSRMatrix<T>::GenerateRandomSPDMatrix(int rows, int cols, int bandwidth, CSRMatrix &output)
{
   // use time as seed to generate different values each time this function is run
   srand(time(0));
   
   // initialise number of non zeroes
   // bandwidth is the number of non-zero rows below the diagonal of the banded matrix
   int nnzs = rows;
   for (int i=0; i<bandwidth; i++)
   {
      nnzs += (rows-(i+1));
   }

   // Cholesky decomposition: real positive definite matrix = L * L.T

   // initialise L
   std::unique_ptr<CSRMatrix<T> > L(new CSRMatrix<T>(rows, cols, nnzs, true));

   // row position
    int row_increment = 0;
    int row_val = 0;
    for (int i=0; i<rows+1; i++)
    {
        L->row_position[i] = row_val;

        if (row_increment < (bandwidth+1))
        {
            row_increment++;
        }
        row_val += row_increment;
    }

   // col index

   int index = 0; // index for col index array
   int extra = 0; // quantity added to index for col index array
   
   // go over row position array
   for (int i=0; i<rows; i++)
   {
      // difference betwen i+1th and ith entry in row position
      // is the number of entries in col index
      int num_entry = L->row_position[i+1] - L->row_position[i];

      for (int j=index; j<index+num_entry; j++)
      {
         L->col_index[j] = j-index+extra;
         
         if (j == (index + num_entry - 1)) 
         {
               // generate larger values at diagonal to ensure 
               // matrix is diagonally dominant
               L->values[j] = rand() % 100 + 100;
         }
         else 
         {
               L->values[j] = rand() % 10 + 1;
         }
      }

      // shift index by the number of entries
      index += num_entry;

      // when number of entry reaches (bandwidth + 1) need to add column index
      // by one over each loop because of the zero values
      if (num_entry==bandwidth+1) extra++;
   }

   // initialise L.T
   std::unique_ptr<CSRMatrix<T> > LT(new CSRMatrix<T>(rows, cols, nnzs, true));
   // transpose L to obtain L.T
   L->transpose(*LT);

   // A = L * L.T
   std::unique_ptr<CSRMatrix<T> > A(new CSRMatrix<T>(rows, cols, nnzs, true));

   LT->matMatMult(*L, *A);

   output = *A;

   return;
}
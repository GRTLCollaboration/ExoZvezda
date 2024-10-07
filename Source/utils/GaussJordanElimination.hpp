/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef GAUSSJORDANELIMINATION_HPP_
#define GAUSSJORDANELIMINATION_HPP_

#include "Matrix.hpp"
#include <stdexcept>
#include <utility> // For std::swap

// Linear equation solution by Gauss-Jordan elimination

class GaussJordan
{
  private:
    Matrix &a; // reference to the input matrix a
    Matrix &b; // reference to the input matrix b
    int num_rows;
    int num_cols;
    std::vector<int> col_index;
    std::vector<int> row_index;
    std::vector<int> pivot_used; // a book-keeping vector

  public:
    // Constructor
    GaussJordan(Matrix &a_matrix, Matrix &b_matrix)
        : a(a_matrix), b(b_matrix), num_rows(a_matrix.nrows()),
          num_cols(b_matrix.ncols()), col_index(num_rows), row_index(num_rows),
          pivot_used(num_rows, 0)
    {
    }

    // Function to perform the Gauss-Jordan elimination
    int solve()
    {
        for (int row = 0; row < num_rows; row++)
        {
            double largest_val = 0.0;
            int temp_row = -1, temp_col = -1;

            // Search for the largest pivot element
            for (int col = 0; col < num_rows; col++)
            {
                if (pivot_used[col] != 1)
                {
                    for (int row2 = 0; row2 < num_rows; row2++)
                    {
                        if (pivot_used[row2] == 0)
                        {
                            double current_val = abs(a.At(col, row2));
                            if (current_val >= largest_val)
                            {
                                largest_val = current_val;
                                temp_row = col;
                                temp_col = row2;
                            }
                        }
                    }
                }
            }

            ++pivot_used[temp_col];

            // Swap rows to move the pivot into place
            if (temp_row != temp_col)
            {
                swap_rows(temp_row, temp_col);
            }

            row_index[row] = temp_row;
            col_index[row] = temp_col;

            // Check if matrix is singular
            if (a.At(temp_col, temp_col) == 0.0)
            {
                MayDay::Error(
                    "Matrix is singular in GaussJordanElimination.hpp!");
                return 0;
            }

            // Normalize the pivot row
            normalise_row(temp_col);

            // Eliminate other rows
            eliminate_rows(temp_col);
        }

        // Unscramble the solution
        unscramble_solution();
        return 1;
    }

  private:
    // Helper function to swap rows in matrices `a` and `b`
    void swap_rows(int row1, int row2)
    {
        for (int col = 0; col < num_rows; col++)
        {
            std::swap(a.At(row1, col), a.At(row2, col));
        }
        for (int col = 0; col < num_cols; col++)
        {
            std::swap(b.At(row1, col), b.At(row2, col));
        }
    }

    // Helper function to normalize a row based on the pivot element
    void normalise_row(int pivot_row)
    {
        double pivot_inv = 1.0 / a.At(pivot_row, pivot_row);
        a.At(pivot_row, pivot_row) = 1.0;

        // Normalize matrix a
        for (int col = 0; col < num_rows; col++)
        {
            a.At(pivot_row, col) *= pivot_inv;
        }

        // Normalize matrix b
        for (int col = 0; col < num_cols; col++)
        {
            b.At(pivot_row, col) *= pivot_inv;
        }
    }

    // Helper function to eliminate rows based on the pivot row
    void eliminate_rows(int pivot_row)
    {
        for (int row = 0; row < num_rows; row++)
        {
            if (row != pivot_row)
            {
                double factor = a.At(row, pivot_row);
                a.At(row, pivot_row) = 0.0;

                // Eliminate in matrix a
                for (int col = 0; col < num_rows; col++)
                {
                    a.At(row, col) -= a.At(pivot_row, col) * factor;
                }

                // Eliminate in matrix b
                for (int col = 0; col < num_cols; col++)
                {
                    b.At(row, col) -= b.At(pivot_row, col) * factor;
                }
            }
        }
    }

    // Helper function to unscramble the solution after elimination
    void unscramble_solution()
    {
        for (int row = num_rows - 1; row >= 0; row--)
        {
            if (row_index[row] != col_index[row])
            {
                for (int swap_row = 0; swap_row < num_rows; swap_row++)
                {
                    std::swap(a.At(swap_row, row_index[row]),
                              a.At(swap_row, col_index[row]));
                }
            }
        }
    }
};

// Overload for the homogeneous system, i.e., b is an empty matrix
class GaussJordanHomogeneous : public GaussJordan
{
  public:
    GaussJordanHomogeneous(Matrix &a_matrix)
        : GaussJordan(a_matrix, *(new Matrix(a_matrix.nrows(), 0)))
    {
    }
};

#endif /* GAUSSJORDANELIMINATION_HPP_ */
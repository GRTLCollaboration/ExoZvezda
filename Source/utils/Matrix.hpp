/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef MATRIX_HPP_
#define MATRIX_HPP_

#include <vector>

// A simple class for a matrix (only used in the non-linear fitting routine)

class Matrix
{
  public:
    Matrix() : data({}) {}

    Matrix(const int &rows, const int &cols) { Reset(rows, cols); }

    void Reset(const int &rows, const int &cols)
    {
        data.resize(rows);
        for (int i = 0; i < rows; ++i)
        {
            data.at(i).resize(cols);
        }
    }

    // Value at given row and column
    double At(const int &row, const int &col) const
    {
        return data.at(row).at(col);
    }

    double &At(const int &row, const int &col) { return data.at(row).at(col); }

    // Number of rows
    int nrows() const { return data.size(); }

    // Number of columns
    int ncols() const
    {
        if (nrows() > 0)
        {
            return data[0].size();
        }

        return 0;
    }

  private:
    std::vector<std::vector<double>> data;
};

#endif /* MATRIX_HPP_ */

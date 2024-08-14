#include <stdexcept>

#include "DimensionDefinitions.hpp"
#include "Matrix.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "VarsTools.hpp"
#include "simd.hpp"

int gaussj(Matrix &a, Matrix &b)
{
    int i, icol, irow, j, k, l, ll, n = a.nrows(), m = b.ncols();
    double big, dum, pivinv;
    std::vector<double> indxc(n), indxr(n), ipiv(n);
    for (j = 0; j < n; j++)
        ipiv[j] = 0;
    for (i = 0; i < n; i++)
    {
        big = 0.0;
        for (j = 0; j < n; j++)
            if (ipiv[j] != 1)
                for (k = 0; k < n; k++)
                {
                    if (ipiv[k] == 0)
                    {
                        if (abs(a.At(j, k)) >= big)
                        {
                            big = abs(a.At(j, k));
                            irow = j;
                            icol = k;
                        }
                    }
                }
        ++(ipiv[icol]);
        if (irow != icol)
        {
            for (l = 0; l < n; l++)
                SWAP(a.At(irow, l), a.At(icol, l));
            for (l = 0; l < m; l++)
                SWAP(b.At(irow, l), b.At(icol, l));
        }
        indxr[i] = irow;
        indxc[i] = icol;
        if (a.At(icol, icol) == 0.0)
            return 0;
        // throw std::runtime_error("gaussj: Singular Matrix");
        pivinv = 1.0 / a.At(icol, icol);
        a.At(icol, icol) = 1.0;
        for (l = 0; l < n; l++)
            a.At(icol, l) *= pivinv;
        for (l = 0; l < m; l++)
            b.At(icol, l) *= pivinv;
        for (ll = 0; ll < n; ll++)
            if (ll != icol)
            {
                dum = a.At(ll, icol);
                a.At(ll, icol) = 0.0;
                for (l = 0; l < n; l++)
                    a.At(ll, l) -= a.At(icol, l) * dum;
                for (l = 0; l < m; l++)
                    b.At(ll, l) -= b.At(icol, l) * dum;
            }
    }
    for (l = n - 1; l >= 0; l--)
    {
        if (indxr[l] != indxc[l])
            for (k = 0; k < n; k++)
                SWAP(a.At(k, indxr[l]), a.At(k, indxc[l]));
    }
    return 1;
}

int gaussj(Matrix &a)
{
    Matrix b(a.nrows(), 0);
    int gsuccess = gaussj(a, b);
    if (gsuccess == 0)
    {
        return 0;
    }
    else
    {
        return 1;
    }
}

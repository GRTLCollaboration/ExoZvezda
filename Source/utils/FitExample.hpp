// #include "DimensionDefinitions.hpp"
// #include "Tensor.hpp"
// #include "TensorAlgebra.hpp"
// #include "VarsTools.hpp"
// #include "simd.hpp"

void fgauss(const double x, std::vector<double> &a, double &y,
            std::vector<double> &dyda)
{
    int i, na = a.size();
    double fac, ex, arg;
    y = 0.;
    for (i = 0; i < na - 1; i += 3)
    {
        arg = (x - a[i + 1]) / a[i + 2];
        ex = exp(-SQR(arg));
        fac = a[i] * ex * 2. * arg;
        y += a[i] * ex + a[i + 3];
        dyda[i] = ex;
        dyda[i + 1] = fac / a[i + 2];
        dyda[i + 2] = fac * arg / a[i + 2];
    }
}

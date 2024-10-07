/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef FITTINGMETHOD_HPP_
#define FITTINGMETHOD_HPP_

#include "GaussJordanElimination.hpp"
#include "Matrix.hpp"
#include <stdexcept>
#include <utility>

// Class for Levenberg-Marquardt nonlinear fitting

// Essentially, this routine aims to apply nonlinear least squares algorithm
// between data and some fitting model of choice. Here we have data arrays "x"
// and "y", along with standard deviation array "error", quantifying the
// error/uncertainty. We model the data with some non-linear function comprised
// of coefficients stored in the "fit_parameters" array; usually you would start
// with some initial guess for these model parameters, that the fitting routine
// will iterate over. The routine (hopefully) returns best-fit parameters in
// "fit_parameters" array.

// The general form of the model may be written as y(x_i) = \sum_k a_k X_k(x_i),
// where X_k are the basis functions and a_k are the coefficients. We want to
// minimise the \chi^2 defined by:
//  \chi^2 = sum_i [(y_i - y(x_i)) / stddev]^2, where y_i defines some 1D NR
//  data.

// The rough algorithm is as follows:
// 1. Linearise the problem around some small perturbations to the model
// parameters (\delta a). So then, the end task is to solve for (\delta a) in
// \alpha (\delta a) = \beta, where \alpha and \beta are as defined below.
// 2. Compute matrix alpha_{kj}, which is the Hessian-like curvature matrix of
// the chi-squared surface with respect to the model parameters. In short, it is
// a second-order partial derivative matrix given by \sum_i (dy(x_i)/da_k *
// dy(x_i)/da_j)/stddev_i^2.
// 3. Compute the vector beta_k, which is the right-hand side of the normal
// equations used to compute the parameter adjustments, it is given by \sum_i
// (y_i - y(x_i) / stddev_i^2 * dy(x_i)/da_k.
// 4. Once we know alpha and beta, compute (\delta a) and update the model
// parameters with a ---> a + (\delta a). 5.

class FittingMethod
{
  public:
    static const int MAX_ITERATIONS = 5000;
    static const int CONVERGENCE_CRITERIA = 4;

    int data_size, num_parameters;
    std::vector<double> &x_values, &y_values, &errors; // x array, y array
    double tolerance;
    void (*model_func)(
        const double, std::vector<double> &, double &,
        std::vector<double>
            &); // fitting model (see e.g. Gaussian model in FittingModel.hpp)
    std::vector<double> fit_parameters; // array for coefficients of the fit
    double chi_square;                  // how well the fitting did

    FittingMethod(std::vector<double> &x_data, std::vector<double> &y_data,
                  std::vector<double> &sigma,
                  std::vector<double> &initial_params,
                  void model_function(const double, std::vector<double> &,
                                      double &, std::vector<double> &),
                  const double tol = 1.e-6)
        : data_size(x_data.size()), num_parameters(initial_params.size()),
          x_values(x_data), y_values(y_data), errors(sigma), tolerance(tol),
          model_func(model_function), fit_parameters(initial_params)
    {
    }

    // Main function to perform the fit
    int fit()
    {
        int iter, done = 0;
        double lambda = 0.001, prev_chi_square; // damping term
        std::vector<double> trial_params(num_parameters),
            beta_vector(num_parameters), delta_params(num_parameters);
        Matrix alpha_matrix(num_parameters, num_parameters),
            covariance_matrix(num_parameters, num_parameters);
        Matrix solution_vector(num_parameters, 1),
            temp_beta_matrix(num_parameters, num_parameters);

        // Compute alpha and beta variables, as discussed in the compute_mrqcof
        // definition
        compute_mrqcof(fit_parameters, alpha_matrix, beta_vector);
        trial_params = fit_parameters;
        prev_chi_square = chi_square;

        for (iter = 0; iter < MAX_ITERATIONS; iter++)
        {
            if (done == CONVERGENCE_CRITERIA)
            {
                lambda = 0.0; // Force final iteration
            }

            // Adjust covariance matrix and solve
            for (int j = 0; j < num_parameters; j++)
            {
                for (int k = 0; k < num_parameters; k++)
                {
                    covariance_matrix.At(j, k) = alpha_matrix.At(j, k);
                }
                covariance_matrix.At(j, j) *=
                    (1.0 +
                     lambda); // to improve stability the Levenberg-Marquardt
                              // algorithm modifies the update step by adding a
                              // damping term (\lambda here) to the diagonal of
                              // the Hessian matrix.
                temp_beta_matrix.At(j, 0) = beta_vector[j];
            }

            // Perform Gauss-Jordan elimination
            GaussJordan gj_solver(covariance_matrix, temp_beta_matrix);
            int success = gj_solver.solve();
            if (success == 0)
                return 0;

            // Write into (\deta a)
            for (int j = 0; j < num_parameters; j++)
            {
                delta_params[j] = temp_beta_matrix.At(j, 0);
            }

            if (done == CONVERGENCE_CRITERIA)
            {
                reorder_alpha_matrix(covariance_matrix);
                reorder_alpha_matrix(alpha_matrix);
                return 1;
            }

            // Update trial parameters
            for (int l = 0, j = 0; l < num_parameters; l++)
            {
                trial_params[l] = fit_parameters[l] + delta_params[j++];
            }

            // Recompute alpha and beta
            compute_mrqcof(trial_params, covariance_matrix, delta_params);
            if (std::abs(chi_square - prev_chi_square) <
                std::max(tolerance, tolerance * chi_square))
            {
                done++;
            }

            if (chi_square < prev_chi_square)
            {
                lambda *= 0.1;
                prev_chi_square = chi_square;
                alpha_matrix = covariance_matrix;
                beta_vector = delta_params;
                fit_parameters = trial_params;
            }
            else
            {
                lambda *= 10.0; // increase the damping
                chi_square = prev_chi_square;
            }
        }
        MayDay::Error(
            "I am out of the maximum number of iterations in the fitting in "
            "FittingMethod.hpp!");
        return 0; // Too many iterations
    }

  private:
    // Matrix alpha_{kj} is the Hessian-like curvature matrix of the chi-squared
    // surface with respect to the model parameters. In short, it is a
    // second-order partial derivative matrix given by \sum_i (dy(x_i)/da_k *
    // dy(x_i)/da_j)/stddev_i^2. The vector beta_k is the right-hand side of the
    // normal equations used to compute the parameter adjustments, it is given
    // by \sum_i (y_i - y(x_i) / stddev_i^2 * dy(x_i)/da_k. Finally, chi-square
    // qunatifies the goodness of the fit performed.
    void compute_mrqcof(std::vector<double> &params, Matrix &alpha,
                        std::vector<double> &beta)
    {
        int i, j, k;
        double model_val, sigma_sq_inv, delta_y;
        std::vector<double> model_derivatives(num_parameters);

        // Initiliase
        for (j = 0; j < num_parameters; j++)
        {
            for (k = 0; k <= j; k++)
                alpha.At(j, k) = 0.0;
            beta[j] = 0.0;
        }

        // Compute alpha and beta
        chi_square = 0.0;
        for (i = 0; i < data_size; i++)
        {
            model_func(x_values[i], params, model_val, model_derivatives);
            sigma_sq_inv = 1.0 / (errors[i] * errors[i]);
            delta_y = y_values[i] - model_val;

            for (j = 0, k = 0; k < num_parameters; k++)
            {
                double weight = model_derivatives[k] * sigma_sq_inv;
                for (int l = 0; l <= k; l++)
                {
                    alpha.At(j, l) += weight * model_derivatives[l];
                }
                beta[j++] += delta_y * weight;
            }
            // Compute chi_square
            chi_square += delta_y * delta_y * sigma_sq_inv;
        }

        // Symmetrise
        for (j = 1; j < num_parameters; j++)
        {
            for (k = 0; k < j; k++)
            {
                alpha.At(k, j) = alpha.At(j, k);
            }
        }
    }

    // Reorder covariance matrix after the fitting
    void reorder_alpha_matrix(Matrix &alpha)
    {
        int i, j, k = num_parameters - 1;
        for (j = num_parameters - 1; j >= 0; j--)
        {
            for (i = 0; i < num_parameters; i++)
            {
                std::swap(alpha.At(i, k), alpha.At(i, j));
                std::swap(alpha.At(k, i), alpha.At(j, i));
            }
            k--;
        }
    }
};

#endif /* FITTINGMETHOD_HPP_ */

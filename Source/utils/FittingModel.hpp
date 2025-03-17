/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef FITTINGMODEL_HPP_
#define FITTINGMODEL_HPP_

// A Gaussian model to be used for fitting. The model is parametrised by
// amplitude A, mean mu and standard deviation sigma (and an off-set term), so
// we have A * exp[-({x - mu} / sigma)^2] + some_offset. In this model we will
// be fitting for these 3 coefficients, or aka model parameters.

// The function that stores 'coefficients' (i.e. the values for model
// parameters) and 'dyda' (i.e. the derivatives of the model with respect to the
// model params)
void gaussian_model(const double x_value, std::vector<double> &coefficients,
                    double &y_model, std::vector<double> &dyda)
{
    int num_coefficients = coefficients.size();
    double factor, exponential_part, argument;
    y_model = 0.0;

    for (int i = 0; i < num_coefficients - 1; i += 3)
    {
        // Parameters are expected to be in sets of 3: amplitude, mean, and
        // standard deviation
        double amplitude = coefficients[i];
        double mean = coefficients[i + 1];
        double stddev = coefficients[i + 2];

        // Calculate the argument for the Gaussian
        argument = (x_value - mean) / stddev;
        exponential_part =
            exp(-argument * argument); // exp(-((x - mean) / stddev)^2)
        factor = amplitude * exponential_part * 2.0 *
                 argument; // Term for derivative

        // Update the model value (y_model)
        y_model += amplitude * exponential_part +
                   coefficients[i + 3]; // coefficients[i+3] is the offset term

        // Compute partial derivatives for each parameter
        dyda[i] = exponential_part;    // Partial derivative w.r.t amplitude
        dyda[i + 1] = factor / stddev; // Partial derivative w.r.t mean
        dyda[i + 2] =
            factor * argument / stddev; // Partial derivative w.r.t stddev
    }
}

#endif /* FITTINGMODELMODEL_HPP_ */
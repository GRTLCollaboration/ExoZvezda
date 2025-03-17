/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef ADMMASSEXTRACTION_HPP_
#define ADMMASSEXTRACTION_HPP_

#include "SphericalExtraction.hpp"
//! This class extracts the ADM mass integrand values, which should be a grid
//! variable, on spherical shells and then integrates them over the spheres
//! It then writes the calculated integrals to a file
class ADMMassExtraction : public SphericalExtraction
{
  public:
    //! The constructor
    ADMMassExtraction(spherical_extraction_params_t &a_params, double a_dt,
                      double a_time, bool a_first_step,
                      double a_restart_time = 0.0)
        : SphericalExtraction(a_params, a_dt, a_time, a_first_step,
                              a_restart_time)
    {
        add_var(c_Madm, VariableType::diagnostic);
    }

    //! Execute the query
    void execute_query(AMRInterpolator<Lagrange<4>> *a_interpolator)
    {
        std::vector<double> adm_integral;

        extract(a_interpolator);

        if (m_params.write_extraction)
        {
            write_extraction("ADMMassExtractionOut_");
        }

        // add integrand
        add_var_integrand(0, adm_integral, IntegrationMethod::simpson);

        // integrate
        integrate();

        // write integrals
        std::string integrals_filename = "ADMmass";
        write_integral(integrals_filename, adm_integral, "ADM mass");
    }
};

#endif /* ADMMASSEXTRACTION_HPP_ */
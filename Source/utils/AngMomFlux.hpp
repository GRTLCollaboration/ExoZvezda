/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef ANGMOMFLUX_HPP_
#define ANGMOMFLUX_HPP_

#include "AMRInterpolator.hpp"
#include "InterpolationQuery.hpp"
#include "Lagrange.hpp"
#include "SmallDataIO.hpp" // for writing data
#include "UserVariables.hpp" // Needs c_mod_phi etc
#include "SimulationParametersBase.hpp"

#include "SphericalExtraction.hpp"

//////
#include <fstream>
#include <string>
#include <vector>
#include <cmath> // for NAN
//////
class AngMomFlux : SphericalExtraction
{
    // variables for AngMomFlux
    double m_dt, m_restart_time;
    string m_filename = "AngMomFlux";
    double m_time;
    bool m_first_step;
    spherical_extraction_params_t m_params;

  public:

      AngMomFlux(spherical_extraction_params_t &a_params, double a_time,
                         double a_dt, double a_restart_time, bool a_first_step)
                         : SphericalExtraction(a_params, a_dt, a_time, a_first_step,
                                               a_restart_time),
                            m_params(a_params),
                            m_time(a_time), m_dt(a_dt), m_restart_time(m_restart_time),
                            m_first_step(a_first_step)
      {
          add_var(c_Fphi_flux, VariableType::diagnostic);
      }

    void run(AMRInterpolator<Lagrange<4>> *a_interpolator)
    {
        std::vector<double> integrals;
        auto integrand = [](std::vector<double> mom_flux_vals, double r,
                                                     double theta, double phi){
            return mom_flux_vals[0];
        };

        // extract the values of the Weyl scalars on the spheres
        extract(a_interpolator);
        add_integrand(integrand, integrals);
        integrate();
        write_to_dat(integrals);
    }

    void write_to_dat(std::vector<double> vals)
    {

        std::vector<string> title_line(m_params.num_extraction_radii);
        string dummy_string;
        for (int i=0; i<m_params.num_extraction_radii; i++)
        {
            dummy_string = "r = " + to_string(m_params.extraction_radii[i]);
            title_line[i] = dummy_string;
        }

        SmallDataIO flux_file(m_filename, m_dt, m_time,
                                      m_restart_time,
                                      SmallDataIO::APPEND,
                                      m_first_step);

        if (m_time > 0) flux_file.remove_duplicate_time_data();

        if (m_time == 0.)
        {
            flux_file.write_header_line(title_line);
        }

        flux_file.write_time_data_line(vals);

    }

    ~AngMomFlux(){;}

    /*void setup_integration_surface();
    void interpolate_over_surface(double r, AMRInterpolator<Lagrange<4>> *a_interpolator);
    void integrate_over_surface(int index);
    void write_to_dat();*/

  public:

    //void run(AMRInterpolator<Lagrange<4>> *a_interpolator);
};

#endif /* ANGMOMFLUX_HPP_ */

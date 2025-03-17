/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef COMPUTEWEIGHTFUNCTION_HPP_
#define COMPUTEWEIGHTFUNCTION_HPP_

#include "BosonStarParams.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "DebuggingTools.hpp"
#include "VarsTools.hpp"
#include "WeightFunction.hpp"
#include "simd.hpp"

/*
 * Class to compute the weight function, useful for visualisation mid run
 */

class ComputeWeightFunction
{
  protected:
    BosonStar_params_t m_params_BosonStar;
    BosonStar_params_t m_params_BosonStar2;
    double m_dx;

  public:
    template <class data_t> struct weightfunc_t
    {
        data_t weight1;
        data_t weight2;
    };

    ComputeWeightFunction(BosonStar_params_t a_params_BosonStar,
                          BosonStar_params_t a_params_BosonStar2,
                          const double a_dx)
        : m_dx(a_dx), m_params_BosonStar(a_params_BosonStar),
          m_params_BosonStar2(a_params_BosonStar2)
    {
    }

    //! Function to compute the value of the initial vars on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        weightfunc_t<data_t> out;

        Coordinates<double> coords(current_cell, m_dx,
                                   m_params_BosonStar.star_centre);

        double separation = m_params_BosonStar.BS_separation;
        double impact_parameter = m_params_BosonStar.BS_impact_parameter;
        double q = m_params_BosonStar.mass_ratio;
        double rapidity = m_params_BosonStar.BS_rapidity;
        double rapidity2 = m_params_BosonStar2.BS_rapidity;
        double radius_width1 = m_params_BosonStar.radius_width1;
        double radius_width2 = m_params_BosonStar.radius_width2;

        WeightFunction weightfunction;

        double profile_func1 = weightfunction.profile_chi(
            (coords.x - q * separation / (q + 1)) * cosh(rapidity),
            coords.y + q * impact_parameter / (q + 1.), coords.z,
            radius_width1);
        double profile_func2 = weightfunction.profile_chi(
            (coords.x + separation / (q + 1)) * cosh(-rapidity2),
            coords.y - impact_parameter / (q + 1.), coords.z, radius_width2);

        current_cell.store_vars(profile_func1, c_profile1);
        current_cell.store_vars(profile_func2, c_profile2);
    }
};

#endif /* COMPUTEWEIGHTFUNCTION_HPP_ */

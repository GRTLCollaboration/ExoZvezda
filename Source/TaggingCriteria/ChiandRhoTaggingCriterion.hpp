/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef CHIANDPHITAGGINGCRITERION_HPP_
#define CHIANDPHITAGGINGCRITERION_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "DiagnosticVariables.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "Tensor.hpp"

class ChiandRhoTaggingCriterion
{
  protected:
    const double m_dx;
    const FourthOrderDerivatives m_deriv;
    const extraction_params_t m_params;
    const int m_level;
    const double m_threshold_rho;
    const double m_threshold_chi;

    template <class data_t>
    using MatterVars = typename ComplexScalarField<>::template Vars<data_t>;

    /// Vars object for chi
    template <class data_t> struct Vars
    {
        data_t chi; //!< Conformal factor
        data_t rho;

        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
            using namespace VarsTools; // define_enum_mapping is part of
                                       // VarsTools
            define_enum_mapping(mapping_function, c_chi, chi);
            define_enum_mapping(mapping_function, c_rho, rho);
        }
    };

  public:
    ChiandRhoTaggingCriterion(const double a_dx, const int a_level,
                              const extraction_params_t a_params,
                              const double a_threshold_rho,
                              const double a_threshold_chi)
        : m_dx(a_dx), m_deriv(a_dx), m_params(a_params), m_level(a_level),
          m_threshold_rho(a_threshold_rho), m_threshold_chi(a_threshold_chi){};

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        const auto d2 = m_deriv.template diff2<MatterVars>(current_cell);
        const auto d2chi = m_deriv.template diff2<Vars>(current_cell);
        const auto d2rho = m_deriv.template diff2<Vars>(current_cell);

        data_t mod_d2_chi = 0;
        data_t mod_d2_rho = 0;

        FOR(idir, jdir)
        {
            mod_d2_chi += d2chi.chi[idir][jdir] * d2chi.chi[idir][jdir];

            mod_d2_rho += d2rho.rho[idir][jdir] * d2rho.rho[idir][jdir];
        }

        data_t criterion_chi = m_dx / m_threshold_chi * sqrt(mod_d2_chi);

        data_t criterion_rho = m_dx / m_threshold_rho * sqrt(mod_d2_rho);

        data_t criterion = simd_max(criterion_chi, criterion_rho);

        for (int iradius = 0; iradius < m_params.num_extraction_radii;
             ++iradius)
        {
            // regrid if within extraction level and not at required refinement
            if (m_level < m_params.extraction_levels[iradius])
            {
                const Coordinates<data_t> coords(current_cell, m_dx,
                                                 m_params.extraction_center);
                const data_t r = coords.get_radius();
                // add a 20% buffer to extraction zone so not too near to
                // boundary
                auto regrid = simd_compare_lt(
                    r, 1.2 * m_params.extraction_radii[iradius]);
                criterion = simd_conditional(regrid, 100.0, criterion);
            }
        }

        // Write back into the flattened Chombo box
        current_cell.store_vars(criterion, 0);
    }
};

#endif /* CHIANDPHITAGGINGCRITERION_HPP_ */

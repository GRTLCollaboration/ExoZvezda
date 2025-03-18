/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef BOSONCHIPUNCTUREEXTRACTIONTAGGINGCRITERION_HPP_
#define BOSONCHIPUNCTUREEXTRACTIONTAGGINGCRITERION_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "SphericalExtraction.hpp"
#include "Tensor.hpp"

/* This class tags cells based on three criteria - the
 value of the second derivs, the extraction regions
 and the puncture 'horizons'. There are no horizon for BSs, but within a
 certain radius of the star, we will aim to place a box around it. Use of
 'horizon' here is not to be taken literally.

 This tagging follows the implementation/idea of Miren Radia, see Section of
 https://arxiv.org/abs/2112.10567

 */
class BosonChiPunctureExtractionTaggingCriterion
{
  protected:
    const double m_dx;
    const int m_level;
    const std::array<int, 2> m_horizon_max_levels;
    const bool m_track_punctures;
    const bool m_activate_extraction;
    const FourthOrderDerivatives m_deriv;
    const spherical_extraction_params_t m_params;
    const std::vector<double> m_puncture_radii;
    const std::vector<double> m_puncture_masses;
    const std::vector<std::array<double, CH_SPACEDIM>> &m_puncture_coords;
    const double m_buffer;

  public:
    template <class data_t> struct Vars
    {
        data_t chi; //!< Conformal factor

        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
            using namespace VarsTools; // define_enum_mapping is part of
                                       // VarsTools
            define_enum_mapping(mapping_function, c_chi, chi);
        }
    };

    // The constructor
    BosonChiPunctureExtractionTaggingCriterion(
        const double dx, const int a_level,
        const std::array<int, 2> a_horizon_max_levels,
        const spherical_extraction_params_t a_params,
        const std::vector<std::array<double, CH_SPACEDIM>> &a_puncture_coords,
        const bool activate_extraction = false,
        const bool track_punctures = false,
        const std::vector<double> a_puncture_radii = {4.0, 4.0},
        const std::vector<double> a_puncture_masses = {1.0, 1.0},
        const double a_buffer = 0.5)
        : m_dx(dx), m_level(a_level),
          m_horizon_max_levels(a_horizon_max_levels),
          m_track_punctures(track_punctures),
          m_activate_extraction(activate_extraction), m_deriv(dx),
          m_params(a_params), m_puncture_radii(a_puncture_radii),
          m_puncture_masses(a_puncture_masses),
          m_puncture_coords(a_puncture_coords), m_buffer(a_buffer){};

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // first test the gradients for regions of high curvature
        const auto d2 = m_deriv.template diff2<Vars>(current_cell);
        data_t mod_d2_chi = 0;
        FOR(idir, jdir)
        {
            mod_d2_chi += d2.chi[idir][jdir] * d2.chi[idir][jdir];
        }
        data_t criterion = m_dx * sqrt(mod_d2_chi);

        // if extracting weyl data at a given radius, enforce a given resolution
        // there
        if (m_activate_extraction)
        {
            for (int iradius = 0; iradius < m_params.num_extraction_radii;
                 ++iradius)
            {
                // regrid if within extraction level and not at required
                // refinement
                if (m_level < m_params.extraction_levels[iradius])
                {
                    const Coordinates<data_t> coords(current_cell, m_dx,
                                                     m_params.center);
                    const data_t r = coords.get_radius();
                    // add a 20% buffer to extraction zone so not too near to
                    // boundary
                    auto regrid = simd_compare_lt(
                        r, 1.2 * m_params.extraction_radii[iradius]);
                    criterion = simd_conditional(regrid, 100.0, criterion);
                }
            }
        }

        if (m_track_punctures == 1)
        {
            double puncture_separation2 = 0.0;

            std::array<double, CH_SPACEDIM> puncture_centre1;
            std::array<double, CH_SPACEDIM> puncture_centre2;

            FOR1(i)
            {
                puncture_centre1[i] = m_puncture_coords[0][i];
                puncture_centre2[i] = m_puncture_coords[1][i];
            }

            FOR1(idir)
            {
                double displacement =
                    puncture_centre1[idir] - puncture_centre2[idir];
                puncture_separation2 += displacement * displacement;
            }

            double puncture_separation = sqrt(puncture_separation2);

            const int merger_horizon_max_level =
                min(m_horizon_max_levels[0], m_horizon_max_levels[1]);

            if (puncture_separation >
                (m_puncture_radii[1] / 2. +
                 m_puncture_radii[0] /
                     2.)) // you may want to also try a puncture_separation >
                          // (m_puncture_radii[1] + m_puncture_radii[0])
            {
                // punctures still far enough apart so tag around each one
                // separately
                for (int ipuncture = 0; ipuncture < 2; ++ipuncture)
                {
                    std::array<double, CH_SPACEDIM> puncture_centre;

                    FOR1(i)
                    {
                        puncture_centre[i] = m_puncture_coords[ipuncture][i];
                    }
                    // where am i?
                    const Coordinates<data_t> coords(current_cell, m_dx,
                                                     puncture_centre);

                    const data_t max_abs_xy =
                        simd_max(abs(coords.x), abs(coords.y));
                    const data_t max_abs_xyz =
                        simd_max(max_abs_xy, abs(coords.z));

                    if (m_level < m_horizon_max_levels[ipuncture])
                    {
                        // we want the 2nd and 3rd levels above
                        // puncture_max_level to be twice the size of the next
                        // finest level
                        const double factor =
                            pow(2.0, min(m_horizon_max_levels[ipuncture] -
                                             m_level - 1,
                                         2));

                        auto regrid = simd_compare_lt(
                            max_abs_xyz,
                            factor *
                                (m_puncture_radii[ipuncture] / 2. + m_buffer));
                        // you can also use
                        // factor * (m_puncture_masses[ipuncture] * 2. +
                        // m_buffer) in the simd_compare_lt call.

                        criterion = simd_conditional(regrid, 100.0, criterion);
                    }
                    else
                    {
                        // remove any finer levels for BHs with
                        auto dont_regrid = simd_compare_lt(
                            max_abs_xyz,
                            m_puncture_radii[ipuncture] / 2. + m_buffer);
                        criterion =
                            simd_conditional(dont_regrid, 0.0, criterion);
                    }
                }
            }

            else
            {
                if (m_level >= merger_horizon_max_level)
                {
                    // drop any finer levels after merger
                    criterion = 0.0;
                }
            }

            double sum_masses = m_puncture_masses[0] + m_puncture_masses[1];

            // if punctures are close enough together tag cells at the
            // center of mass for the remnant
            if (puncture_separation <
                m_puncture_radii[0] / 2. + m_puncture_radii[1] / 2. +
                    m_buffer) // you may want to also try a puncture_separation
                              // > (m_puncture_radii[1] + m_puncture_radii[0])
            {
                std::array<double, CH_SPACEDIM> center_of_mass;
                FOR1(idir)
                {
                    center_of_mass[idir] =
                        (m_puncture_masses[0] * m_puncture_coords[0][idir] +
                         m_puncture_masses[1] * m_puncture_coords[1][idir]) /
                        sum_masses;
                }

                Coordinates<data_t> coords(current_cell, m_dx, center_of_mass);
                const data_t max_abs_xy =
                    simd_max(abs(coords.x), abs(coords.y));
                const data_t max_abs_xyz = simd_max(max_abs_xy, abs(coords.z));

                if (m_level < merger_horizon_max_level)
                {
                    const double factor = pow(
                        2.0, min(merger_horizon_max_level - m_level - 1, 2));
                    auto regrid2 = simd_compare_lt(
                        max_abs_xyz,
                        factor * (m_puncture_radii[0] / 2. +
                                  m_puncture_radii[1] / 2. + m_buffer));
                    criterion = simd_conditional(regrid2, 100.0, criterion);
                }
            }
        }

        // Write back into the flattened Chombo box
        current_cell.store_vars(criterion, 0);
    }
};

#endif /* BOSONCHIPUNCTUREEXTRACTIONTAGGINGCRITERION_HPP_ */

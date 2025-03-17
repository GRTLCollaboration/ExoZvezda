/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef EMTENSORANDFLUXES_HPP
#define EMTENSORANDFLUXES_HPP

#include "CCZ4Geometry.hpp"
#include "Cell.hpp"
#include "FourthOrderDerivatives.hpp"
#include "Interval.H"
#include "MatterCCZ4.hpp"
#include "simd.hpp"

/*
 * Calculates the EM tensor and angular momentum flux
 */

template <class matter_t> class EMTensorAndFluxes
{
  public:
    template <class data_t>
    using Vars = typename MatterCCZ4<matter_t>::template Vars<data_t>;

    //! Constructor
    EMTensorAndFluxes(const matter_t &a_matter, const double dx,
                      const double a_L,
                      const std::array<double, CH_SPACEDIM> a_centre,
                      const int a_c_rho = -1, const int a_c_Fphi_flux = -1,
                      const int a_c_Sphi_source = -1,
                      const int a_c_Qphi_density = -1,
                      const Interval a_c_Si = Interval(),
                      const Interval a_c_Sij = Interval());

    template <class data_t> void compute(Cell<data_t> current_cell) const;

  protected:
    const double m_L, m_dx; // simulation boxsize and cell size
    const matter_t &m_matter;
    FourthOrderDerivatives m_deriv;
    const int m_c_rho;         // var enum for the energy density
    const int m_c_Fphi_flux;   // angular momentum flux
    const int m_c_Sphi_source; // angular source momentum
    const int m_c_Qphi_density;
    const Interval m_c_Si;  // Interval of var enums for the momentum density
    const Interval m_c_Sij; // Interval of var enums for the spatial part
                            // stress-energy density
    const std::array<double, CH_SPACEDIM>
        m_centre; // centre of momentum flux calculation
};

#include "EMTensorAndFluxes.impl.hpp"

#endif /* EMTENSORANDFLUXES_HPP */

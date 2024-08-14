/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(COMPLEXSCALARFIELD_HPP_)
#error "This file should only be included through ComplexScalarField.hpp"
#endif

#ifndef COMPLEXSCALARFIELD_IMPL_HPP_
#define COMPLEXSCALARFIELD_IMPL_HPP_
#include "DimensionDefinitions.hpp"

// Calculate the stress energy tensor elements
template <class potential_t>
template <class data_t, template <typename> class vars_t>
emtensor_t<data_t> ComplexScalarField<potential_t>::compute_emtensor(
    const vars_t<data_t> &vars, const vars_t<Tensor<1, data_t>> &d1,
    const Tensor<2, data_t> &h_UU, const Tensor<3, data_t> &chris_ULL) const
{
    emtensor_t<data_t> out;

    // Copy the field vars into SFObject
    CSFObject<data_t> vars_csf;
    vars_csf.phi_Re = vars.phi_Re;
    vars_csf.phi_Im = vars.phi_Im;
    vars_csf.Pi_Re = vars.Pi_Re;
    vars_csf.Pi_Im = vars.Pi_Im;

    // call the function which computes the em tensor excluding the potential
    emtensor_excl_potential(out, vars, vars_csf, d1.phi_Re, d1.phi_Im, h_UU,
                            chris_ULL);

    // set the default potential values
    data_t V_of_modulus_phi_squared = 0.0;
    data_t dVdmodulus_phi_squared = 0.0;

    // compute potential and add constributions to EM Tensor
    my_potential.compute_potential(V_of_modulus_phi_squared,
                                   dVdmodulus_phi_squared, vars);

    out.rho += 0.5 * V_of_modulus_phi_squared;
    out.S += -1.5 * V_of_modulus_phi_squared;
    FOR2(i, j)
    {
        out.Sij[i][j] +=
            -0.5 * vars.h[i][j] * V_of_modulus_phi_squared / vars.chi;
    }

    return out;
}

// Calculate the stress energy tensor elements
template <class potential_t>
template <class data_t, template <typename> class vars_t>
void ComplexScalarField<potential_t>::emtensor_excl_potential(
    emtensor_t<data_t> &out, const vars_t<data_t> &vars,
    const CSFObject<data_t> &vars_csf, const Tensor<1, data_t> &d1_phi_Re,
    const Tensor<1, data_t> &d1_phi_Im, const Tensor<2, data_t> &h_UU,
    const Tensor<3, data_t> &chris_ULL)
{
    // initialise some useful quantities
    data_t modulus_d1_phi_squared = 0.;
    data_t modulus_Pi_squared =
        vars_csf.Pi_Re * vars_csf.Pi_Re + vars_csf.Pi_Im * vars_csf.Pi_Im;
    FOR2(i, j)
    {
        modulus_d1_phi_squared +=
            vars.chi * h_UU[i][j] *
            (d1_phi_Re[i] * d1_phi_Re[j] + d1_phi_Im[i] * d1_phi_Im[j]);
    }

    // Calculate components of EM Tensor
    // S_ij = T_ij
    FOR2(i, j)
    {
        out.Sij[i][j] =
            d1_phi_Re[i] * d1_phi_Re[j] + d1_phi_Im[i] * d1_phi_Im[j] -
            0.5 * (vars.h[i][j] *
                   (modulus_d1_phi_squared - modulus_Pi_squared) / vars.chi);
    }

    // S = Tr_S_ij
    out.S = vars.chi * TensorAlgebra::compute_trace(out.Sij, h_UU);

    // S_i (note lower index) = n^a T_a0
    FOR1(i)
    {
        out.Si[i] =
            vars_csf.Pi_Re * d1_phi_Re[i] + vars_csf.Pi_Im * d1_phi_Im[i];
    }

    // rho = n^a n^b T_ab
    out.rho = 0.5 * (modulus_Pi_squared + modulus_d1_phi_squared);
}

// Adds in the RHS for the matter vars
template <class potential_t>
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t,
          template <typename> class rhs_vars_t>
void ComplexScalarField<potential_t>::add_matter_rhs(
    rhs_vars_t<data_t> &total_rhs, const vars_t<data_t> &vars,
    const vars_t<Tensor<1, data_t>> &d1,
    const diff2_vars_t<Tensor<2, data_t>> &d2,
    const vars_t<data_t> &advec) const
{
    CSFObject<data_t> rhs_csf;
    // advection terms
    CSFObject<data_t> advec_csf;
    advec_csf.phi_Re = advec.phi_Re;
    advec_csf.phi_Im = advec.phi_Im;
    advec_csf.Pi_Re = advec.Pi_Re;
    advec_csf.Pi_Im = advec.Pi_Im;

    // the vars
    CSFObject<data_t> vars_csf;
    vars_csf.phi_Re = vars.phi_Re;
    vars_csf.phi_Im = vars.phi_Im;
    vars_csf.Pi_Re = vars.Pi_Re;
    vars_csf.Pi_Im = vars.Pi_Im;

    // call the function for the rhs excluding the potential
    matter_rhs_excl_potential(rhs_csf, vars, vars_csf, d1, d1.phi_Re, d1.phi_Im,
                              d2.phi_Re, d2.phi_Im, advec_csf);

    // set the default potential values
    data_t V_of_modulus_phi_squared = 0.0;
    data_t dVdmodulus_phi_squared = 0.0;

    // compute potential and add constributions to EM Tensor
    my_potential.compute_potential(V_of_modulus_phi_squared,
                                   dVdmodulus_phi_squared, vars);

    // adjust RHS for the potential term
    total_rhs.phi_Re = rhs_csf.phi_Re;
    total_rhs.phi_Im = rhs_csf.phi_Im;
    total_rhs.Pi_Re =
        rhs_csf.Pi_Re + vars.lapse * dVdmodulus_phi_squared * vars_csf.phi_Re;
    total_rhs.Pi_Im =
        rhs_csf.Pi_Im + vars.lapse * dVdmodulus_phi_squared * vars_csf.phi_Im;
}

// the RHS excluding the potential terms
template <class potential_t>
template <class data_t, template <typename> class vars_t>
void ComplexScalarField<potential_t>::matter_rhs_excl_potential(
    CSFObject<data_t> &rhs_csf, const vars_t<data_t> &vars,
    const CSFObject<data_t> &vars_csf, const vars_t<Tensor<1, data_t>> &d1,
    const Tensor<1, data_t> &d1_phi_Re, const Tensor<1, data_t> &d1_phi_Im,
    const Tensor<2, data_t> &d2_phi_Re, const Tensor<2, data_t> &d2_phi_Im,
    const CSFObject<data_t> &advec_csf)
{
    using namespace TensorAlgebra;

    const auto h_UU = compute_inverse_sym(vars.h);
    const auto chris = compute_christoffel(d1.h, h_UU);

    // evolution equations for scalar field and (minus) its conjugate momentum
    rhs_csf.phi_Re = advec_csf.phi_Re - vars.lapse * vars_csf.Pi_Re;
    rhs_csf.phi_Im = advec_csf.phi_Im - vars.lapse * vars_csf.Pi_Im;
    rhs_csf.Pi_Re = advec_csf.Pi_Re + vars.lapse * vars_csf.Pi_Re * vars.K;
    rhs_csf.Pi_Im = advec_csf.Pi_Im + vars.lapse * vars_csf.Pi_Im * vars.K;

    FOR1(k)
    {
        rhs_csf.Pi_Re +=
            vars.lapse * vars.chi * chris.contracted[k] * d1_phi_Re[k];
        rhs_csf.Pi_Im +=
            vars.lapse * vars.chi * chris.contracted[k] * d1_phi_Im[k];
    }

    FOR2(k, l)
    {
        rhs_csf.Pi_Re +=
            h_UU[k][l] * (-vars.chi * d1.lapse[k] * d1_phi_Re[l] +
                          vars.lapse * (0.5 * d1.chi[k] * d1_phi_Re[l] -
                                        vars.chi * d2_phi_Re[k][l]));
        rhs_csf.Pi_Im +=
            h_UU[k][l] * (-vars.chi * d1.lapse[k] * d1_phi_Im[l] +
                          vars.lapse * (0.5 * d1.chi[k] * d1_phi_Im[l] -
                                        vars.chi * d2_phi_Im[k][l]));
    }
}

#endif /* COMPLEXSCALARFIELD_IMPL_HPP_ */

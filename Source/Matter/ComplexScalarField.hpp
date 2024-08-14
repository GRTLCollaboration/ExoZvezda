/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef COMPLEXSCALARFIELD_HPP_
#define COMPLEXSCALARFIELD_HPP_

#include "CCZ4Geometry.hpp"
#include "DefaultComplexPotential.hpp"
#include "FourthOrderDerivatives.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS, total num of components
#include "VarsTools.hpp"

//!  Calculates the matter type specific elements such as the EMTensor and
//   matter evolution

template <class potential_t = DefaultComplexPotential> class ComplexScalarField
{
  protected:
    //! The local copy of the potential
    potential_t my_potential;

  public:
    //!  Constructor of class ScalarField, inputs are the matter parameters.
    ComplexScalarField(const potential_t a_potential)
        : my_potential(a_potential)
    {
    }

    //! Structure containing the variables for the matter fields
    template <class data_t> struct CSFObject
    {
        data_t phi_Re;
        data_t phi_Im;
        data_t Pi_Re;
        data_t Pi_Im;
    };

    //! Structure containing the rhs variables for the matter fields
    template <class data_t> struct Vars
    {
        data_t phi_Re;
        data_t phi_Im;
        data_t Pi_Re;
        data_t Pi_Im;

        /// Defines the mapping between members of Vars and Chombo grid
        /// variables (enum in User_Variables)
        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
            VarsTools::define_enum_mapping(mapping_function, c_phi_Re, phi_Re);
            VarsTools::define_enum_mapping(mapping_function, c_phi_Im, phi_Im);
            VarsTools::define_enum_mapping(mapping_function, c_Pi_Re, Pi_Re);
            VarsTools::define_enum_mapping(mapping_function, c_Pi_Im, Pi_Im);
        }
    };

    //! Structure containing the rhs variables for the matter fields requiring
    //!  2nd derivs
    template <class data_t> struct Diff2Vars
    {
        data_t phi_Re;
        data_t phi_Im;

        /// Defines the mapping between members of Vars and Chombo grid
        ///  variables (enum in User_Variables)
        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
            VarsTools::define_enum_mapping(mapping_function, c_phi_Re, phi_Re);
            VarsTools::define_enum_mapping(mapping_function, c_phi_Im, phi_Im);
        }
    };

    //! The function which calculates the EM Tensor, given the vars and
    //! derivatives, including the potential
    template <class data_t, template <typename> class vars_t>
    emtensor_t<data_t> compute_emtensor(
        const vars_t<data_t> &vars,          //!< the value of the variables
        const vars_t<Tensor<1, data_t>> &d1, //!< the value of the 1st derivs
        const Tensor<2, data_t> &h_UU, //!< the inverse metric (raised indices)
        const Tensor<3, data_t> &chris_ULL)
        const; //!< the conformal christoffel symbol

    //! The function which calculates the EM Tensor, given the vars and
    //! derivatives, excluding the potential
    template <class data_t, template <typename> class vars_t>
    static void emtensor_excl_potential(
        emtensor_t<data_t> &out,           //!< the em tensor output
        const vars_t<data_t> &vars,        //!< the value of the variables
        const CSFObject<data_t> &vars_csf, //!< the value of the csf variables
        const Tensor<1, data_t>
            &d1_phi_Re, //!< the value of the first deriv of phi_Re
        const Tensor<1, data_t>
            &d1_phi_Im, //!< the value of the first deriv of phi_Im
        const Tensor<2, data_t> &h_UU, //!< the inverse metric (raised indices).
        const Tensor<3, data_t>
            &chris_ULL); //!< the conformal christoffel symbol

    //! The function which adds in the RHS for the matter field vars,
    //! including the potential
    template <class data_t, template <typename> class vars_t,
              template <typename> class diff2_vars_t,
              template <typename> class rhs_vars_t>
    void add_matter_rhs(
        rhs_vars_t<data_t> &total_rhs,       //!< value of the RHS for all vars
        const vars_t<data_t> &vars,          //!< value of the variables
        const vars_t<Tensor<1, data_t>> &d1, //!< value of the 1st derivs
        const diff2_vars_t<Tensor<2, data_t>> &d2, //!< value of the 2nd derivs
        const vars_t<data_t> &advec) //!< the value of the advection terms
        const;

    //! The function which calculates the RHS for the matter field vars
    //! excluding the potential
    template <class data_t, template <typename> class vars_t>
    static void matter_rhs_excl_potential(
        CSFObject<data_t>
            &rhs_csf, //!< the value of the RHS terms for the sf vars
        const vars_t<data_t> &vars,        //!< the values of all the variables
        const CSFObject<data_t> &vars_csf, //!< the value of the sf variables
        const vars_t<Tensor<1, data_t>> &d1, //!< the value of the 1st derivs
        const Tensor<1, data_t>
            &d1_phi_Re, //!< the value of the 1st derivs of phi_Re
        const Tensor<1, data_t>
            &d1_phi_Im, //!< the value of the 1st derivs of phi_Im
        const Tensor<2, data_t>
            &d2_phi_Re, //!< the value of the 2nd derivs of phi_Re
        const Tensor<2, data_t>
            &d2_phi_Im, //!< the value of the 2nd derivs of phi_Im
        const CSFObject<data_t>
            &advec_csf); //!< advection terms for the csf vars
};

#include "ComplexScalarField.impl.hpp"

#endif /* COMPLEXSCALARFIELD_HPP_ */

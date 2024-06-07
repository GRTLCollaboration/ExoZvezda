/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

 #ifndef DIAGNOSTICVARIABLES_HPP
 #define DIAGNOSTICVARIABLES_HPP

// assign an enum to each variable
enum
{
    c_mod_phi, // the complex scalar field modulus

    c_Madm,
    c_Jadm,

    c_N, // Noether Charge integrand

    c_Ham,

    c_Mom1,
    c_Mom2,
    c_Mom3,

    c_Weyl4_Re,
    c_Weyl4_Im,

    c_rho,   // stress tensor components
    c_s1,
    c_s2,
    c_s3,
    c_s11,
    c_s12,
    c_s13,
    c_s22,
    c_s23,
    c_s33,

    c_Qphi_density,

    c_Fphi_flux,

    c_Sphi_source,

    c_weight1,
    c_weight2,

    c_profile1,
    c_profile2,

    c_testHam,

    NUM_DIAGNOSTIC_VARS
};

namespace DiagnosticVariables
{
static const std::array<std::string, NUM_DIAGNOSTIC_VARS> variable_names = {

    "mod_phi",

    "Madm",   "Jadm",

    "N",

    "Ham",    "Mom1",   "Mom2",   "Mom3",

    "Weyl4_Re",  "Weyl4_Im",

    "rho", "s1", "s2", "s3", "s11", "s12", "s13", "s22", "s23", "s33",

    "Qphi_density", "Fphi_flux", "Sphi_source",
    
    "weight1", "weight2",
        
    "profile1", "profile2",
    
    "testHam"
    
     };
    
    }

#endif /* DIAGNOSTICVARIABLES_HPP */

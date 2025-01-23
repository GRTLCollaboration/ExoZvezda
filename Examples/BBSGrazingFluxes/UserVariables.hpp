/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef USERVARIABLES_HPP
#define USERVARIABLES_HPP

#include "ArrayTools.hpp"
#include "CCZ4UserVariables.hpp"
#include "DiagnosticVariables.hpp"
#include <array>
#include <string>

/// This enum gives the index of every variable stored in the grid
enum
{

    c_phi_Re = NUM_CCZ4_VARS, // real part of scalar field
    c_phi_Im,                 // imaginary part of scalar field
    c_Pi_Re,                  // real part of auxiliary variable Pi = -L_n phi
    c_Pi_Im, // imaginary part of auxiliary variable Pi = -L_n phi
    // Note that it is important that the first enum value is set to 1 more than
    // the last CCZ4 var enum
    NUM_VARS
};

namespace UserVariables
{
static const std::array<std::string, NUM_VARS - NUM_CCZ4_VARS>
    scalarfield_variable_names = {"phi_Re", "phi_Im", "Pi_Re", "Pi_Im"};

static const std::array<std::string, NUM_VARS> variable_names =
    ArrayTools::concatenate(ccz4_variable_names, scalarfield_variable_names);
} // namespace UserVariables

#include "UserVariables.inc.hpp"

#endif /* USERVARIABLES_HPP */

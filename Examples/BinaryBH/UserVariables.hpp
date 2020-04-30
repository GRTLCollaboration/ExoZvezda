/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef USERVARIABLES_HPP
#define USERVARIABLES_HPP

#include "DiagnosticVariables.hpp"

/// This enum gives the index of every variable stored in the grid
enum
{
    // Note that it is important that the first enum value is set to 1 more than
    // the last CCZ4 var enum
    c_Ham = NUM_CCZ4_VARS,

    c_Mom1,
    c_Mom2,
    c_Mom3,

    c_Weyl4_Re,
    c_Weyl4_Im,

    NUM_VARS
};

namespace UserVariables
{
static const std::array<std::string, NUM_VARS - NUM_CCZ4_VARS>
    user_variable_names = {"Ham",

                           "Mom1",     "Mom2",    "Mom3",

                           "Weyl4_Re", "Weyl4_Im"};

static const std::array<std::string, NUM_VARS> variable_names =
    ArrayTools::concatenate(ccz4_variable_names, user_variable_names);
} // namespace UserVariables

#include "UserVariables.inc.hpp"

#endif /* USERVARIABLES_HPP */

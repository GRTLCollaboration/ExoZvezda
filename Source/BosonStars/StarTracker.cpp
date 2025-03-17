/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "StarTracker.hpp"
#include "DimensionDefinitions.hpp"
#include "FittingMethod.hpp"
#include "FittingModel.hpp"
#include "InterpolationQuery.hpp"
#include "SmallDataIO.hpp" // for writing data

void StarTracker::set_up_fitting(int num_star, int fitting_direction)
{
    int success;
    double delta;

    std::vector<double> a_vector(3);        // vector for fitted coefficients
    std::vector<double> vals_chi(m_points); // vector to store chi

    // Define the x/y/z coordinate intervals for fitting. The interval length is
    // determined by 'm_width_A' and 'm_width_B', which can be different for the
    // two BSs. The 'm_points' determines the way the interval is split into
    // parts.
    for (int i = 0; i < m_points; i++)
    {
        if (num_star == 0)
        {
            delta = m_width_A * (double(i) / double(m_points - 1) - 0.5);
        }

        if (num_star == 1)
        {
            delta = m_width_B * (double(i) / double(m_points - 1) - 0.5);
        }

        if (fitting_direction == 0)
        {
            m_x_coords[i] = m_star_coords[num_star][0] + delta;
            m_y_coords[i] = m_star_coords[num_star][1];
            m_z_coords[i] = m_star_coords[num_star][2];
        }

        if (fitting_direction == 1)
        {
            m_x_coords[i] = m_star_coords[num_star][0];
            m_y_coords[i] = m_star_coords[num_star][1] + delta;
            m_z_coords[i] = m_star_coords[num_star][2];
        }

        if (fitting_direction == 2)
        {
            m_x_coords[i] = m_star_coords[num_star][0];
            m_y_coords[i] = m_star_coords[num_star][1];
            m_z_coords[i] = m_star_coords[num_star][2] + delta;
        }

        m_sigma_vector[i] = 1.0; // agnostic about the error vector, so set =
                                 // to 1.
    }

    // Set up interpolator to get chi values along our x/y/z intervals
    bool fill_ghosts = true;
    m_interpolator->refresh(fill_ghosts);

    InterpolationQuery query_st(m_points);
    query_st.setCoords(0, m_x_coords.data())
        .setCoords(1, m_y_coords.data())
        .setCoords(2, m_z_coords.data())
        .addComp(c_chi, vals_chi.data());

    m_interpolator->interp(query_st);

    // We want to fit the Gaussian to (1-chi), since chi asymptotes to 1
    for (int i = 0; i < m_points; i++)
    {
        m_vals_shifted_chi[i] = 1 - vals_chi[i];
    }
}

// Function to find the centres of the BSs using Gaussian Fitting performed via
// FittingMethod routine. Here 'num_star' is the number BSs being simulated, so
// for a binary we have num_star = [0,1], and 'fitting_direction' specifies the
// spatial axis along which we perform the fitting = [0,1,2] correspinding to
// the set of [x,y,z].
double StarTracker::find_centre(int num_star, int fitting_direction)
{
    int success;
    std::vector<double> a_vector(3); // vector with fitted coefficients

    set_up_fitting(num_star, fitting_direction);

    // x-direction
    if (fitting_direction == 0)
    {
        a_vector[0] = m_vals_shifted_chi[(m_points - 1) / 2];
        a_vector[1] = m_star_coords[num_star][0];
        if (num_star == 0)
        {
            a_vector[2] = m_width_A / 2;
        }
        if (num_star == 1)
        {
            a_vector[2] = m_width_B / 2;
        }

        FittingMethod fitting_method(m_x_coords, m_vals_shifted_chi,
                                     m_sigma_vector, a_vector, gaussian_model);

        success = fitting_method.fit();
        if (success == 1)
        {
            return fitting_method.fit_parameters[1];
        }
        else
        {
            MayDay::Error(
                "Oops, help, I cannot fit a Gaussian along x-axis in "
                "StarTracker.cpp! Please check you tracking parameters to "
                "ensure I am fitting something plausible...");
            return 0;
        }
    }

    // y-direction
    if (fitting_direction == 1)
    {
        a_vector[0] = m_vals_shifted_chi[(m_points - 1) / 2];
        a_vector[1] = m_star_coords[num_star][1];
        if (num_star == 0)
        {
            a_vector[2] = m_width_A / 2;
        }
        if (num_star == 1)
        {
            a_vector[2] = m_width_B / 2;
        }

        FittingMethod fitting_method(m_y_coords, m_vals_shifted_chi,
                                     m_sigma_vector, a_vector, gaussian_model);

        success = fitting_method.fit();
        if (success == 1)
        {
            return fitting_method.fit_parameters[1];
        }
        else
        {
            MayDay::Error(
                "Oops, help, I cannot fit a Gaussian along y-axis in "
                "StarTracker.cpp! Please check you tracking parameters to "
                "ensure I am fitting something plausible...");
            return 0;
        }
    }

    // z-direction
    if (fitting_direction == 2)
    {
        a_vector[0] = m_vals_shifted_chi[(m_points - 1) / 2];
        a_vector[1] = m_star_coords[num_star][2];
        if (num_star == 0)
        {
            a_vector[2] = m_width_A / 2;
        }
        if (num_star == 1)
        {
            a_vector[2] = m_width_B / 2;
        }

        FittingMethod fitting_method(m_z_coords, m_vals_shifted_chi,
                                     m_sigma_vector, a_vector, gaussian_model);
        success = fitting_method.fit();
        if (success == 1)
        {
            MayDay::Error(
                "Oops, help, I cannot fit a Gaussian along z-axis in "
                "StarTracker.cpp! Please check you tracking parameters to "
                "ensure I am fitting something plausible...");
            return fitting_method.fit_parameters[1];
        }
        else
        {
            return 0;
        }
    }
    MayDay::Error("I have failed to apply Gaussian fitting!");

    return 0;
}

// Function to find the centres of the BSs near merger.
// We do not update the stars' positions when the coordinate speeds become
// greater than 1 in a given direction (helps to avoid huge jumps around merger,
// where fitting will be harder). If this behaviour happens, we simply switch to
// a "center of mass" version of stars' positions.
void StarTracker::find_centre_merger(int num_star, int fitting_direction)
{
    std::vector<double> a_vector(3);

    set_up_fitting(num_star, fitting_direction);

    if (m_vals_shifted_chi.size() < m_points || m_x_coords.size() < m_points ||
        m_y_coords.size() < m_points || m_z_coords.size() < m_points)
    {
        MayDay::Error("Insufficient data in star tracking!");
    }

    double fmax =
        *max_element(m_vals_shifted_chi.begin(), m_vals_shifted_chi.end());
    double fmin =
        *min_element(m_vals_shifted_chi.begin(), m_vals_shifted_chi.end());

    if (fmax == fmin)
    {
        MayDay::Error("fmax and fmin are equal, and I am dividing by zero!");
    }

    double weight = 0.0;
    double sum1 = 0.0;
    double sum2 = 0.0;

    if (fitting_direction == 0)
    {
        for (int i = 0; i < m_points; i++)
        {
            weight = (m_vals_shifted_chi[i] - fmin) / (fmax - fmin);
            sum1 = sum1 + m_x_coords[i] * weight;
            sum2 = sum2 + weight;
        }

        m_star_coords[num_star][0] = sum1 / sum2;

        if (sum2 == 0.0)
        {
            MayDay::Error("Division by zero detected in find_centre_merger");
        }
    }

    if (fitting_direction == 1)
    {
        for (int i = 0; i < m_points; i++)
        {
            weight = (m_vals_shifted_chi[i] - fmin) / (fmax - fmin);
            sum1 = m_y_coords[i] * weight;
            sum2 = weight;
        }

        m_star_coords[num_star][1] = sum1 / sum2;

        if (sum2 == 0.0)
        {
            MayDay::Error("Division by zero detected in find_centre_merger");
        }
    }

    if (fitting_direction == 2)
    {
        for (int i = 0; i < m_points; i++)
        {
            weight = (m_vals_shifted_chi[i] - fmin) / (fmax - fmin);
            sum1 = m_z_coords[i] * weight;
            sum2 = weight;
        }

        m_star_coords[num_star][2] = sum1 / sum2;

        if (sum2 == 0.0)
        {
            MayDay::Error("Division by zero detected in find_centre_merger");
        }
    }
}

// Finally update the centres either using Gaussian fitting procedure or centre
// of mass calculation.
void StarTracker::update_star_centres(double a_dt)
{
    if (m_fitting_direction == "x")
    {
        double starA_0 = find_centre(0, 0);
        if (abs((starA_0 - m_star_coords[0][0]) / a_dt) < 1.0 && starA_0 != 0)
        {
            m_star_coords[0][0] = starA_0;
        }
        else
        {
            find_centre_merger(0, 0);
        }
        double starB_0 = find_centre(1, 0);
        if ((abs(starB_0 - m_star_coords[1][0]) / a_dt) < 1.0 && starB_0 != 0)
        {
            m_star_coords[1][0] = starB_0;
        }
        else
        {
            find_centre_merger(1, 0);
        }
    }

    if (m_fitting_direction == "xy")
    {
        double starA_0 = find_centre(0, 0);
        if (abs((starA_0 - m_star_coords[0][0]) / a_dt) < 1.0 && starA_0 != 0)
        {
            m_star_coords[0][0] = starA_0;
        }
        else
        {
            find_centre_merger(0, 0);
        }
        double starA_1 = find_centre(0, 1);
        if (abs((starA_1 - m_star_coords[0][1]) / a_dt) < 1.0 && starA_1 != 0)
        {
            m_star_coords[0][1] = starA_1;
        }
        else
        {
            find_centre_merger(0, 1);
        }
        double starB_0 = find_centre(1, 0);
        if (abs((starB_0 - m_star_coords[1][0]) / a_dt) < 1.0 && starB_0 != 0)
        {
            m_star_coords[1][0] = starB_0;
        }
        else
        {
            find_centre_merger(1, 0);
        }
        double starB_1 = find_centre(1, 1);
        if (abs((starB_1 - m_star_coords[1][1]) / a_dt) < 1.0 && starB_1 != 0)
        {
            m_star_coords[1][1] = starB_1;
        }
        else
        {
            find_centre_merger(1, 1);
        }
    }

    if (m_fitting_direction == "xyz")
    {
        double starA_0 = find_centre(0, 0);
        m_star_coords[0][0] = starA_0;
        double starA_1 = find_centre(0, 1);
        m_star_coords[0][1] = starA_1;
        double starA_2 = find_centre(0, 2);
        m_star_coords[0][2] = starA_2;

        double starB_0 = find_centre(1, 0);
        m_star_coords[1][0] = starB_0;
        double starB_1 = find_centre(1, 1);
        m_star_coords[1][1] = starB_1;
        double starB_2 = find_centre(1, 2);
        m_star_coords[1][2] = starB_2;
    }
}

// Read a data line from the previous timestep
void StarTracker::read_in_star_coords(int a_int_step, double a_current_time)
{
    bool first_step = false;
    double dt = (a_current_time / a_int_step);
    SmallDataIO star_file("StarCentres", dt, a_current_time, a_current_time,
                          SmallDataIO::APPEND, first_step);

    // NB need to give the get function an empty vector to fill
    std::vector<double> star_vector;
    star_file.get_specific_data_line(star_vector, a_current_time);

    // check the data returned is the right size
    CH_assert(star_vector.size() % CH_SPACEDIM == 0);

    m_num_stars = star_vector.size() / CH_SPACEDIM;
    m_star_coords.resize(m_num_stars);

    // remove any duplicate data from the file
    const bool keep_m_time_data = true;
    star_file.remove_duplicate_time_data(keep_m_time_data);

    // convert vector to list of coords
    for (int ipuncture = 0; ipuncture < m_num_stars; ipuncture++)
    {
        m_star_coords[ipuncture] = {star_vector[ipuncture * CH_SPACEDIM + 0],
                                    star_vector[ipuncture * CH_SPACEDIM + 1],
                                    star_vector[ipuncture * CH_SPACEDIM + 2]};
    }

    update_star_centres(dt);

    for (int ipuncture = 0; ipuncture < m_num_stars; ipuncture++)
    {
        pout() << "Star " << ipuncture
               << " restarted at : " << m_star_coords[ipuncture][0] << " "
               << m_star_coords[ipuncture][1] << " "
               << m_star_coords[ipuncture][2] << endl;
        pout() << "at time = " << a_current_time << endl;
    }
}

// What to do post restart/how to read star coordinates in
void StarTracker::restart_star_tracking()
{
    int current_step = m_interpolator->getAMR().s_step;

    if (current_step == 0)
    {
        // if it is the first timestep, use the param values
        // rather than look for the output file, e.g. for when
        // restart from IC solver checkpoint
        set_initial_star_coords();
    }
    else
    {
        // look for the current star location in the
        // star output file (it needs to exist!)
        read_in_star_coords(current_step,
                            m_interpolator->getAMR().getCurrentTime());
    }
}

// Set and write initial star locations
void StarTracker::set_initial_star_coords()
{
    CH_assert(m_star_coords.size() > 0); // sanity check

    // now the write out to a new file
    bool first_step = true;
    double dt = 1.; // doesn't matter
    double time = 0.;
    double restart_time = 0.;
    SmallDataIO star_centre_file("StarCentres", dt, time, restart_time,
                                 SmallDataIO::APPEND, first_step);
    std::vector<std::string> header1_strings(CH_SPACEDIM * m_num_stars);
    for (int ipuncture = 0; ipuncture < m_num_stars; ipuncture++)
    {
        std::string idx = std::to_string(ipuncture + 1);
        header1_strings[CH_SPACEDIM * ipuncture + 0] = "x_" + idx;
        header1_strings[CH_SPACEDIM * ipuncture + 1] = "y_" + idx;
        header1_strings[CH_SPACEDIM * ipuncture + 2] = "z_" + idx;
    }
    star_centre_file.write_header_line(header1_strings);

    star_centre_file.write_time_data_line(get_star_vector());
}

// Execute the tracking and write out
void StarTracker::execute_tracking(double a_time, double a_restart_time,
                                   double a_dt, const bool write_data)
{
    CH_assert(m_interpolator != nullptr); // sanity check

    update_star_centres(a_dt);

    // print them out
    if (write_data)
    {
        bool first_step = false;
        SmallDataIO star_file("StarCentres", a_dt, a_time, a_restart_time,
                              SmallDataIO::APPEND, first_step);

        // use a vector for the write out
        star_file.write_time_data_line(get_star_vector());
    }
}

// Get a vector of the star coords - used for write out
std::vector<double> StarTracker::get_star_vector() const
{
    std::vector<double> star_vector;
    star_vector.resize(m_num_stars * CH_SPACEDIM);
    for (int ipuncture = 0; ipuncture < m_num_stars; ipuncture++)
    {
        star_vector[ipuncture * CH_SPACEDIM + 0] = m_star_coords[ipuncture][0];
        star_vector[ipuncture * CH_SPACEDIM + 1] = m_star_coords[ipuncture][1];
        star_vector[ipuncture * CH_SPACEDIM + 2] = m_star_coords[ipuncture][2];
    }
    return star_vector;
}

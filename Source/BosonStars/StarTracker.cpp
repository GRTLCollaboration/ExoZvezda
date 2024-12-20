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
            m_x_coords[i] = m_puncture_coords[num_star][0] + delta;
            m_y_coords[i] = m_puncture_coords[num_star][1];
            m_z_coords[i] = m_puncture_coords[num_star][2];
        }

        if (fitting_direction == 1)
        {
            m_x_coords[i] = m_puncture_coords[num_star][0];
            m_y_coords[i] = m_puncture_coords[num_star][1] + delta;
            m_z_coords[i] = m_puncture_coords[num_star][2];
        }

        if (fitting_direction == 2)
        {
            m_x_coords[i] = m_puncture_coords[num_star][0];
            m_y_coords[i] = m_puncture_coords[num_star][1];
            m_z_coords[i] = m_puncture_coords[num_star][2] + delta;
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
        a_vector[1] = m_puncture_coords[num_star][0];
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
        a_vector[1] = m_puncture_coords[num_star][1];
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
        a_vector[1] = m_puncture_coords[num_star][2];
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
    return 0;
}

// Function to find the centres of the BSs near merger.
// We do not update the stars' positions when the coordinate speeds become
// greater than 1 in a given direction (helps to avoid huge jumps around merger,
// where fitting will be harder). If this behaviour happens, we simply switch to
// a "center of mass" version of stars' positions.
void StarTracker::find_max_min(int num_star, int fitting_direction)
{
    std::vector<double> a_vector(3);

    set_up_fitting(num_star, fitting_direction);

    double fmax =
        *max_element(m_vals_shifted_chi.begin(), m_vals_shifted_chi.end());
    double fmin =
        *min_element(m_vals_shifted_chi.begin(), m_vals_shifted_chi.end());

    double weight;
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

        m_puncture_coords[num_star][0] = sum1 / sum2;
    }

    if (fitting_direction == 1)
    {
        for (int i = 0; i < m_points; i++)
        {
            weight = (m_vals_shifted_chi[i] - fmin) / (fmax - fmin);
            sum1 += m_y_coords[i] * weight;
            sum2 += weight;
        }

        m_puncture_coords[num_star][1] = sum1 / sum2;
    }

    if (fitting_direction == 2)
    {
        for (int i = 0; i < m_points; i++)
        {
            weight = (m_vals_shifted_chi[i] - fmin) / (fmax - fmin);
            sum1 += m_z_coords[i] * weight;
            sum2 += weight;
        }

        m_puncture_coords[num_star][2] = sum1 / sum2;
    }
}

// Finally update the centres either using Gaussian fitting procedure or centre
// of mass calculation.
void StarTracker::update_star_centres(double a_dt)
{
    if (m_fitting_direction == "x")
    {
        double starA_0 = find_centre(0, 0);
        if (abs((starA_0 - m_puncture_coords[0][0]) / a_dt) < 1.0 &&
            starA_0 != 0)
        {
            m_puncture_coords[0][0] = starA_0;
        }
        else
        {
            find_max_min(0, 0);
        }
        double starB_0 = find_centre(1, 0);
        if ((abs(starB_0 - m_puncture_coords[1][0]) / a_dt) < 1.0 &&
            starB_0 != 0)
        {
            m_puncture_coords[1][0] = starB_0;
        }
        else
        {
            find_max_min(1, 0);
        }

        pout() << m_puncture_coords[0][0] << endl;
        pout() << m_puncture_coords[1][0] << endl;
    }

    if (m_fitting_direction == "xy")
    {
        double starA_0 = find_centre(0, 0);
        if (abs((starA_0 - m_puncture_coords[0][0]) / a_dt) < 1.0 &&
            starA_0 != 0)
        {
            m_puncture_coords[0][0] = starA_0;
        }
        else
        {
            find_max_min(0, 0);
        }
        double starA_1 = find_centre(0, 1);
        if (abs((starA_1 - m_puncture_coords[0][1]) / a_dt) < 1.0 &&
            starA_1 != 0)
        {
            m_puncture_coords[0][1] = starA_1;
        }
        else
        {
            find_max_min(0, 1);
        }
        double starB_0 = find_centre(1, 0);
        if (abs((starB_0 - m_puncture_coords[1][0]) / a_dt) < 1.0 &&
            starB_0 != 0)
        {
            m_puncture_coords[1][0] = starB_0;
        }
        else
        {
            find_max_min(1, 0);
        }
        double starB_1 = find_centre(1, 1);
        if (abs((starB_1 - m_puncture_coords[1][1]) / a_dt) < 1.0 &&
            starB_1 != 0)
        {
            m_puncture_coords[1][1] = starB_1;
        }
        else
        {
            find_max_min(1, 1);
        }
    }

    if (m_fitting_direction == "xyz")
    {
        double starA_0 = find_centre(0, 0);
        m_puncture_coords[0][0] = starA_0;
        double starA_1 = find_centre(0, 1);
        m_puncture_coords[0][1] = starA_1;
        double starA_2 = find_centre(0, 2);
        m_puncture_coords[0][2] = starA_2;

        double starB_0 = find_centre(1, 0);
        m_puncture_coords[1][0] = starB_0;
        double starB_1 = find_centre(1, 1);
        m_puncture_coords[1][1] = starB_1;
        double starB_2 = find_centre(1, 2);
        m_puncture_coords[1][2] = starB_2;
    }
}

// Write all data to designated files
void StarTracker::write_to_dat(std::string a_filename, double a_dt,
                               double a_time, double a_restart_time,
                               bool a_first_step)
{
    int size;
    double eps = 10e-8;
    std::vector<double> star_coords;

    SmallDataIO star_centre_file(a_filename, a_dt, a_time, a_restart_time,
                                 SmallDataIO::APPEND, a_first_step);

    if (a_time > a_restart_time + eps)
        star_centre_file.remove_duplicate_time_data();

    std::vector<string> header_line(3. * m_num_stars);

    for (int n = 0; n < m_num_stars; n++)
    {
        header_line[3 * n] = "x" + to_string(n + 1);
        header_line[3 * n + 1] = "y" + to_string(n + 1);
        header_line[3 * n + 2] = "z" + to_string(n + 1);
    }

    if (a_time == 0.)
    {
        star_centre_file.write_header_line(header_line);
    }

    size = CH_SPACEDIM * m_num_stars;
    star_coords.resize(size, 0);

    for (int ipuncture = 0; ipuncture < m_num_stars; ++ipuncture)
    {
        for (int i = 0; i < CH_SPACEDIM; ++i)
        {
            star_coords[CH_SPACEDIM * ipuncture + i] =
                m_puncture_coords[ipuncture][i];
        }
    }
    star_centre_file.write_time_data_line(star_coords);
}

// Read a data line from the previous timestep
void StarTracker::read_old_centre_from_dat(std::string a_filename, double a_dt,
                                           double a_time, double a_restart_time,
                                           bool a_first_step)
{
    int size;
    std::vector<double> data_line;
    std::vector<double> star_coords;

    if (a_time > a_dt / 3.)
    {
        SmallDataIO star_centre_file(a_filename, a_dt, a_time, a_restart_time,
                                     SmallDataIO::READ, a_first_step);
        star_centre_file.get_specific_data_line(data_line, a_time - a_dt);

        bool length_match = data_line.size() == m_num_stars * CH_SPACEDIM;

        size = CH_SPACEDIM * m_num_stars;
        star_coords.resize(size, 0);

        if (length_match)
        {
            for (int i = 0; i < data_line.size(); i++)
            {
                star_coords[i] = data_line[i];
            }

            pout() << "Star A restarted at : " << star_coords[0] << " " << star_coords[1] << " " << star_coords[2] << endl;
            pout() << "Star B restarted at : " << star_coords[3] << " " << star_coords[4] << " " << star_coords[5] << endl;
        }
        else
        {
            for (int i = 0; i < star_coords.size(); i++)
            {
                star_coords[i] = NAN;
            }
            MayDay::Error("Array size mismatch, when loading star positions "
                          "from StarCentres.dat file!");
        }
    }
}

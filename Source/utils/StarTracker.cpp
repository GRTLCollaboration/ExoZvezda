/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "StarTracker.hpp"
#include "ChomboParameters.hpp" // for writing data
#include "DebuggingTools.hpp"
#include "DimensionDefinitions.hpp"
#include "FitExample.hpp"
#include "FitMRQ.hpp"
#include "InterpolationQuery.hpp"
#include "SmallDataIO.hpp"   // for writing data
#include "UserVariables.hpp" // for writing data
#include "nr3.h"
#include <algorithm>

// Function to find the centres of the BSs using Guassian Fitting performed via
// FitMRQ routine 'num_star' is the number of the BS = [0,1] and 'direction' is
// the spatial axis along which we perform the fitting = [0,1,2] correspinding
// to the set of [x,y,z]
double StarTracker::find_centre(int num_star, int direction)
{
    int success;
    double delta;

    std::vector<double> x_coords(m_points);
    std::vector<double> y_coords(m_points);
    std::vector<double> z_coords(m_points);

    std::vector<double> sigma_vector(m_points); // error associated with x/y/z
    std::vector<double> a_vector(3);      // vector with fitted coefficients
    std::vector<double> vals(m_points);   // vector to store chi
    std::vector<double> vals_f(m_points); // vector to store (1-chi)

    // Define the x/y/z coordinate intervals for fitting. The interval length is
    // determined by 'm_width_A' and 'm_width_B', which can be dfferent for the
    // two BSs.
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

        if (direction == 0)
        {
            x_coords[i] = m_star_coords[3 * num_star] + delta;
            y_coords[i] = m_star_coords[3 * num_star + 1];
            z_coords[i] = m_star_coords[3 * num_star + 2];
        }

        if (direction == 1)
        {
            x_coords[i] = m_star_coords[3 * num_star];
            y_coords[i] = m_star_coords[3 * num_star + 1] + delta;
            z_coords[i] = m_star_coords[3 * num_star + 2];
        }

        if (direction == 2)
        {
            x_coords[i] = m_star_coords[3 * num_star];
            y_coords[i] = m_star_coords[3 * num_star + 1];
            z_coords[i] = m_star_coords[3 * num_star + 2] + delta;
        }

        sigma_vector[i] = 1.0; // agnostic about the error vector, so set =
                               // to 1.
    }

    // Set up interpolator to get cho values along our x/y/z intervals
    bool fill_ghosts = true;
    m_interpolator->refresh(fill_ghosts);

    InterpolationQuery query(m_points);
    query.setCoords(0, x_coords.data())
        .setCoords(1, y_coords.data())
        .setCoords(2, z_coords.data())
        .addComp(c_chi, vals.data());

    m_interpolator->interp(query);

    // We want to fit the Gaussian to (1-chi), since chi asymptotes to 1
    for (int i = 0; i < m_points; i++)
    {
        vals_f[i] = 1 - vals[i];
    }

    // Set the initial guesses for all directions and perform the fit
    // individually below:

    if (direction == 0)
    {
        a_vector[0] = vals_f[(m_points - 1) / 2];
        a_vector[1] = m_star_coords[3 * num_star];
        if (num_star == 0)
        {
            a_vector[2] = m_width_A / 2;
        }
        if (num_star == 1)
        {
            a_vector[2] = m_width_B / 2;
        }

        Fitmrq fitmrq1(x_coords, vals_f, sigma_vector, a_vector, fgauss);

        success = fitmrq1.fit();
        if (success == 1)
        {
            return fitmrq1.a[1];
        }
        else
        {
            return 0;
        }
    }

    if (direction == 1)
    {
        a_vector[0] = vals_f[(m_points - 1) / 2];
        a_vector[1] = m_star_coords[3 * num_star + 1];
        if (num_star == 0)
        {
            a_vector[2] = m_width_A / 2;
        }
        if (num_star == 1)
        {
            a_vector[2] = m_width_B / 2;
        }

        Fitmrq fitmrq1(y_coords, vals_f, sigma_vector, a_vector, fgauss);

        success = fitmrq1.fit();
        if (success == 1)
        {
            return fitmrq1.a[1];
        }
        else
        {
            return 0;
        }
    }

    if (direction == 2)
    {
        a_vector[0] = vals_f[(m_points - 1) / 2];
        a_vector[1] = m_star_coords[3 * num_star + 2];
        if (num_star == 0)
        {
            a_vector[2] = m_width_A / 2;
        }
        if (num_star == 1)
        {
            a_vector[2] = m_width_B / 2;
        }

        Fitmrq fitmrq1(z_coords, vals_f, sigma_vector, a_vector, fgauss);
        success = fitmrq1.fit();
        if (success == 1)
        {
            return fitmrq1.a[1];
        }
        else
        {
            return 0;
        }
    }

    return 0;
}

// Function to find the centres of the BSs near merger.
// We do not update the position if it implies a coordinate speed greater than 1
// in this direction. This is to avoid jumps around merger, where we temporarily
// may develop multi-peaked profiles. If this happens, we switch to a "center of
// mass" version of the position.
void StarTracker::find_max_min(int num_star, int direction)
{
    double delta;

    std::vector<double> x_coords(m_points);
    std::vector<double> y_coords(m_points);
    std::vector<double> z_coords(m_points);

    std::vector<double> sigma_vector(m_points);
    std::vector<double> a_vector(3);
    std::vector<double> vals(m_points);
    std::vector<double> vals_f(m_points);

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

        if (direction == 0)
        {
            x_coords[i] = m_star_coords[3 * num_star] + delta;
            y_coords[i] = m_star_coords[3 * num_star + 1];
            z_coords[i] = m_star_coords[3 * num_star + 2];
        }

        if (direction == 1)
        {
            x_coords[i] = m_star_coords[3 * num_star];
            y_coords[i] = m_star_coords[3 * num_star + 1] + delta;
            z_coords[i] = m_star_coords[3 * num_star + 2];
        }

        if (direction == 2)
        {
            x_coords[i] = m_star_coords[3 * num_star];
            y_coords[i] = m_star_coords[3 * num_star + 1];
            z_coords[i] = m_star_coords[3 * num_star + 2] + delta;
        }

        sigma_vector[i] = 1.0;
    }

    bool fill_ghosts = true;
    m_interpolator->refresh(fill_ghosts);

    InterpolationQuery query(m_points);
    query.setCoords(0, x_coords.data())
        .setCoords(1, y_coords.data())
        .setCoords(2, z_coords.data())
        .addComp(c_chi, vals.data());

    m_interpolator->interp(query);

    for (int i = 0; i < m_points; i++)
    {
        vals_f[i] = 1 - vals[i];
    }

    double fmax = *max_element(vals_f.begin(), vals_f.end());
    double fmin = *min_element(vals_f.begin(), vals_f.end());

    double weight;
    double sum1 = 0.0;
    double sum2 = 0.0;

    if (direction == 0)
    {
        for (int i = 0; i < m_points; i++)
        {
            weight = (vals_f[i] - fmin) / (fmax - fmin);
            sum1 = sum1 + x_coords[i] * weight;
            sum2 = sum2 + weight;
        }

        m_star_coords[3 * num_star] = sum1 / sum2;
    }

    if (direction == 1)
    {
        for (int i = 0; i < m_points; i++)
        {
            weight = (vals_f[i] - fmin) / (fmax - fmin);
            sum1 += y_coords[i] * weight;
            sum2 += weight;
        }

        m_star_coords[3 * num_star + 1] = sum1 / sum2;
    }

    if (direction == 2)
    {
        for (int i = 0; i < m_points; i++)
        {
            weight = (vals_f[i] - fmin) / (fmax - fmin);
            sum1 += z_coords[i] * weight;
            sum2 += weight;
        }

        m_star_coords[3 * num_star + 2] = sum1 / sum2;
    }
}

// Finally update the centres either using Gaussian fitting procedure or centre
// of mass calculation.
void StarTracker::update_star_centres(double a_dt)
{
    if (m_direction == "x")
    {
        double starA_0 = find_centre(0, 0);
        if (abs((starA_0 - m_star_coords[0]) / a_dt) < 1.0 && starA_0 != 0)
        {
            m_star_coords[0] = starA_0;
        }
        else
        {
            find_max_min(0, 0);
        }
        double starB_0 = find_centre(1, 0);
        if ((abs(starB_0 - m_star_coords[3]) / a_dt) < 1.0 && starB_0 != 0)
        {
            m_star_coords[3] = starB_0;
        }
        else
        {
            find_max_min(1, 0);
        }
    }

    if (m_direction == "xy")
    {
        double starA_0 = find_centre(0, 0);
        if (abs((starA_0 - m_star_coords[0]) / a_dt) < 1.0 && starA_0 != 0)
        {
            m_star_coords[0] = starA_0;
        }
        else
        {
            find_max_min(0, 0);
        }
        double starA_1 = find_centre(0, 1);
        if (abs((starA_1 - m_star_coords[1]) / a_dt) < 1.0 && starA_1 != 0)
        {
            m_star_coords[1] = starA_1;
        }
        else
        {
            find_max_min(0, 1);
        }
        double starB_0 = find_centre(1, 0);
        if (abs((starB_0 - m_star_coords[3]) / a_dt) < 1.0 && starB_0 != 0)
        {
            m_star_coords[3] = starB_0;
        }
        else
        {
            find_max_min(1, 0);
        }
        double starB_1 = find_centre(1, 1);
        if (abs((starB_1 - m_star_coords[4]) / a_dt) < 1.0 && starB_1 != 0)
        {
            m_star_coords[4] = starB_1;
        }
        else
        {
            find_max_min(1, 1);
        }
    }

    if (m_direction == "xyz")
    {
        double starA_0 = find_centre(0, 0);
        m_star_coords[0] = starA_0;
        double starA_1 = find_centre(0, 1);
        m_star_coords[1] = starA_1;
        double starA_2 = find_centre(0, 2);
        m_star_coords[2] = starA_2;

        double starB_0 = find_centre(1, 0);
        m_star_coords[3] = starB_0;
        double starB_1 = find_centre(1, 1);
        m_star_coords[4] = starB_1;
        double starB_2 = find_centre(1, 2);
        m_star_coords[5] = starB_2;
    }
}

// Write all data to designated files
void StarTracker::write_to_dat(std::string a_filename, double a_dt,
                               double a_time, double a_restart_time,
                               bool a_first_step)
{
    double eps = 10e-8;
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

    star_centre_file.write_time_data_line(m_star_coords);
}

// Read a data line from the previous timestep
void StarTracker::read_old_centre_from_dat(std::string a_filename, double a_dt,
                                           double a_time, double a_restart_time,
                                           bool a_first_step)
{
    std::vector<double> data_line;
    if (a_time > a_dt / 3.)
    {
        SmallDataIO star_centre_file(a_filename, a_dt, a_time, a_restart_time,
                                     SmallDataIO::READ, a_first_step);
        star_centre_file.get_specific_data_line(data_line, a_time - a_dt);

        bool length_match = data_line.size() == m_num_stars * CH_SPACEDIM;

        if (length_match)
        {
            for (int i = 0; i < data_line.size(); i++)
            {
                m_star_coords[i] = data_line[i];
            }
        }
        else
        {
            for (int i = 0; i < m_star_coords.size(); i++)
            {
                m_star_coords[i] = NAN;
            }
            std::cout
                << "Array Size Mismatch While Loading Star Centre From File ! "
                << std::endl;
        }
    }
}

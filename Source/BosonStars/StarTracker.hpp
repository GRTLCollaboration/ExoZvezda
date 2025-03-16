/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef STARTRACKER_HPP_
#define STARTRACKER_HPP_

#include "AMRInterpolator.hpp"
#include "AlwaysInline.hpp"
#include "Lagrange.hpp"

/*
 * Class for tracking boson star positions by fitting a Gaussian to (1-conformal factor).
 */

class StarTracker
{
  private:
    int m_num_stars; // number of stars
    std::vector<std::array<double, CH_SPACEDIM>> m_puncture_coords; // here we will store the star positions
    std::array<double, CH_SPACEDIM> m_centre;
    int m_tracking_level; // levels to execute tracking on
    int m_points;         // number of points used for tracking
    std::vector<double> m_x_coords; // array to store x-coords
    std::vector<double> m_y_coords; // array to store y-coords
    std::vector<double> m_z_coords; // array to store z-coords
    std::vector<double> m_sigma_vector;     // vector to store the error
    std::vector<double> m_vals_shifted_chi; // vector to store (1-chi)
    double m_width_A;                       // width for fitting around star A
    double m_width_B;                       // width for fitting around star B
    std::string
        m_fitting_direction; // along which direction to fit (x or y or z)

    // saved pointer to external interpolator
    AMRInterpolator<Lagrange<4>> *m_interpolator;

  public:
    //! The constructor
    StarTracker() : m_interpolator(nullptr) {}

    //! set star locations on start (or restart)
    //! this needs to be done before 'setupAMRObject'
    //! if the stars' locations are required for Tagging Criteria
    void
    initialise_star_tracking(bool a_do_star_track, int a_number_of_stars,
                             const std::vector<std::array<double, CH_SPACEDIM>>
                                 &a_initial_star_centres,
                             int a_star_points, double a_star_track_width_A,
                             double a_star_track_width_B,
                             std::string a_fitting_direction)
    {
        m_num_stars = a_number_of_stars;
        m_points = a_star_points;

        m_x_coords.resize(m_points, 0);
        m_y_coords.resize(m_points, 0);
        m_z_coords.resize(m_points, 0);
        m_sigma_vector.resize(m_points, 0);
        m_vals_shifted_chi.resize(m_points, 0);

        m_width_A = a_star_track_width_A;
        m_width_B = a_star_track_width_B;
        m_fitting_direction = a_fitting_direction;

        m_puncture_coords = a_initial_star_centres;
    }

    ALWAYS_INLINE void
    set_interpolator(AMRInterpolator<Lagrange<4>> *a_interpolator)
    {
        m_interpolator = a_interpolator;
    }

    void set_up_fitting(int num_star, int fitting_direction);

    // General function for fitting the stars' centres
    double find_centre(int num_star, int fitting_direction);

    // Function for finding star positions near merger 
    void find_centre_merger(int num_star, int fitting_direction);

    void update_star_centres(double a_dt);

    void write_to_dat(std::string a_filename, double a_dt, double a_time,
                      double a_restart_time, bool a_first_step);

    void read_old_centre_from_dat(std::string a_filename, double a_dt,
                                  double a_time, double a_restart_time,
                                  bool a_first_step);

     // Execute the tracking and write out
     void execute_tracking(double a_time, double a_restart_time, double a_dt,
        const bool write_punctures = true);

    ALWAYS_INLINE const std::vector<std::array<double, CH_SPACEDIM>> &
    get_puncture_coords() const
    {
        return m_puncture_coords;
    }
};

#endif /* STARTRACKER_HPP_ */

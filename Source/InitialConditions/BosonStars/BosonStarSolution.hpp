
#ifndef BOSONSTARSOLUTION_HPP_
#define BOSONSTARSOLUTION_HPP_

#include "cmath"
#include <vector>

class BosonStarSolution
{

  private:         // private member variables/arrays
    double MM, A0; // Klein Gordon mass squared, KG scalr field central
                   // amplitude
    double PSC, OMC;
    // double PSC=2.52749283231717, OMC=0.160056080542207; // central density of
    // scalar field (0.272 for kaup)  PSC and OMC are central values of
    // conformal factor and lapse, not important as long as they are sensible
    // (i.e. order 1) PSC=2.0, OMC=0.5
    double lambda;  // phi 4 coupling in Klein gordon potential
    double sigma;   // 0.2 works with PC = 0.05 // parameter for solitonic stars
    bool solitonic; // false fro mini/lambda star. true for solitonic star
    bool BS_verbosity; // outputs more messages whilst finding the solution
    double EIGEN = 0;  // the desired eigenstate, 0 for ground
    int gridsize, adaptive_buffer;                // anywhere from 2k-200k is ok
    const int adaptive_stepsize_repetitions = 20; // 50; // 0 for no adaptive
    double L, dx, omega_ansatz, omega_true;       // L, length of domain, dx.
    double OM_INF, PSI_INF; // asymptotics of lapse and cpnformal factpr
    int matching_index;     // integer where growing mode becomes relevant
    double radius,
        compactness_value; // for storing the radius and compactness of the BS
    double eps = 10e-20;
    double omega_tolerance = 1e-20;

    std::vector<double> A;                  // scalar field modulus
    std::vector<double> dA;                 // scalar field modulus gradient
    std::vector<double> psi;                // conformal factor
    std::vector<double> dpsi;               // conformal factor gradient
    std::vector<double> omega;              // lapse
    std::vector<double> radius_array;       // isotropic radius
    std::vector<double> areal_radius_array; // areal radius
    std::vector<double> boson_mass;         // mass
    std::vector<double> adm_mass;           // mass
    std::vector<double> boson_radius;       // radius of the BS
    std::vector<double> compactness;        // compactness array

  private: // private member fucntions functions
    void rk4(const double ww_);
    void rk4_asymp(const int iter, const bool adaptive, const double ww_);
    void rk4_match(const int iter, const bool adaptive, const double ww_);
    double small_A_RHS(const double x, const double A, const double DA,
                       const double PSI, const double DPSI, const double OM,
                       const double ww_);
    double A_RHS(const double x, const double A, const double DA,
                 const double PSI, const double DPSI, const double OM,
                 const double ww_);
    double DA_RHS(const double x, const double A, const double DA,
                  const double PSI, const double DPSI, const double OM,
                  const double ww_);
    double PSI_RHS(const double x, const double A, const double DA,
                   const double PSI, const double DPSI, const double OM,
                   const double ww_);
    double DPSI_RHS(const double x, const double A, const double DA,
                    const double PSI, const double DPSI, const double OM,
                    const double ww_);
    double OMEGA_RHS(const double x, const double A, const double DA,
                     const double PSI, const double DPSI, const double OM,
                     const double ww_);
    void initialise();
    double amplitude_criterion(double omega, double index);
    int crossings();
    int zero_crossings();
    void truncate_solution();
    void force_flat(const int iter_crit);
    double bisect_omega(double omega_min, double omega_max);
    int found_zero_crossing();
    double find_upper_omega();
    int find_matching_index();
    int iofr(double rtarget);
    double V(const double A);
    double DV(const double A);
    void calculate_aspect_mass();
    void calculate_adm_mass();
    double calculate_radius();

  public:
    BosonStarSolution();
    void set_initialcondition_params(BosonStar_params_t m_params_BosonStar,
                                     Potential::params_t m_params_potential,
                                     const double max_r);
    double get_A_interp(const double r) const;
    double get_lapse_interp(const double r) const;
    double get_psi_interp(const double r) const;
    double get_dpsi_interp(const double r) const;
    double get_dA_interp(const double r) const;
    double get_dlapse_interp(const double r) const;
    double get_mass() const;
    double get_w() const;
    void shout() const;
    void output_csv();
    void main();
};

#include "BosonStarSolution.impl.hpp"

#endif /* BOSONSTARSOLUTION_HPP_ */

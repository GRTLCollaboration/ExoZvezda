
#ifndef BOSONSTARSOLUTION_HPP_
#define BOSONSTARSOLUTION_HPP_

#include "cmath"
#include <vector>

class BosonStarSolution
{

  private:
    double mu, A0; // Klein Gordon mass squared, central
                   // amplitude of the scalar field
    double PSC,
        OMC; // central values of the conformal factor and lapse (usually we set
             // PSC=2.0, OMC=0.5, but for more complex BS solutions, these
             // initial guesses may need to be more carefully iterated over)
    double lambda;     // phi^4 coupling in the potential
    double sigma;      // self-interaction term for solitonic stars
    bool solitonic;    // false for mini/repulsive star, true for solitonic star
    bool BS_verbosity; // outputs more messages whilst finding the solution
    double EIGEN = 0;  // the desired eigenstate, 0 for ground
    int gridsize;      // number of grid points (10^6 is very good)
    int adaptive_buffer; // number of gridpoints to intergate more carefully
    const int adaptive_stepsize_repetitions = 20; // 50; // 0 for no adaptive
    double L, dx;                                 // L, length of domain, dx.
    double omega_ansatz,
        omega_true; // ansatz of omega, square of BS frequency (note! here we
                    // work with squares, so don't forget to take a square root
                    // when finding the right value)
    double OM_INF, PSI_INF; // asymptotics of the lapse and conformal factpr
    int matching_index;     // integer where growing mode becomes relevant
    double eps = 10e-20;    // some very small number
    double omega_tolerance = 1e-20; // omega threshold used in the bisection
    int niter; // number of iterations for finding the solution

    void rk4(const double ww_);
    void rk4_asymp(const int iter, const bool adaptive, const double ww_);
    void rk4_match(const int iter, const bool adaptive, const double ww_);
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
    void truncate_solution();
    void force_to_zero(const int iter_crit);
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

    std::vector<double> A;            // scalar field modulus
    std::vector<double> dA;           // scalar field modulus gradient
    std::vector<double> psi;          // conformal factor
    std::vector<double> dpsi;         // conformal factor gradient
    std::vector<double> omega;        // lapse
    std::vector<double> radius_array; // isotropic radius
    std::vector<double> boson_mass;   // mass
    std::vector<double> adm_mass;     // mass
    double radius,
        compactness_value; // for storing the radius and compactness of the BS

    void set_initialcondition_params(BosonStar_params_t m_params_BosonStar,
                                     ComplexPotential::params_t m_params_potential,
                                     const double max_r);
    double get_A_interp(const double r) const;
    double get_lapse_interp(const double r) const;
    double get_psi_interp(const double r) const;
    double get_dpsi_interp(const double r) const;
    double get_dA_interp(const double r) const;
    double get_dlapse_interp(const double r) const;
    double get_BSfrequency() const;
    void output_csv();
    void main();
};

#include "BosonStarSolution.impl.hpp"

#endif /* BOSONSTARSOLUTION_HPP_ */

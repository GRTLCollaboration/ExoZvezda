#ifndef BOSONSTARSOLUTION_HPP_
#define BOSONSTARSOLUTION_HPP_

class BosonStarSolution
{

  private:         // private member variables/arrays
    double MM, PC; // Klein Gordon mass squared, KG scalr field central
                   // amplitude
    double G;      // 1./(M_PI*4.); // rescales KG field in ODE's
    double PSC = 2.,
           OMC =
               0.5; // central density of scalar field (0.272 for kaup)  PSC and
                    // OMC are central values of conformal factor and lapse, not
                    // important as long as they are sensible (i.e. order 1)
    double lambda;  // phi 4 coupling in Klein gordon potential
    double sigma;   // 0.2 works with PC = 0.05 // parameter for solitonic stars
    bool solitonic; // false fro mini/lambda star. true for solitonic star
    double EIGEN;   // the desired eigenstate, 0 for ground
    int gridsize, adaptive_buffer;                // anywhere from 2k-200k is ok
    const int adaptive_stepsize_repetitions = 20; // 50; // 0 for no adaptive
    double L, dx, WW, ww;                         // L, length of domain, dx.
    double OM_INF, PSI_INF; // asymptotics of lapse and cpnformal factpr
    int mid_int;            // integer where growing mode becomes relevant
    double upper_ww, lower_ww, middle_ww;
    double adm_mass, aspect_mass;
    double eps = 10e-20;
    double ww_tolerance = 10e-20;

    std::vector<double> p;            // scalar field modulus
    std::vector<double> dp;           // scalar field modulus gradient
    std::vector<double> psi;          // conformal factor
    std::vector<double> dpsi;         // conformal factor gradient
    std::vector<double> omega;        // lapse
    std::vector<double> radius_array; // radius

  private: // private member fucntions functions
    void rk4(const double ww_);
    void rk4_asymp(const int iter, const bool adaptive, const double ww_);
    double small_P_RHS(const double x, const double P, const double DP,
                       const double PSI, const double DPSI, const double OM,
                       const double ww_);
    double P_RHS(const double x, const double P, const double DP,
                 const double PSI, const double DPSI, const double OM,
                 const double ww_);
    double DP_RHS(const double x, const double P, const double DP,
                  const double PSI, const double DPSI, const double OM,
                  const double ww_);
    double PSI_RHS(const double x, const double P, const double DP,
                   const double PSI, const double DPSI, const double OM,
                   const double ww_);
    double DPSI_RHS(const double x, const double P, const double DP,
                    const double PSI, const double DPSI, const double OM,
                    const double ww_);
    double OMEGA_RHS(const double x, const double P, const double DP,
                     const double PSI, const double DPSI, const double OM,
                     const double ww_);
    void initialise();
    int crossings();
    void fix();
    void force_flat(const int iter_crit);
    double ww_min(const double WW_);
    double ww_max(const double WW_, const double lower_ww_);
    double ww_IB(double lower_ww_, double upper_ww_);
    double ww_IB_soliton(double lower_ww_, double upper_ww_);
    bool soliton_eigen();
    double find_WW();
    double find_WW_soliton();
    int find_midint();
    double V(const double P);
    double DV(const double P);

  public:
    BosonStarSolution();
    void set_initialcondition_params(BosonStar_params_t m_params_BosonStar,
                                     Potential::params_t m_params_potential,
                                     const double max_r);
    double get_p_interp(const double r) const;
    double get_lapse_interp(const double r) const;
    double get_psi_interp(const double r) const;
    double get_dpsi_interp(const double r) const;
    double get_dp_interp(const double r) const;
    double get_dlapse_interp(const double r) const;
    double get_mass() const;
    double get_w() const;
    double get_r(const double frac) const;
    void shout() const;
    void main();
};

#include "BosonStarSolution.impl.hpp"

#endif /* BOSONSTARSOLUTION_HPP_ */
/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

 #ifndef BOSONSTARSOLUTION_HPP_
 #define BOSONSTARSOLUTION_HPP_
 
 #include "BosonStarParams.hpp"
 #include "BosonStarSolver.hpp"
 #include "ComplexPotential.hpp"
 #include "Coordinates.hpp"
 #include "MatterCCZ4.hpp"
 #include "Tensor.hpp"
 #include "UserVariables.hpp" 
 #include "VarsTools.hpp"
 #include "parstream.H" 
 #include "simd.hpp"
 
 /*
  * Class for constructing initial data for equal-mass BS binary using equal-mass fix, as in https://arxiv.org/abs/2108.11995. This construction is recommended over plain superposition!
  */
 
 class BosonStarSolution
 {
 
   public:
     //! The constructor
     BosonStarSolution(BosonStar_params_t a_params_BosonStar,
                        ComplexPotential::params_t a_params_potential) : m_params_BosonStar(a_params_BosonStar),
                        m_params_potential(a_params_potential)
                  {
                  }

    void compute_3d_solution(const BosonStarSolver &sol, BosonStar_params_t a_params_BosonStar, Tensor<2, double> &star_metricLL, Tensor<2, double> &star_KLL,
    Tensor<1, double> &star_shift, double &star_lapse, double &star_phi_Re, double &star_phi_Im,
    double &star_Pi_Re, double &star_Pi_Im, const Tensor<1, double> &coords, double t) const
    {   
        std::cout << "In compute_3d_solution" << std::endl;
        //Initialise 
        FOR(i, j) { star_metricLL[i][j] = 0.0;
                    star_KLL[i][j] = 0.0; }
        
        double x = coords[0];
        double y = coords[1];
        double z = coords[2];

        double r = sqrt(x * x + y * y + z * z);

        std::cout << "Coords" << std::endl;
    
        // Import BS parameters 
        double rapidity = m_params_BosonStar.BS_rapidity;
        bool antiboson = m_params_BosonStar.antiboson;
        double separation = m_params_BosonStar.BS_separation;
        double impact_parameter = m_params_BosonStar.BS_impact_parameter;
    
        // First star positioning
        double c_ = cosh(rapidity);
        double s_ = sinh(rapidity);
        double v_ = tanh(rapidity);

        std::cout << "Interpolate" << std::endl; 

        double p_, dp_, omega_, omega_prime_, psi_, psi_prime_, pc_os;
        
       // First star physical variables
        p_ = sol.get_A_interp(r); // BS amplitude
        dp_ = sol.get_dA_interp(r); // derivative of BS amplitude
        omega_ = sol.get_lapse_interp(r); // lapse
         omega_prime_ = sol.get_dlapse_interp(r); // derivative of the lapse
         psi_ = sol.get_psi_interp(r); // conformal factor
         psi_prime_ = sol.get_dpsi_interp(r); // derivative of the conformal factor
         pc_os = psi_ * psi_ * c_ * c_ - omega_ * omega_ * s_ * s_;
       
        std::cout << "Gauge" << std::endl;

        // Get lapse, shift, BS frequency and phase
        star_lapse = omega_ * psi_ / (sqrt(pc_os));
        star_shift[0] = s_ * c_ * (psi_ * psi_ - omega_ * omega_) / (pc_os);
        star_shift[1] = 0.;
        star_shift[2] = 0.;
        double w_;
        if (antiboson)
        {
            w_ = -sol.get_BSfrequency();
        }
        else
        {
            w_ = sol.get_BSfrequency();
        }
        double phase_ = m_params_BosonStar.phase * M_PI + w_ * t;

        std::cout << "Gamma" << std::endl;

    
        // Write metric componnets of star 1
        star_metricLL[0][0] = psi_ * psi_;
        star_metricLL[1][1] = psi_ * psi_;
        star_metricLL[2][2] = pc_os;
    
        // Compute real and imaginary parts of the sclaar field and its conjugate momentum. Write into the evolution varaibles.
        star_phi_Re = p_ * cos(phase_);
        star_phi_Im = p_ * sin(phase_);
        star_Pi_Re =
            -(1. / star_lapse) * ((x / r) * (s_ - star_shift[0] * c_) * dp_ * cos(phase_) -
                               w_ * (c_ - star_shift[0] * s_) * p_ * sin(phase_));
        star_Pi_Re =
            -(1. / star_lapse) * ((x / r) * (s_ - star_shift[0] * c_) * dp_ * sin(phase_) +
                               w_ * (c_ - star_shift[0] * s_) * p_ * cos(phase_));
    
        // Fill extrinsic curvature components  
        star_KLL[2][2] = -star_lapse * s_ * x * psi_prime_ / (r * psi_);
        star_KLL[1][1] = star_KLL[2][2];
        star_KLL[0][1] = star_lapse * c_ * s_ * (y / r) *
                      (psi_prime_ / psi_ - omega_prime_ / omega_);
        star_KLL[0][2] = star_lapse * c_ * s_ * (z / r) *
                      (psi_prime_ / psi_ - omega_prime_ / omega_);
        star_KLL[1][0] = star_KLL[0][1];
        star_KLL[2][0] = star_KLL[0][2];
        star_KLL[2][1] = 0.;
        star_KLL[1][2] = 0.;
        star_KLL[0][0] = star_lapse * (x / r) * s_ * c_ * c_ *
                      (psi_prime_ / psi_ - 2. * omega_prime_ / omega_ +
                       v_ * v_ * omega_ * omega_prime_ * pow(psi_, -2));
        std::cout << "Finished compute_3d_solution" << std::endl;
    }

    void find_helfer_correction(const BosonStarSolver &sol, Tensor<2, double> &helferLL)
    {
    //Initialise 
    FOR(i, j) { helferLL[i][j] = 0.0;}
    // We now calculate the effect of star 1 on star 2, which is equal to the the effect of star 2 on star 1 (for equal-mass binaries!). This will then represent the value
    // to be substracted in the initial data from the position of object 2
    double c_ = cosh(m_params_BosonStar.BS_rapidity);
    double s_ = sinh(m_params_BosonStar.BS_rapidity);
    
    double t_p = (-m_params_BosonStar.BS_separation) * s_; // set /tilde{t} to zero
    double x_p = (-m_params_BosonStar.BS_separation) * c_;
    double z_p = 0.; // set /tilde{t} to zero
    double y_p = m_params_BosonStar.BS_impact_parameter;
    double r_p = sqrt(x_p * x_p + y_p * y_p + z_p * z_p);
    double p_p = sol.get_A_interp(r_p);
    double dp_p = sol.get_dA_interp(r_p);
    double omega_p = sol.get_lapse_interp(r_p);
    double omega_prime_p = sol.get_dlapse_interp(r_p);
    double psi_p = sol.get_psi_interp(r_p);
    double psi_prime_p = sol.get_dpsi_interp(r_p);
    double pc_os_p = psi_p * psi_p * c_ * c_ - omega_p * omega_p * s_ * s_;

    helferLL[1][1] = psi_p * psi_p;
    helferLL[2][2] = psi_p * psi_p;
    helferLL[0][0] = pc_os_p;
    }
 
     double central_amplitude1, central_amplitude2;
     double mass1, mass2;
     double frequency1, frequency2;
     double radius1, radius2;
     double compactness1, compactness2;
 
   protected:
        BosonStar_params_t m_params_BosonStar;
        ComplexPotential::params_t m_params_potential; 
 };
 
 
 #endif /* BOSONSTARSOLUTION_HPP_ */

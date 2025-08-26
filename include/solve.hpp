#pragma one
#include "discharge.hpp"
#include <Eigen/Dense>
#include <fstream>

class Solver: protected Discharge{
    private:
        double n_i_max_, T_e_anode_, T_e_cathode_;

    public:
        Solver();

        void set_init_conditions(double n_i_min, double n_i_max, double T_e_anode_eV, double T_e_cathode_eV);

        Plasma_state right_boundary_condition(std::size_t time_ind);

        Plasma_state left_boundary_condition(std::size_t time_ind);
        
        Electrons_state calc_electrons_state(const std::size_t time_ind, const std::size_t spatial_ind) const;

        Eigen::MatrixXd Omega_L(const std::size_t time_ind, const std::size_t spatial_ind) const;

        Eigen::MatrixXd Omega_R(const std::size_t time_ind, const std::size_t spatial_ind) const;

        Eigen::MatrixXd Lamda_plus(const std::size_t time_ind, const std::size_t spatial_ind) const;

        Eigen::MatrixXd Lamda_minus(const std::size_t time_ind, const std::size_t spatial_ind) const;

        Eigen::MatrixXd A_plus(const std::size_t time_ind, const std::size_t spatial_ind) const;

        Eigen::MatrixXd A_minus(const std::size_t time_ind, const std::size_t spatial_ind) const;

        Eigen::VectorXd right_hand_side(const std::size_t time_ind, const std::size_t spatial_ind) const;

        std::vector<std::vector<Plasma_state>> solve();

};
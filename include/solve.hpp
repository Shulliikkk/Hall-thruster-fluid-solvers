#pragma one
#include "discharge.hpp"
#include <Eigen/Dense>
#include <fstream>
#include <utility>

class Solver: protected Discharge{
    private:
        double n_i_max_, T_e_anode_, T_e_cathode_;

    public:
        Solver();

        void set_init_conditions(double n_i_min, double n_i_max, double T_e_anode_eV, double T_e_cathode_eV);

        Plasma_state right_boundary_condition(const std::size_t time_ind) const;

        Plasma_state left_boundary_condition(const std::size_t time_ind) const;

        Eigen::VectorXd hyperbolic_system_flux(const std::size_t time_ind, const Plasma_state curr_state) const;

        Eigen::VectorXd hyperbolic_numerical_flux(const std::size_t time_ind, const Plasma_state left_state, const Plasma_state right_state) const;

        Eigen::VectorXd hyperbolic_right_hand_side(const std::size_t time_ind, const std::size_t spatial_ind) const;

        Eigen::VectorXd tridiagonal_matrix_algorithm(const Eigen::VectorXd& left_secondary_diagonal, const Eigen::VectorXd& main_diagonal, const Eigen::VectorXd& right_secondary_diagonal, const Eigen::VectorXd& column) const;

        double electron_cross_field_mobility(const Plasma_state state, const double spatial_ind) const;

        std::pair<double, double> heat_flux_coefficient_plus(const std::size_t time_ind, const std::size_t spatial_ind, const double T_e_curr, const double T_e_next) const;

        std::pair<double, double> heat_flux_coefficient_minus(const std::size_t time_ind, const std::size_t spatial_ind, const double T_e_curr, const double T_e_prev) const;

        double flux_coefficient(const std::size_t time_ind, const std::size_t spatial_ind) const;

        double upwind_coefficient_plus(const std::size_t time_ind, const std::size_t spatial_ind) const;

        double upwind_coefficient_minus(const std::size_t time_ind, const std::size_t spatial_ind) const;

        std::pair<std::array<Eigen::VectorXd, 3>, Eigen::MatrixXd> Jacobian_matrix(const std::size_t time_ind, const Eigen::VectorXd vec_T_e) const;

        Eigen::VectorXd nonlin_equation_function(const std::size_t time_ind, const Eigen::VectorXd T_e_vec) const;

        std::vector<std::vector<Plasma_state>> solve(const double eps = 0.5, const size_t max_iters = 100);
};
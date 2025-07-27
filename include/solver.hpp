#pragma one
#include "walls.hpp"

class Solver: protected Walls {
    public:
        Solver(double discharge_voltage, double v_n, double input_mass_flow, std::string dataset_ioniz_name, std::string dataset_neutral_name);

        Electrons_state calc_electrons_state(const std::size_t time_step, const std::size_t space_step) const;

        Eigen::MatrixX<double> Omega_L(const std::size_t time_step, const std::size_t space_step) const;

        Eigen::MatrixX<double> Omega_R(const std::size_t time_step, const std::size_t space_step) const;

        Eigen::MatrixX<double> Lamda(const std::size_t time_step, const std::size_t space_step) const;

        Eigen::VectorX<double> flux(const std::size_t time_step, const std::size_t space_step) const;

        Eigen::VectorX<double> right_hand_side(const std::size_t time_step, const std::size_t space_step) const;

        Eigen::MatrixX<double> matrix_halph_plus(const std::size_t time_step, const std::size_t space_step) const; // Holodov A. S. 1978

        Eigen::MatrixX<double> matrix_halph_minus(const std::size_t time_step, const std::size_t space_step) const; // Holodov A. S. 1978

        Eigen::VectorX<double> num_flux_halph_plus(const std::size_t time_step, const std::size_t space_step) const;

        Eigen::VectorX<double> num_flux_halph_minus(const std::size_t time_step, const std::size_t space_step) const;

        std::array<std::array<Plasma_state, grid::M>, grid::N> solve();
};
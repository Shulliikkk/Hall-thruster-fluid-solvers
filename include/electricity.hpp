#pragma one
#include <array>
#include <vector>
#include <string>
#include <Eigen/Dense>
#include "rapidcsv.h"
#include "interpolator.hpp"
#include "states.hpp"
#include "constants.hpp"

class Electricity {
    protected:
        const double discharge_voltage;
        const double v_n, input_mass_flow;
        std::vector<double> Electron_energy_data;
        std::vector<double> beta_n_data;
        std::vector<double> beta_i_data;
        std::array<std::array<Plasma_state, grid::M>, grid::N> states;

    public:
        Electricity(double discharge_voltage, double v_n, double input_mass_flow, std::string dataset_ioniz_name, std::string dataset_neutral_name);

        template <typename func>
        double trapez_numeric_int(func integrand_func, const std::size_t time_step, const std::size_t space_step) const {
            double sum = 0;
            for (std::size_t space_curr_step = 0; space_curr_step < grid::M - 1; space_curr_step++) {
                sum += (integrand_func(time_step, space_step + 1) + integrand_func(time_step, space_step)) / 2 * grid::h;
            }
            return sum;
        }

        double magnetic_field(const std::size_t space_step) const;

        double cicl_friq(const std::size_t space_step) const;

        double electron_neutral_coll_friq(const std::size_t time_step, const std::size_t space_step) const;

        double anomal_coll_friq(const std::size_t space_step) const;

        double moment_exch_coll_friq(const std::size_t time_step, const std::size_t space_step) const;

        double ioniz_coll_friq(const std::size_t time_step, const std::size_t space_step) const;

        double electron_cross_field_mobility(const std::size_t time_step, const std::size_t space_step) const;

        double numerat_integr_func_disch_curr(const std::size_t time_step, const std::size_t space_step) const;

        double denominat_integr_func_disch_curr(const std::size_t time_step, const std::size_t space_step) const;

        double discharge_curr(const std::size_t time_step, const std::size_t space_step) const;

        double electric_field(const std::size_t time_step, const std::size_t space_step) const;
};
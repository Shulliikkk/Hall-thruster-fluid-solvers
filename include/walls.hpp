#pragma one
#include <utility>
#include "electricity.hpp"


class Walls: protected Electricity {
    public:
        Walls(double discharge_voltage, double v_n, double input_mass_flow, std::string dataset_ioniz_name, std::string dataset_neutral_name);

        double SEE_coef(const std::size_t time_step, const std::size_t space_step) const;

        std::pair<double, double> wall_coll_friq(const std::size_t time_step, const std::size_t space_step) const;

        double wall_potential(const std::size_t time_step, const std::size_t space_step) const;

        double energy_lose_rate(const std::size_t time_step, const std::size_t space_step) const;
};
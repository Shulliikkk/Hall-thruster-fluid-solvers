#include "walls.hpp"

Walls::Walls(double discharge_voltage, double v_n, double input_mass_flow, std::string dataset_ioniz_name, std::string dataset_neutral_name): Electricity(discharge_voltage, v_n, input_mass_flow, dataset_ioniz_name, dataset_neutral_name) {}

double Walls::SEE_coef(const size_t time_step, const size_t space_step) const {
    Plasma_state curr_state = states[time_step][space_step];
    return BNSiO2::gamma_func * BNSiO2::a * std::pow((phys::k_B * curr_state.T_e), BNSiO2::b);
}

std::pair<double, double> Walls::wall_coll_friq(const size_t time_step, const size_t space_step) const {
    Plasma_state curr_state = states[time_step][space_step];
    double nu_iw = 1.2 * std::sqrt(phys::k_B * curr_state.T_e / phys::m_i) / (geom::R_out - geom::R_in);
    double nu_ew = nu_iw / (1 - SEE_coef(time_step, space_step));
    return std::make_pair(nu_iw, nu_ew);
}

double Walls::wall_potential(const size_t time_step, const size_t space_step) const {
    Plasma_state curr_state = states[time_step][space_step];
    return phys::k_B * curr_state.T_e / phys::e * std::log((1 - SEE_coef(time_step, space_step)) * std::sqrt(phys::m_i / (2 * M_PI * phys::m_e)));
}

double Walls::energy_lose_rate(const size_t time_step, const size_t space_step) const {
    Plasma_state curr_state = states[time_step][space_step];
    double nu_ew = wall_coll_friq(time_step, space_step).second;
    double phi_w = wall_potential(time_step, space_step);
    return nu_ew * (2 * phys::k_B * curr_state.T_e + (1 - SEE_coef(time_step, space_step)) * phys::e * phi_w);
}
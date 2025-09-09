#include "walls.hpp"

double Walls::SEE_coef(const Plasma_state curr_state) const {
    return std::min(wall::gamma_func * wall::a * std::pow(curr_state.T_e / phys::e, wall::b), 0.983);
}

double Walls::ions_wall_coll_friq(const std::size_t time_ind, const std::size_t spatial_ind) const {
    Plasma_state curr_state = states[time_ind][spatial_ind];
    return (spatial_ind * sim::spatial_step <= geom::channel_length) ? 1.2 * std::sqrt(curr_state.T_e / phys::m_i) / (geom::channel_radius_out - geom::channel_radius_in) : 0;
}

double Walls::electrons_wall_coll_friq(const std::size_t time_ind, const std::size_t spatial_ind) const {
    Plasma_state curr_state = states[time_ind][spatial_ind];
    return ions_wall_coll_friq(time_ind, spatial_ind) / (1 - SEE_coef(curr_state));
}

double Walls::wall_potential(const std::size_t time_ind, const std::size_t spatial_ind) const {
    Plasma_state curr_state = states[time_ind][spatial_ind];
    double SEE = SEE_coef(curr_state);
    return curr_state.T_e * std::log((1 - SEE) * std::sqrt(phys::m_i / (2 * M_PI * phys::m_e)));
    return 0;
}

double Walls::energy_lose_rate(const std::size_t time_ind, const std::size_t spatial_ind) const {
    Plasma_state curr_state = states[time_ind][spatial_ind];
    double nu_ew = electrons_wall_coll_friq(time_ind, spatial_ind);
    double phi_w = wall_potential(time_ind, spatial_ind);
    return nu_ew * (2 * curr_state.T_e + (1 - SEE_coef(curr_state)) * phi_w);
}
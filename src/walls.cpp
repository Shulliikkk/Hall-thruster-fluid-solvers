#include "walls.hpp"

double Walls::SEE_coef(const std::size_t time_ind, const std::size_t spatial_ind) const {
    Plasma_state curr_state = states[time_ind][spatial_ind];
    return wall::gamma_func * wall::a * std::pow(curr_state.T_e / phys::e, wall::b);
}

double Walls::ions_wall_coll_friq(const std::size_t time_ind, const std::size_t spatial_ind) const {
    Plasma_state curr_state = states[time_ind][spatial_ind];
    return (spatial_ind * sim::spatial_step <= geom::channel_length) ? 1.2 * std::sqrt(curr_state.T_e / phys::m_i) / (geom::channel_radius_out - geom::channel_radius_in) : 0;
}

double Walls::electrons_wall_coll_friq(const std::size_t time_ind, const std::size_t spatial_ind) const {
    return (spatial_ind * sim::spatial_step <= geom::channel_length) ? ions_wall_coll_friq(time_ind, spatial_ind) / (1 - SEE_coef(time_ind, spatial_ind)) : 0;
}

double Walls::wall_potential(const std::size_t time_ind, const std::size_t spatial_ind) const {
    Plasma_state curr_state = states[time_ind][spatial_ind];
    double SEE = SEE_coef(time_ind, spatial_ind);
    if (SEE < 1) {
        return curr_state.T_e / phys::e * std::log((1 - SEE_coef(time_ind, spatial_ind)) * std::sqrt(phys::m_i / (2 * M_PI * phys::m_e)));
    }
    else {
        return -1.02 * curr_state.T_e / phys::e;
    }
    return 0;
}

double Walls::energy_lose_rate(const std::size_t time_ind, const std::size_t spatial_ind) const {
    Plasma_state curr_state = states[time_ind][spatial_ind];
    double nu_ew = electrons_wall_coll_friq(time_ind, spatial_ind);
    double phi_w = wall_potential(time_ind, spatial_ind);
    return nu_ew * (2 * curr_state.T_e + (1 - SEE_coef(time_ind, spatial_ind)) * phys::e * phi_w);
}
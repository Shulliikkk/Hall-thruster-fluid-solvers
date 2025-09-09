#include "conduction.hpp"

Conduction::Conduction() {
    states.resize(sim::time_steps, std::vector<Plasma_state>(sim::spatial_steps));
}

double Conduction::cicl_friq(const std::size_t spatial_ind) const {
    double B = magfld::magnetic_field(spatial_ind);
    return phys::e * B / phys::m_e;
}

double Conduction::electron_neutral_coll_friq(const std::size_t time_ind, const std::size_t spatial_ind) const {
    Plasma_state curr_state = states[time_ind][spatial_ind];
    JustInterp::LinearInterpolator<double> interpolator(electron_energy_data, beta_n_data);
    return curr_state.n_n * interpolator(curr_state.T_e / phys::e);
}

double Conduction::anomal_coll_friq(const std::size_t spatial_ind) const {
    double anomal_coll_friq_coef = (spatial_ind * sim::spatial_step <= geom::channel_length) ? 0.1 : 1;
    return  anomal_coll_friq_coef / 16. * cicl_friq(spatial_ind);
}

double Conduction::moment_exch_coll_friq(const std::size_t time_ind, const std::size_t spatial_ind) const {
    return electron_neutral_coll_friq(time_ind, spatial_ind) + anomal_coll_friq(spatial_ind) + electrons_wall_coll_friq(time_ind, spatial_ind);
}

double Conduction::ioniz_coll_friq(const std::size_t time_ind, const std::size_t spatial_ind) const {
    Plasma_state curr_state = states[time_ind][spatial_ind];
    JustInterp::LinearInterpolator<double> interpolator(electron_energy_data, beta_i_data);
    return curr_state.n_n * interpolator(curr_state.T_e / phys::e);
}

double Conduction::electron_cross_field_mobility(const std::size_t time_ind, const std::size_t spatial_ind) const {
    Plasma_state curr_state = states[time_ind][spatial_ind];
    double nu_m = moment_exch_coll_friq(time_ind, spatial_ind);
    double omega_c = cicl_friq(spatial_ind);
    return phys::e / (phys::m_e * nu_m) * 1. / (1. + omega_c * omega_c / (nu_m * nu_m));
}
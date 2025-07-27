#include "electricity.hpp"

Electricity::Electricity(double discharge_voltage, double v_n, double input_mass_flow, std::string dataset_ioniz_name, std::string dataset_neutral_name): discharge_voltage(discharge_voltage), v_n(v_n), input_mass_flow(input_mass_flow) {
    rapidcsv::Document ioniz_rate_dataset(dataset_ioniz_name);
    rapidcsv::Document neutral_rate_dataset(dataset_neutral_name);
    Electron_energy_data = neutral_rate_dataset.GetColumn<double>("EnergyEv");
    beta_n_data = neutral_rate_dataset.GetColumn<double>("betan");
    beta_i_data =  ioniz_rate_dataset.GetColumn<double>("betai");
    // initial condition
    for (size_t space_step = 1; space_step < grid::M - 1; space_step++) {
        states[0][space_step].n_i = 0;
        states[0][space_step].v_i = 0;
        states[0][space_step].T_e = 0;
        states[0][space_step].n_n = 0;
    }
}

double Electricity::magnetic_field(const size_t space_step) const {
    return magfield::B_0 * std::exp(-(space_step * grid::h - grid::L) * (space_step * grid::h - grid::L) / (2 * magfield::delta_B * magfield::delta_B));
}

double Electricity::cicl_friq(const size_t space_step) const {
    double B = magnetic_field(space_step);
    return phys::e * B / phys::m_e;
}

double Electricity::electron_neutral_coll_friq(const size_t time_step, const size_t space_step) const {
    Plasma_state curr_state = states[time_step][space_step];
    interp::Cubic_Spline interpolator(Electron_energy_data, beta_n_data);
    return curr_state.n_n * interpolator.interpolate(phys::k_B * curr_state.T_e / phys::e);
}

double Electricity::anomal_coll_friq(const size_t space_step) const {
    return 1. / 16. * cicl_friq(space_step);
}

double Electricity::moment_exch_coll_friq(const size_t time_step, const size_t space_step) const {
    return electron_neutral_coll_friq(time_step, space_step) + anomal_coll_friq(space_step);
}

double Electricity::ioniz_coll_friq(const size_t time_step, const size_t space_step) const {
    Plasma_state curr_state = states[time_step][space_step];
    interp::Cubic_Spline interpolator(Electron_energy_data, beta_i_data);
    return curr_state.n_n * interpolator.interpolate(phys::k_B * curr_state.T_e / phys::e);
}

double Electricity::electron_cross_field_mobility(const size_t time_step, const size_t space_step) const {
    double nu_m = moment_exch_coll_friq(time_step, space_step);
    double omega_c = cicl_friq(space_step);
    return phys::e / (phys::m_e * nu_m) * 1 / (1 + omega_c * omega_c / (nu_m * nu_m));
}

double Electricity::numerat_integr_func_disch_curr(const size_t time_step, const size_t space_step) const {
    Plasma_state curr_state = states[time_step][space_step];
    Plasma_state next_state = states[time_step][space_step + 1];
    double mu_perp = electron_cross_field_mobility(time_step, space_step);
    return curr_state.v_i / mu_perp + 1 / (phys::e * curr_state.n_i) * phys::k_B * (next_state.n_i * next_state.T_e - curr_state.n_i * curr_state.T_e) / grid::h;
}

double Electricity::denominat_integr_func_disch_curr(const size_t time_step, const size_t space_step) const {
    double mu_perp = electron_cross_field_mobility(time_step, space_step);
    Plasma_state curr_state = states[time_step][space_step];
    return 1 / (phys::e * curr_state.n_i * mu_perp);
}

double Electricity::discharge_curr(const size_t time_step, const size_t space_step) const {
    auto numerator = [this](size_t t, size_t s) { return numerat_integr_func_disch_curr(t, s); };
    auto denominator = [this](size_t t, size_t s) { return denominat_integr_func_disch_curr(t, s); };
    return -(discharge_voltage - trapez_numeric_int(numerator, time_step, space_step)) /
            trapez_numeric_int(denominator, time_step, space_step);
}

double Electricity::electric_field(const size_t time_step, const size_t space_step) const {
    Plasma_state curr_state = states[time_step][space_step];
    Plasma_state next_state = states[time_step][space_step + 1];
    double j = discharge_curr(time_step, space_step);
    double mu_perp = electron_cross_field_mobility(time_step, space_step);
    return -j / (phys::e * curr_state.n_i * mu_perp) +
            curr_state.v_i / mu_perp + 1 / (phys::e * curr_state.n_i) * (next_state.n_i * next_state.T_e -
            curr_state.n_i * curr_state.T_e) / grid::h;

}
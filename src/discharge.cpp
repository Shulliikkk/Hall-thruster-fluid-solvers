#include "discharge.hpp"

Discharge::Discharge() : Conduction() {}

double Discharge::derive_dniTe_dz(const std::size_t time_ind, const std::size_t spatial_ind) const {
    Plasma_state curr_state = states[time_ind][spatial_ind];
    auto n_iT_e = [](Plasma_state state) {return state.n_i * state.T_e;};
    if (spatial_ind == 0) {
        Plasma_state next_state = states[time_ind][spatial_ind + 1];
        Plasma_state next_next_state = states[time_ind][spatial_ind + 2];
        return (-3 * n_iT_e(curr_state) + 4 * n_iT_e(next_state) - n_iT_e(next_next_state)) / (2 * sim::spatial_step);  
    }
    else if (spatial_ind == sim::spatial_steps - 1) {
        Plasma_state prev_prev_state = states[time_ind][spatial_ind - 2];
        Plasma_state prev_state = states[time_ind][spatial_ind - 1];
        return (3 * n_iT_e(curr_state) - 4 * n_iT_e(prev_state) + n_iT_e(prev_prev_state)) / (2 * sim::spatial_step);  
    }
    else {
        Plasma_state prev_state = states[time_ind][spatial_ind - 1];
        Plasma_state next_state = states[time_ind][spatial_ind + 1];
        return (n_iT_e(next_state) - n_iT_e(prev_state)) / (2 * sim::spatial_step);
    }
}

double Discharge::numerat_integr_func_disch_curr(const std::size_t time_ind, const std::size_t spatial_ind) const {
    Plasma_state curr_state = states[time_ind][spatial_ind];
    double mu_perp = electron_cross_field_mobility(time_ind, spatial_ind);
    //if (time_ind == 122) std::cout << "num" << ' ' << time_ind << ' ' << spatial_ind << ' ' << mu_perp << ' '  <<  derive_dniTe_dz(time_ind, spatial_ind) << ' ' << curr_state.v_i / mu_perp + derive_dniTe_dz(time_ind, spatial_ind) / (curr_state.n_i * phys::e) << std::endl;
    return curr_state.v_i / mu_perp + derive_dniTe_dz(time_ind, spatial_ind) / (curr_state.n_i * phys::e);
}

double Discharge::denominat_integr_func_disch_curr(const std::size_t time_ind, const std::size_t spatial_ind) const {
    double mu_perp = electron_cross_field_mobility(time_ind, spatial_ind);
    Plasma_state curr_state = states[time_ind][spatial_ind];
    //std::cout << "num" << 1. / (phys::e * curr_state.n_i * mu_perp) << std::endl;
    return 1. / (phys::e * curr_state.n_i * mu_perp);
}

double Discharge::discharge_curr(const std::size_t time_ind) const {
    auto numerator = [this](std::size_t t, std::size_t s) { return numerat_integr_func_disch_curr(t, s); };
    auto denominator = [this](std::size_t t, std::size_t s) { return denominat_integr_func_disch_curr(t, s); };
    //std::cout << trapez_numeric_int(numerator, time_ind) << ' ' << trapez_numeric_int(denominator, time_ind) << std::endl;
    return (disch::discharge_voltage + trapez_numeric_int(numerator, time_ind)) /
        trapez_numeric_int(denominator, time_ind); // Нужен ли минус в знаменателе и после разрядного напряжения???
}

double Discharge::electric_field(const std::size_t time_ind, const std::size_t spatial_ind) const {
    Plasma_state curr_state = states[time_ind][spatial_ind];
    double j = discharge_curr(time_ind);
    double mu_perp = electron_cross_field_mobility(time_ind, spatial_ind);
    return j / (phys::e * curr_state.n_i * mu_perp) -
            curr_state.v_i / mu_perp - derive_dniTe_dz(time_ind, spatial_ind) / (curr_state.n_i * phys::e); // Знак???

}
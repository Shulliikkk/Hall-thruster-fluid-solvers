#include "solve.hpp"
#include <omp.h>

Solver::Solver(): Discharge() {}

void Solver::set_init_conditions(double n_i_min, double n_i_max, double T_e_anode_eV, double T_e_cathode_eV) {
    double T_e_max = disch::discharge_voltage * phys::e / 100;
    double T_e_anode = T_e_anode_eV * phys::e, T_e_cathode = T_e_cathode_eV * phys::e;
    n_i_max_ = n_i_max;
    T_e_anode_ = T_e_anode;
    T_e_cathode_ = T_e_cathode;
    for (std::size_t spatial_ind = 0; spatial_ind < sim::spatial_steps; spatial_ind++) {
        double v_bohm_anode = -std::sqrt(T_e_anode / phys::m_i);
        double v_max = std::sqrt(2 * phys::e * disch::discharge_voltage / phys::m_i);
        double axial_coordianat = spatial_ind * sim::spatial_step;
        states[0][spatial_ind].n_i = std::sqrt(disch::discharge_voltage / 300) * disch::mass_flow / 5e-6 * (n_i_min + (n_i_max - n_i_min) *
                                      std::exp(-(axial_coordianat - geom::channel_length / 2) / (geom::channel_length / 3) * (axial_coordianat - geom::channel_length / 2) / (geom::channel_length / 3)));

        states[0][spatial_ind].v_i = (spatial_ind * sim::spatial_step < geom::channel_length) ? v_bohm_anode +
                                      2. / 3. * (v_max - v_bohm_anode) * (axial_coordianat / geom::channel_length) * (axial_coordianat / geom::channel_length) :
                                      1. / 3. * (v_bohm_anode + v_max) * (1 - (axial_coordianat - geom::channel_length) / (sim::simulation_area_length -
                                      geom::channel_length)) + v_max * ((axial_coordianat - geom::channel_length) / (sim::simulation_area_length - geom::channel_length));

        states[0][spatial_ind].T_e = (1 - axial_coordianat / sim::simulation_area_length) * T_e_anode + axial_coordianat / sim::simulation_area_length * T_e_cathode +
                                      T_e_max * std::exp(-(axial_coordianat - geom::channel_length) / (geom::channel_length / 3) *
                                      (axial_coordianat - geom::channel_length) / (geom::channel_length / 3));
        
        double n_n_anode =  disch::mass_flow / (phys::m_i * geom::channel_area * disch::v_n) - states[0][0].n_i * states[0][0].v_i / disch::v_n;
        double n_n_cathode = 0.01 * n_n_anode;
        states[0][spatial_ind].n_n = 0.5 * (n_n_anode + n_n_cathode + (n_n_cathode - n_n_anode) * std::tanh((axial_coordianat - geom::channel_length / 2) /
                                      (geom::channel_length / 24)));
    }
}

Plasma_state Solver::right_boundary_condition(const std::size_t time_ind) const{
    Plasma_state rbc;
    Plasma_state prev_state = states[time_ind][sim::spatial_steps - 2];
    double T_e_max = disch::discharge_voltage * phys::e / 100;
    rbc.n_i = prev_state.n_i;
    rbc.v_i = prev_state.v_i;
    rbc.T_e = (1 - (sim::spatial_steps - 1) * sim::spatial_step / sim::simulation_area_length) * T_e_anode_ + (sim::spatial_steps - 1) * sim::spatial_step / sim::simulation_area_length * T_e_cathode_ +
        T_e_max * std::exp(-((sim::spatial_steps - 1) * sim::spatial_step - geom::channel_length) / (geom::channel_length / 3) *
        ((sim::spatial_steps - 1) * sim::spatial_step - geom::channel_length) / (geom::channel_length / 3));
    rbc.n_n = prev_state.n_n;
    return rbc;
}

Plasma_state Solver::left_boundary_condition(const std::size_t time_ind) const {
    Plasma_state lbc;
    double T_e_max = disch::discharge_voltage * phys::e / 100;
    double v_bohm_anode = -std::sqrt(T_e_anode_ / phys::m_i);
    lbc.n_i = std::sqrt(disch::discharge_voltage / 300) * disch::mass_flow / 5e-6 * n_i_max_ * std::exp(-9. / 4.);
    lbc.T_e = T_e_anode_ + T_e_max * std::exp(-9);
    lbc.v_i = v_bohm_anode;
    double n_n_anode =  disch::mass_flow / (phys::m_i * geom::channel_area * disch::v_n) - lbc.n_i * lbc.v_i / disch::v_n;
    double n_n_cathode = 0.01 * n_n_anode;
    lbc.n_n = 0.5 * (n_n_anode + n_n_cathode + (n_n_cathode - n_n_anode) * std::tanh(-12));
    return lbc;
}

Eigen::VectorXd Solver::hyperbolic_system_flux(const std::size_t time_ind, const Plasma_state curr_state) const {
    double j = discharge_curr(time_ind);
    double v_e = curr_state.v_i - j / (phys::e * curr_state.n_i);
    Eigen::VectorXd flux(sim::Dim - 1);
    flux << curr_state.n_i * curr_state.v_i, curr_state.n_i * curr_state.v_i * curr_state.v_i, curr_state.n_n * disch::v_n;
    return flux;
}

Eigen::VectorXd Solver::hyperbolic_numerical_flux(const std::size_t time_ind, const Plasma_state left_state, const Plasma_state right_state) const {
    double s_right_max = std::max(std::abs(right_state.v_i - std::sqrt(right_state.T_e / phys::m_i)), std::abs(right_state.v_i + std::sqrt(right_state.T_e / phys::m_i)));
    double s_left_max = std::max(std::abs(left_state.v_i - std::sqrt(left_state.T_e / phys::m_i)), std::abs(left_state.v_i + std::sqrt(left_state.T_e / phys::m_i)));
    double s_max = std::max(s_left_max, s_right_max);
    Eigen::VectorXd right(sim::Dim - 1);
    Eigen::VectorXd left(sim::Dim - 1);
    right << right_state.n_i, right_state.n_i * right_state.v_i, right_state.n_n;
    left << left_state.n_i, left_state.n_i * left_state.v_i, left_state.n_n;
    return 0.5 * (hyperbolic_system_flux(time_ind, right_state) + hyperbolic_system_flux(time_ind, left_state)) - 0.5 * s_max * (right - left);
}

Eigen::VectorXd Solver::hyperbolic_right_hand_side(const std::size_t time_ind, const std::size_t spatial_ind) const {
    Plasma_state curr_state = states[time_ind][spatial_ind];
    double nu_i = ioniz_coll_friq(time_ind, spatial_ind);
    double nu_iw = ions_wall_coll_friq(time_ind, spatial_ind);
    double j = discharge_curr(time_ind);
    double E =  electric_field(time_ind, spatial_ind);
    double delta_eps_w = energy_lose_rate(time_ind, spatial_ind);
    Eigen::VectorXd rhs(sim::Dim - 1);
    rhs << curr_state.n_i * (nu_i - nu_iw),
            phys::e * curr_state.n_i * E / phys::m_i + curr_state.n_i * (nu_i * disch::v_n - nu_iw * curr_state.v_i),
            - curr_state.n_i * (nu_i - nu_iw); 
    return rhs;
}

Eigen::VectorXd Solver::tridiagonal_matrix_algorithm(const Eigen::VectorXd& left_secondary_diagonal, 
                                                    const Eigen::VectorXd& main_diagonal,
                                                    const Eigen::VectorXd& right_secondary_diagonal,
                                                    const Eigen::VectorXd& column) const {
    std::size_t matrix_size = main_diagonal.size();
    if (matrix_size == 0) {
        return Eigen::VectorXd();
    }
    Eigen::VectorXd p(matrix_size - 1);
    Eigen::VectorXd q(matrix_size - 1);
    Eigen::VectorXd x(matrix_size);
    p(0) = -right_secondary_diagonal(0) / main_diagonal(0);
    q(0) = column(0) / main_diagonal(0);
    for (std::size_t i = 1; i < matrix_size - 1; i++) {
        double denominator = left_secondary_diagonal(i - 1) * p(i - 1) + main_diagonal(i);
        p(i) = -right_secondary_diagonal(i) / denominator;
        q(i) = (column(i) - left_secondary_diagonal(i - 1) * q(i - 1)) / denominator;
    }
    if (matrix_size > 1) {
        std::size_t last = matrix_size - 1;
        double denominator = left_secondary_diagonal(last - 1) * p(last - 1) + main_diagonal(last);
        x(last) = (column(last) - left_secondary_diagonal(last - 1) * q(last - 1)) / denominator;
    } else {
        x(0) = q(0);
    }
    for (int i = static_cast<int>(matrix_size) - 2; i >= 0; i--) {
        x(i) = p(i) * x(i + 1) + q(i);
    }
    return x;
}

double Solver::electron_cross_field_mobility(const Plasma_state state, const double spatial_ind) const {
    JustInterp::LinearInterpolator<double> interpolator(electron_energy_data, beta_n_data);
    double B = magfld::magnetic_field(spatial_ind);
    double omega_c = phys::e * B / phys::m_e;
    double nu_n = state.n_n * interpolator(state.T_e / phys::e);
    double anomal_coll_friq_coef = (spatial_ind * sim::spatial_step <= geom::channel_length) ? 0.1 : 1;
    double nu_ano = anomal_coll_friq_coef / 16. * cicl_friq(spatial_ind);
    double nu_wall = (spatial_ind * sim::spatial_step <= geom::channel_length) ? 1.2 * std::sqrt(state.T_e / phys::m_i) / (geom::channel_radius_out - geom::channel_radius_in) / (1 - SEE_coef(state)) : 0;
    double nu_m = nu_n + nu_ano + nu_wall;
    return  phys::e / (phys::m_e * nu_m) * 1. / (1. + omega_c * omega_c / (nu_m * nu_m));
}

std::pair<double, double> Solver::heat_flux_coefficient_plus(const std::size_t time_ind, const std::size_t spatial_ind, const double T_e_curr, const double T_e_next) const {
    Plasma_state curr_state = states[time_ind + 1][spatial_ind];
    Plasma_state next_state = states[time_ind + 1][spatial_ind + 1];
    curr_state.T_e = T_e_curr;
    next_state.T_e = T_e_next;
    double mu_perp_next = electron_cross_field_mobility(next_state, spatial_ind + 1);
    double mu_perp_curr = electron_cross_field_mobility(curr_state, spatial_ind);
    std::pair<double, double> HFC_plus(-5. / (3. * phys::e) * (mu_perp_next + mu_perp_curr) / 2 * (curr_state.n_i + next_state.n_i) / 2 * (T_e_curr + T_e_next) / 2, 
        -5. / (3. * phys::e) * (mu_perp_next + mu_perp_curr) / 2 * (curr_state.n_i + next_state.n_i) / 2); // HFC, d/dT_e(HFC)
    return HFC_plus;
}

std::pair<double, double> Solver::heat_flux_coefficient_minus(const std::size_t time_ind, const std::size_t spatial_ind, const double T_e_curr, const double T_e_prev) const {
    Plasma_state curr_state = states[time_ind + 1][spatial_ind];
    Plasma_state prev_state = states[time_ind + 1][spatial_ind - 1];
    curr_state.T_e = T_e_curr;
    prev_state.T_e = T_e_prev;
    double mu_perp_prev = electron_cross_field_mobility(prev_state, spatial_ind - 1);
    double mu_perp_curr = electron_cross_field_mobility(curr_state, spatial_ind);
    std::pair<double, double> HFC_minus(-5. / (3. * phys::e) * (mu_perp_prev + mu_perp_curr) / 2 * (curr_state.n_i + prev_state.n_i) / 2 * (T_e_curr + T_e_prev) / 2, 
        -5. / (3. * phys::e) * (mu_perp_prev + mu_perp_curr) / 2 * (curr_state.n_i + prev_state.n_i) / 2); // HFC, d/dT_e(HFC)
    return HFC_minus;
}

double Solver::flux_coefficient(const std::size_t time_ind, const std::size_t spatial_ind) const {
    Plasma_state curr_state = states[time_ind][spatial_ind];
    Plasma_state next_state = states[time_ind + 1][spatial_ind];
    double j = discharge_curr(time_ind);
    double v_e = curr_state.v_i - j / (phys::e * curr_state.n_i);
    return 5. / 3. * next_state.n_i * v_e;
}

double Solver::upwind_coefficient_plus(const std::size_t time_ind, const std::size_t spatial_ind) const {
    return (std::abs(flux_coefficient(time_ind, spatial_ind)) + flux_coefficient(time_ind, spatial_ind)) / (2 * std::abs(flux_coefficient(time_ind, spatial_ind)));
}

double Solver::upwind_coefficient_minus(const std::size_t time_ind, const std::size_t spatial_ind) const {
     return (std::abs(flux_coefficient(time_ind, spatial_ind)) - flux_coefficient(time_ind, spatial_ind)) / (2 * std::abs(flux_coefficient(time_ind, spatial_ind)));
}

std::pair<std::array<Eigen::VectorXd, 3>, Eigen::MatrixXd> Solver::Jacobian_matrix(const std::size_t time_ind, const Eigen::VectorXd T_e_vec) const {
    Eigen::VectorXd left_secondary_diagonal(sim::spatial_steps - 3);
    Eigen::VectorXd right_secondary_diagonal(sim::spatial_steps - 3);
    Eigen::VectorXd main_diagonal(sim::spatial_steps - 2);
    Eigen::MatrixXd Jacobian_matrix(sim::spatial_steps - 2, sim::spatial_steps - 2);
    double T_e_prev, T_e, T_e_next, flux_coefficient_plus, flux_coefficient_minus;
    std::pair<double, double> HFC_p, HFC_m;
    for (std::size_t spatial_ind = 1; spatial_ind < sim::spatial_steps - 1; spatial_ind++) {
        T_e = T_e_vec(spatial_ind - 1);
        if (spatial_ind == 1) {
            Plasma_state lbc = left_boundary_condition(time_ind);
            T_e_prev = lbc.T_e;
            T_e_next = T_e_vec(spatial_ind);
            flux_coefficient_minus = flux_coefficient(time_ind, 0);
            flux_coefficient_plus = flux_coefficient(time_ind, spatial_ind + 1);
            HFC_m =  heat_flux_coefficient_minus(time_ind, spatial_ind + 1, T_e, T_e);
            HFC_p = heat_flux_coefficient_plus(time_ind, spatial_ind, T_e, T_e_next);
        }
        else if (spatial_ind == sim::spatial_steps - 2) {
            Plasma_state rbc = right_boundary_condition(time_ind);
            T_e_next = rbc.T_e;
            T_e_prev = T_e_vec(spatial_ind - 2);
            flux_coefficient_plus = flux_coefficient(time_ind, sim::spatial_steps - 1);
            flux_coefficient_minus = flux_coefficient(time_ind, spatial_ind - 1);
            HFC_m =  heat_flux_coefficient_minus(time_ind, spatial_ind, T_e, T_e_prev);
            HFC_p = heat_flux_coefficient_plus(time_ind, spatial_ind, T_e, T_e_next);
        }
        else {
            T_e_prev = T_e_vec(spatial_ind - 2);
            T_e_next = T_e_vec(spatial_ind);
            flux_coefficient_minus = flux_coefficient(time_ind, spatial_ind - 1);
            flux_coefficient_plus = flux_coefficient(time_ind, spatial_ind + 1);
            HFC_m =  heat_flux_coefficient_minus(time_ind, spatial_ind, T_e, T_e_prev);
            HFC_p = heat_flux_coefficient_plus(time_ind, spatial_ind, T_e, T_e_next);
        }
        double HFC_minus = HFC_m.first;
        double derivative_HFC_minus = HFC_m.second;
        double HFC_plus = HFC_p.first;
        double derivative_HFC_plus = HFC_p.second;
        Plasma_state next_time_state = states[time_ind + 1][spatial_ind];
        if (spatial_ind != 1) {
            left_secondary_diagonal(spatial_ind - 2) = -derivative_HFC_minus * 0.5 * (T_e - T_e_prev) / (sim::spatial_step * sim::spatial_step) +  
            HFC_minus / (sim::spatial_step * sim::spatial_step) - 
            upwind_coefficient_plus(time_ind, spatial_ind) * flux_coefficient_minus / sim::spatial_step;
            Jacobian_matrix(spatial_ind - 1, spatial_ind - 2) = left_secondary_diagonal(spatial_ind - 2);
        }
        main_diagonal(spatial_ind - 1) = next_time_state.n_i / sim::time_step + 
            derivative_HFC_plus * 0.5 * (T_e_next - T_e) / (sim::spatial_step * sim::spatial_step) -
            HFC_plus / (sim::spatial_step * sim::spatial_step) -
            derivative_HFC_minus * 0.5 * (T_e - T_e_prev) / (sim::spatial_step * sim::spatial_step) -
            HFC_minus / (sim::spatial_step * sim::spatial_step) + 
            (upwind_coefficient_plus(time_ind, spatial_ind) - upwind_coefficient_minus(time_ind, spatial_ind)) * flux_coefficient(time_ind, spatial_ind) / sim::spatial_step;
        Jacobian_matrix(spatial_ind - 1, spatial_ind - 1) = main_diagonal(spatial_ind - 1);
        if (spatial_ind != sim::spatial_steps - 2) {
            right_secondary_diagonal(spatial_ind - 1) = derivative_HFC_plus * 0.5 * (T_e_next - T_e) / (sim::spatial_step * sim::spatial_step) -  
            HFC_plus / (sim::spatial_step * sim::spatial_step) + 
            upwind_coefficient_minus(time_ind, spatial_ind) * flux_coefficient_plus / sim::spatial_step;
            Jacobian_matrix(spatial_ind - 1, spatial_ind) = right_secondary_diagonal(spatial_ind - 1);
        }
    }
    std::array<Eigen::VectorXd, 3> diagonals = {left_secondary_diagonal,  main_diagonal, right_secondary_diagonal};
    std::pair<std::array<Eigen::VectorXd, 3>, Eigen::MatrixXd> result(diagonals, Jacobian_matrix);
    return result;
}

Eigen::VectorXd Solver::nonlin_equation_function(const std::size_t time_ind, const Eigen::VectorXd T_e_vec) const {
    Eigen::VectorXd f(sim::spatial_steps - 2);
    double T_e_prev, T_e, T_e_next, flux_coefficient_plus, flux_coefficient_minus;
    std::pair<double, double> HFC_p, HFC_m;
    double j = discharge_curr(time_ind);
    for (std::size_t spatial_ind = 1; spatial_ind < sim::spatial_steps - 1; spatial_ind++) {
        T_e = T_e_vec(spatial_ind - 1);
        if (spatial_ind == 1) {
            Plasma_state lbc = left_boundary_condition(time_ind);
            T_e_prev = lbc.T_e;
            T_e_next = T_e_vec(spatial_ind);
            flux_coefficient_minus = flux_coefficient(time_ind, 0);
            flux_coefficient_plus = flux_coefficient(time_ind, spatial_ind + 1);
            HFC_m =  heat_flux_coefficient_minus(time_ind, spatial_ind + 1, T_e, T_e);
            HFC_p = heat_flux_coefficient_plus(time_ind, spatial_ind, T_e, T_e_next);
        }
        else if (spatial_ind == sim::spatial_steps - 2) {
            Plasma_state rbc = right_boundary_condition(time_ind);
            T_e_next = rbc.T_e;
            T_e_prev = T_e_vec(spatial_ind - 2);
            flux_coefficient_plus = flux_coefficient(time_ind, sim::spatial_steps - 1);
            flux_coefficient_minus = flux_coefficient(time_ind, spatial_ind - 1);
            HFC_m =  heat_flux_coefficient_minus(time_ind, spatial_ind, T_e, T_e_prev);
            HFC_p = heat_flux_coefficient_plus(time_ind, spatial_ind, T_e, T_e_next);
        }
        else {
            T_e_prev = T_e_vec(spatial_ind - 2);
            T_e_next = T_e_vec(spatial_ind);
            flux_coefficient_minus = flux_coefficient(time_ind, spatial_ind - 1);
            flux_coefficient_plus = flux_coefficient(time_ind, spatial_ind + 1);
            HFC_m =  heat_flux_coefficient_minus(time_ind, spatial_ind, T_e, T_e_prev);
            HFC_p = heat_flux_coefficient_plus(time_ind, spatial_ind, T_e, T_e_next);
        }
        double HFC_minus = HFC_m.first;
        double HFC_plus = HFC_p.first;

        Plasma_state curr_state = states[time_ind][spatial_ind];
        Plasma_state next_time_state = states[time_ind + 1][spatial_ind];
        double nu_i = ioniz_coll_friq(time_ind, spatial_ind);
        double nu_iw = ions_wall_coll_friq(time_ind, spatial_ind);
        double v_e = curr_state.v_i - j / (phys::e * curr_state.n_i);
        double E =  electric_field(time_ind, spatial_ind);
        double delta_eps_w = energy_lose_rate(time_ind, spatial_ind);
        
        f(spatial_ind - 1) = (next_time_state.n_i * T_e - curr_state.n_i * curr_state.T_e) / sim::time_step + 
            HFC_plus * (T_e_next - T_e) / (sim::spatial_step * sim::spatial_step) -
            HFC_minus * (T_e - T_e_prev) / (sim::spatial_step * sim::spatial_step) + 
            upwind_coefficient_plus(time_ind, spatial_ind) * (flux_coefficient(time_ind, spatial_ind) * T_e - flux_coefficient_minus * T_e_prev) / sim::spatial_step +
            upwind_coefficient_minus(time_ind, spatial_ind) * (flux_coefficient_plus * T_e_next - flux_coefficient(time_ind, spatial_ind) * T_e) / sim::spatial_step - 
            2. / 3. * (-phys::e * curr_state.n_i * v_e * E - curr_state.n_i * (nu_i * colliz::Eps_i + delta_eps_w)); 
    }
    return f; 
}

std::vector<std::vector<Plasma_state>> Solver::solve(const double eps, const size_t max_iters) {
    Eigen::VectorXd vec_curr_ion_neutral_state(sim::Dim - 1);
    Eigen::VectorXd vec_prev_ion_neutral_state(sim::Dim - 1);
    Eigen::VectorXd vec_next_ion_neutral_state(sim::Dim - 1);
    Eigen::VectorXd vec_next_time_ion_neutral_state(sim::Dim);
    double steps_ratio = sim::time_step / sim::spatial_step;
    std::ofstream file("results.txt");
    file << "time_ind" << ',' << "spatial_ind" << ',' << "n_i" << ',' << "v_i" << ',' << "T_e" << ',' << "n_n" << ',' << "I" << std::endl;
    for (std::size_t time_ind = 0; time_ind < sim::time_steps - 1; time_ind++) {
        Plasma_state lbc = left_boundary_condition(time_ind);
        Plasma_state rbc = right_boundary_condition(time_ind);
        states[time_ind + 1][0].n_i = lbc.n_i;
        states[time_ind + 1][0].v_i = lbc.v_i;
        states[time_ind + 1][0].n_n = lbc.n_n;
        //#pragma omp parallel
        for (std::size_t spatial_ind = 1;  spatial_ind < sim::spatial_steps - 1; spatial_ind++) {
            Plasma_state curr_state = states[time_ind][spatial_ind];
            Plasma_state next_state = (spatial_ind == sim::spatial_steps - 2) ? rbc : states[time_ind][spatial_ind + 1];
            Plasma_state prev_state = (spatial_ind == 1) ? lbc : states[time_ind][spatial_ind - 1];
            vec_curr_ion_neutral_state << curr_state.n_i, curr_state.n_i * curr_state.v_i, curr_state.n_n;
            vec_next_ion_neutral_state << next_state.n_i, next_state.n_i * next_state.v_i, next_state.n_n;
            vec_prev_ion_neutral_state << prev_state.n_i, prev_state.n_i * prev_state.v_i, prev_state.n_n;
            vec_next_time_ion_neutral_state = vec_curr_ion_neutral_state - steps_ratio * (hyperbolic_numerical_flux(time_ind, curr_state, next_state) - hyperbolic_numerical_flux(time_ind, prev_state, curr_state)) + sim::time_step * hyperbolic_right_hand_side(time_ind, spatial_ind);        
            states[time_ind + 1][spatial_ind].n_i = vec_next_time_ion_neutral_state(0);
            states[time_ind + 1][spatial_ind].v_i = vec_next_time_ion_neutral_state(1) / vec_next_time_ion_neutral_state(0);
            states[time_ind + 1][spatial_ind].n_n = vec_next_time_ion_neutral_state(2);
        }
        states[time_ind + 1][sim::spatial_steps - 1].n_i = rbc.n_i;
        states[time_ind + 1][sim::spatial_steps - 1].v_i = rbc.v_i;
        states[time_ind + 1][sim::spatial_steps - 1].n_n = rbc.n_n;
        
        Eigen::VectorXd prev_T_e_layer(sim::spatial_steps - 2);
        Eigen::VectorXd next_T_e_layer(sim::spatial_steps - 2);
        Eigen::VectorXd column(sim::spatial_steps - 2);
        for (std::size_t spatial_ind = 1; spatial_ind < sim::spatial_steps - 1; spatial_ind++) {
            prev_T_e_layer(spatial_ind - 1) = states[time_ind][spatial_ind].T_e;
        }
        std::size_t i;
        double error;
        for (i = 0; i < max_iters; i++) {
            auto [diagonals, Jacobian] = Jacobian_matrix(time_ind, prev_T_e_layer);
            column = Jacobian * prev_T_e_layer - nonlin_equation_function(time_ind, prev_T_e_layer);
            next_T_e_layer = tridiagonal_matrix_algorithm(diagonals[0], diagonals[1], diagonals[2], column);
            error = ((next_T_e_layer - prev_T_e_layer) / phys::e).cwiseAbs().maxCoeff();
            if (error < eps) {
                for (std::size_t spatial_ind = 1; spatial_ind < sim::spatial_steps - 1; spatial_ind++) {
                    states[time_ind + 1][spatial_ind].T_e = next_T_e_layer(spatial_ind - 1);
                }
                Plasma_state lbc = left_boundary_condition(time_ind);
                Plasma_state rbc = right_boundary_condition(time_ind);
                states[time_ind + 1][0].T_e = lbc.T_e;
                states[time_ind + 1][sim::spatial_steps - 1].T_e = rbc.T_e;
                break;
            }
            next_T_e_layer.swap(prev_T_e_layer);
        }
        if (i == max_iters) std::cout << "limit iterations!" << std::endl;

        if (time_ind % 10 == 0) {
            for (std::size_t spatial_ind = 0; spatial_ind < sim::spatial_steps; spatial_ind++) {
                file << time_ind << ',' << spatial_ind << ',' << states[time_ind + 1][spatial_ind].n_i << ',' << states[time_ind + 1][spatial_ind].v_i << ',' << states[time_ind + 1][spatial_ind].T_e / phys::e << ',' << states[time_ind + 1][spatial_ind].n_n << ',' << discharge_curr(time_ind) * geom::channel_area << std::endl;
            }
        }
    }
    return states;
}

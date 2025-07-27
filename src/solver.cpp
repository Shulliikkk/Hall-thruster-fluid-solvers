#include "solver.hpp"

Solver::Solver(double discharge_voltage, double v_n, double input_mass_flow, std::string dataset_ioniz_name, std::string dataset_neutral_name): Walls(discharge_voltage, v_n, input_mass_flow, dataset_ioniz_name, dataset_neutral_name) {}

Electrons_state Solver::calc_electrons_state(const size_t time_step, const size_t space_step) const {
    Electrons_state electrons_state;
    Plasma_state curr_state = states[time_step][space_step];
    double j = discharge_curr(time_step, space_step);
    electrons_state.v_e = curr_state.v_i - j / (phys::e * curr_state.n_i);
    electrons_state.c_s = std::sqrt(5 * curr_state.T_e / (3 * phys::m_e));
    electrons_state.u = curr_state.v_i - electrons_state.v_e;
    return electrons_state;
}

Eigen::MatrixX<double> Solver::Omega_L(const size_t time_step, const size_t space_step) const {
    Electrons_state electrons_state = calc_electrons_state(time_step, space_step);
    Plasma_state curr_state = states[time_step][space_step];
    Eigen::MatrixX<double> Omega_L(grid::Dim, grid::Dim);
    Omega_L << - 2. / 3. * curr_state.T_e, 0, curr_state.n_i, 0,
                curr_state.T_e - phys::m_i * electrons_state.c_s * electrons_state.u, phys::m_i * curr_state.n_i * (electrons_state.u - electrons_state.c_s), curr_state.n_i, 0,
                curr_state.T_e + phys::m_i * electrons_state.c_s * electrons_state.u, phys::m_i * curr_state.n_i * (electrons_state.u + electrons_state.c_s), curr_state.n_i, 0,
                0, 0, 0, 1;
    return Omega_L;
}

Eigen::MatrixX<double> Solver::Omega_R(const size_t time_step, const size_t space_step) const {
    Electrons_state electrons_state = calc_electrons_state(time_step, space_step);
    Plasma_state curr_state = states[time_step][space_step];
    double frac = 1 / (5 * curr_state.T_e + 3 * phys::m_i * electrons_state.u * electrons_state.u);
    Eigen::MatrixX<double> Omega_R(grid::Dim, grid::Dim);
    Omega_R << - 3 * frac, 3 * (electrons_state.c_s + electrons_state.u) * frac / (2 * electrons_state.c_s), 3 * (electrons_state.c_s + electrons_state.u) * frac / (2 * electrons_state.c_s), 0,
                3 * electrons_state.u * frac / curr_state.n_i, (5 * curr_state.T_e + 3 * electrons_state.c_s * phys::m_i * electrons_state.u) * (-frac) / (2 * electrons_state.c_s * phys::m_i * curr_state.n_i), (5 * curr_state.T_e - 3 * electrons_state.c_s * phys::m_i * electrons_state.u) * frac / (2 * electrons_state.c_s * phys::m_i * curr_state.n_i), 0,
                3 * (curr_state.T_e - phys::m_i * electrons_state.u * electrons_state.u) * frac / curr_state.n_i, curr_state.T_e * (electrons_state.c_s + electrons_state.u) * frac / (electrons_state.c_s * curr_state.n_i), curr_state.T_e * (electrons_state.c_s - electrons_state.u) * frac / (electrons_state.c_s * curr_state.n_i), 0,
                0, 0, 0, 1;
    return Omega_R;
}

Eigen::MatrixX<double> Solver::Lamda(const size_t time_step, const size_t space_step) const {
    Electrons_state electrons_state = calc_electrons_state(time_step, space_step);
    Plasma_state curr_state = states[time_step][space_step];
    Eigen::MatrixX<double> Lamda(grid::Dim, grid::Dim);
    Lamda << electrons_state.v_e, 0, 0, 0,
                0, curr_state.v_i - electrons_state.c_s, 0, 0,
                0, 0, curr_state.v_i + electrons_state.c_s, 0,
                0, 0, 0, 1;
    return Lamda;
}

Eigen::VectorX<double> Solver::flux(const size_t time_step, const size_t space_step) const {
    Plasma_state curr_state = states[time_step][space_step];
    double j = discharge_curr(time_step, space_step);
    double v_e = curr_state.v_i - j / (phys::e * curr_state.n_i);
    Eigen::VectorX<double> flux(grid::Dim);
    flux << curr_state.n_i * curr_state.v_i,
            phys::m_i * curr_state.n_i * curr_state.v_i * curr_state.v_i,
            5. / 2. * curr_state.T_e * curr_state.n_i * v_e,
            curr_state.n_n * v_n;
    return flux;
}

Eigen::VectorX<double> Solver::right_hand_side(const size_t time_step, const size_t space_step) const {
    Plasma_state curr_state = states[time_step][space_step];
    double j = discharge_curr(time_step, space_step);
    double v_e = curr_state.v_i - j / (phys::e * curr_state.n_i);
    double delta = SEE_coef(time_step, space_step);
    double nu_iw = wall_coll_friq(time_step, space_step).first;
    double nu_ew = wall_coll_friq(time_step, space_step).second;
    double delta_eps_w = energy_lose_rate(time_step, space_step);
    double nu_i = ioniz_coll_friq(time_step, space_step);
    double E =  electric_field(time_step, space_step);
    Eigen::VectorX<double> rhs(grid::Dim);
    rhs << curr_state.n_i * (nu_i - nu_iw),
            - phys::e * curr_state.n_i * E + phys::m_i * curr_state.n_i * (nu_i * v_n - nu_iw * curr_state.v_i),
            curr_state.n_i * (phys::e * v_e * E - nu_i * ioniz::alpha_i * ioniz::Eps_i - delta_eps_w),
            - curr_state.n_i * (nu_i - nu_iw);
    return rhs;
}

// Holodov A. S. 1978
Eigen::MatrixX<double> Solver::matrix_halph_plus(const size_t time_step, const size_t space_step) const {
    return 0.5 * (Omega_R(time_step, space_step) * Lamda(time_step, space_step) * Omega_L(time_step, space_step) +
                    Omega_R(time_step, space_step + 1) * Lamda(time_step, space_step + 1) * Omega_L(time_step, space_step + 1));
}

// Holodov A. S. 1978
Eigen::MatrixX<double> Solver::matrix_halph_minus(const size_t time_step, const size_t space_step) const {
    return 0.5 * (Omega_R(time_step, space_step) * Lamda(time_step, space_step) * Omega_L(time_step, space_step) +
                    Omega_R(time_step, space_step - 1) * Lamda(time_step, space_step - 1) * Omega_L(time_step, space_step - 1));
}

Eigen::VectorX<double> Solver::num_flux_halph_plus(const size_t time_step, const size_t space_step) const {
    Plasma_state curr_state = states[time_step][space_step];
    Plasma_state next_state;
    Eigen::VectorX<double> vec_curr_state(grid::Dim);
    Eigen::VectorX<double> vec_next_state(grid::Dim);
    vec_curr_state << curr_state.n_i,  curr_state.v_i,  curr_state.T_e,  curr_state.n_n;
    vec_next_state << next_state.n_i,  next_state.v_i,  next_state.T_e,  next_state.n_n;
    if (space_step == grid::M - 2) {
        next_state.n_i = curr_state.n_i;
        next_state.v_i = curr_state.v_i;
        next_state.T_e = 2 / (phys::k_B * phys::e);
        next_state.n_n = curr_state.n_n;
    }
    else {
        next_state = states[time_step][space_step + 1];
    }
    return 0.5 * (flux(time_step, space_step + 1) + flux(time_step, space_step)) - 0.5 * matrix_halph_plus(time_step, space_step) * (vec_next_state - vec_curr_state);
}

Eigen::VectorX<double> Solver::num_flux_halph_minus(const size_t time_step, const size_t space_step) const {
    Plasma_state curr_state = states[time_step][space_step];
    Plasma_state prev_state;
    Eigen::VectorX<double> vec_curr_state(grid::Dim);
    Eigen::VectorX<double> vec_prev_state(grid::Dim);
    vec_curr_state << curr_state.n_i,  curr_state.v_i,  curr_state.T_e,  curr_state.n_n;
    vec_prev_state << prev_state.n_i,  prev_state.v_i,  prev_state.T_e,  prev_state.n_n;
    if (space_step == 1) {
        prev_state.n_i = curr_state.n_i;
        prev_state.T_e = 2 * phys::e / phys::k_B; // 2 eV
        prev_state.v_i = std::sqrt(phys::k_B * prev_state.T_e); // v_Bohm
        prev_state.n_n = input_mass_flow / (phys::m_i * geom::A * v_n) - prev_state.n_i *  prev_state.v_i / v_n;
    }
    else {
        prev_state = states[time_step][space_step - 1];
    }
    return 0.5 * (flux(time_step, space_step) + flux(time_step, space_step - 1)) - 0.5 * matrix_halph_minus(time_step, space_step) * (vec_curr_state - vec_prev_state);
}

std::array<std::array<Plasma_state, grid::M>, grid::N> Solver::solve() {
    for (size_t time_step = 0; time_step < grid::N - 1; time_step++) {
        for (size_t space_step = 1; space_step < grid::M - 1; space_step++) {
            Eigen::VectorX<double> plus = num_flux_halph_plus(time_step, space_step);
            Eigen::VectorX<double> minus = num_flux_halph_minus(time_step, space_step);
            Eigen::VectorX<double> rhs = right_hand_side(time_step, space_step);
            states[time_step + 1][space_step].n_i = states[time_step][space_step].n_i - grid::tau / grid::h * (plus(0) - minus(0)) + grid::tau * rhs(0);
            states[time_step + 1][space_step].v_i = states[time_step][space_step].v_i - grid::tau / grid::h * (plus(1) - minus(1)) + grid::tau * rhs(1);
            states[time_step + 1][space_step].T_e = states[time_step][space_step].T_e - grid::tau / grid::h * (plus(2) - minus(2)) + grid::tau * rhs(2);
            states[time_step + 1][space_step].n_n = states[time_step][space_step].n_n - grid::tau / grid::h * (plus(3) - minus(3)) + grid::tau * rhs(3);
        }
    }
    return states;
}
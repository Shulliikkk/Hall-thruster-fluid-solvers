#include "solve.hpp"
#include <iostream>

Solver::Solver(): Discharge() {}

void Solver::set_init_conditions(double n_i_min, double n_i_max, double T_e_anode_eV, double T_e_cathode_eV) {
    double T_e_max = disch::discharge_voltage * phys::e / 50;
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
    //std::cout << discharge_curr(0) * geom::channel_area << std::endl;
    /*for (std::size_t spatial_ind = 0; spatial_ind < sim::spatial_steps; spatial_ind++) {
        //std::cout << states[0][spatial_ind].n_i << ' ' << states[0][spatial_ind].v_i << ' ' << states[0][spatial_ind].T_e << ' ' << states[0][spatial_ind].n_n << std::endl;
        //std::cout << SEE_coef(0, spatial_ind) << ' ' << ions_wall_coll_friq(0, spatial_ind) << ' ' << electrons_wall_coll_friq(0, spatial_ind) << std::endl;
        //std::cout << magfld::magnetic_field(spatial_ind) << ' ' << phys::e / (moment_exch_coll_friq(0, spatial_ind) * phys::m_e) << ' ' << states[0][spatial_ind].n_i << ' ' << cicl_friq(spatial_ind) / moment_exch_coll_friq(0, spatial_ind) << ' ' << electron_cross_field_mobility(0, spatial_ind) << std::endl;
    }*/
}

Plasma_state Solver::right_boundary_condition(std::size_t time_ind) {
    Plasma_state rbc;
    Plasma_state prev_state = states[time_ind][sim::spatial_steps - 2];
    double T_e_max = disch::discharge_voltage * phys::e / 50;
    rbc.n_i = prev_state.n_i;
    rbc.v_i = prev_state.v_i;
    rbc.T_e = (1 - (sim::spatial_steps - 1) * sim::spatial_step / sim::simulation_area_length) * T_e_anode_ + (sim::spatial_steps - 1) * sim::spatial_step / sim::simulation_area_length * T_e_cathode_ +
        T_e_max * std::exp(-((sim::spatial_steps - 1) * sim::spatial_step - geom::channel_length) / (geom::channel_length / 3) *
        ((sim::spatial_steps - 1) * sim::spatial_step - geom::channel_length) / (geom::channel_length / 3));
    rbc.n_n = prev_state.n_n;
    return rbc;
}

Plasma_state Solver::left_boundary_condition(std::size_t time_ind) {
    Plasma_state lbc;
    double T_e_max = disch::discharge_voltage * phys::e / 50;
    double v_bohm_anode = -std::sqrt(T_e_anode_ / phys::m_i);
    lbc.n_i = std::sqrt(disch::discharge_voltage / 300) * disch::mass_flow / 5e-6 * n_i_max_ * std::exp(-9. / 4.);
    lbc.T_e = T_e_anode_ + T_e_max * std::exp(-9);
    lbc.v_i = v_bohm_anode;
    double n_n_anode =  disch::mass_flow / (phys::m_i * geom::channel_area * disch::v_n) - lbc.n_i * lbc.v_i / disch::v_n;
    double n_n_cathode = 0.01 * n_n_anode;
    lbc.n_n = 0.5 * (n_n_anode + n_n_cathode + (n_n_cathode - n_n_anode) * std::tanh(-12));
    return lbc;
}

Electrons_state Solver::calc_electrons_state(const std::size_t time_ind, const std::size_t spatial_ind) const {
    Electrons_state electrons_state;
    Plasma_state curr_state = states[time_ind][spatial_ind];
    double j = discharge_curr(time_ind);
    //if (time_ind > 110) std::cout << 'j' <<  time_ind << ' ' << j << std::endl;
    electrons_state.v_e = curr_state.v_i - j / (phys::e * curr_state.n_i);
    electrons_state.c_s = std::sqrt(5 * curr_state.T_e / (3 * phys::m_i));
    //if (time_ind == 122 && spatial_ind == 36) std::cout << "disch" << ' ' << time_ind << ' ' << spatial_ind <<  ' ' << electrons_state.c_s << std::endl;
    electrons_state.u = curr_state.v_i - electrons_state.v_e;
    return electrons_state;
}

Eigen::MatrixXd Solver::Omega_L(const std::size_t time_ind, const std::size_t spatial_ind) const {
    Electrons_state electrons_state = calc_electrons_state(time_ind, spatial_ind);
    Plasma_state curr_state = states[time_ind][spatial_ind];
    Eigen::MatrixX<double> Omega_L(sim::Dim, sim::Dim);
    Omega_L << - 2. / 3. * curr_state.T_e, 0, curr_state.n_i, 0,
                curr_state.T_e - phys::m_i * electrons_state.c_s * electrons_state.u, phys::m_i * curr_state.n_i * (electrons_state.u - electrons_state.c_s), curr_state.n_i, 0,
                curr_state.T_e + phys::m_i * electrons_state.c_s * electrons_state.u, phys::m_i * curr_state.n_i * (electrons_state.u + electrons_state.c_s), curr_state.n_i, 0,
                0, 0, 0, 1;
    //if (time_ind > 120) std::cout << "Omega_L" << ' ' << time_ind << ' ' << spatial_ind << ' ' << electrons_state.u << ' ' << std::endl << Omega_L << std::endl;
    return Omega_L;
}

Eigen::MatrixXd Solver::Omega_R(const std::size_t time_ind, const std::size_t spatial_ind) const {
    Electrons_state electrons_state = calc_electrons_state(time_ind, spatial_ind);
    Plasma_state curr_state = states[time_ind][spatial_ind];
    double frac = 1. / (5 * curr_state.T_e - 3 * phys::m_i * electrons_state.u * electrons_state.u);
    Eigen::MatrixXd Omega_R(sim::Dim, sim::Dim);
    Omega_R << - 3 * frac, 3 * (electrons_state.c_s + electrons_state.u) * frac / (2 * electrons_state.c_s), 3 * (electrons_state.c_s - electrons_state.u) * frac / (2 * electrons_state.c_s), 0,
                3 * electrons_state.u * frac / curr_state.n_i, (5 * curr_state.T_e + 3 * electrons_state.c_s * phys::m_i * electrons_state.u) * (-frac) / (2 * electrons_state.c_s * phys::m_i * curr_state.n_i), (5  * curr_state.T_e - 3 * electrons_state.c_s * phys::m_i * electrons_state.u) * frac / (2 * electrons_state.c_s * phys::m_i * curr_state.n_i), 0,
                3 * (curr_state.T_e - phys::m_i * electrons_state.u * electrons_state.u) * frac / curr_state.n_i, curr_state.T_e * (electrons_state.c_s + electrons_state.u) * frac / (electrons_state.c_s * curr_state.n_i), curr_state.T_e * (electrons_state.c_s - electrons_state.u) * frac / (electrons_state.c_s * curr_state.n_i), 0,
                0, 0, 0, 1;
    //if (time_ind > 120) std::cout << "Omega_R" << ' ' << time_ind << ' ' << spatial_ind << ' ' << electrons_state.u << ' ' << std::endl << Omega_R << std::endl;
    return Omega_R;
}

Eigen::MatrixXd Solver::Lamda_plus(const std::size_t time_ind, const std::size_t spatial_ind) const {
    Electrons_state electrons_state = calc_electrons_state(time_ind, spatial_ind);
    Plasma_state curr_state = states[time_ind][spatial_ind];
    Eigen::MatrixXd Lamda_plus(sim::Dim, sim::Dim);
    Lamda_plus << (std::abs(electrons_state.v_e) + electrons_state.v_e) / 2, 0, 0, 0,
                0, (std::abs(curr_state.v_i - electrons_state.c_s) + (curr_state.v_i - electrons_state.c_s)) / 2, 0, 0,
                0, 0, (std::abs(curr_state.v_i + electrons_state.c_s) + (curr_state.v_i + electrons_state.c_s)) / 2, 0,
                0, 0, 0, (std::abs(disch::v_n) + disch::v_n) / 2;
    //std::cout << "Lamda_plus" << ' ' << time_ind << ' ' << spatial_ind << std::endl << Lamda_plus << std::endl;
    return Lamda_plus;
}

Eigen::MatrixXd Solver::Lamda_minus(const std::size_t time_ind, const std::size_t spatial_ind) const {
    Electrons_state electrons_state = calc_electrons_state(time_ind, spatial_ind);
    Plasma_state curr_state = states[time_ind][spatial_ind];
    Eigen::MatrixXd Lamda_minus(sim::Dim, sim::Dim);
    Lamda_minus << (electrons_state.v_e - std::abs(electrons_state.v_e)) / 2, 0, 0, 0,
                0, ((curr_state.v_i - electrons_state.c_s) - std::abs(curr_state.v_i - electrons_state.c_s)) / 2, 0, 0,
                0, 0, ((curr_state.v_i + electrons_state.c_s) - std::abs(curr_state.v_i + electrons_state.c_s)) / 2, 0,
                0, 0, 0, (disch::v_n - std::abs(disch::v_n)) / 2;
    //std::cout << "Lamda_minus" << ' ' << time_ind << ' ' << spatial_ind << std::endl << Lamda_minus << std::endl;
    return Lamda_minus;
}

Eigen::VectorXd Solver::right_hand_side(const std::size_t time_ind, const std::size_t spatial_ind) const {
    Plasma_state curr_state = states[time_ind][spatial_ind];
    double omega_c = cicl_friq(spatial_ind);
    double nu_iw = ions_wall_coll_friq(time_ind, spatial_ind);
    double nu_ew = electrons_wall_coll_friq(time_ind, spatial_ind);
    double nu_i = ioniz_coll_friq(time_ind, spatial_ind);
    //std::cout << "nui" << time_ind << ' ' << spatial_ind << ' ' << curr_state.T_e / phys::e << ' ' << nu_i << std::endl;
    double nu_m = moment_exch_coll_friq(time_ind, spatial_ind);
    double nu_d = phys::m_e * omega_c * omega_c / (phys::m_i * nu_m);
    double j = discharge_curr(time_ind);
    //double E =  electric_field(time_ind, spatial_ind);
    double v_e = curr_state.v_i - j / (phys::e * curr_state.n_i);
    double delta_eps_w = energy_lose_rate(time_ind, spatial_ind);
    //std::cout << omega_c << ' ' << nu_iw << ' ' << nu_ew << ' ' << nu_i << ' ' << nu_m << ' ' << nu_d << ' ' << j << ' ' << E << ' ' << v_e << std::endl;
    Eigen::VectorXd rhs(sim::Dim);
    rhs << curr_state.n_i * (nu_i - nu_iw),
            nu_i * (disch::v_n - curr_state.v_i) - nu_iw * curr_state.v_i - nu_d * v_e,
            2. / 3. * nu_d * phys::m_i * v_e * v_e - nu_i * (curr_state.T_e + 2. / 3. * colliz::alpha_i * phys::e * colliz::Eps_i) - 2. / 3. * delta_eps_w,
            - curr_state.n_i * (nu_i - nu_iw);
    //std::cout << "rhs" << ' ' << time_ind << ' ' << spatial_ind << std::endl << nu_d << ' ' << v_e << ' ' << nu_i << ' ' << delta_eps_w << ' ' << SEE_coef(time_ind, spatial_ind) << std::endl;
    return rhs;
}

Eigen::MatrixXd Solver::A_plus(const std::size_t time_ind, const std::size_t spatial_ind) const {
    //std::cout << "A_plus" << ' ' << time_ind << ' ' << spatial_ind << std::endl << Omega_R(time_ind, spatial_ind) * Lamda_plus(time_ind, spatial_ind) * Omega_L(time_ind, spatial_ind) << std::endl;
    return Omega_R(time_ind, spatial_ind) * Lamda_plus(time_ind, spatial_ind) * Omega_L(time_ind, spatial_ind);
}

Eigen::MatrixXd Solver::A_minus(const std::size_t time_ind, const std::size_t spatial_ind) const {
    //std::cout << "A_minus" << ' ' << time_ind << ' ' << spatial_ind << std::endl << Omega_R(time_ind, spatial_ind) * Lamda_minus(time_ind, spatial_ind) * Omega_L(time_ind, spatial_ind) << std::endl;
    return Omega_R(time_ind, spatial_ind) * Lamda_minus(time_ind, spatial_ind) * Omega_L(time_ind, spatial_ind);
}

std::vector<std::vector<Plasma_state>> Solver::solve() {
    Eigen::VectorXd vec_curr_state(sim::Dim);
    Eigen::VectorXd vec_prev_state(sim::Dim);
    Eigen::VectorXd vec_next_state(sim::Dim);
    Eigen::VectorXd vec_next_time_state(sim::Dim);
    double steps_ratio = sim::time_step / sim::spatial_step;
    std::ofstream file("results.txt");
    file << "time_ind" << ',' << "spatial_ind" << ',' << "n_i" << ',' << "v_i" << ',' << "T_e" << ',' << "n_n" << "I" << std::endl;
    for (std::size_t time_ind = 0; time_ind < sim::time_steps - 1; time_ind++) {
        Plasma_state lbc = left_boundary_condition(time_ind);
        states[time_ind + 1][0].n_i = lbc.n_i;
        states[time_ind + 1][0].v_i = lbc.v_i;
        states[time_ind + 1][0].T_e = lbc.T_e;
        states[time_ind + 1][0].n_n = lbc.n_n;
        for (std::size_t spatial_ind = 1;  spatial_ind < sim::spatial_steps - 1; spatial_ind++) {
            Plasma_state curr_state = states[time_ind][spatial_ind];
            Plasma_state next_state = (spatial_ind == sim::spatial_steps - 2) ? right_boundary_condition(time_ind) : states[time_ind][spatial_ind + 1];
            Plasma_state prev_state = (spatial_ind == 1) ? left_boundary_condition(time_ind) : states[time_ind][spatial_ind - 1];
            vec_curr_state << curr_state.n_i, curr_state.v_i, curr_state.T_e, curr_state.n_n;
            vec_next_state << next_state.n_i, next_state.v_i, next_state.T_e, next_state.n_n;
            vec_prev_state << prev_state.n_i, prev_state.v_i, prev_state.T_e, prev_state.n_n;
            /*if (time_ind == 16) {
                //std::cout << "lol";
                std::cout << Omega_L(time_ind, spatial_ind) << std::endl;
                std::cout << Omega_R(time_ind, spatial_ind) << std::endl;
                std::cout << Lamda_plus(time_ind, spatial_ind) << std::endl;
                std::cout << Lamda_minus(time_ind, spatial_ind) << std::endl;
                //std::cout << A_minus(time_ind, spatial_ind) << std::endl;
            }*/
            vec_next_time_state = vec_curr_state - steps_ratio * A_plus(time_ind, spatial_ind) * (vec_curr_state - vec_prev_state) - steps_ratio * A_minus(time_ind, spatial_ind) * (vec_next_state - vec_curr_state) + sim::time_step * right_hand_side(time_ind, spatial_ind);        
            states[time_ind + 1][spatial_ind].n_i = vec_next_time_state(0);
            states[time_ind + 1][spatial_ind].v_i = vec_next_time_state(1);
            states[time_ind + 1][spatial_ind].T_e = vec_next_time_state(2);
            states[time_ind + 1][spatial_ind].n_n = vec_next_time_state(3);
            //Electrons_state electrons_state = calc_electrons_state(time_ind, spatial_ind);
            //Plasma_state curr_st = states[time_ind][spatial_ind];
            //std::array<double, 4> Co = {std::abs(electrons_state.v_e), std::abs(curr_st.v_i - electrons_state.c_s), std::abs(curr_st.v_i + electrons_state.c_s), disch::v_n};
            //std::cout << sim::time_step / sim::spatial_step * *std::max_element(Co.begin(), Co.end()) << std::endl;
            //std::cout << states[time_ind + 1][spatial_ind].n_i << ' ' << states[time_ind + 1][spatial_ind].v_i << ' ' << states[time_ind + 1][spatial_ind].T_e << ' ' << states[time_ind + 1][spatial_ind].n_n << std::endl;
        }
        Plasma_state rbc = right_boundary_condition(time_ind);
        states[time_ind + 1][sim::spatial_steps - 1].n_i = rbc.n_i;
        states[time_ind + 1][sim::spatial_steps - 1].v_i = rbc.v_i;
        states[time_ind + 1][sim::spatial_steps - 1].T_e = rbc.T_e;
        states[time_ind + 1][sim::spatial_steps - 1].n_n = rbc.n_n;
        for (std::size_t spatial_ind = 0; spatial_ind < sim::spatial_steps; spatial_ind++) {
            file << time_ind << ',' << spatial_ind << ',' << states[time_ind + 1][spatial_ind].n_i << ',' << states[time_ind + 1][spatial_ind].v_i << ',' << states[time_ind + 1][spatial_ind].T_e / phys::e << ',' << states[time_ind + 1][spatial_ind].n_n << ',' << discharge_curr(time_ind) * geom::channel_area << std::endl;
        }
    }
    return states;
}
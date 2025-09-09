#include "params.hpp"

namespace geom {
    double channel_length = 0.0, channel_radius_in = 0.0, channel_radius_out = 0.0, channel_area = 0.0;
}
void geom::set_geom(double chann_len, double chann_rad_in, double chann_rad_out) {
    geom::channel_length = chann_len;
    geom::channel_radius_in = chann_rad_in;
    geom::channel_radius_out = chann_rad_out;
    geom::channel_area = M_PI * (geom::channel_radius_out * geom::channel_radius_out - geom::channel_radius_in * geom::channel_radius_in);
}

namespace wall {
    double a = 0.0, b = 0.0, gamma_func = 0.0;
}
void wall::set_wall_material(std::string material_name) try {
    if (material_name == "BNSiO2") {
        wall::a = 0.123;
        wall::b = 0.528;
        wall::gamma_func = 1.36;
    }
    else {
        throw std::invalid_argument("Wrong type wall material");
    }
}
catch(std::invalid_argument& e) {
    std::cout << "Some exception of type ExceptionType: " << e.what() << std::endl;
}

namespace magfld {
    double B_0 = 0.0, delta_B_in = 0.0, delta_B_out = 0.0;
    std::vector<double> magnetic_field_distribution(0);
}
void magfld::set_magnetic_field(std::string case_name) try {
    magnetic_field_distribution.resize(sim::spatial_steps);
    if (case_name == "default gaussian-distributed") {
        magfld::B_0 = 150 * 1e-4; 
        magfld::delta_B_in = 1.1e-2;
        magfld::delta_B_out = 1.8e-2;
        double delta_B;
        for (std::size_t spatial_ind = 0; spatial_ind < sim::spatial_steps; spatial_ind++) {
            delta_B = (spatial_ind * sim::spatial_step <= geom::channel_length) ? magfld::delta_B_in : magfld::delta_B_out;
            magfld::magnetic_field_distribution[spatial_ind] = magfld::B_0 * std::exp(-(spatial_ind * sim::spatial_step - geom::channel_length) * (spatial_ind * sim::spatial_step - geom::channel_length) / (2 * delta_B * delta_B));
        }
    }
    else if (case_name == "tabulated-distributed") {
        // ДОПИСАТЬ
    }
    else {
        throw std::invalid_argument("Wrong type magnetic field");
    }
}
catch(std::invalid_argument& e) {
    std::cout << "Some exception of type ExceptionType: " << e.what() << std::endl;
}

void magfld::set_magnetic_field(std::string case_name, double amplitude_mag_field, double delta_mag_field_in, double delta_mag_field_out) try {
    magnetic_field_distribution.resize(sim::spatial_steps);
    if (case_name == "gaussian-distributed") {
        magfld::B_0 = amplitude_mag_field;
        magfld::delta_B_in = delta_mag_field_in;
        magfld::delta_B_out = delta_mag_field_out;
        double delta_B;
        for (std::size_t spatial_ind = 0; spatial_ind < sim::spatial_steps; spatial_ind++) {
            delta_B = (spatial_ind * sim::spatial_step <= geom::channel_length) ? magfld::delta_B_in : magfld::delta_B_out;
            magfld::magnetic_field_distribution[spatial_ind] = magfld::B_0 * std::exp(-(spatial_ind * sim::spatial_step - geom::channel_length) * (spatial_ind * sim::spatial_step - geom::channel_length) / (2 * delta_B * delta_B));
        }
    }
    else {
        throw std::invalid_argument("Wrong type magnetic field");
    }
}
catch(std::invalid_argument& e) {
    std::cout << "Some exception of type ExceptionType: " << e.what() << std::endl;
}

double magfld::magnetic_field(std::size_t spatial_ind) {
    return magfld::magnetic_field_distribution[spatial_ind];
}

namespace disch {
    double discharge_voltage = 0.0, mass_flow = 0.0, v_n = 0.0;
}
void disch::set_disch_params(double disch_volt, double vn, double ms_fl) {
    disch::discharge_voltage = disch_volt;
    disch::v_n = vn;
    disch::mass_flow = ms_fl;
}

namespace sim {
    double simulation_area_length = 0.0, simulation_time = 0.0, spatial_step = 0.0, time_step = 0.0;
    std::size_t time_steps = 0, spatial_steps = 0, simulation_area_steps = 0;
}
void sim::set_sim_params(double sim_area_len, double sim_time, double sp_step, double t_step) {
    sim::simulation_area_length = sim_area_len;
    sim::simulation_time = sim_time;
    sim::spatial_step = sp_step;
    sim::time_step = t_step;
    sim::time_steps = static_cast<std::size_t>(simulation_time / time_step);
    sim::spatial_steps = static_cast<std::size_t>(sim::simulation_area_length / spatial_step);
    sim::simulation_area_steps =  static_cast<std::size_t>(geom::channel_length / sim::spatial_step);
}
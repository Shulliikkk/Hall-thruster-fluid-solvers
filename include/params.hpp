// params
#pragma one
#include <cmath>
#include <stdexcept>
#include <cstddef>
#include <vector>
#include <string>
#include <iostream>
#include "rapidcsv.h"

namespace phys {
    constexpr double k_B = 1.380649e-23;
    constexpr double e = 1.602176634e-19;
    constexpr double m_e = 9.10938356e-31;
    constexpr double m_i = 131.293 * 1.6605402e-27;
}

namespace geom {
    extern double channel_length, channel_radius_in, channel_radius_out, channel_area;
    void set_geom(double chann_len, double chann_rad_in, double chann_rad_out);
}

namespace wall {
    extern double a, b, gamma_func;
    void set_wall_material(std::string material_name);
}

namespace magfld {
    extern double B_0, delta_B_in, delta_B_out;
    extern std::vector<double> magnetic_field_distribution;
    void set_magnetic_field(std::string case_name);
    void set_magnetic_field(std::string case_name, double amplitude_mag_field, double delta_mag_field_in, double delta_mag_field_out);
    double magnetic_field(std::size_t spatial_ind);
}

namespace colliz {
    constexpr double alpha_i = 1;
    constexpr double Eps_i = 12.13 * phys::e;
    void load_collisions_data(std::string dataset_ioniz_name, std::string dataset_neutral_name);
}

namespace disch {
    extern double discharge_voltage, mass_flow, v_n;
    void set_disch_params(double disch_volt, double vn, double ms_fl);
}

namespace sim {
    extern double simulation_area_length, simulation_time, spatial_step, time_step;
    extern std::size_t time_steps, spatial_steps, simulation_area_steps;
    constexpr std::size_t Dim = 4;
    void set_sim_params(double sim_area_len, double sim_time, double sp_step, double t_step);
}

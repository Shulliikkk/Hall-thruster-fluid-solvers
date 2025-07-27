#pragma one
#include <cmath>

namespace phys {
    constexpr double k_B = 1.380649e-23;
    constexpr double e = 1.602176634e-19;
    constexpr double m_e = 9.10938356e-31;
    constexpr double m_i = 131.293 * 1.6605402e-27;
}

namespace grid {
    constexpr double L = 2.5e-2, T = 10e-6; // sm, s
    constexpr double h = 1e-3, tau = 1e-6;
    constexpr std::size_t M = static_cast<std::size_t>(L / h);
    constexpr std::size_t N = static_cast<std::size_t>(T / tau);
    constexpr std::size_t Dim = 4;
}

namespace geom {
    constexpr double R_in = 3.5e-2, R_out = 5e-2; // sm
    constexpr double A = M_PI * (R_out * R_out - R_in * R_in); // sm
}

namespace BNSiO2 {
    constexpr double a = 0.123, b = 0.528;
    constexpr double gamma_func = 1.36;
}

namespace magfield {
    constexpr double B_0 = 150 * 1e-4; // Gauss
    constexpr double delta_B = 1.1e-2; // sm
}

namespace ioniz {
    constexpr double alpha_i = 1;
    constexpr double Eps_i = 12.13 * phys::e;
}
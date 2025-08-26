#pragma one
#include "conduction.hpp"

class Discharge: protected Conduction{
    public:
        template <typename func>
        double trapez_numeric_int(func integrand_func, const std::size_t time_ind) const {
            double sum = 0;
            for (std::size_t spatial_ind = 0; spatial_ind < sim::simulation_area_steps - 1; spatial_ind++) {
                sum += (integrand_func(time_ind, spatial_ind + 1) + integrand_func(time_ind, spatial_ind)) / 2 * sim::spatial_step;
            }
            return sum;
        }

        Discharge();

        double derive_dniTe_dz(const std::size_t time_ind, const std::size_t spatial_ind) const;

        double numerat_integr_func_disch_curr(const std::size_t time_ind, const std::size_t spatial_ind) const;

        double denominat_integr_func_disch_curr(const std::size_t time_ind, const std::size_t spatial_ind) const;

        double discharge_curr(const std::size_t time_ind) const;

        double electric_field(const std::size_t time_ind, const std::size_t spatial_ind) const;
};
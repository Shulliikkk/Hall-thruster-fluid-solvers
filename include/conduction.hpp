#pragma one
#include "JustInterp/JustInterp.hpp"
#include "walls.hpp"

class Conduction: protected Walls{
    protected:
        const std::vector<double> electron_energy_data = {
            0.1000, 0.1831, 0.4325, 0.8480, 1.430, 2.178, 3.092, 4.173, 5.419, 6.832,
            8.412, 10.16, 12.07, 14.15, 16.39, 18.80, 21.38, 24.12, 27.03, 30.10,
            33.35, 36.75, 40.33, 44.07, 47.97, 52.05, 56.29, 60.69, 65.26, 70.00
        }; 
        const std::vector<double> beta_n_data = {
            4.311e-14, 3.192e-14, 1.804e-14, 2.358e-14, 5.624e-14, 1.104e-13,
            1.689e-13, 2.174e-13, 2.506e-13, 2.695e-13, 2.774e-13, 2.778e-13,
            2.734e-13, 2.662e-13, 2.575e-13, 2.482e-13, 2.387e-13, 2.295e-13,
            2.206e-13, 2.121e-13, 2.043e-13, 1.970e-13, 1.903e-13, 1.841e-13,
            1.785e-13, 1.734e-13, 1.688e-13, 1.647e-13, 1.609e-13, 1.575e-13
        }; 
        const std::vector<double> beta_i_data = {
            0.0, 0.0, 0.0, 
            8.783e-24, 1.061e-19, 1.082e-17, 1.577e-16, 8.640e-16,
            2.745e-15, 6.283e-15, 1.166e-14, 1.879e-14, 2.744e-14, 3.730e-14, 4.808e-14,
            5.952e-14, 7.141e-14, 8.358e-14, 9.588e-14, 1.082e-13, 1.205e-13, 1.326e-13,
            1.446e-13, 1.563e-13, 1.678e-13, 1.789e-13, 1.897e-13, 2.001e-13, 2.102e-13, 2.198e-13
        };
        
    public:
        Conduction();

        double cicl_friq(const std::size_t spatial_ind) const;

        double electron_neutral_coll_friq(const std::size_t time_ind, const std::size_t spatial_ind) const;

        double anomal_coll_friq(const std::size_t spatial_ind) const;

        double moment_exch_coll_friq(const std::size_t time_ind, const std::size_t spatial_ind) const;

        double ioniz_coll_friq(const std::size_t time_ind, const std::size_t spatial_ind) const;

        double electron_cross_field_mobility(const std::size_t time_ind, const std::size_t spatial_ind) const;
};
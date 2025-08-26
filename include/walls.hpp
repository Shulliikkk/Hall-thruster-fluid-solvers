#pragma one
#include <utility>
#include "states.hpp"
#include "params.hpp"

class Walls
 {
    protected:
        std::vector<std::vector<Plasma_state>> states;
    public:
        double SEE_coef(const std::size_t time_ind, const std::size_t spatial_ind) const;

        double ions_wall_coll_friq(const std::size_t time_ind, const size_t spatial_ind) const;

        double electrons_wall_coll_friq(const size_t time_ind, const size_t spatial_ind) const;

        double wall_potential(const std::size_t time_ind, const std::size_t spatial_ind) const;

        double energy_lose_rate(const std::size_t time_ind, const std::size_t spatial_ind) const;
};
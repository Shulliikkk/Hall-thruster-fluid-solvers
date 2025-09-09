#include "solve.hpp"

int main() {
    geom::set_geom(2.5e-2, 3.5e-2, 5e-2);
    disch::set_disch_params(300, 100, 3e-6);
    sim::set_sim_params(5e-2, 5e-4, 5e-4, 5e-9);
    wall::set_wall_material("BNSiO2");
    magfld::set_magnetic_field("default gaussian-distributed");

    Solver solver;
    solver.set_init_conditions(2e17, 2e18, 3, 5);
    auto results = solver.solve();
}

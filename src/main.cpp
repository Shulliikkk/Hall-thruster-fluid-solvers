#include "solver.hpp"

int main() {
    Solver solver(400, 100, 0.005, "ionization_rate_data.csv", "neutral_rate_data.csv");
    auto results = solver.solve();
    for (size_t space_step = 0; space_step < grid::M; space_step++) {
        std::cout << results[grid::N - 1][space_step].T_e << '\n';
    }
}
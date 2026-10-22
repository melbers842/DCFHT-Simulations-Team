#include "Mesh.hpp"
#include "PlasmaState.hpp"
#include "Vec.hpp"

#include <cassert>

int main() {
    Mesh<2> mesh(
        {4, 4},
        {0.0, 0.0},
        {0.001, 0.001}
    );

    PlasmaState<2> state(mesh);

    assert(state.size() == mesh.numNodes());

    state.n_e[0] = 1.0e17;
    state.n_i[0] = 1.0e17;
    state.n_n[0] = 1.0e20;
    state.T_e[0] = 3.0;

    assert(state.isPhysical());
    assert(state.maxChargeSeparation() == 0.0);

    return 0;
}
from sage.all import (
    is_prime,
    GF,
    set_random_seed,
    randint,
)
from richelot_rm.product_point import ProductPoint
from richelot_rm.richelot_product_isogenies import get_arbitrary_square_example
from richelot_rm_utils import (
    RMVertex,
    get_computable_isogeny,
    golden_ratio_action_on_symplectic_torsion,
)

### Fixed Parameters ###
# e = 43
set_random_seed(3)
e = 11
f = 3
p = 2**e * f - 1
assert is_prime(p)
assert p % 4 == 3
E1, E2 = get_arbitrary_square_example(p)
P2e_1, Q2e_1 = E1.torsion_basis(2**e)
P2e_2, Q2e_2 = E2(P2e_1), E2(Q2e_1)
torsion_generators = [
    ProductPoint(P2e_1, E2(0)),
    ProductPoint(E1(0), P2e_2),
    ProductPoint(Q2e_1, E2(0)),
    ProductPoint(E1(0), Q2e_2),
]
initial_vertex = RMVertex(
    (E1, E2), e, torsion_generators, golden_ratio_action_on_symplectic_torsion(2, e)
)

print(type(E1))
raise NotImplementedError("This file is not ready to run yet.")

def hash_message():
    ### Hash a random message ###
    m = [randint(0, 4) for _ in range(e - 2)]
    current_vertex = initial_vertex
    for i in range(len(m)):
        print(f"Step {i+1}/{len(m)}: current vertex is \n {current_vertex}\n")
        kernels, subspaces = current_vertex.generate_RM_kernels()
        assert len(kernels) == 5, f"Expected 5 kernels, got {len(kernels)}"

        phi_kernel = kernels[m[i]]
        phi_subspace = subspaces[m[i]]
        av, phi = get_computable_isogeny(current_vertex, phi_kernel)
        assert (
            phi(phi_kernel[0]) == 0 and phi(phi_kernel[1]) == 0
        ), f"Kernel not mapped to 0: {phi(phi_kernel[0])}, {phi(phi_kernel[1])}"

        next_vertex = current_vertex.get_neighbor(av, phi, phi_subspace)

        current_vertex = next_vertex


# hash_message()

for i in range(100):
    print(f"Hashing iteration {i}/100")
    hash_message()

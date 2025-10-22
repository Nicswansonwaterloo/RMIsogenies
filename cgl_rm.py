from sage.all import is_prime, GF, Integers, identity_matrix, vector, HyperellipticCurve, set_random_seed, randint
from dependencies.Theta_SageMath.theta_structures.couple_point import CouplePoint
from richelot_products import get_arbitrary_square_example
from richelot_rm_utils import RMVertex, get_computable_isogeny, golden_ratio_action_on_symplectic_torsion

### Fixed Parameters ###
set_random_seed(31)
e = 11
f = 3
p = 2**e * f - 1
assert is_prime(p)
assert p % 4 == 3
E1, E2 = get_arbitrary_square_example(p)
P2e_1, Q2e_1 = E1.torsion_basis(2 ** e)
P2e_2, Q2e_2 = E2(P2e_1), E2(Q2e_1)
torsion_generators = [CouplePoint(P2e_1, E2(0)), CouplePoint(E1(0), P2e_2), CouplePoint(Q2e_1, E2(0)), CouplePoint(E1(0), Q2e_2)]
initial_vertex = RMVertex((E1, E2), e, torsion_generators, golden_ratio_action_on_symplectic_torsion(2, e))

### Hash a random message ###
m = [randint(0, 4) for _ in range(e - 1)]
current_vertex = initial_vertex
for i in range(len(m)):
    kernels, subspaces = current_vertex.generate_RM_kernels()
    assert len(kernels) == 5, f"Expected 5 kernels, got {len(kernels)}"
    phi_kernel = kernels[m[i]]
    phi_subspace = subspaces[m[i]]
    av, phi = get_computable_isogeny(current_vertex, phi_kernel)
    next_vertex = current_vertex.get_neighbor(av, phi, phi_subspace)
    current_vertex = next_vertex


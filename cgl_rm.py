from sage.all import is_prime, GF, Integers, identity_matrix, vector
# from sage.schemes.elliptic_curves.weierstrass_morphism import WeierstrassIsomorphism
from sage.schemes.hyperelliptic_curves.invariants import absolute_igusa_invariants_kohel
from dependencies.Theta_SageMath.theta_structures.couple_point import CouplePoint
from richelot_products import get_arbitrary_square_example, get_isogeny_from_product_two_kernel, get_maximal_isotropic_subspaces
from richelot_rm_utils import RMVertex, get_isogeny_from_jacobian_two_kernel, golden_ratio_action_on_symplectic_torsion

e = 11
f = 3
p = 2**e * f - 1
assert is_prime(p)
assert p % 4 == 3

E1, E2 = get_arbitrary_square_example(p)
P2e_1, Q2e_1 = E1.torsion_basis(2 ** e)
P2e_2, Q2e_2 = E2(P2e_1), E2(Q2e_1)


torsion_basis = [CouplePoint(P2e_1, E2(0)), CouplePoint(E1(0), P2e_2), CouplePoint(Q2e_1, E2(0)), CouplePoint(E1(0), Q2e_2)]
initial_vertex = RMVertex((E1, E2), e, golden_ratio_action_on_symplectic_torsion(2, e), torsion_basis)

def is_subspace_invariant(subspace, action):
    """
    Check if the subspace is invariant under the action.
    subspace: a matrix whose rows are a basis for the subspace (over GF(2))
    action: a matrix representing the action (over GF(2))
    """
    for v in subspace.rows():
        w = action * v
        try:
            if not subspace.transpose().solve_right(w):
                return False
        except ValueError: #This is the case the augmented matrix is not full rank, i.e. no solution.
            return False
    return True


def compute_neighbor_vertices(vertex):
    r = vertex.r
    torsion_basis = vertex.torsion_basis
    two_torsion_basis = [2**(r - 1)*P for P in torsion_basis]
    new_action = vertex.action.change_ring(Integers(2**(r - 1)))
    print(type(two_torsion_basis[0]))
    
    def vec_to_point(vec):
                components = [int(vec[i]) * two_torsion_basis[i] for i in range(4)]
                return components[0] + components[1] + components[2] + components[3]
    
    isotropic_subspaces = get_maximal_isotropic_subspaces(2)
    two_torsion_action = vertex.action.change_ring(GF(2))
    neighboring_vertices = []
 
    for subspace in isotropic_subspaces:
        if is_subspace_invariant(subspace, two_torsion_action):
            v1 = subspace.row(0)
            v2 = subspace.row(1)
            kernel_generators = (vec_to_point(v1), vec_to_point(v2))
            if vertex.is_product():
                av, phi = get_isogeny_from_product_two_kernel(kernel_generators)
            else:
                av, phi = get_isogeny_from_jacobian_two_kernel(kernel_generators)
            
            # phi kills "half" of the 2-torsion, meaning we have to adjust so that the RM action is now on 2^(r-1)-torsion.
            new_torsion_basis = []
            for i in range(4):
                if two_torsion_basis[i] not in kernel_generators:
                    new_torsion_basis.append(phi(2 * torsion_basis[i]))
                else:
                    new_torsion_basis.append(phi(torsion_basis[i]))

            neighboring_vertices.append(RMVertex(av, r - 1, new_action, new_torsion_basis))

    return neighboring_vertices


for neighbor in compute_neighbor_vertices(initial_vertex):
    if str(neighbor)[0] == "J":
        jacobian_neighbor = neighbor
        print(compute_neighbor_vertices(jacobian_neighbor))
        break
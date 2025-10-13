from sage.all import is_prime, GF, ZZ, identity_matrix
# from sage.schemes.elliptic_curves.weierstrass_morphism import WeierstrassIsomorphism
from sage.schemes.hyperelliptic_curves.invariants import absolute_igusa_invariants_kohel
from dependencies.Theta_SageMath.theta_structures.couple_point import CouplePoint
from richelot_products import get_arbitrary_square_example, get_isogeny_from_product_two_kernel, get_maximal_isotropic_subgroups_of_N_torsion, get_maximal_isotropic_subspaces
from richelot_rm_utils import RMVertex, golden_ratio_action_on_symplectic_torsion

# from dependencies.Theta_SageMath.utilities.supersingular import torsion_basis, torsion_basis_with_pairing


    # def get_neighbors(self):
    #     if isinstance(self.variety, tuple):
    #         E1, E2 = self.variety
    #         P2_1, Q2_1 = isotropic_torsion_basis(E1, 2)
    #         P2_2, Q2_2 = isotropic_torsion_basis(E2, 2)
    #         basis = [CouplePoint(P2_1, E2(0)), CouplePoint(Q2_1, E2(0)), CouplePoint(E1(0), P2_2), CouplePoint(E1(0), Q2_2)]
    #         isotropic_subspaces = get_maximal_isotropic_subspaces(2) # Given canonically
    #         two_torsion_action = self.RM_action_on_2e_torsion.change_ring(GF(2))
    #         neighboring_vertices = []
    #         # Check which subspaces are invariant under the RM action
    #         for subspace in isotropic_subspaces:
    #             if (subspace * two_torsion_action) * subspace.transpose() == identity_matrix(GF(2), 2):
    #                 v1 = subspace.row(0)
    #                 v2 = subspace.row(1)
    #                 stable_points = (vector_to_point(v1, basis), vector_to_point(v2, basis))
    #                 av, phi = get_isogeny_from_product_two_kernel(stable_points)
    #                 print(av)
    #                 print(phi(basis[0]))
                    
    #         return None
    #     else:
    #         J = self.variety
    #         H = isotropic_torsion_basis(J, 2)
    #         phi = FromJacToJac(H[0], H[1])
    #         E1, E2 = phi.codomain()
    #         RM_action_on_2e_torsion = self.RM_action_on_2e_torsion.change_ring(Integers(2**(e-1)))
    #         return [RM_graph_richelot_vertex((E1, E2), RM_action_on_2e_torsion)]




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

def compute_neighbor_vertices(vertex):
    if vertex.is_product():
        torsion_basis = vertex.torsion_basis
        r = vertex.r
        two_torsion_basis = [2**(r - 1)*P for P in torsion_basis]
        
        def vec_to_point(vec):
            components = [vec[i] * two_torsion_basis[i] for i in range(4)]
            return components[0] + components[1] + components[2] + components[3]

        isotropic_subspaces = get_maximal_isotropic_subspaces(2) # Returns canonical output.
        two_torsion_action = vertex.action.change_ring(GF(2))

        # Check which subspaces are invariant under the RM action
        for subspace in isotropic_subspaces:
            if (subspace * two_torsion_action) * subspace.transpose() == identity_matrix(GF(2), 2):
                v1 = subspace.row(0)
                v2 = subspace.row(1)
                kernel_generators = (vec_to_point(v1), vec_to_point(v2))
                av, phi = get_isogeny_from_product_two_kernel(kernel_generators)

print(compute_neighbor_vertices(initial_vertex))
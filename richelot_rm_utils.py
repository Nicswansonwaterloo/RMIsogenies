import copy
from sage.all import Matrix, Integers, GF, VectorSpace, Integer, matrix, identity_matrix

from sage.schemes.hyperelliptic_curves.invariants import absolute_igusa_invariants_kohel
from couple_point import CouplePoint
from richelot_products import get_isogeny_from_product_two_kernel
from richelot_jacobians import get_isogeny_from_jacobian_two_kernel


def golden_ratio_endomorphism(cp_pt):
    """ The endomorphism of E x E given by (P, Q) |-> (Q, P + Q). """
    P, Q = cp_pt
    return CouplePoint(Q, P + Q)

def golden_ratio_action_on_symplectic_torsion(ell=2, e=1):
    """ Returns the matrix of the action of the golden ratio endomorphism on the ell^e-torsion,
        with respect to the basis [ (P2_1, 0), (0, P2_2), (Q2_1, 0), (0, Q2_2) ].
    """
    Zle = Integers(ell**e)
    return Matrix(Zle, [[0, 1, 0, 0],
                        [1, 1, 0, 0],
                        [0, 0, 0, 1],
                        [0, 0, 1, 1]])

# The issue is here: I need library that computes arbitrary 2,2 isogenies! This has proven to be much harder than expected.
def get_computable_isogeny(domain_vertex, kernel):
    if domain_vertex.is_product():
        av, phi = get_isogeny_from_product_two_kernel(kernel)
    else:
        av, phi = get_isogeny_from_jacobian_two_kernel(kernel, domain_vertex.variety)
    return av, phi

class RMVertex:
    """
        A vertex in the (2,2)-isogeny graph of supersingular abelian surfaces with RM.
        
        - The vertex can either be a product of two elliptic curves or the Jacobian of a genus 2 curve. This is stored in the attribute 'variety'.
        - The RM is represented by its action on 2^r-torsion, for some r. This integer is stored in the attribute 'r'.
        - The action on the 2^r-torsion is stored in the attribute 'action', which is a 4x4 matrix with coefficients in Z/(2^r)Z.
        - The weil pairing on the 2-torsion is stored in the attribute 'weil_pairing', which is a 4x4 matrix with coefficients in GF(2).
    """
    def __init__(self, variety, r, torsion_generators, rm_action):
        self.variety = variety
        self.r = r
        self.torsion_generators = torsion_generators
        self.two_torsion_basis = [2**(r - 1) * P for P in torsion_generators]
        self.rm_action = rm_action
        self.weil_pairing = self._compute_weil_pairing()


    def _compute_weil_pairing(self):
        Me = Matrix(GF(2), 4, 4)
        for i, P in enumerate(self.two_torsion_basis):
            for j, Q in enumerate(self.two_torsion_basis):
                if self.is_jacobian():
                    if i == j:
                        entry = 1
                    else:
                        G1 = P[0]
                        G2 = Q[0]
                        entry = 1 if G1.gcd(G2) == 1 else -1
                else:
                    entry = P.weil_pairing(Q, 2)
                # Apply a reduction map mu_{2^r} -> GF(2)
                if entry == 1:
                    Me[i, j] = 0
                else:
                    Me[i, j] = 1
        return Me

    def is_product(self):
        return isinstance(self.variety, tuple)

    def is_jacobian(self):
        return not self.is_product()

    def __repr__(self):
        if isinstance(self.variety, tuple):
            return f"Product: ({self.variety[0].j_invariant()}, {self.variety[1].j_invariant()})"
        return f"Jacobian: {absolute_igusa_invariants_kohel(self.variety)}"
    
    def __eq__(self, other):
        return repr(self.variety) == repr(other.variety)
    
    def _vector_to_point(self, vec, two_torsion = False):
        generators = self.two_torsion_basis if two_torsion else self.torsion_generators
        modulus = 2 if two_torsion else 2 ** self.r
        vec = [int(vec[i]) % modulus for i in range(4)]
        components = [int(vec[i]) * generators[i] for i in range(4)]
        return components[0] + components[1] + components[2] + components[3]
    
    def _get_maximal_isotropic_subspaces(self):
        V = VectorSpace(GF(2), 4)
        isotropic_subspaces = []
        for W in V.subspaces(2):
            basis_matrix = W.basis_matrix().transpose()
            if (basis_matrix.transpose() * self.weil_pairing * basis_matrix).is_zero():
                isotropic_subspaces.append(basis_matrix)

        return isotropic_subspaces

    def generate_RM_kernels(self):
        maximal_isotropic_subspaces = self._get_maximal_isotropic_subspaces()
        assert len(maximal_isotropic_subspaces) == 15, f"Expected 15 maximal isotropic subspaces, got {len(maximal_isotropic_subspaces)}"
        Mphi = self.rm_action.change_ring(GF(2))
        kernels = []
        subspaces = []
        for subspace in maximal_isotropic_subspaces:
            phi_subspace = Mphi * subspace
            P = subspace.augment(phi_subspace)
            if P.rank() == 2:
                kernel = [self._vector_to_point(subspace.column(i), two_torsion=True) for i in range(2)]
                kernels.append(kernel)
                subspaces.append(subspace)

        return kernels, subspaces

    def get_neighbor(self, codomain, phi, phi_subspace):
        W = phi_subspace
        id_2 = identity_matrix(GF(2), 2)
        A = W.transpose() * self.weil_pairing
        # W_perp = A.solve_right(id_2) # The Lagrangian complement of W
        V = A.solve_right(id_2)
        P = V.transpose() * self.weil_pairing * V
        c = P[0,1]
        A_corr = matrix(GF(2), [[0, c], [0, 0]])
        W_perp = V + W * A_corr # Possible correction to ensure W_perp is Lagrangian

        # Check that 2 torsion basis is well formed.
        assert W_perp.transpose() * self.weil_pairing * W_perp == 0, f"W_perp is not isotropic:\n {W_perp.transpose() * self.weil_pairing * W_perp }"

        # Form the change of basis matrix C from the old torsion basis to the new one
        C = W.augment(W_perp)
        assert C.is_invertible(), f"{C} \n is not invertible."
        C_inv = C.inverse()

        # Compute the new torsion generators on the domain and codomain
        new_torsion_gens = [self._vector_to_point(col, two_torsion=False) for col in C.columns()]
        codomain_torsion_gens = [phi(new_torsion_gens[0]),phi(new_torsion_gens[1]),phi(2 * new_torsion_gens[2]),phi(2 * new_torsion_gens[3])]

        # Check that the codomain gens are linearly independent
        print(f"Codomain: {codomain}")
        for i in range(4):
            print(f"Codomain torsion gen {i}: {codomain_torsion_gens[i]}")
            if isinstance(codomain, tuple):
                print(f" - order: {codomain_torsion_gens[i][0].order()}, {codomain_torsion_gens[i][1].order()}")
                

        # Check that the new torsion generators have correct orders
        should_be_zero = [Integer(2**(self.r - 1)) * P for P in codomain_torsion_gens]
        two_torsion_orders = [Integer(2**(self.r - 2)) * P for P in codomain_torsion_gens]        
        assert all(P == 0 for P in should_be_zero), f"Should be zero check failed:\n {should_be_zero}"
        assert all(P != 0 for P in two_torsion_orders), f"New torsion generators do not have correct orders:\n {two_torsion_orders}"

        # Compute the RM action on the new_torsion_gens basis
        Mphi = self.rm_action.change_ring(GF(2))
        rm_action_on_C = C_inv * Mphi * C # This behaves as expected mod 2

        # This can probabily be done lifting to 2** r or 2**(r-1)
        C_lifted = C.change_ring(Integers())
        rm_action_lifted = self.rm_action.change_ring(Integers())
        assert C_lifted.is_invertible(), f"C_lifted is not invertible:\n {C_lifted}"
        C_inv_lifted = C_lifted.inverse()
        rm_action_on_C_lifted = C_inv_lifted * rm_action_lifted * C_lifted # Have checked to be correctly representing the action of RM on the changed basis on the domain.

        # Compute the RM action on the codomain with respect to codomain_torsion_gens
        assert rm_action_on_C_lifted[2:4, 0:2] % 2 == 0, f"lower-left block not even:\n{rm_action_on_C}"
        rm_action_prime = copy.deepcopy(rm_action_on_C_lifted)
        rm_action_prime[0:2, 2:4] *= 2; rm_action_prime[2:4, 0:2] /= 2
        rm_action_prime = rm_action_prime.change_ring(Integers(2**(self.r - 1)))

        next_vertex = RMVertex(codomain, self.r - 1, codomain_torsion_gens, rm_action_prime)
        print(f"rm_action on domain:\n{rm_action_on_C_lifted} \n")
        print(f"rm_action on domain mod 2:\n{rm_action_on_C} \n")
        print(f"rm_action on codomain:\n{next_vertex.rm_action} \n")
        print(f"rm_action on codomain mod 2:\n{next_vertex.rm_action.change_ring(GF(2))} \n")

        # Check that rm' \circ phi = phi \circ rm on the reduced generators of 2^r
        # for i in range(4):
        #     print(f"Checking commutation on generator {i} \n \n \n")
        #     pi_four_torsion = [Integer(2**(self.r - 2)) * P for P in new_torsion_gens]
        #     def Pi_basis_to_point(col):
        #         components = [int(col[j]) * pi_four_torsion[j] for j in range(4)]
        #         print(f" ----- components: {components} ----- ")
        #         return components[0] + components[1] + components[2] + components[3]
        #     rm_Pi_col = rm_action_on_C.column(i) # This is the image of rm(P_i) in the P_i basis
        #     print(f"rm(P_{i}) in P_i basis (mod 2): {rm_Pi_col}")
        #     rm_Pi = Pi_basis_to_point(rm_Pi_col)
        #     print(f"rm(P_{i}) as point (mod 2): {rm_Pi}")
        #     phi_rm_Pi = phi(rm_Pi)
        #     print(f"phi(rm(P_{i})) in codomain (mod 2): {phi_rm_Pi}")
        #     if i > 1:
        #         phi_rm_Pi = 2 * phi_rm_Pi # This is rm(T_i') in the codomain
        #     assert phi_rm_Pi * Integer(2) == 0, f"phi(rm(P_{i})) does not have correct order:\n {phi_rm_Pi}"

        #     Ti_prime = two_torsion_orders[i]
        #     print(f"T_{i}' in codomain: {Ti_prime}")
        #     rm_Ti_prime_col = next_vertex.rm_action.column(i)
        #     print(f"rm(T_{i}') in codomain basis (mod 2): {rm_Ti_prime_col}")
        #     rm_Ti_prime = next_vertex._vector_to_point(rm_Ti_prime_col, two_torsion=True)
        #     assert rm_Ti_prime * Integer(2) == 0, f"rm(T_{i}') does not have correct order:\n {rm_Ti_prime}"

        #     assert phi_rm_Pi - rm_Ti_prime == 0, f"RM action does not commute with isogeny on generator {i}:\n phi(rm(P_{i})) = \n {phi_rm_Pi}\n rm(T_{i}') = \n{rm_Ti_prime} \n diff: \n {phi_rm_Pi - rm_Ti_prime}"
            
        #     # print("mod 2 ----")

        #     print("Lifted ----")
        #     def Pi_basis_to_point(col):
        #         components = [(int(col[j]) % 2 ** self.r) * new_torsion_gens[j] for j in range(4)]
        #         return components[0] + components[1] + components[2] + components[3]
        #     rm_Pi_col = rm_action_on_C_lifted.column(i) # This is the image of rm(P_i) in the P_i basis
        #     print(f"rm(P_{i}) in P_i basis: {rm_Pi_col}")
        #     rm_Pi = Pi_basis_to_point(rm_Pi_col)
        #     print(f"rm(P_{i}) as point: {rm_Pi}")
        #     phi_rm_Pi = phi(rm_Pi)
        #     print(f"phi(rm(P_{i})) in codomain: {phi_rm_Pi}")
        #     if i > 1:
        #         phi_rm_Pi = 2 * phi_rm_Pi # This is rm(T_i') in the codomain

        #     assert phi_rm_Pi * Integer(2 ** (self.r - 1)) == 0, f"phi(rm(P_{i})) does not have correct order:\n {phi_rm_Pi}"
        #     assert phi_rm_Pi * Integer(2 ** (self.r - 2)) != 0, f"phi(rm(P_{i})) does not have correct order:\n {phi_rm_Pi}"

        #     Ti_prime = codomain_torsion_gens[i]
        #     rm_Ti_prime_col = next_vertex.rm_action.column(i)
        #     print(f"T_{i}' in codomain: {Ti_prime}")
        #     rm_Ti_prime = next_vertex._vector_to_point(rm_Ti_prime_col, two_torsion=False)
        #     assert phi_rm_Pi - rm_Ti_prime == 0, f"RM action does not commute with isogeny on generator {i}:\n phi(rm(P_{i})) = \n {phi_rm_Pi}\n rm(T_{i}') = \n{rm_Ti_prime} \n diff: \n {phi_rm_Pi - rm_Ti_prime}"

        return next_vertex
                
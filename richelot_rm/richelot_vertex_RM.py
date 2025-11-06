from sage.all import GF, identity_matrix, Integers, Matrix
import copy
from richelot_rm.richelot_vertex import RichelotVertex
from richelot_rm.genus_two_structures import GenusTwoStructure


class RMVertex(RichelotVertex):
    def __init__(
        self,
        g2_structure: GenusTwoStructure,
        r,
        two_r_torsion_generators,
        rm_action_on_two_r_torsion,
    ):
        self.r = r  # int
        self.two_r_torsion_generators = (
            two_r_torsion_generators  # List of 2^r-torsion points generating A[2^r]
        )
        self.rm_action_on_two_r_torsion = rm_action_on_two_r_torsion  # Matrix representing the RM action on the 2^r-torsion points (over Z/2^rZ)
        two_torsion_generators = [2 ** (r - 1) * P for P in two_r_torsion_generators]
        super().__init__(g2_structure, two_torsion_generators)

    def _vector_to_large_torsion_point(self, vec):
        components = [int(vec[i]) * self.two_r_torsion_generators[i] for i in range(4)]
        return components[0] + components[1] + components[2] + components[3]

    # The way that the torsion is passed in get_neighbors, we can do something clever to make sure that the dual is always first.
    # In particular, the dual kernel is always the one spanned by the last two kernel vectors:
    # [ 0 0 ]
    # [ 0 0 ]
    # [ 1 0 ]
    # [ 0 1 ]
    # After this is added, the rest are ordered in a deterministic way.
    def _get_all_two_kernels(self):
        maximal_isotropic_subspaces = self._get_maximal_isotropic_subspaces()
        M_rm = self.rm_action_on_two_r_torsion.change_ring(GF(2))
        kernels = []
        subspaces = []
        for subspace in maximal_isotropic_subspaces:
            phi_subspace = M_rm * subspace
            P = subspace.augment(phi_subspace)
            if P.rank() == 2:
                subspaces.append(subspace)

        W_dual = Matrix(GF(2), [[0, 0], [0, 0], [1, 0], [0, 1]])
        assert W_dual in subspaces, "Dual kernel does not seem to perserve RM." # Should probably be an equality check, possibly does not happen on initial vertex depending...
        # The way that the torsion is passed in get_neighbors, we can do something clever to make sure that the dual is always first.
        # In particular, the dual kernel is always the one spanned by the last two kernel vectors:
        # [ 0 0 ]
        # [ 0 0 ]
        # [ 1 0 ]
        # [ 0 1 ]
        # After this is added, the rest are ordered in a deterministic way.
        subspaces.remove(W_dual)
        subspaces.sort(key=lambda M: M.list())
        subspaces.insert(0, W_dual)

        for subspace in subspaces:
            kernel_gens = [self._vector_to_point(col) for col in subspace.columns()]
            kernels.append(kernel_gens)

        return kernels, subspaces

    def _compute_neighboring_isogenies(self):
        if self.computed_neighbors is None:
            kernels, subspaces = self._get_all_two_kernels()
            neighbors_with_edges = []
            for kernel, subspace in zip(kernels, subspaces):
                codomain, isogeny = self._compute_isogeny(kernel)
                neighbors_with_edges.append((codomain, isogeny, subspace))
            self.computed_neighbors = neighbors_with_edges

        return self.computed_neighbors

    def get_neighbors_with_multiplicities(self):
        neighbors_with_edges = self._compute_neighboring_isogenies()
        neighbors = {}
        for neighbor, phi, W in neighbors_with_edges:
            id_2 = identity_matrix(GF(2), 2)
            A = W.transpose() * self.weil_pairing_two_torsion_action
            W_c = A.solve_right(id_2)
            C = W.augment(W_c)
            assert C.is_invertible(), f"{C} \n is not invertible."
            C_inv = C.inverse()

            changed_torsion_gens = [
                self._vector_to_large_torsion_point(col) for col in C.columns()
            ]
            codomain_torsion_gens = [
                phi(changed_torsion_gens[0]),
                phi(changed_torsion_gens[1]),
                phi(2 * changed_torsion_gens[2]),
                phi(2 * changed_torsion_gens[3]),
            ]

            # Check that the new torsion generators have correct orders
            should_be_zero = [(2 ** (self.r - 1)) * P for P in codomain_torsion_gens]
            two_torsion_orders = [
                (2 ** (self.r - 2)) * P for P in codomain_torsion_gens
            ]
            assert all(
                P == 0 for P in should_be_zero
            ), f"Should be zero check failed:\n {should_be_zero}"
            assert all(
                P != 0 for P in two_torsion_orders
            ), f"New torsion generators do not have correct orders:\n {two_torsion_orders}"

            # Compute the RM action on the new_torsion_gens basis
            Mphi = self.rm_action_on_two_r_torsion.change_ring(GF(2))
            rm_action_on_C = C_inv * Mphi * C  # This behaves as expected mod 2

            # This can probabily be done lifting to 2** r or 2**(r-1)
            C_lifted = C.change_ring(Integers())
            rm_action_lifted = self.rm_action_on_two_r_torsion.change_ring(Integers())
            assert C_lifted.is_invertible(), f"C_lifted is not invertible:\n {C_lifted}"
            C_inv_lifted = C_lifted.inverse()
            rm_action_on_C_lifted = (
                C_inv_lifted * rm_action_lifted * C_lifted
            )  # Have checked to be correctly representing the action of RM on the changed basis on the domain.

            # Compute the RM action on the codomain with respect to codomain_torsion_gens
            assert (
                rm_action_on_C_lifted[2:4, 0:2] % 2 == 0
            ), f"lower-left block not even:\n{rm_action_on_C}"
            rm_action_prime = copy.deepcopy(rm_action_on_C_lifted)
            rm_action_prime[0:2, 2:4] *= 2
            rm_action_prime[2:4, 0:2] /= 2
            rm_action_prime = rm_action_prime.change_ring(Integers(2 ** (self.r - 1)))

            neighbor = RMVertex(
                neighbor, self.r - 1, codomain_torsion_gens, rm_action_prime
            )

            if neighbor in neighbors:
                neighbors[neighbor] += 1
            else:
                neighbors[neighbor] = 1

        return neighbors

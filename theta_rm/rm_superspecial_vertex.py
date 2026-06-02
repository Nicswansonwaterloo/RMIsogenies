from sage.all import (
    GF,
    identity_matrix,
    Integers,
    Matrix,
    cached_method,
    matrix,
    copy,
    ZZ,
)
from theta_rm.superspecial_vertex import OMEGA, SYMPLECTIC_GL2_SUBSPACES, PPASVertex
from theta_rm.genus_two_structures import ThetaStructure
from theta_rm.torsion_wrappers import ThetaTorsion

class RMVertex(PPASVertex):
    """Similar to PPASVertex, but keep track of RM action on large torsion."""

    def __init__(
        self, g2_structure, r, two_r_theta_torsion, M_rm
    ):
        if not isinstance(g2_structure, ThetaStructure):
            raise TypeError("g2_structure must be a genus-2 structure.")
        if r < 3:
            raise ValueError("r must be at least 3.")
        if not isinstance(two_r_theta_torsion, ThetaTorsion):
            raise ValueError("two_r_theta_torsion must be a ThetaTorsion instance.")
        if (
            M_rm.nrows() != 4
            or M_rm.ncols() != 4
        ):
            raise ValueError("M_rm must be a 4x4 matrix.")

        self.r = r
        self.two_r_theta_torsion = two_r_theta_torsion
        self.M_rm = M_rm

        super().__init__(g2_structure)

    def _get_all_kernels(self, force_deterministic=True):
        # Compute matrix action on the 2-torsion basis:
        # TODO: it may be a problem that this is over GF(2) now...
        M_rm = self.M_rm.change_ring(GF(2))
        kernels = []
        subspaces = []

        # Check for maximally isotropic subspaces that are preserved by the RM action, i.e. M_rm * W is contained in W.
        for subspace in SYMPLECTIC_GL2_SUBSPACES:
            phi_subspace = M_rm * subspace
            P = subspace.augment(phi_subspace)
            if P.rank() == 2:
                subspaces.append(subspace)

        if force_deterministic:
            # This naturally puts the dual kernel first.
            subspaces.sort(key=lambda M: M.list())
            # W_dual = Matrix(GF(2), [[0, 0], [0, 0], [1, 0], [0, 1]])
            # assert W_dual == subspaces[0], "Dual kernel is not first after sorting."

        if self.eight_torsion is None:
            self.eight_torsion = self.two_r_theta_torsion.get_projection(3)
        for subspace in subspaces:
            kernels.append(ThetaTorsion(self.eight_torsion.matrix_to_kernel(subspace), self.r))

        return kernels, subspaces

    # def _compute_neighboring_isogenies(self):
    #     kernels, subspaces = self._get_all_kernels()
    #     neighbors_with_edges = []
    #     for kernel, subspace in zip(kernels, subspaces):
    #         codomain, isogeny = self._compute_isogeny(kernel)
    #         neighbors_with_edges.append((codomain, isogeny, subspace))

    #     return neighbors_with_edges

    def _lift_symplectic(self, C):
        """
        Lifts a mod 2 integer symplectic matrix C to a mod 2^r symplectic matrix.
        Requires C^T * Omega * C == Omega (mod 2) and C invertible mod 2.
        """
        n = C.nrows()
        r = self.r
        assert (
            C.transpose() * matrix(GF(2), OMEGA) * C == OMEGA
        ), "C does not satisfy the symplectic condition mod 2."
        assert C.is_invertible(), "C is not invertible mod 2."

        C = Matrix(ZZ, C)
        O_inv = Matrix(GF(2), OMEGA).inverse()

        # Now we hensel lift C to satisfy C^T * Omega * C == Omega (mod 2^r) iteratively.
        for k in range(1, r):
            E = (C.T * OMEGA * C - OMEGA) / (2**k)
            E_mod2 = Matrix(GF(2), E)

            # Form S as the strictly upper triangular part of E (mod 2)
            S = Matrix(GF(2), n, n)
            for i in range(n):
                for j in range(i + 1, n):
                    S[i, j] = E_mod2[i, j]
                    
            W_mod2 = O_inv * S
            W = Matrix(ZZ, W_mod2)

            C = C + (2**k) * C * W

        return C.apply_map(lambda x: x % (2**r))

    def _form_special_basis(self, W):
        """Given a subspace W over GF(2) representing ker phi, returns a symplectic basis of A[2^r] that decends down to ker phi."""
        id_2 = identity_matrix(GF(2), 2)

        # First solve for V in W^T * M_e * V_0 = I_2. Solution may not be maximally isotropic.
        A = W.transpose() * OMEGA
        V_0 = A.solve_right(id_2)
        S = V_0.transpose() * OMEGA * V_0

        # Correct V such that V^T * M_e * V = 0, which ensures the new basis for the kernel is symplectic.
        if S == 0:
            V = V_0
        else:
            symplectic_correction = matrix(GF(2), [[0, 1], [0, 0]])
            V = V_0 + (W * symplectic_correction)

        assert (
            V.transpose() * OMEGA * V == 0
        ), f"V is not isotropic:\n {V.transpose() * OMEGA * V}"

        # C0, the change-of-basis matrix from the original basis to the new basis, but only on the 2-torsion.
        C0 = W.augment(V)
        assert C0.is_invertible(), f"C is not invertible: \n{C0}"

        # Now, we hensel lift C0 to C:
        C = self._lift_symplectic(C0)
        C = C.change_ring(Integers(2 ** (self.r - 1)))
        assert C.is_invertible(), f"C_lifted is not invertible:\n {C}"
        C_inv = C.inverse()

        # Compute the RM action on the codomain in the new basis.
        M_rm = copy(self.M_rm)
        M_rm_new = C_inv * M_rm * C
        M_rm_new = M_rm_new.change_ring(ZZ)

        # Correct issues coming from non-coprime torsion.
        assert M_rm_new[2:4, 0:2] % 2 == 0, f"Not of expected form.\n{M_rm_new}"
        M_rm_new[0:2, 2:4] *= 2
        M_rm_new[2:4, 0:2] /= 2
        M_rm_new = M_rm_new.change_ring(Integers(2 ** (self.r - 1)))

        # Compute the new well-formed generators from the columns of C after pushing phi.
        special_basis = [
            self._vector_to_point(col, self.two_r_theta_torsion)
            for col in C.columns()
        ]
        return special_basis, M_rm_new

    @cached_method
    def get_neighbors_with_multiplicities(self):
        neighbors_with_edges = self._compute_neighboring_isogenies()
        neighbors = {}

        for neighbor, phi, W in neighbors_with_edges:

            special_basis, M_rm_new = self._form_special_basis(W)
            codomain_torsion_gens = [
                phi(special_basis[0]),
                phi(special_basis[1]),
                phi(2 * special_basis[2]),
                phi(2 * special_basis[3]),
            ]
            # assert all(P.has_order(2, self.r - 1) for P in codomain_torsion_gens)

            neighbor = RMVertex(neighbor, self.r - 1, codomain_torsion_gens, M_rm_new)

            if neighbor in neighbors:
                neighbors[neighbor] += 1
            else:
                neighbors[neighbor] = 1

        return neighbors



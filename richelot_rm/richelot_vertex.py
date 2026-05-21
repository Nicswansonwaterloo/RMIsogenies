from sage.all import Matrix, GF, VectorSpace, cached_method, ZZ
from richelot_rm.genus_two_structures import GenusTwoStructureAbstract
from richelot_rm.richelot_product_isogenies import (
    compute_2_isogeny_from_product,
    get_symplectic_two_torsion_prod,
)
from richelot_rm.richelot_jacobian_isogeny import (
    compute_2_isogeny_from_jacobian,
    get_symplectic_two_torsion_jac,
)

# The action of Omega with respect to a symplectic basis.
OMEGA = Matrix(ZZ, [[0, 0, 1, 0], [0, 0, 0, 1], [-1, 0, 0, 0], [0, -1, 0, 0]])
# All 15 maximally isotropic subspaces of (Z/2Z)^4 with respect to OMEGA, represented as the row space of a matrix in reduced row echelon form.
SYMPLECTIC_GL2_SUBSPACES = [
    Matrix(GF(2), [
        [1, 0],
        [0, 1],
        [0, 0],
        [0, 0],
    ]),
    Matrix(GF(2), [
        [1, 0],
        [0, 1],
        [0, 0],
        [0, 1],
    ]),
    Matrix(GF(2), [
        [1, 0],
        [0, 1],
        [0, 1],
        [1, 0],
    ]),
    Matrix(GF(2), [
        [1, 0],
        [0, 1],
        [0, 1],
        [1, 1],
    ]),
    Matrix(GF(2), [
        [1, 0],
        [0, 1],
        [1, 0],
        [0, 0],
    ]),
    Matrix(GF(2), [
        [1, 0],
        [0, 1],
        [1, 0],
        [0, 1],
    ]),
    Matrix(GF(2), [
        [1, 0],
        [0, 1],
        [1, 1],
        [1, 0],
    ]),
    Matrix(GF(2), [
        [1, 0],
        [0, 1],
        [1, 1],
        [1, 1],
    ]),
    Matrix(GF(2), [
        [1, 0],
        [1, 0],
        [0, 1],
        [0, 1],
    ]),
    Matrix(GF(2), [
        [1, 0],
        [1, 0],
        [0, 1],
        [1, 1],
    ]),
    Matrix(GF(2), [
        [1, 0],
        [0, 0],
        [0, 0],
        [0, 1],
    ]),
    Matrix(GF(2), [
        [1, 0],
        [0, 0],
        [1, 0],
        [0, 1],
    ]),
    Matrix(GF(2), [
        [0, 0],
        [1, 0],
        [0, 1],
        [0, 0],
    ]),
    Matrix(GF(2), [
        [0, 0],
        [1, 0],
        [0, 1],
        [0, 1],
    ]),
    Matrix(GF(2), [
        [0, 0],
        [0, 0],
        [1, 0],
        [0, 1],
    ]),
]


class RichelotVertex:
    """Vertex wrapper for a genus-2 structure with cached 2-torsion data. The two_torsion_generators must form a symplectic basis of the 2-torsion."""

    def __init__(
        self, g2_structure: GenusTwoStructureAbstract, two_torsion_generators=None
    ):
        self.g2_structure = g2_structure
        self.invariants = g2_structure.get_isomorphism_class_invariants()
        self.two_torsion_generators = two_torsion_generators

    def _compute_two_torsion_generators(self):
        """Return a symplectic 2-torsion basis for the current structure."""
        if self.g2_structure.is_product:
            return get_symplectic_two_torsion_prod(self.g2_structure)
        return get_symplectic_two_torsion_jac(self.g2_structure)

    def _compute_two_torsion(self):
        """Return a cached symplectic 2-torsion basis for the current structure."""
        if self.two_torsion_generators is None:
            self.two_torsion_generators = self._compute_two_torsion_generators()
        return self.two_torsion_generators

    def __repr__(self):
        if self.g2_structure.is_product:
            return f"Product: {self.invariants}"
        return f"Jacobian: {self.invariants}"

    def __eq__(self, other):
        """Equality is determined by the isomorphism class invariants."""
        if not isinstance(other, RichelotVertex):
            return False
        return self.invariants == other.invariants

    def __hash__(self):
        return hash(self.invariants)

    # The types given by Florian and Smith: https://eprint.iacr.org/2021/013.pdf
    def get_type(self):
        """Return the Florian-Smith vertex type label."""
        if self.g2_structure.is_product:
            if self.invariants[0] == 1728 or self.invariants[1] == 1728:
                if self.invariants[0] == self.invariants[1]:
                    return R"S_1728"
                elif self.invariants[0] == 0 or self.invariants[1] == 0:
                    return R"P_0_1728"
                else:
                    return R"P_1728"
            elif self.invariants[0] == 0 or self.invariants[1] == 0:
                if self.invariants[0] == self.invariants[1]:
                    return R"S_0"
                else:
                    return R"P_0"
            elif self.invariants[0] == self.invariants[1]:
                return R"S"
            return R"P"
        elif self.g2_structure.is_jacobian:
            return R"J"

        raise ValueError("Unknown genus 2 structure type.")

    def get_type_latex(self):
        """Return the vertex type as a LaTeX string."""
        regular_string = self.get_type()
        conversion = {
            "S_1728": R"\sum_{1728}",
            "P_0_1728": R"\prod_{0,1728}",
            "P_1728": R"\prod_{1728}",
            "S_0": R"\sum_{0}",
            "P_0": R"\prod_{0}",
            "S": R"\sum",
            "P": R"\prod",
            "VI": R"VI",
            "V": R"V",
            "IV": R"IV",
            "III": R"III",
            "II": R"II",
            "I": R"I",
            "A": R"A",
            "J": R"J",
        }
        return conversion[regular_string]

    def _vector_to_point(self, vec, generators=None):
        """Convert a vector to a torsion point using the provided generators."""
        basis = self._compute_two_torsion() if generators is None else generators
        components = [vec[i].lift_centered() * basis[i] for i in range(4)]
        return components[0] + components[1] + components[2] + components[3]


    def _get_all_RM_two_kernels(self):
        """Return all rank-2 kernels as 2-torsion points."""
        kernels = []
        for subspace in SYMPLECTIC_GL2_SUBSPACES:
            kernel = [self._vector_to_point(subspace.column(i)) for i in range(2)]
            kernels.append(kernel)

        return kernels

    def _compute_isogeny(self, kernel):
        """Compute the 2-isogeny with the given kernel."""
        if self.g2_structure.is_jacobian:
            codomain, isogeny = compute_2_isogeny_from_jacobian(kernel)
        else:
            codomain, isogeny = compute_2_isogeny_from_product(kernel)
        return codomain, isogeny

    @cached_method
    def _compute_neighboring_isogenies(self):
        """Compute all neighboring 2-isogenies."""
        kernels = self._get_all_RM_two_kernels()
        neighbors_with_edges = []
        for kernel in kernels:
            codomain, isogeny = self._compute_isogeny(kernel)
            neighbors_with_edges.append((codomain, isogeny))

        return neighbors_with_edges

    # This must be structured so that the first neighbor returned is the one corresponding to the dual.
    def get_neighbors(self):
        """Return neighboring vertices, deduplicated."""
        neighbors_with_multiplicities = self.get_neighbors_with_multiplicities()
        return list(neighbors_with_multiplicities.keys())

    def get_neighbors_with_multiplicities(self):
        """Return neighboring vertices with multiplicities."""
        neighbors_with_edges = self._compute_neighboring_isogenies()
        codomain_counts = {}
        for codomain, isogeny in neighbors_with_edges:
            vertex = RichelotVertex(codomain)
            if vertex in codomain_counts:
                codomain_counts[vertex] += 1
            else:
                codomain_counts[vertex] = 1
        return codomain_counts

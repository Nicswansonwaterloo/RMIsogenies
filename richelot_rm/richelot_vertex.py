from sage.all import Matrix, GF, VectorSpace

from richelot_rm.genus_two_structures import GenusTwoStructure
from richelot_rm.richelot_product_isogenies import (
    compute_2_isogeny_from_product,
    get_symplectic_two_torsion_prod,
)
from richelot_rm.richelot_jacobian_isogeny import (
    compute_2_isogeny_from_jacobian,
    get_symplectic_two_torsion_jac,
)


class RichelotVertex:
    def __init__(self, g2_structure: GenusTwoStructure, two_torsion_generators=None):
        self.g2_structure = g2_structure
        self.two_torsion_generators = two_torsion_generators
        self.weil_pairing_two_torsion_action = None
        self.invariants = g2_structure.get_isomorphism_class_invariants()
        self.computed_neighbors = None
        

    def _initialize_two_torsion_generators(self):
        if self.two_torsion_generators is None:
            if self.g2_structure.is_product:
                self.two_torsion_generators = get_symplectic_two_torsion_prod(
                    self.g2_structure
                )
            else:
                self.two_torsion_generators = get_symplectic_two_torsion_jac(
                    self.g2_structure
                )
        return self.two_torsion_generators

    def _compute_weil_pairing(self):
        if self.two_torsion_generators is None:
            self._initialize_two_torsion_generators()
        
        Me = Matrix(GF(2), 4, 4)
        for i, P in enumerate(self.two_torsion_generators):
            for j, Q in enumerate(self.two_torsion_generators):
                if i == j:
                    Me[i, j] = 0
                    continue
                entry = P.weil_pairing(Q, 2)
                if entry == 1:
                    Me[i, j] = 0
                else:
                    Me[i, j] = 1
        return Me

    def __repr__(self):
        if self.g2_structure.is_product:
            return f"Product: {self.invariants}"
        return f"Jacobian: {self.invariants}"

    def __eq__(self, other):
        if not isinstance(other, RichelotVertex):
            return False
        return self.invariants == other.invariants

    def __hash__(self):
        return hash(self.invariants)

    # The types given by the Florian and Smith paper: https://eprint.iacr.org/2021/013.pdf
    def get_type(self):
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

    def _vector_to_point(self, vec):
        if self.two_torsion_generators is None:
            self._initialize_two_torsion_generators()
        
        generators = self.two_torsion_generators
        components = [int(vec[i]) * generators[i] for i in range(4)]
        return components[0] + components[1] + components[2] + components[3]

    def _get_maximal_isotropic_subspaces(self):
        if self.weil_pairing_two_torsion_action is None:
            self.weil_pairing_two_torsion_action = self._compute_weil_pairing()
        
        V = VectorSpace(GF(2), 4)
        Me = self.weil_pairing_two_torsion_action
        isotropic_subspaces = []
        for W in V.subspaces(2):
            basis_matrix = W.basis_matrix().transpose()
            if (basis_matrix.transpose() * Me * basis_matrix).is_zero():
                isotropic_subspaces.append(basis_matrix)

        return isotropic_subspaces

    def _get_all_two_kernels(self):
        maximal_isotropic_subspaces = self._get_maximal_isotropic_subspaces()
        kernels = []
        for subspace in maximal_isotropic_subspaces:
            kernel = [self._vector_to_point(subspace.column(i)) for i in range(2)]
            kernels.append(kernel)

        return kernels

    def _compute_isogeny(self, kernel):
        if self.g2_structure.is_jacobian:
            codomain, isogeny = compute_2_isogeny_from_jacobian(kernel)
        else:
            codomain, isogeny = compute_2_isogeny_from_product(kernel)
        return codomain, isogeny

    def _compute_neighboring_isogenies(self):
        if self.computed_neighbors is None:
            kernels = self._get_all_two_kernels()
            neighbors_with_edges = []
            for kernel in kernels:
                codomain, isogeny = self._compute_isogeny(kernel)
                neighbors_with_edges.append((codomain, isogeny))
            self.computed_neighbors = neighbors_with_edges

        return self.computed_neighbors

    # This must be structured so that the first neighbor returned is the one corresponding to the dual
    def get_neighbors(self):
        neighbors_with_multiplicities = self.get_neighbors_with_multiplicities()
        return list(neighbors_with_multiplicities.keys())
    
    def get_neighbors_with_multiplicities(self):
        neighbors_with_edges = self._compute_neighboring_isogenies()
        codomain_counts = {}
        for codomain, isogeny in neighbors_with_edges:
            vertex = RichelotVertex(codomain)
            if vertex in codomain_counts:
                codomain_counts[vertex] += 1
            else:
                codomain_counts[vertex] = 1
        return codomain_counts

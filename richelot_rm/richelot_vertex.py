from sage.all import Matrix, GF, VectorSpace

from richelot_rm.genus_two_structures import GenusTwoStructure
from richelot_rm.richelot_product_isogenies import compute_2_isogeny_from_product, get_symplectic_two_torsion_prod
from richelot_rm.richelot_jacobian_isogeny import compute_2_isogeny_from_jacobian, get_symplectic_two_torsion_jac

class RichelotVertex:
    def __init__(self, g2_structure: GenusTwoStructure, two_torsion_generators = None):
        self.g2_structure = g2_structure
        self.invariants = g2_structure.get_isomorphism_class_invariants()
        
        if two_torsion_generators is None:
            if g2_structure.is_product:
                self.two_torsion_generators = get_symplectic_two_torsion_prod(g2_structure)
            else:
                self.two_torsion_generators = get_symplectic_two_torsion_jac(g2_structure)
        else:
            self.two_torsion_generators = two_torsion_generators

        self.weil_pairing_two_torsion_action = self._compute_weil_pairing()
        self.computed_neighbors = None

    def _compute_weil_pairing(self):
        Me = Matrix(GF(2), 4, 4)
        for i, P in enumerate(self.two_torsion_generators):
            for j, Q in enumerate(self.two_torsion_generators):
                entry = P.weil_pairing(Q, 2)
                if entry == 1:
                    Me[i, j] = 0
                else:
                    Me[i, j] = 1
        return Me

    def __repr__(self):
        if self.g2_structure.is_product:
            return f"Product: {self.invariants }"
        return f"Jacobian: {self.invariants}"

    def __eq__(self, other):
        if not isinstance(other, RichelotVertex):
            return False
        return self.invariants == other.invariants

    def __hash__(self):
        return hash(self.invariants)
    
    def get_type(self):
        if self.g2_structure.is_product:
            return R"\Pi" if self.invariants[0] != self.invariants[1] else R"\Sigma"
        elif self.g2_structure.is_jacobian:
            return R"J"
        return "Error"

    def _vector_to_point(self, vec):
        generators = self.two_torsion_generators
        components = [int(vec[i]) * generators[i] for i in range(4)]
        return components[0] + components[1] + components[2] + components[3]
    
    def _get_maximal_isotropic_subspaces(self):
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
        if self.g2_structure.is_jacobian():
            codomain, isogeny = compute_2_isogeny_from_jacobian(kernel)
        else:
            codomain, isogeny = compute_2_isogeny_from_product(kernel)
        return codomain, isogeny

    def _compute_all_neighbors(self):
        if self.computed_neighbors is None:
            kernels = self._get_all_two_kernels()
            neighbors_with_edges = []
            for kernel in kernels:
                codomain, isogeny = self._compute_isogeny(kernel)
                neighbors_with_edges.append((codomain, isogeny))
            self.computed_neighbors = neighbors_with_edges
        
        return self.computed_neighbors
    
    def get_neighbors(self):
        neighbors_with_edges = self._compute_all_neighbors()
        codomains = [neighbor for neighbor, edge in neighbors_with_edges]
        neighbors = [RichelotVertex(neighbor) for neighbor in codomains]
        neighbors = list(set(neighbors))  # Remove duplicates
        return neighbors
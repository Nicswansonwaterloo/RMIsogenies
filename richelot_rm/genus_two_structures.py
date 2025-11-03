from sage.all import HyperellipticCurve, GF
from sage.schemes.hyperelliptic_curves.invariants import absolute_igusa_invariants_kohel

class GenusTwoStructure:
    """
    A base class for representing an abelian surface.
    """

    def __init__(self, is_product):
        self.is_product = is_product

    def is_product(self):
        return self.is_product

    def is_jacobian(self):
        return not self.is_product

class GenusTwoProductStructure(GenusTwoStructure):
    """
    A base class for representing an abelian surface which is a product of two elliptic curves.
    """

    def __init__(self, E1, E2):
        self.E1 = E1
        self.E2 = E2
        super().__init__(is_product=True)

    def __repr__(self):
        return f"Product: ( {self.E1.a_invariants()}, {self.E2.a_invariants()} )"

    def __eq__(self, other):
        if not isinstance(other, GenusTwoProductStructure):
            return False
        return self.E1 == other.E1 and self.E2 == other.E2
    
    def get_isomorphism_class_invariants(self):
        j1 = self.E1.j_invariant()
        j2 = self.E2.j_invariant()
        return tuple(sorted((j1, j2)))
    
    def is_isomorphic_to(self, other):
        if not isinstance(other, GenusTwoProductStructure):
            return False
        j1 = self.E1.j_invariant()
        j2 = self.E2.j_invariant()
        j1_other = other.E1.j_invariant()
        j2_other = other.E2.j_invariant()
        if j1 == j1_other and j2 == j2_other:
            return True
        if j1 == j2_other and j2 == j1_other:
            return True
        return False

    def __getitem__(self, i):
        # Operator to get self[i].
        if i == 0:
            return self.E1
        elif i == 1:
            return self.E2
        else:
            raise IndexError("Index {} is out of range.".format(i))

class GenusTwoJacobianStructure(GenusTwoStructure):
    """
    A base class for representing an abelian surface which is a jacobian of a genus two curve.
    The representation uses richelot isogenies, so this requires a hyperelliptic curve.
    """
    def __init__(self, h):
        self.Rx = h.parent()
        self.x = self.Rx.gen()
        self.h = h
        # h will be monic
        self.h = self.h.monic()
        self.H = HyperellipticCurve(self.h)
        self.jac = self.H.jacobian()
        super().__init__(is_product=False)

    def __repr__(self):
        return f"{self.h}"
    
    def __eq__(self, other):
        if not isinstance(other, GenusTwoJacobianStructure):
            return False
        return self.h == other.h
    
    def get_isomorphism_class_invariants(self):
        return absolute_igusa_invariants_kohel(self.h)
    
    def is_isomorphic_to(self, other):
        if not isinstance(other, GenusTwoJacobianStructure):
            return False
        return self.get_isomorphism_class_invariants() == other.get_isomorphism_class_invariants()
    
    # Points will be divisors on the jacobian. Handled by Sage internally.
    def __call__(self, *args):
        if len(args) != 1:
            raise ValueError(f"Invalid arguments {args} to create a point on {self}.")
        D = self.jac(args[0])
        return D
    
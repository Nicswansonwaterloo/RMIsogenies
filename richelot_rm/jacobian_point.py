from sage.schemes.hyperelliptic_curves.jacobian_morphism import JacobianMorphism_divisor_class_field
from richelot_rm.genus_two_structures import GenusTwoJacobianStructure

"""
A Wrapper for points that live on a genus_two_jacobian_structure
"""
class JacobianPoint:
    """
    A helper class which represents an element on the Jacobian of a genus two curve.
    This will mostly be handled internally by sage, but we wrap it here to have a consistent interface with ProductPoint.
    """
    def __init__(self, D: JacobianMorphism_divisor_class_field):
        self.D = D
        self.D1 = D[0]
        self.D2 = D[1]

    def __repr__(self):
        return repr(self.D)

    def parent(self):
        h, f = self.D.scheme().curve().hyperelliptic_polynomials()
        if f != 0:
            raise ValueError("The divisor's parent is not of the correct form.")
        return GenusTwoJacobianStructure(h)

    def order(self):
        return self.D.order()

    def __getitem__(self, i):
        # Operator to get self[i].
        if i == 0:
            return self.D1
        elif i == 1:
            return self.D2
        else:
            raise IndexError("Index {} is out of range.".format(i))

    def __setitem__(self, i, Px):
        # Operator to set self[i]=P.
        if i == 0:
            self.D1 = Px
        elif i == 1:
            self.D2 = Px
        else:
            raise IndexError("Index {} is out of range.".format(i))

    def __eq__(self, other):
        if other == 0:
            return self.D == 0
        if not isinstance(other, JacobianPoint):
            return False
        return self.D == other.D

    def __add__(self, other):
        if not isinstance(other, JacobianPoint):
            raise TypeError(f"Cannot add Jacobian point to non-Jacobian point. {type(other)}")
        return JacobianPoint(self.D + other.D)

    def __sub__(self, other):
        if not isinstance(other, JacobianPoint):
            raise TypeError(f"Cannot subtract Jacobian point from non-Jacobian point. {type(other)}")
        return JacobianPoint(self.D - other.D)

    def __neg__(self):
        return JacobianPoint(-self.D)

    def __mul__(self, m):
        return JacobianPoint(m * self.D)

    def __rmul__(self, m):
        return JacobianPoint(m * self.D)
    
    def __hash__(self):
        return hash((hash(self.D1), hash(self.D2)))

    def weil_pairing(self, other, n):
        if not isinstance(other, JacobianPoint):
            raise TypeError("Both inputs must be jacobian points")
        if n != 2:
            raise NotImplementedError("Weil pairing is only implemented for n=2")
        if self.D2 != 0 or other.D2 != 0:
            raise ValueError("Trying to call Weil pairing n = 2 on non 2-torsion points")
        
        D1_prime = other.D1
        if self.D1.gcd(D1_prime) == 1:
            return 1
        else:
            return -1
from sage.all import ZZ, lcm
from sage.schemes.elliptic_curves.ell_point import EllipticCurvePoint
from richelot_rm.genus_two_structures import GenusTwoProductStructure

"""
Borrowed and modified from https://github.com/ThetaIsogenies/two-isogenies
"""
class ProductPoint:
    """
    A helper class which represents an element P = (P1, P2) in E1 x E2
    and allows us to compute certain useful functions, such as adding,
    doubling or computing the Weil pairing of e(P,Q) for P,Q in E1 x E2
    """

    def __init__(self, P1, P2):
        if not isinstance(P1, EllipticCurvePoint) or not isinstance(P2, EllipticCurvePoint):
            raise ValueError(f"Points must be points on an elliptic curve: {P1}, {P2}")
        self.P1 = P1
        self.P2 = P2

    def __repr__(self):
        return "[ {} , {} ]".format(self.P1, self.P2)

    def parent(self):
        return GenusTwoProductStructure(*self.curves())

    def curves(self):
        return (self.P1.curve(), self.P2.curve())

    def points(self):
        return self.P1, self.P2

    def order(self):
        return lcm(self.P1.order(), self.P2.order())

    def double(self):
        """
        Computes [2] P = ([2] P1, [2] P2)
        """
        return ZZ(2) * self

    def __getitem__(self, i):
        # Operator to get self[i].
        if i == 0:
            return self.P1
        elif i == 1:
            return self.P2
        else:
            raise IndexError("Index {} is out of range.".format(i))

    def __setitem__(self, i, P):
        # Operator to set self[i]=P.
        if i == 0:
            self.P1 = P
        elif i == 1:
            self.P2 = P
        else:
            raise IndexError("Index {} is out of range.".format(i))

    def __eq__(self, other):
        if other == 0:
            E1, E2 = self.parent()
            return self.P1 == E1(0) and self.P2 == E2(0)
        
        return self.P1 == other.P1 and self.P2 == other.P2

    def __add__(self, other):
        return ProductPoint(self.P1 + other.P1, self.P2 + other.P2)

    def __sub__(self, other):
        return ProductPoint(self.P1 - other.P1, self.P2 - other.P2)

    def __neg__(self):
        return ProductPoint(-self.P1, -self.P2)

    def __mul__(self, m):
        """
        Compute [m] P = ([m] P1, [m] P2)
        """
        # When the scalar is a python int, then
        # sagemath does multiplication naively, when
        # the scalar in a Sage type, it instead calls
        # _acted_upon_, which calls pari, which is fast
        m = ZZ(m)
        return ProductPoint(m * self.P1, m * self.P2)

    def __rmul__(self, m):
        return self * m
    
    def __hash__(self):
        return hash((hash(self.P1), hash(self.P2)))

    def weil_pairing(self, other, n):
        """
        The Weil pairing e_n(P, Q) for P = (P1, P2) and Q = (Q1, Q2)
        is defined as

            e_n(P, Q) = e_n(P1, Q1) * e_n(P2, Q2)
        """
        if not isinstance(other, ProductPoint):
            raise TypeError("Both inputs must be product points")

        P1, P2 = self.points()
        Q1, Q2 = other.points()


        ePQ1 = P1.weil_pairing(Q1, n)
        ePQ2 = P2.weil_pairing(Q2, n)
        
        Fp2 = P1.base_ring()
        return Fp2(ePQ1 * ePQ2)

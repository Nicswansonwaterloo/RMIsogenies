from sage.all import (
    HyperellipticCurve,
    get_coercion_model,
    PolynomialRing,
    cached_method,
)
from sage.schemes.hyperelliptic_curves.invariants import absolute_igusa_invariants_kohel
from theta_rm.couple_point import CouplePoint
from theta_rm.theta_point import ThetaPoint


from vendors.Theta_SageMath.theta_structures.dimension_one import (
    montgomery_curve_to_theta_null_point,
    montgomery_point_to_theta_point,
    montgomery_torsion_to_theta_null_point,
)

cm = get_coercion_model()


class ThetaStructure:
    """
    Class for the ThetaStructure, defined by it's theta null point. This type
    represents the generic domain/codomain of the (2,2)-isogeny in the theta
    model.
    """

    def __init__(self, null_point):
        if not len(null_point) == 4:
            raise ValueError

        self._base_ring = cm.common_parent(*(c.parent() for c in null_point))
        self._point = ThetaPoint
        self._precomputation = None

        self._null_point = self._point(self, null_point)

    def null_point(self):
        """
        Return the null point of the given theta structure
        """
        return self._null_point

    def is_product(self):
        return isinstance(self, ProductThetaStructure)

    def base_ring(self):
        """
        Return the base ring of the common parent of the coordinates of the null point
        """
        return self._base_ring

    def zero(self):
        """
        The additive identity is the theta null point
        """
        return self.null_point()

    def __repr__(self):
        return f"Theta structure over {self.base_ring()} with null point: {self.null_point()}"

    def coords(self):
        """
        Return the coordinates of the theta null point of the theta structure
        """
        return self.null_point().coords()

    def hadamard(self):
        """
        Compute the Hadamard transformation of the theta null point of the theta structure
        """
        return self.null_point().hadamard()

    def squared_theta(self):
        """
        Square the coefficients and then compute the Hadamard transformation of
        the theta null point of the theta structure
        """
        return self.null_point().squared_theta()

    def _arithmetic_precomputation(self):
        """
        Precompute 6 field elements used in arithmetic and isogeny computations
        """
        if self._precomputation is None:
            a, b, c, d = self.null_point().coords()

            # Technically this computes 4A^2, 4B^2, ...
            # but as we take quotients this doesnt matter
            # Cost: 4S
            AA, BB, CC, DD = self.squared_theta()

            # Precomputed constants for addition and doubling
            b_inv, c_inv, d_inv, BB_inv, CC_inv, DD_inv = batched_inversion(
                b, c, d, BB, CC, DD
            )

            y0 = a * b_inv
            z0 = a * c_inv
            t0 = a * d_inv

            Y0 = AA * BB_inv
            Z0 = AA * CC_inv
            T0 = AA * DD_inv

            self._precomputation = (y0, z0, t0, Y0, Z0, T0)
        return self._precomputation

    @cached_method
    def rosenhain_from_theta(self):
        """
        From a theta null point structure, return the Rosenhain invariants

        NOTE: used for debugging
        """
        # Extract out the hadamard transform from the point class
        to_hadamard = self.zero().to_hadamard

        a, b, c, d = self.coords()
        A, C, B, D = to_hadamard(a, d, b, c)  # beware the weird order
        if A * B * C * D == 0:
            a, b, c, d = to_hadamard(a, b, c, d)
            A, C, B, D = to_hadamard(a, d, b, c)

        try:
            al = a * a + b * b + c * c + d * d
            be = 2 * (a * d + b * c)
            ga = 2 * (a * b + c * d)
            de = 2 * (a * c + b * d)
            CD_AB = C * D / (A * B)
            ep_phi = (1 + CD_AB) / (1 - CD_AB)
            lam = al * ga / (de * be)
            mu = ga / de * ep_phi
            nu = al / be * ep_phi

        except Exception as e:
            raise ValueError(f"Conversion to rosenhain failed because of: {e}")

        return lam, mu, nu

    @cached_method
    def hyperelliptic_from_theta(self):
        """
        Convert a theta null point structure to an hyperelliptic curve
        """

        lam, mu, nu = self.rosenhain_from_theta()
        if lam is None:
            raise ValueError("Could not compute Rosenhain roots from the null point")

        R = PolynomialRing(self.base_ring(), name="x")
        x = R.gens()[0]

        f_poly = x * (x - 1) * (x - lam) * (x - mu) * (x - nu)

        try:
            _ = HyperellipticCurve(f_poly)
        except:
            raise ValueError("Converted curve is not a hyperelliptic curve")
        return f_poly

    def __call__(self, coords):
        return self._point(self, coords)

    def __eq__(self, other):
        if not isinstance(other, ThetaStructure):
            return False
        return self.null_point() == other.null_point()

    def get_isomorphism_class_invariants(self):
        return absolute_igusa_invariants_kohel(self.hyperelliptic_from_theta())

    def is_isomorphic_to(self, other):
        if not isinstance(other, ThetaStructure):
            return False
        if self.is_product() != other.is_product():
            return False
        return (
            self.get_isomorphism_class_invariants()
            == other.get_isomorphism_class_invariants()
        )



class ProductThetaStructure(ThetaStructure):
    """
    A base class for representing an abelian surface which is a product of two elliptic curves.
    """

    def __init__(self, E1, E2, B1=None, B2=None):
        self.E1 = E1
        self.E2 = E2
        if B1 is None or B2 is None:
            null_point, O01, O02 = self.montgomery_curves_to_theta_null_point(E1, E2)
        else:
            O01 = montgomery_torsion_to_theta_null_point(B1[0])
            O02 = montgomery_torsion_to_theta_null_point(B2[0])
            null_point = self.pairwise_product(O01, O02)

        self.O01 = O01
        self.O02 = O02
        super().__init__(null_point)

    def __repr__(self):
        return f"Pair of Elliptic curves with a invariants: ({self.E1.a_invariants()}, {self.E2.a_invariants()} )"

    def __eq__(self, other):
        if not isinstance(other, ProductThetaStructure):
            return False
        return self.E1 == other.E1 and self.E2 == other.E2

    def __call__(self, *args):
        """
        Create a theta point on the product from either a pair of points,
        a CouplePoint or coefficients of a theta point.
        """
        if len(args) == 2:
            P1, P2 = args
            coords = self.montgomery_points_to_theta_point(P1, P2)
            return self._point(self, coords)

        if len(args) == 1:
            input_value = args[0]
            if isinstance(input_value, CouplePoint):
                P1, P2 = input_value.points()
                coords = self.montgomery_points_to_theta_point(P1, P2)
                return self._point(self, coords)

        raise ValueError(f"Invalid arguments {args} to create a point on {self}.")

    def get_isomorphism_class_invariants(self):
        j1 = self.E1.j_invariant()
        j2 = self.E2.j_invariant()
        return tuple(sorted((j1, j2)))

    def is_isomorphic_to(self, other):
        if not isinstance(other, ProductThetaStructure):
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

    def theta_null_point(self):
        return self.null_point()

    def base_ring(self):
        return self.base_ring()

    def pairwise_product(self, O1, O2):
        """
        Given two dim-1 theta points, compute the level two theta point on the
        product.

        NOTE: The order used is implicitly 00, 10, 01, 11
        """
        a1, b1 = O1
        a2, b2 = O2

        return (a1 * a2, b1 * a2, a1 * b2, b1 * b2)

    def montgomery_curves_to_theta_null_point(self, E1, E2):
        """
        Given a pair of elliptic curves, compute the corresponding
        dim-1 theta points as well as the coordinates of the
        level-2 null point.
        """
        O1 = montgomery_curve_to_theta_null_point(E1)
        O2 = montgomery_curve_to_theta_null_point(E2)

        null_coords = self.pairwise_product(O1, O2)

        return null_coords, O1, O2

    def montgomery_points_to_theta_point(self, P1, P2):
        """
        Given two points P1, P2 on curves E1, E2 where it is
        assumed we have computed the corresponding dim-1
        theta null points, compute the product theta point
        from P1 x P2.
        """
        O1 = montgomery_point_to_theta_point(self.O01, P1)
        O2 = montgomery_point_to_theta_point(self.O02, P2)

        return self.pairwise_product(O1, O2)


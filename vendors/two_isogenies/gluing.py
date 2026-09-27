from sage.all import Matrix

from theta_rm.product import ProductPoint as CouplePoint
from theta_rm.theta import INDICES, ThetaSurface, hadamard, normalise

from .batched_inversion import batched_inversion


class GluingThetaIsogeny:
    """
    Compute the gluing isogeny from E1 x E2 (Elliptic Product) -> A (Theta Model)

    Expected input:

    - (K1_8, K2_8) The 8-torsion above the kernel generating the isogeny
    - M (Optional) a base change matrix, if this is not including, it can
      be derived from [2](K1_8, K2_8)
    """

    def __init__(self, K1_8, K2_8, M=None):
        # Double points to get four-torsion, we always need one of these, used
        # for the image computations but we'll need both if we wish to derived
        # the base change matrix as well
        K1_4 = K1_8.double()
        self.S_shift = K2_8.double()

        # If M is not included, compute the matrix on the fly from the four
        # torsion.
        if M is None:
            K2_4 = K2_8.double()
            M = self.get_base_change_matrix(K1_4, K2_4)

        # Initalise self
        self._base_change_matrix = M
        self.T_shift = K1_4
        self._precomputation = None
        self._zero_idx = INDICES[0]

        # Map points from elliptic product onto the product theta structure
        # using the base change matrix
        T1_8 = self.base_change(K1_8)
        T2_8 = self.base_change(K2_8)

        # Compute the codomain of the gluing isogeny
        self._codomain = self._special_compute_codomain(T1_8, T2_8)

    @staticmethod
    def get_base_change_matrix(T1, T2):
        """
        Given the four torsion above the kernel generating the gluing isogeny,
        compute the matrix M which allows us to map points on an elliptic
        product to the compatible theta structure.
        """

        def get_matrix(T):
            """
            Compute the matrix [a, b]
                               [c, d]
            From the (X : Z) coordinates of a point
            T and its double TT = (T + T)

            NOTE: Assumes that Z = 1 or 0
            """
            TT = T + T
            x = T[0]
            u = TT[0]

            # Our Z-coordinates are always 1 or 0
            # as Sage uses affine addition
            assert T[2] == 1
            assert TT[2] == 1

            # Compute the matrix coefficents
            det = x - u
            inv_det = 1 / det
            a = -u * inv_det
            b = -inv_det
            c = x * (u - det) * inv_det
            d = -a
            return Matrix(2, 2, [a, b, c, d])

        # Extract elliptic curve points from CouplePoints
        P1, P2 = T1.points()
        Q1, Q2 = T2.points()

        # Compute matrices from points
        # TODO: totally overkill, but if these were all computed together
        #       you could batch all 4 inversions into one.
        g1 = get_matrix(P1)
        g2 = get_matrix(P2)
        h1 = get_matrix(Q1)
        h2 = get_matrix(Q2)

        # the matrices gi, hi does not commute, but g1 \tens g2 should commute with h1 \tens h2
        gh1 = g1 * h1
        gh2 = g2 * h2

        # Access Coefficients once
        # Notice some coeffs are never used
        g00_1, g01_1, g10_1, g11_1 = g1.list()
        g00_2, _, g10_2, _ = g2.list()
        h00_1, _, h10_1, _ = h1.list()
        h00_2, h01_2, h10_2, h11_2 = h2.list()

        # First row from the product of gi * hi
        gh00_1, _, gh10_1, _ = gh1.list()
        gh00_2, _, gh10_2, _ = gh2.list()

        # start the trace with id
        a = 1
        b = 0
        c = 0
        d = 0

        # T1
        a += g00_1 * g00_2
        b += g00_1 * g10_2
        c += g10_1 * g00_2
        d += g10_1 * g10_2

        # T2
        a += h00_1 * h00_2
        b += h00_1 * h10_2
        c += h10_1 * h00_2
        d += h10_1 * h10_2

        # T1+T2
        a += gh00_1 * gh00_2
        b += gh00_1 * gh10_2
        c += gh10_1 * gh00_2
        d += gh10_1 * gh10_2

        if not any((a, b, c, d)):
            from sage.all import identity_matrix

            projector = (
                identity_matrix(g1.base_ring(), 4)
                + g1.tensor_product(g2)
                + h1.tensor_product(h2)
                + gh1.tensor_product(gh2)
            )
            a, b, c, d = next(column for column in projector.columns() if any(column))

        # Now we act by (0, Q2)
        a1 = h00_2 * a + h01_2 * b
        b1 = h10_2 * a + h11_2 * b
        c1 = h00_2 * c + h01_2 * d
        d1 = h10_2 * c + h11_2 * d

        # Now we act by (P1, 0)
        a2 = g00_1 * a + g01_1 * c
        b2 = g00_1 * b + g01_1 * d
        c2 = g10_1 * a + g11_1 * c
        d2 = g10_1 * b + g11_1 * d

        # Now we act by (P1, Q2)
        a3 = g00_1 * a1 + g01_1 * c1
        b3 = g00_1 * b1 + g01_1 * d1
        c3 = g10_1 * a1 + g11_1 * c1
        d3 = g10_1 * b1 + g11_1 * d1

        return Matrix(
            [[a, b, c, d], [a1, b1, c1, d1], [a2, b2, c2, d2], [a3, b3, c3, d3]]
        )

    def apply_base_change(self, coords):
        """
        Apply the basis change by acting with matrix multiplication, treating
        the coordinates as a vector
        """
        N = self._base_change_matrix
        x, y, z, t = coords
        X = N[0, 0] * x + N[0, 1] * y + N[0, 2] * z + N[0, 3] * t
        Y = N[1, 0] * x + N[1, 1] * y + N[1, 2] * z + N[1, 3] * t
        Z = N[2, 0] * x + N[2, 1] * y + N[2, 2] * z + N[2, 3] * t
        T = N[3, 0] * x + N[3, 1] * y + N[3, 2] * z + N[3, 3] * t

        return (X, Y, Z, T)

    def base_change(self, P):
        """
        Compute the basis change on a CouplePoint to recover a ThetaPoint of
        compatible form
        """
        if not isinstance(P, CouplePoint):
            raise TypeError("Function assumes that the input is of type `CouplePoint`")

        # extract X,Z coordinates on pairs of points
        P1, P2 = P.points()
        X1, Z1 = P1[0], P1[2]
        X2, Z2 = P2[0], P2[2]

        # Correct in the case of (0 : 0)
        if X1 == 0 and Z1 == 0:
            X1 = 1
            Z1 = 0
        if X2 == 0 and Z2 == 0:
            X2 = 1
            Z2 = 0

        # Apply the basis transformation on the product
        coords = self.apply_base_change([X1 * X2, X1 * Z2, Z1 * X2, Z1 * Z2])
        return coords

    def _special_compute_codomain(self, T1, T2):
        """
        Given two isotropic points of 8-torsion T1 and T2, compatible with
        the theta null point, compute the level two theta null point A/K_2
        """
        xAxByCyD = hadamard([x * x for x in T1])
        zAtBzYtD = hadamard([x * x for x in T2])

        # Find the value of the non-zero index
        zero_idx = next((i for i, x in enumerate(xAxByCyD) if x == 0), None)
        self._zero_idx = INDICES[zero_idx]

        # Dumb check to make sure everything is OK
        assert xAxByCyD[zero_idx] == zAtBzYtD[zero_idx] == 0

        # Initialize lists
        # The zero index described the permutation
        ABCD = [0 for _ in range(4)]
        precomp = [0 for _ in range(4)]

        i0, i1, i2, i3 = (INDICES.index(k + self._zero_idx) for k in INDICES)

        # Compute non-trivial numerators (Others are either 1 or 0)
        num_1 = zAtBzYtD[i1]
        num_2 = xAxByCyD[i2]
        num_3 = zAtBzYtD[i3]
        num_4 = xAxByCyD[i3]

        # Compute and invert non-trivial denominators
        den_1, den_2, den_3, den_4 = batched_inversion(num_1, num_2, num_3, num_4)

        # Compute A, B, C, D
        ABCD[i0] = 0
        ABCD[i1] = num_1 * den_3
        ABCD[i2] = num_2 * den_4
        ABCD[i3] = 1

        # Compute precomputation for isogeny images
        precomp[i0] = 0
        precomp[i1] = den_1 * num_3
        precomp[i2] = den_2 * num_4
        precomp[i3] = 1
        self._precomputation = precomp

        # Final Hadamard of the above coordinates
        a, b, c, d = hadamard(ABCD)

        return ThetaSurface([a, b, c, d])

    def special_image(self, P, translate):
        """
        When the domain is a non product theta structure on a product of
        elliptic curves, we will have one of A,B,C,D=0, so the image is more
        difficult. We need to give the coordinates of P but also of
        P+Ti, Ti one of the point of 4-torsion used in the isogeny
        normalisation
        """
        AxByCzDt = hadamard([x * x for x in P])

        # We are in the case where at most one of A, B, C, D is
        # zero, so we need to account for this
        #
        # To recover values, we use the translated point to get
        AyBxCtDz = hadamard([x * x for x in translate])

        i0, i1, i2, i3 = (INDICES.index(k + self._zero_idx) for k in INDICES)

        # Directly compute y,z,t
        y = AxByCzDt[i1] * self._precomputation[i1]
        z = AxByCzDt[i2] * self._precomputation[i2]
        t = AxByCzDt[i3]

        if z == 0 and t == 0:
            return None

        # We can compute x from the translation
        # First we need a normalisation
        if z != 0:
            zb = AyBxCtDz[i3]
            lam = z / zb
        else:
            tb = AyBxCtDz[i2] * self._precomputation[i2]
            lam = t / tb

        # Finally we recover x
        xb = AyBxCtDz[i1] * self._precomputation[i1]
        x = xb * lam

        xyzt = [0 for _ in range(4)]
        xyzt[i0] = x
        xyzt[i1] = y
        xyzt[i2] = z
        xyzt[i3] = t

        image = hadamard(xyzt)
        return normalise(image)

    def _image(self, P):
        """
        Take into input the theta null point of A/K_2, and return the image
        of the point by the isogeny
        """
        if not isinstance(P, CouplePoint):
            raise TypeError(
                "Isogeny image for the gluing isogeny is defined to act on CouplePoints"
            )

        # Compute sum of points on elliptic curve
        P_sum_T = P + self.T_shift

        # Push both the point and the translation through the
        # completion
        iso_P = self.base_change(P)
        iso_P_sum_T = self.base_change(P_sum_T)

        return self.special_image(iso_P, iso_P_sum_T)

    def __call__(self, P):
        image = self._image(P)
        if image is not None:
            return image
        node = self._image(self.S_shift)
        a, b = next(v for v, T in self._codomain.two_torsion().items() if T == node)
        image = self._image(P + self.S_shift)
        assert image is not None
        return normalise(
            [
                (-1) ** int(a.dot_product(k)) * image[INDICES.index(k + b)]
                for k in INDICES
            ]
        )

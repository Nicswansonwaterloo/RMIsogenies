from sage.all import (
    Integer,
)
from sage.structure.element import get_coercion_model, RingElement
# from theta_rm.genus_two_structures import ThetaStructure # Removed due to circular import, but we can add back if needed

cm = get_coercion_model()


class ThetaPoint:
    """
    A Theta Point in the level-2 Theta Structure is defined with four projective
    coordinates

    We cannot perform arbitrary arithmetic, but we can compute doubles and
    differential addition, which like x-only points on the Kummer line, allows
    for scalar multiplication
    """

    def __init__(self, parent, coords):
        # if not isinstance(parent, ThetaStructure):
        #     raise ValueError
        # Removed due to circular import, but we can add back if needed

        self._parent = parent
        self._coords = tuple(coords)

        self._hadamard = None
        self._squared_theta = None

    def parent(self):
        """
        Return the parent of the element, of type ThetaStructure
        """
        return self._parent

    def theta(self):
        """
        Return the parent theta structure of this ThetaPoint"""
        return self.parent()

    def coords(self):
        """
        Return the projective coordinates of the ThetaPoint
        """
        return self._coords

    def is_zero(self):
        """
        An element is zero if it is equivalent to the null point of the parent
        ThetaStrcuture
        """
        return self == self.parent().zero()

    @staticmethod
    def to_hadamard(x_00, x_10, x_01, x_11):
        """
        Compute the Hadamard transformation of four coordinates, using recursive
        formula.
        """
        x_00, x_10 = (x_00 + x_10, x_00 - x_10)
        x_01, x_11 = (x_01 + x_11, x_01 - x_11)
        return x_00 + x_01, x_10 + x_11, x_00 - x_01, x_10 - x_11

    def hadamard(self):
        """
        Compute the Hadamard transformation of this element
        """
        if self._hadamard is None:
            self._hadamard = self.to_hadamard(*self.coords())
        return self._hadamard

    @staticmethod
    def to_squared_theta(x, y, z, t):
        """
        Square the coordinates and then compute the Hadamard transform of the
        input
        """
        return ThetaPoint.to_hadamard(x * x, y * y, z * z, t * t)

    def squared_theta(self):
        """
        Compute the Squared Theta transformation of this element
        which is the square operator followed by Hadamard.
        """
        if self._squared_theta is None:
            self._squared_theta = self.to_squared_theta(*self.coords())
        return self._squared_theta

    def double(self):
        """
        Computes [2]*self

        NOTE: Assumes that no coordinate is zero

        Cost: 8S 6M
        """
        # If a,b,c,d = 0, then the codomain of A->A/K_2 is a product of
        # elliptic curves with a non product theta structure.
        # Unless we are very unlucky, A/K_1 will not be in this case, so we
        # just need to Hadamard, double, and Hadamard inverse
        # If A,B,C,D=0 then the domain itself is a product of elliptic
        # curves with a non product theta structure. The Hadamard transform
        # will not change this, we need a symplectic change of variable
        # that puts us back in a product theta structure
        y0, z0, t0, Y0, Z0, T0 = self.parent()._arithmetic_precomputation()

        # Temp coordinates
        # Cost 8S 3M
        xp, yp, zp, tp = self.squared_theta()
        xp = xp**2
        yp = Y0 * yp**2
        zp = Z0 * zp**2
        tp = T0 * tp**2

        # Final coordinates
        # Cost 3M
        X, Y, Z, T = self.to_hadamard(xp, yp, zp, tp)
        X = X
        Y = y0 * Y
        Z = z0 * Z
        T = t0 * T

        coords = (X, Y, Z, T)
        return self._parent(coords)

    def diff_addition(P, Q, PQ):
        """
        Given the theta points of P, Q and P-Q computes the theta point of
        P + Q.

        NOTE: Assumes that no coordinate is zero

        Cost: 8S 17M
        """
        # Extract out the precomputations
        Y0, Z0, T0 = P.parent()._arithmetic_precomputation()[-3:]

        # Transform with the Hadamard matrix and multiply
        # Cost: 8S 7M
        p1, p2, p3, p4 = P.squared_theta()
        q1, q2, q3, q4 = Q.squared_theta()

        xp = p1 * q1
        yp = Y0 * p2 * q2
        zp = Z0 * p3 * q3
        tp = T0 * p4 * q4

        # Final coordinates
        PQx, PQy, PQz, PQt = PQ.coords()

        # Note:
        # We replace the four divisions by
        # PQx, PQy, PQz, PQt by 10 multiplications
        # Cost: 10M
        PQxy = PQx * PQy
        PQzt = PQz * PQt

        X, Y, Z, T = P.to_hadamard(xp, yp, zp, tp)
        X = X * PQzt * PQy
        Y = Y * PQzt * PQx
        Z = Z * PQxy * PQt
        T = T * PQxy * PQz

        coords = (X, Y, Z, T)
        return P.parent()(coords)
    
    def diff_add_ladder(P, Q, PQ, m):
            """
            Given the theta points of Q and self - Q (PQ), computes the theta point of
            self + [m]Q using a 3-point Montgomery ladder.

            NOTE: Assumes that no coordinate is zero at any point during the doubling.

            Cost: |m|_bits * (1 diff_addition + 1 double)
            """
            # Ensure m is an integer
            if not isinstance(m, (int, Integer)):
                try:
                    m = Integer(m)
                except Exception:
                    raise TypeError(f"Cannot coerce input scalar {m = } to an integer")

            if m == 0:
                return P

            # Handle negative scalars:
            # self - |m|Q is equivalent to self + |m|(-Q). 
            # The required difference point is self - (-Q) = self + Q.
            # We can compute this initial difference using one standard diff_addition.
            if m < 0:
                PQ = P.diff_addition(Q, PQ)
                m = abs(m)

            A = Q
            B = P
            C = PQ

            # Process bits from Least Significant Bit (LSB) to Most Significant Bit (MSB)
            # bin(m)[2:] gives the binary string; reversed() lets us read LSB first.
            for bit in reversed(bin(m)[2:]):
                if bit == '1':
                    # B_new = A + B. The difference between A and B is maintained in C.
                    B = A.diff_addition(B, C)
                else:
                    # C_new = A - C. By passing B (which is A + C), diff_addition returns the difference.
                    C = A.diff_addition(C, B)
                    
                A = A.double()

            return B
        
    def extended_addition(P, Q, R, PQ, QR, PR):
            """
            Given the theta points of Q, R, P+Q (PQ), Q+R (QR), and P+R (PR),
            computes the theta point of P + Q + R using 3-way addition Riemann relations.
            
            NOTE: Assumes self is P.
            """
            # Extract precomputations
            Y0, Z0, T0 = P.parent()._arithmetic_precomputation()[-3:]

            # Step 1: Compute the squared theta coordinates of the partial sums
            # We need the Hadamard transforms of the pairs to build the Riemann system
            pq1, pq2, pq3, pq4 = PQ.squared_theta()
            r1, r2, r3, r4 = R.squared_theta()
            
            pr1, pr2, pr3, pr4 = PR.squared_theta()
            q1, q2, q3, q4 = Q.squared_theta()
            
            qr1, qr2, qr3, qr4 = QR.squared_theta()
            p1, p2, p3, p4 = P.squared_theta()

            # Step 2: Cross-multiply the squared coordinates
            # According to the Riemann relations for 3-way addition, the products of the 
            # theta constants of the sums and the base points are algebraically linked to P+Q+R.
            # We construct the unscaled coordinates of the target point.
            
            # Note: We use the symmetric combinations of (P+Q, R), (P+R, Q), and (Q+R, P)
            X_sym = (pq1 * r1) + (pr1 * q1) + (qr1 * p1)
            Y_sym = Y0 * ((pq2 * r2) + (pr2 * q2) + (qr2 * p2))
            Z_sym = Z0 * ((pq3 * r3) + (pr3 * q3) + (qr3 * p3))
            T_sym = T0 * ((pq4 * r4) + (pr4 * q4) + (qr4 * p4))

            # Step 3: Apply the Hadamard transform to map back to standard coordinates
            X, Y, Z, T = P.to_hadamard(X_sym, Y_sym, Z_sym, T_sym)

            # Step 4: Scale by the base point components to normalize the projective equivalence
            # (This avoids a division by zero if we were to solve the quadratic exactly)
            Px, Py, Pz, Pt = P.coords()
            Qx, Qy, Qz, Qt = Q.coords()
            Rx, Ry, Rz, Rt = R.coords()

            # Scale by the product of the base coordinates to balance the projective weights
            scale_x = Px * Qx * Rx
            scale_y = Py * Qy * Ry
            scale_z = Pz * Qz * Rz
            scale_t = Pt * Qt * Rt

            X = X * scale_x
            Y = Y * scale_y
            Z = Z * scale_z
            T = T * scale_t

            coords = (X, Y, Z, T)
            return P.parent()(coords)

    def scale(self, n):
        """
        Scale all coordinates of the ThetaPoint by `n`
        """
        x, y, z, t = self.coords()
        if not isinstance(n, RingElement):
            raise ValueError(f"Cannot scale by element {n} of type {type(n)}")
        scaled_coords = (n * x, n * y, n * z, n * t)
        return self._parent(scaled_coords)

    def order(self, ell):
        """Requires ell to be small"""
        # ell should be small
        P = self
        e = 1
        for e in range(1, 10000):
            if P.is_zero():
                break
            P *= ell
            e += 1
            
        if not P.is_zero():
            raise ValueError("Theta point does not have order a power of 2")
        return 2**e
    
    def has_order(self, ell, e):
        """checks if self (assumed to be in E[ell^e]) has order ell^e"""
        return (ell**(e - 1)) * self != 0

    def double_iter(self, m):
        """
        Compute [2^n] Self

        NOTE: Assumes that no coordinate is zero at any point during the doubling
        """
        if not isinstance(m, Integer):
            try:
                m = Integer(m)
            except:
                raise TypeError(f"Cannot coerce input scalar {m = } to an integer")

        if m.is_zero():
            return self.parent().zero()

        P1 = self
        for _ in range(m):
            P1 = P1.double()
        return P1

    def __mul__(self, m):
        """
        Uses Montgomery ladder to compute [m] Self

        NOTE: Assumes that no coordinate is zero at any point during the doubling
        """
        # Make sure we're multiplying by something value
        if not isinstance(m, (int, Integer)):
            try:
                m = Integer(m)
            except:
                raise TypeError(f"Cannot coerce input scalar {m = } to an integer")

        # If m is zero, return the null point
        if not m:
            return self.parent().zero()

        # We are with ±1 identified, so we take the absolute value of m
        m = abs(m)

        P0, P1 = self, self
        P2 = P1.double()
        # If we are multiplying by two, the chain stops here
        if m == 2:
            return P2

        # Montgomery double and add.
        for bit in bin(m)[3:]:
            Q = P2.diff_addition(P1, P0)
            if bit == "1":
                P2 = P2.double()
                P1 = Q
            else:
                P1 = P1.double()
                P2 = Q

        return P1

    def __rmul__(self, m):
        return self * m

    def __imul__(self, m):
        self = self * m
        return self

    def __eq__(self, other):
        """
        Check the quality of two ThetaPoints. Note that as this is a
        projective equality, we must be careful for when certain coefficients may
        be zero.
        """
        if other.is_zero():
            return self.is_zero()
            
        if not isinstance(other, ThetaPoint):
            return False

        a1, b1, c1, d1 = self.coords()
        a2, b2, c2, d2 = other.coords()

        if d1 != 0 or d2 != 0:
            return all([a1 * d2 == a2 * d1, b1 * d2 == b2 * d1, c1 * d2 == c2 * d1])
        elif c1 != 0 or c2 != 0:
            return all([a1 * c2 == a2 * c1, b1 * c2 == b2 * c1])
        elif b1 != 0 or b2 != 0:
            return a1 * b2 == a2 * b1
        else:
            return True

    def normalize(self):
        a, b, c, d = self.coords()
        if d != 0:
            inv = 1 / d
            return (a * inv, b * inv, c * inv, 1)
        elif c != 0:
            inv = 1 / c
            return (a * inv, b * inv, 1, 0)
        elif b != 0:
            inv = 1 / b
            return (a * inv, 1, 0, 0)
        else:
            return (1, 0, 0, 0)
        
    def __hash__(self):
        return hash(self.normalize())

    def __repr__(self):
        return f"Theta point: {self.coords()}"

    # def weil_pairing(self, other, n):
    #     if not isinstance(other, JacobianPoint):
    #         raise TypeError("Both inputs must be jacobian points")
    #     if n != 2:
    #         raise NotImplementedError("Weil pairing is only implemented for n=2")
    #     if self.D2 != 0 or other.D2 != 0:
    #         raise ValueError(
    #             "Trying to call Weil pairing n = 2 on non 2-torsion points"
    #         )

    #     D1_prime = other.D1
    #     if self.D1.gcd(D1_prime) == 1:
    #         return 1
    #     else:
    #         return -1



from functools import cached_property
from itertools import product

from sage.all import (
    GF,
    Matrix,
    PolynomialRing,
    cached_function,
    cached_method,
    vector,
)

from vendors.two_isogenies.batched_inversion import batched_inversion


# Coordinate order: (00, 10, 01, 11).
INDICES = tuple(
    vector(GF(2), v, immutable=True) for v in ((0, 0), (1, 0), (0, 1), (1, 1))
)


ZERO_LABEL = (INDICES[0], INDICES[0])
TWO_TORSION_BASIS = (
    (INDICES[1], INDICES[0]),
    (INDICES[2], INDICES[0]),
    (INDICES[0], INDICES[1]),
    (INDICES[0], INDICES[2]),
)


def torsion_label(a, b):
    return (vector(GF(2), a, immutable=True), vector(GF(2), b, immutable=True))


def normalise(P):
    """Select the affine representative whose first nonzero coordinate is 1.

    This does not enforce compatibility between affine lifts.
    """
    cinv = 1 / next(c for c in P if c)
    return tuple(t * cinv for t in P)


def hadamard(P):
    x_00, x_10, x_01, x_11 = P
    x_00, x_10 = (x_00 + x_10, x_00 - x_10)
    x_01, x_11 = (x_01 + x_11, x_01 - x_11)
    return x_00 + x_01, x_10 + x_11, x_00 - x_01, x_10 - x_11


@cached_function
def quadratic_basis(Fpp):
    i = Fpp(-1).sqrt(extend=False)
    a, b = i.list()
    return i, a, 1 / b


def sqrt_fp2(x):
    Fpp = x.parent()
    p = Fpp.characteristic()
    if Fpp.degree() != 2 or p % 4 != 3:
        r = x.sqrt(extend=False)
        return min(r, -r)
    i, a, binv = quadratic_basis(Fpp)
    u, v = x.list()
    v *= binv
    u -= a * v
    e = (p + 1) // 4
    if not v:
        r = u**e
        r = Fpp(r) if r * r == u else i * (-u) ** e
    else:
        norm = u * u + v * v
        d = norm**e
        if d * d != norm:
            raise ValueError('Not a square in the base field.')
        s = (u + d) / 2
        r = s**e
        if r * r != s:
            s -= d
            r = s**e
        r = Fpp(r) + i * (v / (2 * r))
    if r * r != x:
        raise ValueError('Not a square in the base field.')
    return min(r, -r)


def normalizer_generators(Fpp):
    z = sqrt_fp2(Fpp(-1))
    return (
        Matrix(Fpp, [[1, 1, 0, 0], [1, -1, 0, 0], [0, 0, 1, 1], [0, 0, 1, -1]]),
        Matrix(Fpp, [[1, 0, 1, 0], [0, 1, 0, 1], [1, 0, -1, 0], [0, 1, 0, -1]]),
        Matrix.diagonal(Fpp, [1, z, 1, z]),
        Matrix.diagonal(Fpp, [1, 1, z, z]),
        Matrix(Fpp, [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]]),
    )


def theta_changes(Fpp):
    generators = normalizer_generators(Fpp)
    I = Matrix.identity(Fpp, 4)
    todo = [I]
    seen = {tuple(I.list())}
    for M in todo:
        yield M
        for G in generators:
            entries = normalise((G * M).list())
            if entries not in seen:
                seen.add(entries)
                todo.append(Matrix(Fpp, 4, entries))


class ThetaSurface:
    """A theta surface over the initialized quadratic field Fpp."""

    def __init__(self, O):
        self.O = normalise(O)
        self.Fpp = self.O[0].parent()
        self.H = hadamard([c * c for c in self.O])
        self.Hinv = tuple(batched_inversion(*self.H)) if all(self.H) else None
        self.Oinv = tuple(batched_inversion(*self.O)) if all(self.O) else None

    @cached_property
    def is_product(self):
        return any(
            not sum(
                (-1) ** int(a.dot_product(k))
                * self.O[INDICES.index(k + b)]
                * self.O[INDICES.index(k)]
                for k in INDICES
            )
            for a in INDICES
            for b in INDICES
            if a.dot_product(b) == 0
        )

    @cached_method
    def quartic(self):
        """The quartic equation in this theta structure (Jacobian case)."""
        if self.is_product:
            raise ValueError('A decomposable theta model has a quadric image.')
        a, b, c, d = self.O
        # Invariance under the two-torsion operators leaves sum(x_i^4) and
        # the four terms below. The theta null is singular: grad(q)(O) = 0.
        # The determinant is 16(a^2*d^2-b^2*c^2)(a^2*c^2-b^2*d^2)
        #                         (a^2*b^2-c^2*d^2), nonzero on Jacobians.
        derivatives = Matrix(
            self.Fpp,
            [
                [b * c * d, 2 * a * d * d, 2 * a * c * c, 2 * a * b * b],
                [a * c * d, 2 * b * c * c, 2 * b * d * d, 2 * b * a * a],
                [a * b * d, 2 * c * b * b, 2 * c * a * a, 2 * c * d * d],
                [a * b * c, 2 * d * a * a, 2 * d * b * b, 2 * d * c * c],
            ],
        )
        coefficients = derivatives.solve_right(vector([-4 * a**3 for a in self.O]))
        x, y, z, t = PolynomialRing(self.Fpp, 'x,y,z,t').gens()
        terms = (
            x * y * z * t,
            x * x * t * t + y * y * z * z,
            x * x * z * z + y * y * t * t,
            x * x * y * y + z * z * t * t,
        )
        return x**4 + y**4 + z**4 + t**4 + sum(c * f for c, f in zip(coefficients, terms))

    def random_point(self):
        """Sample rational theta points on the quartic."""
        q = self.quartic()
        R = PolynomialRing(self.Fpp, 't')
        while True:
            x, y = self.Fpp.random_element(), self.Fpp.random_element()
            roots = R(q(x, y, 1, R.gen())).roots(multiplicities=False)
            if roots:
                return normalise((x, y, self.Fpp(1), roots[0]))

    @cached_method
    def biquadratic_coefficients(self):
        """Reciprocals 1/U_{d,chi}(O) for the factored Riemann biquadratics.

        U_{d,chi}(O) = sum_t (-1)^(chi.t) O_t O_(t+d).
        See Lubicz--Robert, Fast change of level, Theorem 3.3.

        Entries with chi.d != 0 are None and contribute zero.
        """
        O = self.O
        return tuple(
            tuple(
                1 / x if chi.dot_product(d) == 0 else None
                for chi, x in zip(
                    INDICES,
                    hadamard(
                        [O[INDICES.index(t)] * O[INDICES.index(t + d)] for t in INDICES]
                    ),
                )
            )
            for d in INDICES
        )

    def biquadratic(self, P, Q):
        """Riemann biquadratics, evaluated before expanding their quartics."""
        B = Matrix(self.Fpp, 4)
        for d, coefficients in zip(INDICES, self.biquadratic_coefficients()):
            U = hadamard(
                [P[INDICES.index(t)] * P[INDICES.index(t + d)] for t in INDICES]
            )
            V = hadamard(
                [Q[INDICES.index(t)] * Q[INDICES.index(t + d)] for t in INDICES]
            )
            C = hadamard(
                [
                    u * v * c if c is not None else self.Fpp(0)
                    for u, v, c in zip(U, V, coefficients)
                ]
            )
            for i in range(4):
                j = INDICES.index(INDICES[i] + d)
                if i <= j:
                    B[i, j] = B[j, i] = C[i] / 2
        return B

    def sums(self, P, Q):
        """The unordered pair {theta(P+Q), theta(P-Q)} on a Jacobian."""
        if self.is_product:
            raise ValueError(
                'Product theta coordinates identify (P, Q) with '
                '(+/-P, +/-Q), either sign separately.'
            )
        P, Q = normalise(P), normalise(Q)
        if P == Q:
            return (normalise(self.double(P)), self.O)
        two_torsion_points = {T: v for v, T in self.two_torsion().items()}
        if P in two_torsion_points:
            P, Q = Q, P
        if Q in two_torsion_points:
            a, b = two_torsion_points[Q]
            T = normalise(
                [
                    (-1) ** int(a.dot_product(k)) * P[INDICES.index(k + b)]
                    for k in INDICES
                ]
            )
            return (T, T)

        B = self.biquadratic(P, Q)
        i = next((i for i in range(4) if B[i, i]), None)
        if i is None:
            # S_k D_k = 0 for every k. If B_ij != 0, one of S_i D_j
            # and D_i S_j is nonzero, so rows i and j recover D and S.
            i, j = next((i, j) for i in range(4) for j in range(i + 1, 4) if B[i, j])
            return normalise(B.row(i)), normalise(B.row(j))
        j = next((j for j in range(4) if B[i, j] ** 2 != B[i, i] * B[j, j]), None)
        if j is None:
            T = normalise(B.row(i))
            return T, T

        a, b, c = B[i, i] / 2, B[i, j], B[j, j] / 2
        d = (b + sqrt_fp2(b * b - 4 * a * c)) / (2 * a)
        s = b - a * d
        D = [(B[k, j] - B[k, i] * d) / (s - a * d) for k in range(4)]
        S = [B[k, i] - a * D[k] for k in range(4)]
        return normalise(S), normalise(D)

    def _diff_add(self, P, Q, Dinv):
        p1, p2, p3, p4 = hadamard([c * c for c in P])
        q1, q2, q3, q4 = hadamard([c * c for c in Q])
        AA_inv, BB_inv, CC_inv, DD_inv = self.Hinv

        xp = AA_inv * p1 * q1
        yp = BB_inv * p2 * q2
        zp = CC_inv * p3 * q3
        tp = DD_inv * p4 * q4

        X, Y, Z, T = hadamard((xp, yp, zp, tp))
        return tuple(c * d / 4 for c, d in zip((X, Y, Z, T), Dinv))

    def diff_add(self, P, Q, D):
        """An affine lift of theta(P+Q), given theta(P-Q)=D.

        The scaling is lambda_P^2 lambda_Q^2 / lambda_D. In particular,
        do not normalize intermediate values for compatible affine lifts
        (Lubicz--Robert, Fast change of level, Section 3 and Theorem 3.8).
        """
        if all(D) and self.Hinv is not None:
            return self._diff_add(P, Q, tuple(1 / d for d in D))
        # Choose i with D_i != 0, even when another coordinate of theta(P-Q) vanishes.
        B = self.biquadratic(P, Q)
        i = next(i for i in range(4) if D[i])
        s = B[i, i] / (2 * D[i])
        return tuple((B[i, j] - s * D[j]) / D[i] if i != j else s for j in range(4))

    def double(self, P):
        if self.Oinv is not None and self.Hinv is not None:
            return self._diff_add(P, P, self.Oinv)
        return self.diff_add(P, P, self.O)

    def mul(self, n, P):
        """Montgomery ladder; negative scalars have the same image."""
        U, V = self.O, normalise(P)
        Dinv = tuple(1 / d for d in P) if all(P) and self.Hinv is not None else None
        n = abs(int(n))
        if not n:
            return self.O
        e = (n & -n).bit_length() - 1
        for bit in bin(n >> e)[2:]:
            S = (
                self._diff_add(U, V, Dinv)
                if Dinv is not None
                else self.diff_add(U, V, P)
            )
            if bit == '0':
                U, V = self.double(U), S
            else:
                U, V = S, self.double(V)
        # [2^e m]P: use the ladder only for odd m, then e doublings.
        for _ in range(e):
            U = self.double(U)
        return normalise(U)

    @cached_method
    def two_torsion(self):
        """Two-torsion coordinates (a, b), pairing (-1)^(a.d+b.c).

        The labels are pairs of vectors in F_2^2. Their ordered basis is
        ((1,0),0), ((0,1),0), (0,(1,0)), (0,(0,1)),
        with J = [[0, I_2], [I_2, 0]], relative to this level-two structure.
        """
        return {
            (a, b): normalise(
                [
                    (-1) ** int(a.dot_product(k)) * self.O[INDICES.index(k + b)]
                    for k in INDICES
                ]
            )
            for a in INDICES
            for b in INDICES
        }

    def four_torsion_basis(self):
        """Halving two-torsion points: return Q_i with 2Q_i = P_i for a basis of A[2]."""
        if self.is_product:
            raise ValueError('Use elliptic sampling on a product.')
        basis = []
        for v in TWO_TORSION_BASIS:
            T = self.two_torsion()[v]
            # H(P^2)^2 is projectively H(T*O)*H(O^2). Normalize before
            # each square root: the scalar multiplying the coordinate tuple
            # need not be a square.
            squares = normalise(
                tuple(
                    x * h
                    for x, h in zip(
                        hadamard([t * o for t, o in zip(T, self.O)]), self.H
                    )
                )
            )
            U = tuple(sqrt_fp2(x) for x in squares)
            found = None

            for signs in product((1, -1), repeat=3):
                squares = normalise(hadamard([u * s for u, s in zip(U, (1,) + signs)]))
                try:
                    roots = tuple(sqrt_fp2(x) for x in squares)
                except ValueError:
                    continue
                for signs in product((1, -1), repeat=3):
                    P = tuple(r * s for r, s in zip(roots, (1,) + signs))
                    if not self.quartic()(*P) and normalise(self.double(P)) == T:
                        found = normalise(P)
                        break
                if found is not None:
                    break

            assert found is not None
            basis.append(found)
        return tuple(basis)

"""Split theta structures and elliptic curve point pairs."""

from itertools import product

from sage.all import (
    EllipticCurve,
    Matrix,
    PolynomialRing,
    ZZ,
    cached_function,
    cached_method,
    vector,
)

from theta_rm.theta import INDICES, ThetaSurface, normalise, sqrt_fp2
from vendors.two_isogenies.splitting import splitting_matrix


class ProductPoint(tuple):
    def __new__(cls, P, Q):
        return tuple.__new__(cls, (P, Q))

    def points(self):
        return tuple(self)

    def double(self):
        return self + self

    def __add__(self, Q):
        return ProductPoint(self[0] + Q[0], self[1] + Q[1])

    def __sub__(self, Q):
        return ProductPoint(self[0] - Q[0], self[1] - Q[1])

    def __neg__(self):
        return ProductPoint(-self[0], -self[1])

    def __mul__(self, n):
        return ProductPoint(ZZ(n) * self[0], ZZ(n) * self[1])

    __rmul__ = __mul__

    def weil_pairing(self, Q, n):
        return self[0].weil_pairing(Q[0], n) * self[1].weil_pairing(Q[1], n)


class EllipticTheta:
    def __init__(self, E):
        self.E = E
        Fpp = E.base_ring()
        x = PolynomialRing(Fpp, 'x').gen()
        self.h = (
            x**3
            + (E.a2() + E.a1() ** 2 / 4) * x * x
            + (E.a4() + E.a1() * E.a3() / 2) * x
            + E.a6()
            + E.a3() ** 2 / 4
        )
        e0, e1, e2 = sorted(self.h.roots(multiplicities=False))
        r = sqrt_fp2((e2 - e0) / (e1 - e0))
        self.t = sqrt_fp2((1 - r) / (1 + r))
        self.e0 = e0
        self.c = (e1 - e0) * r
        self.O = (Fpp(1), self.t)

    @classmethod
    def from_null(cls, O):
        t = O[1] / O[0]
        Fpp = t.parent()
        e1, e2 = (1 + t * t) ** 2, (1 - t * t) ** 2
        E = EllipticCurve(Fpp, [0, -e1 - e2, 0, e1 * e2, 0])
        chart = cls.__new__(cls)
        chart.E = E
        chart.t = t
        chart.e0 = Fpp(0)
        chart.c = 1 - t**4
        x = PolynomialRing(Fpp, 'x').gen()
        chart.h = x * (x - e1) * (x - e2)
        chart.O = (Fpp(1), t)
        return chart

    def __call__(self, P):
        if P == 0:
            return self.O
        z = P[0] - self.e0
        return normalise([z - self.c, self.t * (z + self.c)])

    def points_from_theta(self, T):
        """Return the elliptic points with theta coordinates T."""
        u, v = T
        if v == self.t * u:
            return [self.E(0)]
        x = self.e0 + self.c * (v + self.t * u) / (v - self.t * u)
        try:
            y = sqrt_fp2(self.h(x)) - (self.E.a1() * x + self.E.a3()) / 2
        except ValueError:
            return []
        P = self.E(x, y)
        return [P, -P] if P != -P else [P]


class ProductSurface(ThetaSurface):
    def __init__(self, O=None, curves=None):
        if curves is not None:
            self.charts = [EllipticTheta(E) for E in curves]
            U, V = [C.O for C in self.charts]
            O = [U[i] * V[j] for j in range(2) for i in range(2)]
            self.M = Matrix.identity(O[0].parent(), 4)
        else:
            Fpp = O[0].parent()
            zero = next(
                (a, b)
                for a in INDICES
                for b in INDICES
                if a.dot_product(b) == 0
                and not sum(
                    (-1) ** int(a.dot_product(k))
                    * O[INDICES.index(k + b)]
                    * O[INDICES.index(k)]
                    for k in INDICES
                )
            )
            self.M = splitting_matrix(O, zero, sqrt_fp2(Fpp(-1)))
            T = self.M * vector(O)
            self.charts = [
                EllipticTheta.from_null(T[:2]),
                EllipticTheta.from_null((T[0], T[2])),
            ]

        self.Minv = self.M.inverse()
        self.E1, self.E2 = (C.E for C in self.charts)
        self.zero = ProductPoint(self.E1(0), self.E2(0))
        super().__init__(O)
        assert self.theta_representative(self.zero) == self.O and self.is_product

    def theta_representative(self, P):
        """Theta coordinates of P in this surface's theta structure.

        The split coordinates (u_i v_j) are ordered with i varying first.
        """
        U, V = [C(Q) for C, Q in zip(self.charts, P)]
        T = vector([U[i] * V[j] for j in range(2) for i in range(2)])
        return normalise(self.Minv * T)

    def points_from_theta(self, P):
        """Return the product points with theta coordinates P.

        These are (U, V) -> (+/-U, +/-V), with either sign chosen separately.
        """
        T = self.M * vector(P)
        # Choose a nonzero row and column of (u_i v_j); these also cover infinity.
        i = next(i for i in range(2) if T[i] or T[i + 2])
        j = next(j for j in range(2) if T[2 * j] or T[2 * j + 1])
        return [
            ProductPoint(P, Q)
            for P in self.charts[0].points_from_theta(T[2 * j : 2 * j + 2])
            for Q in self.charts[1].points_from_theta((T[i], T[i + 2]))
        ]

    def sign_compatible_points(self, points, sums):
        """Return tuples B with theta(B[j]) = points[j] and
        theta(B[i]+B[j]) = sums[i, j] for every i < j.

        Product theta coordinates identify (P, Q) with (+/-P, +/-Q),
        with either sign chosen separately; impose the sums on elliptic curve point pairs."""
        todo = [()]
        for j, P in enumerate(points):
            todo = [
                B + (Q,)
                for B in todo
                for Q in self.points_from_theta(P)
                if all(
                    self.theta_representative(B[i] + Q) == sums[i, j] for i in range(j)
                )
            ]
        return todo

    def random_point(self):
        """Return an elliptic curve point pair on E1 x E2."""
        return ProductPoint(self.E1.random_point(), self.E2.random_point())

    @cached_method
    def torsion_basis(self, n):
        """Return an ordered basis of elliptic curve point pairs in (E1 x E2)[n]."""
        P, Q = self.E1.torsion_basis(n)
        R, S = (P, Q) if self.E1 == self.E2 else self.E2.torsion_basis(n)
        return (
            ProductPoint(P, self.E2(0)),
            ProductPoint(self.E1(0), R),
            ProductPoint(Q, self.E2(0)),
            ProductPoint(self.E1(0), S),
        )

    def four_torsion_basis(self):
        return tuple(self.theta_representative(P) for P in self.torsion_basis(4))


class ProductTwoIsogeny:
    def __init__(self, A, kernel):
        self.domain = A
        P, Q = [A.points_from_theta(T)[0] for T in kernel]
        H = (P, Q, P + Q)
        U = next((T[0] for T in H if T[1] == 0), None)
        # A product kernel has a nonzero point on each elliptic factor.
        if U is not None:
            V = next(T[1] for T in H if T[0] == 0)
            f, g = A.E1.isogeny(U), A.E2.isogeny(V)
            B = ProductSurface(curves=(f.codomain(), g.codomain()))
            self.evaluate = lambda T: B.theta_representative(
                ProductPoint(f(T[0]), g(T[1]))
            )
        else:
            # The kernel is the graph of E1[2] -> E2[2]. Check whether it
            # extends to an isomorphism E1 -> E2 before using gluing.
            f = next(
                (f for f in A.E1.isomorphisms(A.E2) if all(f(T[0]) == T[1] for T in H)),
                None,
            )
            if f is not None:
                B = ProductSurface(curves=(A.E2, A.E2))
                self.evaluate = lambda T: B.theta_representative(
                    ProductPoint(f(T[0]) + T[1], f(T[0]) - T[1])
                )
            else:
                from vendors.two_isogenies.gluing import GluingThetaIsogeny

                basis = A.torsion_basis(8)
                fourth_roots = []
                for T in (P, Q):
                    factors = [
                        [
                            a * basis[j][j] + b * basis[j + 2][j]
                            for a, b in product(range(8), repeat=2)
                            if 4 * (a * basis[j][j] + b * basis[j + 2][j]) == T[j]
                        ]
                        for j in range(2)
                    ]
                    fourth_roots.append(
                        # Fix the lifts independently of the random torsion basis.
                        sorted(
                            (ProductPoint(U, V) for U, V in product(*factors)),
                            key=lambda T: tuple(tuple(P) for P in T),
                        )
                    )
                U = fourth_roots[0][0]
                V = next(V for V in fourth_roots[1] if U.weil_pairing(V, 8) == 1)
                glue = GluingThetaIsogeny(U, V)
                B = glue._codomain
                self.evaluate = glue

        self.codomain = B
        assert all(self.evaluate(T) == B.O for T in (A.zero, *H))

    def __call__(self, P):
        return self.evaluate(P)


@cached_function
def product_isomorphism_maps(E1, E2, D1, D2):
    result = []
    for swap in (False, True):
        E, D = (E2, E1) if swap else (E1, E2)
        for f, g in product(E.isomorphisms(D1), D.isomorphisms(D2)):
            result.append(
                lambda P, f=f, g=g, i=int(swap): ProductPoint(f(P[i]), g(P[1 - i]))
            )
    return tuple(result)


def product_isomorphisms(A, B):
    return iter(product_isomorphism_maps(A.E1, A.E2, B.E1, B.E2))

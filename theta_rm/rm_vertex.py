"""RM-data on abelian surfaces in level-two theta coordinates."""

from itertools import combinations, product

from sage.all import GF, Matrix, cached_method, prod, vector
from sage.rings.finite_rings.element_base import FinitePolyExtElement

from theta_rm.theta import INDICES, torsion_label, ZERO_LABEL, hadamard, normalise
from theta_rm.five_isogeny import ThetaFiveIsogeny
from theta_rm.isomorphisms import normalizer_matrices, theta_isomorphisms
from theta_rm.two_isogeny import ThetaTwoIsogeny
from theta_rm.product import (
    ProductPoint,
    ProductSurface,
    ProductTwoIsogeny,
    product_isomorphisms,
)


class RMVertex:
    def __init__(self, A, five_kernel, three_torsion):
        self.surface = (
            ProductSurface(A.O)
            if A.is_product and not isinstance(A, ProductSurface)
            else A
        )
        self.five_kernel = tuple(map(normalise, five_kernel))
        # theta(P_1), ..., theta(P_4), followed by theta(P_i+P_j)
        # for (i,j) = (1,2), (1,3), (1,4), (2,3), (2,4), (3,4).
        # The six sums retain relative signs, up to simultaneous negation.
        self.three_torsion = tuple(map(normalise, three_torsion))
        if len(self.three_torsion) != 10:
            raise ValueError('Store four basis points and their six pairwise sums.')
        assert len(self.five_kernel) == 3

    @cached_method
    def partial_ppas_invariant(self):
        """Return Igusa i_3, or the unordered elliptic j-pair for products, ignoring RM."""
        # i_3 is only a cheap equality prefilter: matching values do not imply equality.
        # __eq__ still checks for an isomorphism respecting the RM data.
        A = self.surface
        if A.is_product:
            return ('P', *sorted((A.E1.j_invariant(), A.E2.j_invariant())))
        a, b, c, d = A.O
        # Squares of the ten even theta constants.
        t = list(hadamard((a * a, b * b, c * c, d * d)))
        t += [2 * (a * b + c * d), 2 * (a * b - c * d)]
        t += [2 * (a * c + b * d), 2 * (a * c - b * d)]
        t += [2 * (a * d + b * c), 2 * (a * d - b * c)]

        h4 = sum(x**4 for x in t)
        h10 = prod(t)
        return ('J', h4**5 / h10**2)  # Streng's absolute Igusa invariant i_3.

    def get_type(self):
        """For debugging. Returns a string describing the type of underlying PPAS."""
        A = self.surface
        if not A.is_product:
            return 'J'
        j1, j2 = self.partial_ppas_invariant()[1:]
        if j1 == j2:
            return 'S_0' if j1 == 0 else 'S_1728' if j1 == A.Fpp(1728) else 'S'
        special = [s for s in (0, 1728) if A.Fpp(s) in (j1, j2)]
        return 'P' + ''.join('_' + str(s) for s in special)

    def __str__(self):
        return f'{self.get_type()}: {self.partial_ppas_invariant()[1:]}'

    @cached_method
    def three_torsion_basis(self):
        """Return the stored basis: elliptic curve point pairs on products, theta points otherwise."""
        A = self.surface
        if not A.is_product:
            return self.three_torsion[:4]
        points = self.three_torsion
        pairwise_sums = dict(zip(combinations(range(4), 2), points[4:]))
        bases = A.sign_compatible_points(points[:4], pairwise_sums)
        if not bases:
            raise ValueError(
                'The supplied theta points and pairwise sums are inconsistent on the product.'
            )
        basis = bases[0]
        assert all(3 * P == A.zero for P in basis)
        assert (
            len(
                {
                    sum((a * P for a, P in zip(v, basis)), A.zero)
                    for v in product(range(3), repeat=4)
                }
            )
            == 81
        )

        return basis

    @cached_method
    def three_torsion_point(self, v):
        """Return theta(sum(v_i P_i)) for v in F_3^4 from the ten stored points."""
        A = self.surface
        v = tuple(int(c) % 3 for c in v)
        support = [i for i, c in enumerate(v) if c]
        # theta(v) = theta(-v); reconstruct only one representative.
        if support and v[support[0]] == 2:
            return self.three_torsion_point(tuple(-c % 3 for c in v))
        if not support:
            return A.O
        if len(support) == 1:
            return self.three_torsion[support[0]]
        if len(support) == 2:
            i, j = support
            pairs = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
            S = self.three_torsion[4 + pairs.index((i, j))]
            return (
                S
                if v[i] % 3 == v[j] % 3
                else normalise(A.diff_add(self.three_torsion[i], self.three_torsion[j], S))
            )
        # Intersect the possible sums from different decompositions. Independence
        # of the basis leaves exactly one theta point in the intersection.
        candidates = None
        for i in support:
            w = tuple(0 if j == i else c for j, c in enumerate(v))
            choices = set(A.sums(self.three_torsion_point(w), self.three_torsion[i]))
            candidates = choices if candidates is None else candidates & choices
            if len(candidates) == 1:
                break
        assert len(candidates) == 1
        return candidates.pop()

    @cached_method
    def three_torsion_coordinates(self):
        # Choose a vector representative of each pair {v, -v}.
        return {
            self.three_torsion_point(v): vector(GF(3), v)
            for v in product(range(3), repeat=4)
        }

    @cached_method
    def three_torsion_action(self):
        """The action of tau on A[3]: elliptic curve point pairs, theta points otherwise."""
        A = self.surface
        if not A.is_product:
            assert all(P != A.O and A.mul(3, P) == A.O for P in self.three_torsion[:4])
            assert all(
                S in A.sums(self.three_torsion[i], self.three_torsion[j])
                for (i, j), S in zip(combinations(range(4), 2), self.three_torsion[4:])
            )
            result = {}
            for v in product(range(3), repeat=4):
                a, b, c, d = v
                w = tuple(t % 3 for t in (-a + 2 * b, 2 * a + b, -c + 2 * d, 2 * c + d))
                result[self.three_torsion_point(v)] = self.three_torsion_point(w)
            assert len(result) == 41
            return result
        basis = self.three_torsion_basis()
        P, Q, R, S = basis
        images = (-P + 2 * Q, 2 * P + Q, -R + 2 * S, 2 * R + S)
        zero = self.surface.zero
        return {
            sum((a * P for a, P in zip(v, basis)), zero): sum(
                (a * P for a, P in zip(v, images)), zero
            )
            for v in product(range(3), repeat=4)
        }

    @cached_method
    def _quotient(self):
        return ThetaFiveIsogeny(self.surface, self.five_kernel)

    @cached_method
    def tau(self):
        """The endomorphism tau = 2*sigma - 1 on this surface."""
        A = self.surface
        if A.is_product:
            basis = self.three_torsion_basis()
            P, Q, R, S = basis
            expected = (-P + 2 * Q, 2 * P + Q, -R + 2 * S, 2 * R + S)
            # Rosati symmetry and sigma^2-sigma=1 force sigma = [[0, f^-1], [f, 1]]
            # or the conjugate embedding, with f an elliptic isomorphism.
            # In particular (1+tau)/2 is an endomorphism.
            for f in A.E1.isomorphisms(A.E2):
                for sign in (1, -1):

                    def tau(T, f=f, sign=sign):
                        P, Q = T
                        return sign * ProductPoint(-P + 2 * (~f)(Q), 2 * f(P) + Q)

                    if any(tau(P) != Q for P, Q in zip(basis, expected)):
                        continue
                    return tau
            raise ValueError('No product endomorphism matches the three-torsion action.')
        phi = self._quotient()
        # The kernel sum transports the level-two theta structure
        # (Lubicz--Robert, Fast change of level, Proposition 2.15 and Corollary 4.5).
        # tau = 2*sigma-1 acts as the identity on A[2]. Finding the polarized
        # isomorphism therefore uses the 16 translations determined by this structure:
        # H_{a,b}(x)_k = (-1)^(a.k) x_{k+b}.
        O = phi.codomain.O
        a, b = next(
            (a, b)
            for a, b in product(INDICES, repeat=2)
            if normalise(
                [
                    (-1) ** int(a.dot_product(k)) * O[INDICES.index(k + b)]
                    for k in INDICES
                ]
            )
            == phi.source.O
        )
        N = phi.M.inverse()

        def tau(P):
            Q = phi(P)
            return normalise(
                N
                * vector(
                    [
                        (-1) ** int(a.dot_product(k)) * Q[INDICES.index(k + b)]
                        for k in INDICES
                    ]
                )
            )

        return tau

    @cached_method
    def five_kernel_basis(self):
        """Return kernel generators, choosing product lifts annihilated by tau."""
        A = self.surface
        if not A.is_product:
            return self.five_kernel[:2]
        tau = self.tau()
        for P, Q in A.sign_compatible_points(
            self.five_kernel[:2], {(0, 1): self.five_kernel[2]}
        ):
            if tau(P) == A.zero and tau(Q) == A.zero:
                assert len({a * P + b * Q for a, b in product(range(5), repeat=2)}) == 25
                return P, Q
        raise ValueError('No product lifts of the five-kernel lie in ker(tau).')

    @cached_method
    def five_kernel_set(self):
        """Return the kernel as actual product points or Jacobian theta points."""
        if self.surface.is_product:
            P, Q = self.five_kernel_basis()
            return frozenset(a * P + b * Q for a, b in product(range(5), repeat=2))
        return frozenset(self._quotient().kernel)

    five_kernel_points = five_kernel_set  # Keep the existing public name.

    @cached_method
    def four_torsion_basis(self):
        A = self.surface
        return A.torsion_basis(4) if A.is_product else A.four_torsion_basis()

    def kernel_set(self, P, Q):
        A = self.surface
        two_torsion_points = A.two_torsion()
        labels = {T: v for v, T in two_torsion_points.items()}
        a, b = labels[P]
        c, d = labels[Q]
        return frozenset((A.O, P, Q, two_torsion_points[torsion_label(a + c, b + d)]))

    @cached_method
    def two_kernels(self):
        A = self.surface
        tau = self.tau()
        B = self.four_torsion_basis()
        # For 2Q = P, sigma(P) = Q+tau(Q). Products use componentwise elliptic curve addition;
        # on Jacobians, either sign of tau(Q) gives the same kernel-membership test.
        if A.is_product:
            S = [Q + tau(Q) for Q in B]
            assert all(2 * P == A.zero for P in S)
            images = {}
            for v in product(range(2), repeat=4):
                if any(v):
                    P = sum((2 * a * Q for a, Q in zip(v, B)), A.zero)
                    sigma_P = sum((a * Q for a, Q in zip(v, S)), A.zero)
                    images[A.theta_representative(P)] = (
                        A.theta_representative(sigma_P),
                    ) * 2
        else:
            # Either sum of halves doubles to the sum of their two-torsion points.
            halves = [A.O]
            for Q in B:
                halves += [Q if R == A.O else A.sums(R, Q)[0] for R in halves]
            # For 2Q=P, Q+tau(Q) and Q-tau(Q) differ by tau(P)=P.
            # Thus both belong to K, or neither does, whenever P belongs to K.
            images = {
                normalise(A.double(Q)): A.sums(Q, tau(Q)) for Q in halves[1:]
            }
            assert len(images) == 15

        two_torsion_points = A.two_torsion()
        assert all(
            T in two_torsion_points.values() for pair in images.values() for T in pair
        )
        kernels = []
        seen = set()
        for (a, b), (c, d) in combinations(list(two_torsion_points)[1:], 2):
            if a.dot_product(d) + b.dot_product(c):
                continue
            P, Q = two_torsion_points[a, b], two_torsion_points[c, d]
            K = frozenset((A.O, P, Q, two_torsion_points[torsion_label(a + c, b + d)]))
            if K not in seen and all(T in K for R in (P, Q) for T in images[R]):
                kernels.append((P, Q))
                seen.add(K)
        assert len(kernels) == 5
        return tuple(kernels)

    def dual_kernel_set(self, phi):
        A = self.surface
        images = [
            phi(2 * Q) if A.is_product else phi(normalise(A.double(Q)))
            for Q in self.four_torsion_basis()
        ]
        two_torsion_points = phi.codomain.two_torsion()
        labels = {P: v for v, P in two_torsion_points.items()}
        span = {ZERO_LABEL}
        for P in images:
            a, b = labels[P]
            span |= {torsion_label(c + a, d + b) for c, d in span}
        assert len(span) == 4
        return frozenset(two_torsion_points[v] for v in span)

    dual_kernel = dual_kernel_set  # Keep the existing public name.

    def push(self, phi):
        if self.surface.is_product:
            basis = self.three_torsion_basis()
            P, Q = self.five_kernel_basis()
            kernel = [P, Q, P + Q]
            three_torsion = list(basis) + [
                P + Q for P, Q in combinations(basis, 2)
            ]
        else:
            kernel = self.five_kernel
            three_torsion = self.three_torsion
        return RMVertex(
            phi.codomain,
            [phi(P) for P in kernel],
            [phi(P) for P in three_torsion],
        )

    @cached_method
    def get_neighbor_edges(self, dual_kernel=None):
        """Return outgoing edges, placing the supplied incoming dual kernel first."""
        if dual_kernel is not None:
            edges = self.get_neighbor_edges()
            i = next(
                i for i, (_, _, K) in enumerate(edges)
                if self.kernel_set(*K) == dual_kernel
            )
            return [edges[i], *edges[:i], *edges[i + 1:]]
        A = self.surface
        result = []
        for K in self.two_kernels():
            phi = ProductTwoIsogeny(A, K) if A.is_product else ThetaTwoIsogeny(A, K)
            result.append((self.push(phi), phi, K))
        return result

    @cached_method
    def get_neighbors(self):
        """Return a dictionary mapping neighbors to edge multiplicities."""
        result = {}
        for V, _, _ in self.get_neighbor_edges():
            result[V] = result.get(V, 0) + 1
        return result

    def is_K_linear(self, other, f):
        """Check f tau_A = tau_B f on A[3], allowing simultaneously changing all signs."""
        # In particular, a dual composition is [2], hence [-1] on A[3].
        # The theta coordinates of the four points and all six sums then agree.
        if all(f(P) == Q for P, Q in zip(self.three_torsion, other.three_torsion)):
            return True
        coordinates = other.three_torsion_coordinates()
        columns = [coordinates.get(f(P)) for P in self.three_torsion[:4]]
        if any(v is None for v in columns):
            return False
        images = {
            pair: f(P)
            for pair, P in zip(combinations(range(4), 2), self.three_torsion[4:])
        }
        T = Matrix(GF(3), [[-1, 2, 0, 0], [2, 1, 0, 0], [0, 0, -1, 2], [0, 0, 2, 1]])
        # Fix the first column: simultaneously changing all signs replaces f by -f.
        # The condition f tau_A = tau_B f is unchanged.
        for signs in product((1, -1), repeat=3):
            C = [columns[0]] + [s * v for s, v in zip(signs, columns[1:])]
            if not all(
                other.three_torsion_point(tuple(C[i] + C[j])) == Q for (i, j), Q in images.items()
            ):
                continue
            M = Matrix(C).transpose()
            return M.is_invertible() and M * T == T * M
        return False

    def __eq__(self, other):
        if self is other:
            return True
        if (
            not isinstance(other, RMVertex)
            or self.partial_ppas_invariant() != other.partial_ppas_invariant()
        ):
            return False
        A, B = self.surface, other.surface
        if A.is_product:
            basis = self.three_torsion_basis()
            tau = self.tau()
            other_tau = other.tau()
            for f in product_isomorphisms(A, B):
                if all(f(tau(P)) == other_tau(f(P)) for P in basis) and all(
                    f(P) in other.five_kernel_set() for P in self.five_kernel_basis()
                ):
                    return True
        else:
            for f in theta_isomorphisms(A, B):
                if self.is_K_linear(other, f) and all(
                    f(P) in other.five_kernel_set() for P in self.five_kernel_basis()
                ):
                    return True
        return False

    def __hash__(self):
        """Equal RM vertices have equal underlying PPAS invariants."""
        return hash(self.partial_ppas_invariant())

    @cached_method
    def rm_invariants(self):
        r"""Return a canonical representative of the K-equivalence class.

        Choose O_can as the lexicographically minimal normalized theta-null
        representative among all changes of theta structure.
        This fixes the theta model B; for products, reconstruct the ordered
        elliptic models B.E1 x B.E2.
        Transport the kernel and three-torsion action together through each
        polarized isomorphism f: A -> B. Choose the lexicographically minimal
        pair (G, R) representing f(ker(tau)) and
        (f tau f^{-1})|_{B[3]}, respectively.

        Return (type, O_can, G, R), where:

        - type is 'P' for products or 'J' for Jacobians.
        - On products, G is the sorted tuple of the 25 points of f(ker(tau)).
          R is the sorted tuple (f(P), f(tau(P))) for all 81 points P in A[3].
          Each product point is an ordered pair of elliptic projective
          coordinate tuples on B.E1 x B.E2.
        - On Jacobians, G is the sorted tuple of the 13 distinct theta points
          of f(ker(tau)). R is the sorted tuple of the 41 triples
          (theta(f(P)), theta(f(tau(P))), theta(f(P + tau(P))))
          for P in A[3]/{+/-1}, including zero.

        All theta coordinates are normalized. Sorting and minimization use
        Sage's ordering over the fixed common field.
        """
        A = self.surface
        Fpp = A.Fpp

        def encode(x):
            if isinstance(x, FinitePolyExtElement):
                return x
            if isinstance(x, dict):
                return tuple(sorted(encode(pair) for pair in x.items()))
            if isinstance(x, (set, frozenset)):
                return tuple(sorted(encode(P) for P in x))
            return tuple(encode(P) for P in x)

        # Keep every minimizer, including automorphisms at special theta nulls.
        O_can = None
        maps = []
        for M in normalizer_matrices(Fpp):
            O = normalise(M * vector(A.O))
            if O_can is None or O < O_can:
                O_can, maps = O, [M]
            elif O == O_can:
                maps.append(M)

        if A.is_product:
            B = ProductSurface(O_can)
            P, Q = self.five_kernel_basis()
            kernel = [a * P + b * Q for a, b in product(range(5), repeat=2)]
            action = self.three_torsion_action()

            candidates = []
            for f in product_isomorphisms(A, B):
                G = encode({f(P) for P in kernel})
                R = encode({f(P): f(Q) for P, Q in action.items()})
                candidates.append((G, R))
        else:
            kernel = self._quotient().kernel
            T = Matrix(
                GF(3), [[-1, 2, 0, 0], [2, 1, 0, 0], [0, 0, -1, 2], [0, 0, 2, 1]]
            )
            triples = [
                (
                    self.three_torsion_point(tuple(v)),
                    self.three_torsion_point(tuple(T * v)),
                    self.three_torsion_point(tuple(v + T * v)),
                )
                for v in self.three_torsion_coordinates().values()
            ]
            candidates = []
            for M in maps:

                def f(P):
                    return normalise(M * vector(P))

                G = encode({f(P) for P in kernel})
                R = encode({tuple(f(P) for P in triple) for triple in triples})
                candidates.append((G, R))

        return ('P' if A.is_product else 'J', O_can, *min(candidates))

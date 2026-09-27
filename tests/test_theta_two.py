"""All maximal isotropic subgroups: kernels, phi(P+R)=phi(P), and polarized duals."""

from itertools import combinations, product

from sage.all import EllipticCurve, GF, Matrix, set_random_seed

from theta_rm.theta import torsion_label, TWO_TORSION_BASIS, normalise
from theta_rm.initialization import gen_rm_hash_prime, get_initial_vertex
from theta_rm.product import ProductSurface, ProductTwoIsogeny
from theta_rm.two_isogeny import ThetaTwoIsogeny
from theta_rm.isomorphisms import theta_isomorphisms


# Check nine product kernels and six graphs of isomorphisms E1[2] -> E2[2].
set_random_seed(239)
Fpp = GF(239**2, 'i', modulus=[1, 0, 1])
E, J = EllipticCurve(Fpp, [1, 0]), EllipticCurve(Fpp, [0, 1])
for E1, E2, products in ((E, E, 11), (J, J, 12), (E, J, 9)):
    A = ProductSurface(curves=(E1, E2))
    B = A.torsion_basis(2)
    points = [
        sum((a * P for a, P in zip(v, B)), A.zero) for v in product(range(2), repeat=4)
    ]
    samples = [A.zero, *A.torsion_basis(8), A.random_point()]
    seen = set()
    split = 0

    for P, Q in combinations(points[1:], 2):
        H = frozenset((A.zero, P, Q, P + Q))
        if H in seen or P.weil_pairing(Q, 2) != 1:
            continue
        seen.add(H)
        phi = ProductTwoIsogeny(
            A, [A.theta_representative(P), A.theta_representative(Q)]
        )
        C = phi.codomain
        split += C.is_product
        assert {T for T in points if phi(T) == C.O} == H
        assert all(phi(T + R) == phi(T) for T in samples for R in H)
    assert len(seen) == 15 and split == products
    print(
        f'j={E1.j_invariant(), E2.j_invariant()}: {split} products, {15 - split} Jacobians.'
    )

# Check all 15 Jacobian quotients at small and cryptographic exponents.
for e in (3, 256):
    set_random_seed(31000 + e)
    p, _ = gen_rm_hash_prime(e)
    V = get_initial_vertex(p, e)
    W, rho, _ = next(
        edge for edge in V.get_neighbor_edges() if not edge[0].surface.is_product
    )
    A, B = V.surface, W.surface
    two_torsion_points = B.two_torsion()
    # Halving two-torsion points gives 2Q_i = P_i for a basis of B[2].
    halves = B.four_torsion_basis()
    assert tuple(normalise(B.double(P)) for P in halves) == tuple(
        two_torsion_points[v] for v in TWO_TORSION_BASIS
    )
    assert all(B.quartic()(*P) == 0 and B.mul(4, P) == B.O for P in halves)
    seen = set()
    samples = [A.zero, *A.torsion_basis(8), *A.torsion_basis(3), *A.torsion_basis(5)]
    samples += [A.random_point() for _ in range(3)]

    # Check differential addition and {theta(P+Q), theta(P-Q)}.
    for P, Q in zip(samples[1:], samples[2:]):
        assert set(B.sums(rho(P), rho(Q))) == {rho(P + Q), rho(P - Q)}
        assert normalise(B.diff_add(rho(P), rho(Q), rho(P - Q))) == rho(P + Q)

    # Compare theta arithmetic with elliptic arithmetic before gluing.
    for P in samples:
        for n in (0, -1, 2, 3, 12, (p + 1) // 4, p + 1):
            assert B.mul(n, rho(P)) == rho(n * P)
            assert B.mul(n, tuple(3 * x for x in rho(P))) == rho(n * P)
    for P, Q in zip(samples[1:], samples[2:]):
        X, Y, S, D = (rho(T) for T in (P, Q, P + Q, P - Q))
        C = B.biquadratic(X, Y)
        expected = Matrix(B.Fpp, 4, lambda i, j: S[i] * D[j] + D[i] * S[j])
        assert normalise(C.list()) == normalise(expected.list())
        assert C.rank() <= 2

    # Each quotient kills precisely the supplied subgroup of B[2].
    for (a, b), (c, d) in combinations(list(two_torsion_points)[1:], 2):
        if a.dot_product(d) + b.dot_product(c):
            continue
        H = frozenset(
            (
                B.O,
                two_torsion_points[a, b],
                two_torsion_points[c, d],
                two_torsion_points[torsion_label(a + c, b + d)],
            )
        )
        if H in seen:
            continue
        seen.add(H)
        phi = ThetaTwoIsogeny(B, [two_torsion_points[a, b], two_torsion_points[c, d]])
        C = phi.codomain
        assert {T for T in two_torsion_points.values() if phi(T) == C.O} == H
        # The dual kernel is phi(B[2]), obtained without RMVertex.dual_kernel.
        dual = set(map(phi, two_torsion_points.values()))
        assert len(dual) == 4
        K = list(dual - {C.O})
        if C.is_product:
            D = ProductSurface(C.O)
            psi = ProductTwoIsogeny(D, K[:2])
            points = [phi(rho(P)) for P in samples]
            sums = {
                (i, j): phi(rho(P + Q))
                for i, P in enumerate(samples)
                for j, Q in enumerate(samples)
                if i < j
            }
            # (P, Q) -> (+/-P, +/-Q), either sign separately: the supplied
            # pairwise sums select product points before composing with psi.
            bases = D.sign_compatible_points(points, sums)
            assert bases
            images = [psi(P) for P in bases[0]]
        else:
            psi = ThetaTwoIsogeny(C, K[:2])
            images = [psi(phi(rho(P))) for P in samples]
        expected = [B.mul(2, rho(P)) for P in samples]
        assert any(
            all(f(X) == Y for X, Y in zip(images, expected))
            for f in theta_isomorphisms(psi.codomain, B)
        )
    assert len(seen) == 15

    print(f'e={e}: passed: theta arithmetic and all 15 (2,2)-isogenies; composing with the dual gives [2].')

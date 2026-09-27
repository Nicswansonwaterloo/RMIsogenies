"""Degree-25 isogenies and RM checked against the explicit endomorphism on a product."""

from itertools import combinations, product

from sage.all import set_random_seed

from theta_rm.theta import torsion_label
from theta_rm.initialization import gen_rm_hash_prime, get_initial_vertex
from theta_rm.five_isogeny import ThetaFiveIsogeny
from theta_rm.isomorphisms import theta_isomorphisms
from theta_rm.product import ProductPoint


def torsion(A, n):
    B = A.torsion_basis(n)
    return [
        sum((a * P for a, P in zip(v, B)), A.zero) for v in product(range(n), repeat=4)
    ]


# Check the quotient and RM transport at small and cryptographic exponents.
for e in (3, 256):
    set_random_seed(31000 + e)
    p, _ = gen_rm_hash_prime(e)
    V = get_initial_vertex(p, e)
    W, rho, _ = next(
        edge for edge in V.get_neighbor_edges() if not edge[0].surface.is_product
    )
    A, B = V.surface, W.surface
    tau = lambda P: ProductPoint(-P[0] + 2 * P[1], 2 * P[0] + P[1])
    sigma = lambda P: ProductPoint(P[1], P[0] + P[1])

    # Pushing the basis reconstructs all 81 three-torsion combinations.
    basis = A.torsion_basis(3)
    for v in product(range(3), repeat=4):
        P = sum((a * Q for a, Q in zip(v, basis)), A.zero)
        assert W.three_torsion_point(v) == rho(P)

    # The kernel has 25 points; theta coordinates give 13 pairs including zero.
    five = torsion(A, 5)
    kernel = [P for P in five if tau(P) == A.zero]
    assert len(kernel) == 25
    assert all(P.weil_pairing(Q, 5) == 1 for P, Q in combinations(kernel, 2))
    phi = ThetaFiveIsogeny(B, W.five_kernel)
    assert set(phi.kernel) == {rho(P) for P in kernel}
    images = {rho(P): phi(rho(P)) for P in five}
    assert len(images) == 313
    assert sum(X == phi.codomain.O for X in images.values()) == 13
    assert len(set(images.values())) == 13

    # Check phi rho = f^-1 rho tau and translation by kernel points.
    samples = [A.zero, *A.torsion_basis(8), *A.torsion_basis(3), *A.torsion_basis(5)]
    samples += [A.random_point() for _ in range(4)]
    # Compare with the product matrix using one polarized isomorphism to the target.
    f = next(
        (
            f
            for f in theta_isomorphisms(phi.codomain, B)
            if all(f(phi(rho(P))) == rho(tau(P)) for P in samples)
        ),
        None,
    )
    assert f is not None
    assert all(f(X) == rho(tau(P)) for P in five for X in [images[rho(P)]])
    for P in samples:
        assert f(phi(f(phi(rho(P))))) == rho(5 * P)
        assert sigma(sigma(P)) - sigma(P) == P
        assert W.tau()(rho(P)) == rho(tau(P))
        assert all(phi(rho(P + R)) == phi(rho(P)) for R in kernel[:4])

    # All of B[2] lifts to A[4]. Transport sigma itself, not the scalar
    # tau|B[2], to predict the exact five Jacobian kernels independently.
    sigma_images = {}
    for P in torsion(A, 4):
        if rho(2 * P) == B.O:
            X, Y = rho(P), rho(sigma(P))
            assert X not in sigma_images or sigma_images[X] == Y
            sigma_images[X] = Y
    assert len(sigma_images) == 16
    two_torsion_points = B.two_torsion()
    expected = set()
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
        if all(sigma_images[X] in H for X in H):
            expected.add(H)
    assert len(expected) == 5
    assert {W.kernel_set(*K) for K in W.two_kernels()} == expected

    # Changing generators and affine representatives gives the same quotient,
    # up to a single polarized isomorphism.
    P, Q = V.five_kernel_basis()
    R, S = P + Q, Q
    K = [
        tuple(B.Fpp(c) * x for x in rho(T)) for c, T in zip((2, 3, 4), (R, S, R + S))
    ]
    psi = ThetaFiveIsogeny(B, K)
    assert set(psi.kernel) == set(phi.kernel)
    assert any(
        all(g(psi(rho(T))) == phi(rho(T)) for T in samples)
        for g in theta_isomorphisms(psi.codomain, phi.codomain)
    )

    print(f'e={e}: passed: degree-25 isogenies, their kernels, the induced RM and five RM-stable subgroups of B[2].')

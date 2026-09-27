"""Norm-two loop prediction and edge multiplicities: w_W m(V,W) = w_V m(W,V)."""

from itertools import combinations, product

from sage.all import EllipticCurve, set_random_seed

from theta_rm.initialization import gen_rm_hash_prime, get_initial_vertex
from theta_rm.product import ProductPoint, ProductSurface, product_isomorphisms
from theta_rm.rm_vertex import RMVertex
from theta_rm.isomorphisms import theta_isomorphisms


e = 3
for jcase in (1728, 0):
    set_random_seed(31000 + e)
    p, _ = gen_rm_hash_prime(e)
    initial = get_initial_vertex(p, e)
    Fpp = initial.surface.Fpp
    E = EllipticCurve(Fpp, [1, 0] if jcase == 1728 else [0, 1])
    A = ProductSurface(curves=(E, E))
    P, Q = E.torsion_basis(5)
    P, Q = ProductPoint(2 * P, P), ProductPoint(2 * Q, Q)
    B = A.torsion_basis(3)
    V = RMVertex(
        A,
        [A.theta_representative(T) for T in (P, Q, P + Q)],
        [A.theta_representative(T) for T in B]
        + [A.theta_representative(P + Q) for P, Q in combinations(B, 2)],
    )

    tau = lambda P: ProductPoint(-P[0] + 2 * P[1], 2 * P[0] + P[1])
    basis = A.torsion_basis(8)
    automorphisms = [
        f
        for f in product_isomorphisms(A, A)
        if all(f(tau(P)) == tau(f(P)) for P in basis)
    ]
    assert len(automorphisms) == (4 if jcase == 1728 else 6)
    w = len(automorphisms) // 2

    # alpha=[[a,b],[b,a+b]] commutes with sigma. alpha^dagger alpha=2 iff
    # deg(a)+deg(b)=2 and a^dagger b+b^dagger a=-deg(b).
    endos = [(0, lambda P: E(0), lambda P: E(0))]
    endos += [(1, f, ~f) for f in E.automorphisms()]
    for P in E(0).division_points(2):
        if P == E(0):
            continue
        q = E.isogeny(P)
        for f in q.codomain().isomorphisms(E):
            endos.append(
                (2, lambda P, q=q, f=f: f(q(P)), lambda P, q=q, f=f: q.dual()((~f)(P)))
            )
    norm_two = []
    for (da, a, ad), (db, b, bd) in product(endos, repeat=2):
        if da + db != 2:
            continue
        # The difference has elliptic degree <=9; E[8] detects zero.
        if all(ad(b(P)) + bd(a(P)) == -db * P for P in E.torsion_basis(8)):
            norm_two.append(
                lambda P, a=a, b=b: ProductPoint(
                    a(P[0]) + b(P[1]), b(P[0]) + a(P[1]) + b(P[1])
                )
            )
    predicted = len(norm_two) // len(automorphisms)
    assert len(norm_two) == predicted * len(automorphisms) and predicted == (
        1 if jcase == 1728 else 2
    )

    # Compare the predicted loop kernels with the outgoing RM-subgroups.
    B2 = A.torsion_basis(2)
    two = [
        sum((a * P for a, P in zip(v, B2)), A.zero) for v in product(range(2), repeat=4)
    ]
    loop_kernels = {
        frozenset(A.theta_representative(P) for P in two if f(P) == A.zero)
        for f in norm_two
    }
    edges = V.get_neighbor_edges()
    assert len(edges) == 5
    assert len({V.kernel_set(*K) for _, _, K in edges}) == 5
    actual = [K for W, phi, K in edges if W == V]
    assert len(actual) == predicted
    assert {V.kernel_set(*K) for K in actual} == loop_kernels

    # Duality weights multiplicities by the automorphism groups.
    unequal = False
    for W in V.get_neighbors():
        C = W.surface
        if C.is_product:
            B = W.three_torsion_basis()
            action = W.three_torsion_action()
            wy = (
                sum(
                    all(f(action[P]) == action[f(P)] for P in B)
                    for f in product_isomorphisms(C, C)
                )
                // 2
            )
        else:
            wy = sum(
                W.is_K_linear(W, f)
                for f in theta_isomorphisms(C, C)
            )
        forward = sum(X == W for X, _, _ in edges)
        outgoing = W.get_neighbor_edges()
        assert len(outgoing) == 5
        reverse = sum(X == V for X, _, _ in outgoing)
        assert wy * forward == w * reverse
        unequal |= wy != w
    assert unequal
    print(f'j={jcase}: {predicted} loops and w_W m(V,W) = w_V m(W,V).')

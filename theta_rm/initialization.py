"""The square with sigma(P, Q) = (Q, P+Q)."""

from itertools import combinations

from sage.all import EllipticCurve, GF, is_prime

from theta_rm.product import ProductSurface, ProductPoint
from theta_rm.rm_vertex import RMVertex


def gen_rm_hash_prime(e, d=5):
    assert e >= 3 and d == 5
    f = 1
    while not is_prime(2**e * 15 * f - 1):
        f += 2
    return 2**e * 15 * f - 1, f


def get_initial_vertex(p, e=4):
    assert e >= 3 and (p + 1) % (2**e * 15) == 0
    Fpp = GF(p * p, "i", modulus=[1, 0, 1])
    E = EllipticCurve(Fpp, [1, 0])
    A = ProductSurface(curves=(E, E))
    P, Q = E.torsion_basis(5)
    P, Q = ProductPoint(2 * P, P), ProductPoint(2 * Q, Q)
    basis = A.torsion_basis(3)
    five_kernel = [A.theta_representative(T) for T in (P, Q, P + Q)]
    points = [A.theta_representative(T) for T in basis]
    pairwise_sums = [A.theta_representative(P + Q) for P, Q in combinations(basis, 2)]

    return RMVertex(A, five_kernel, points + pairwise_sums)

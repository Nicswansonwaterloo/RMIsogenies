"""Check RM-stable kernels and isogeny evaluation by summing over the kernel."""

from itertools import combinations, product

from sage.all import set_random_seed, vector

from theta_rm.initialization import gen_rm_hash_prime, get_initial_vertex
from theta_rm.product import ProductPoint
from theta_rm.rm_vertex import RMVertex
from theta_rm.theta import ThetaSurface, normalise, normalizer_generators


e = 3
set_random_seed(61000 + e)
V = get_initial_vertex(gen_rm_hash_prime(e)[0], e)
W, rho, _ = next(
    edge for edge in V.get_neighbor_edges() if not edge[0].surface.is_product
)
samples = (
    *V.surface.torsion_basis(8),
    *V.surface.torsion_basis(3),
    V.surface.random_point(),
    V.surface.random_point(),
)
M = normalizer_generators(W.surface.Fpp)[0]
f = lambda P: normalise(M * vector(P))
U = RMVertex(
    ThetaSurface(f(W.surface.O)),
    [f(P) for P in W.five_kernel],
    [f(P) for P in W.three_torsion],
)
tau = U.tau()
phi = U._quotient()
lines = [
    {phi.source.mul(n, R) for n in (1, 2)} for R in phi.generators
]
assert len(lines) == 6 and all(len(line) == 2 for line in lines)
assert set.union(*lines) == set(phi.delta)
for P in samples:
    X = f(rho(P))
    expected = phi(X)
    assert phi(tuple(3 * x for x in X)) == expected
    # Independently evaluate the kernel sum in the source theta structure.
    R = normalise(phi.M * vector(X))
    total = vector([r**5 for r in R])
    for T, c in phi.delta.items():
        for S in phi.source.sums(R, T):
            total += vector([s**5 for s in S]) / (phi.cycle(R, T, S) * c**4)
    assert normalise(total) == expected
for P in samples:
    T = ProductPoint(-P[0] + 2 * P[1], 2 * P[0] + P[1])
    assert tau(f(rho(P))) == f(rho(T))
assert len(U.two_kernels()) == 5
# Independently transport sigma(P,Q)=(Q,P+Q) through rho and the theta change.
# Preimages in the product's four-torsion cover the target's two-torsion.
two = U.surface.two_torsion()
basis = V.surface.torsion_basis(4)
images = {}
for coefficients in product(range(4), repeat=4):
    P = sum((a * Q for a, Q in zip(coefficients, basis)), V.surface.zero)
    X = f(rho(P))
    if X in two.values():
        Y = f(rho(ProductPoint(P[1], P[0] + P[1])))
        assert X not in images or images[X] == Y
        images[X] = Y
assert len(images) == 16
expected = set()
for (a, b), (c, d) in combinations(list(two)[1:], 2):
    if a.dot_product(d) + b.dot_product(c):
        continue
    K = U.kernel_set(two[a, b], two[c, d])
    if all(images[P] in K for P in K):
        expected.add(K)
assert len(expected) == 5
assert {U.kernel_set(*K) for K in U.two_kernels()} == expected
print('Passed: RM-stable kernels agree with the explicit action of sigma; isogeny evaluation agrees with the sum over the kernel.')

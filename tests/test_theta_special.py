"""Zero coordinates and the sign ambiguity on products."""

from itertools import product

from sage.all import GF, vector, set_random_seed

from theta_rm.theta import (
    ThetaSurface,
    normalise,
    theta_changes,
)
from theta_rm.initialization import get_initial_vertex
from theta_rm.product import ProductPoint
from theta_rm.five_isogeny import ThetaFiveIsogeny
from theta_rm.isomorphisms import theta_isomorphisms


def torsion(A, n):
    B = A.torsion_basis(n)
    return [
        sum((a * P for a, P in zip(v, B)), A.zero) for v in product(range(n), repeat=4)
    ]


# Check differential addition when a coordinate of theta(P-Q) vanishes.
set_random_seed(31004)
V = get_initial_vertex(239, 4)
W, rho, _ = next(
    edge for edge in V.get_neighbor_edges() if not edge[0].surface.is_product
)
A, B = V.surface, W.surface
points = torsion(A, 4)
found = False
for j, M in enumerate(theta_changes(B.Fpp)):
    changed = ThetaSurface(M * vector(B.O))
    for R in points:
        D = normalise(M * vector(rho(R)))
        if all(D):
            continue
        Q = A.random_point()
        P = Q + R
        X, Y, S = [normalise(M * vector(rho(T))) for T in (P, Q, P + Q)]
        Z = changed.diff_add(X, Y, D)
        assert normalise(Z) == S
        # Affine scaling must be 2^2*3^2/5; compare tuples before normalising.
        assert changed.diff_add(
            tuple(2 * x for x in X), tuple(3 * y for y in Y), tuple(5 * d for d in D)
        ) == tuple(B.Fpp(36) / 5 * z for z in Z)
        assert set(changed.sums(X, Y)) == {S, D}
        found = True
        break
    if found:
        break
    assert j < 100, (
        'No vanishing coordinate of theta(P-Q) in the first 101 theta structures.'
    )
assert found

# The first gluing chart has z=0 and t=0; check the second chart and translation.
exceptional = next((P for P in torsion(A, 8) if rho.evaluate._image(P) is None), None)
assert exceptional is not None
X = rho(exceptional)
assert B.quartic()(*X) == 0
assert B.mul(3, X) == rho(3 * exceptional)
assert all(rho(exceptional + P) == X for P in torsion(A, 2) if rho(P) == B.O)

# A zero target tuple must trigger a change of theta structure and recomputation.
Fpp = GF(239**2, 'i', modulus=[1, 0, 1])
i = Fpp.gen()
C = ThetaSurface((Fpp(1), 238 * i, 224 * i + 224, 137 * i + 137))
K = [
    (Fpp(1), 157 * i + 41, 155 * i + 32, 118 * i + 43),
    (Fpp(1), 43 * i + 65, 39 * i + 228, 12 * i + 161),
    (Fpp(1), 147 * i + 97, 231 * i + 40, 120 * i + 84),
]
probe = ThetaFiveIsogeny.__new__(ThetaFiveIsogeny)
probe.source = C
assert not any(probe._compute(C, K))
phi = ThetaFiveIsogeny(C, K)
assert not phi.M.is_one()
assert len(phi.kernel) == 13
assert all(phi(P) == phi.codomain.O for P in phi.kernel)
samples = [C.O, *K, *C.four_torsion_basis()]
# A separately changed source gives the same quotient map up to one
# polarized target isomorphism, including evaluations outside the kernel.
for N in theta_changes(Fpp):
    D = ThetaSurface(N * vector(C.O))
    L = [normalise(N * vector(P)) for P in K]
    probe.source = D
    if any(probe._compute(D, L)) and N != phi.M:
        break
psi = ThetaFiveIsogeny(D, L)
assert any(
    all(f(phi(P)) == psi(normalise(N * vector(P))) for P in samples)
    for f in theta_isomorphisms(phi.codomain, psi.codomain)
)

# On products, matching tau modulo 3 does not imply (1+tau)/2 is integral.
A = V.surface
correct = lambda P: ProductPoint(-P[0] + 2 * P[1], 2 * P[0] + P[1])
wrong = lambda P: ProductPoint(2 * P[0] - P[1], -P[0] - 2 * P[1])
assert all(wrong(P) == correct(P) for P in A.torsion_basis(3))
assert all(wrong(wrong(P)) == 5 * P for P in A.torsion_basis(8))
assert any(wrong(P) != P for P in A.torsion_basis(2))
tau = V.tau()
assert all(tau(P) == correct(P) for P in A.torsion_basis(8))
assert all(tau(P) == P for P in A.torsion_basis(2))

# A product theta point loses independent signs, including points at infinity.
P, Q = A.torsion_basis(3)[:2]
samples = (A.zero, P, Q, P + Q)
for R in samples:
    X = A.theta_representative(R)
    assert R in A.points_from_theta(X)
    assert A.theta_representative(ProductPoint(-R[0], R[1])) == X
pairs = {
    (i, j): A.theta_representative(R + S)
    for i, R in enumerate(samples)
    for j, S in enumerate(samples)
    if i < j
}
assert tuple(samples) in A.sign_compatible_points(
    [A.theta_representative(R) for R in samples], pairs
)

# Here theta(P+Q)_i * theta(P-Q)_i = 0 for every i, so every B_ii vanishes.
B8 = A.torsion_basis(8)
# Choose disjoint supports for S and D, then halve their sum and difference.
points = []
for v in product(range(8), repeat=4):
    T = sum((a * P for a, P in zip(v, B8)), A.zero)
    image = rho(T)
    if sum(bool(c) for c in image) <= 2:
        points.append((v, image))
u, v = next(
    (u, v) for u, S in points for v, D in points
    if all((a - b) % 2 == 0 for a, b in zip(u, v))
    and all(s * d == 0 for s, d in zip(S, D))
)
P = sum(((a + b) // 2 * T for a, b, T in zip(u, v, B8)), A.zero)
Q = sum(((a - b) // 2 * T for a, b, T in zip(u, v, B8)), A.zero)
X, Y, S, D = (rho(T) for T in (P, Q, P + Q, P - Q))
assert all(s * d == 0 for s, d in zip(S, D))
C = B.biquadratic(X, Y)
assert all(C[i, i] == 0 for i in range(4)) and C.rank() == 2
assert set(B.sums(X, Y)) == {S, D}

# A nonzero double excludes the quadratic twist: its exponent p-1 is 2 mod 4.
p = B.Fpp.characteristic()
basis4 = B.four_torsion_basis()
assert all(B.mul(p + 1, Q) == B.O and B.mul(p - 1, Q) != B.O for Q in basis4)
print('Passed: vanishing theta coordinates, exceptional gluing points, signs on products and rational four-torsion.')

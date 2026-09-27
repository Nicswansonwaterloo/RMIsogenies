"""Independent product enumeration and change-of-polarization checks.

Run from the repository root: python -m tests.test_enumerate_rm
"""

from itertools import product
from math import isqrt

from sage.all import QQ, ZZ, QuaternionAlgebra, identity_matrix, matrix

from enumerate_rm import enumerate_rm, norm_vectors, star, trd


def frozen(X):
    X.set_immutable()
    return X


# Compare the recursive ellipsoid enumeration with a complete coordinate box.
G = matrix(QQ, [[3, QQ(1)/2], [QQ(1)/2, 2]])
for N in range(10):
    expected = {(x, y) for x, y in product(range(-4, 5), repeat=2)
                if 3*x*x + x*y + 2*y*y == N}
    assert {tuple(c) for c in norm_vectors(G, N)} == expected

for p in (3, 7, 19, 31):
    B = QuaternionAlgebra(QQ, -1, -p)
    i, j, k = B.gens()
    O = B.quaternion_order([B(1), i, (i+j)/2, (1+k)/2])
    assert O.is_maximal()
    I = identity_matrix(B, 2)
    for d in (2, 3, 5, 29):
        t, n = (1, (1-d)//4) if d % 4 == 1 else (0, -d)
        bound = QQ(t*t)/4 - n
        # q=(a+b*i+c*j+e*k)/2 lies in O iff a=e, b=c mod 2.
        # Its norm is (a^2+b^2+p*c^2+p*e^2)/4.
        R, S = isqrt(ZZ((4*bound).floor())), isqrt(ZZ((4*bound/p).floor()))
        quaternions = []
        for a, b, c, e in product(range(-R, R+1), range(-R, R+1),
                                 range(-S, S+1), range(-S, S+1)):
            if (a-e) % 2 == 0 and (b-c) % 2 == 0:
                q = (a+b*i+c*j+e*k)/2
                if q.reduced_norm() <= bound:
                    quaternions.append(q)
        expected = set()
        for m in range(-d, d+1):
            for q in quaternions:
                if q.reduced_norm() == t*m - m*m - n:
                    expected.add(frozen(matrix(B, [[m, q], [q.conjugate(), t-m]])))
        rms, units, classes = enumerate_rm(d, O)
        assert set(rms) == expected
        norm_one = [q for q in quaternions if q.reduced_norm() == 1]
        expected_units = {frozen(matrix(B, entries)) for a, b in product(norm_one, repeat=2)
                          for entries in ([[a, 0], [0, b]], [[0, a], [b, 0]])}
        assert set(units) == expected_units
        if p == 31 and d in (3, 5, 29):
            assert len(classes) == {3: 1, 5: 1, 29: 3}[d]
            representatives = {
                3: [[[1, -1-i], [-1+i, -1]]],
                5: [[[0, 1], [1, 1]]],
                29: [[[2, -1-2*i], [-1+2*i, -1]],
                     [[2, 1-2*i], [1+2*i, -1]], [[3, -i], [i, -2]]],
            }
            for orbit in classes:
                assert sum(matrix(B, entries) in orbit
                           for entries in representatives[d]) == 1
        if d == 29 and p in (7, 19):
            assert len(classes) == {7: 6, 19: 5}[p]
        print(f"p={p}, d={d}: {len(rms)} matrices, {len(units)} automorphisms, {len(classes)} classes")

# A noncentral change of basis checks the twist and the reduced trace.
h = matrix(B, [[1, j], [0, 1]])
hi = matrix(B, [[1, -j], [0, 1]])
g = star(h)*h
rms, units, classes = enumerate_rm(5, O)
twisted, twisted_units, twisted_classes = enumerate_rm(5, O, g)
assert set(twisted) == {frozen(hi*a*h) for a in rms}
assert set(twisted_units) == {frozen(hi*u*h) for u in units}
assert sorted(map(len, classes)) == sorted(map(len, twisted_classes))
a = hi*matrix(B, [[0, i], [-i, 1]])*h
assert a.trace() != B(1) and trd(a) == 2

# A further principal polarization, not supplied as a change of basis.
B = QuaternionAlgebra(QQ, -1, -19)
i, j, k = B.gens()
O = B.quaternion_order([B(1), i, (i+j)/2, (1+k)/2])
z = (i+j)/2
g = matrix(B, [[2, z], [z.conjugate(), 3]])
rms, units, classes = enumerate_rm(5, O, g)
assert units and sum(map(len, classes)) == len(rms)
for a in rms:
    assert a*a-a == identity_matrix(B, 2) and g*a == star(a)*g
print(f"p=19, g=[[2,z],[conjugate(z),3]], d=5: {len(rms)} matrices, {len(units)} automorphisms, {len(classes)} classes")

print("All enumeration checks passed.")

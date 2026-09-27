"""RM equality and conjugate embeddings on a depth-three theta BFS."""

from contextlib import redirect_stdout
from io import StringIO
from itertools import combinations

from sage.all import GF, Matrix, set_random_seed, vector

from bfs_traversal import bfs_rm
from theta_rm.initialization import get_initial_vertex
from theta_rm.isomorphisms import theta_isomorphisms
from theta_rm.rm_vertex import RMVertex
from theta_rm.theta import ThetaSurface, normalise, normalizer_generators


# Check equality, hashing, and five outgoing RM-subgroups along a depth-three BFS.
set_random_seed(421)
with redirect_stdout(StringIO()):
    vertices, edges = bfs_rm(get_initial_vertex(359, 3), 3)
assert all(len(neighbors) == 5 for neighbors in edges.values())

# The new basis anticommutes with T, so it represents the conjugate embedding.
generic = False
conjugate_isomorphic = False
B = ((0, 1, 0, 0), (2, 0, 0, 0), (0, 0, 0, 1), (0, 0, 2, 0))
for W in vertices:
    copy = RMVertex(W.surface, W.five_kernel, W.three_torsion)
    assert copy == W and W == copy and hash(copy) == hash(W)
    U = RMVertex(
        W.surface,
        W.five_kernel,
        [W.three_torsion_point(v) for v in B]
        + [
            W.three_torsion_point(tuple((a + b) % 3 for a, b in zip(v, w)))
            for v, w in combinations(B, 2)
        ],
    )
    equal = W == U
    assert (W.rm_invariants() == U.rm_invariants()) == equal
    assert (U == W) == equal
    assert len({W, U}) == (1 if equal else 2)
    assert {W.kernel_set(*K) for K in W.two_kernels()} == {
        U.kernel_set(*K) for K in U.two_kernels()
    }
    if not W.surface.is_product:
        if sum(1 for _ in theta_isomorphisms(W.surface, W.surface)) == 1:
            assert not equal
            assert W.three_torsion_action() == U.three_torsion_action()
            generic = True
        elif equal:
            conjugate_isomorphic = True
assert generic and conjugate_isomorphic

# A change of theta structure preserves the surface and its RM-data.
V = next(W for W in vertices if not W.surface.is_product)
for M in normalizer_generators(V.surface.Fpp):
    f = lambda P, M=M: normalise(M * vector(P))
    W = RMVertex(
        ThetaSurface(f(V.surface.O)),
        [f(P) for P in V.five_kernel],
        [f(P) for P in V.three_torsion],
    )
    assert W == V and hash(W) == hash(V)
    assert W.rm_invariants() == V.rm_invariants()


# Commuting changes of the three-torsion basis preserve the endomorphism.
A = V.surface
T = Matrix(GF(3), [[-1, 2, 0, 0], [2, 1, 0, 0], [0, 0, -1, 2], [0, 0, 2, 1]])
for C in (T, Matrix(GF(3), [[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]])):
    assert C * T == T * C and C.is_invertible()
    columns = C.columns()
    W = RMVertex(
        A,
        V.five_kernel,
        [V.three_torsion_point(tuple(v)) for v in columns]
        + [V.three_torsion_point(tuple(v + w)) for v, w in combinations(columns, 2)],
    )
    assert W == V and W.rm_invariants() == V.rm_invariants()
scales = range(2, 15)
points = V.five_kernel + V.three_torsion
scaled = [tuple(A.Fpp(c) * x for x in P) for c, P in zip(scales, points)]
assert RMVertex(ThetaSurface(tuple(7 * x for x in A.O)), scaled[:3], scaled[3:]) == V

print('Passed: RM equality and invariants under conjugation, basis changes and theta structure changes.')

"""BFS multiplicities and CGL walks with the dual subgroup excluded."""

from collections import Counter
from contextlib import redirect_stdout
from io import StringIO

from sage.all import set_random_seed

from bfs_traversal import bfs_rm
from theta_rm.initialization import get_initial_vertex


set_random_seed(51003)
V = get_initial_vertex(359, 3)

# Compute successive spheres from individual edges, before BFS groups targets.
distance = {V: 0}
sphere = {V}
for depth in range(4):
    sphere = {W for U in sphere for W, _, _ in U.get_neighbor_edges()} - set(distance)
    distance.update((W, depth + 1) for W in sphere)
for depth in (0, 1, 4):
    with redirect_stdout(StringIO()):
        visited, adjacency = bfs_rm(V, depth)
    assert visited == {W for W, d in distance.items() if d <= depth}
    assert set(adjacency) == {W for W, d in distance.items() if d < depth}
    for W, targets in adjacency.items():
        assert Counter(targets) == Counter(U for U, _, _ in W.get_neighbor_edges())
        assert sum(W.get_neighbors().values()) == 5

# A loop is an edge; a return along a parallel edge is not its dual.
W, phi, K = next(edge for edge in V.get_neighbor_edges() if edge[0] == V)
dual = V.dual_kernel(phi)
reverse, _, L = W.get_neighbor_edges(dual)[0]
assert reverse == V and W.kernel_set(*L) == dual
assert len(W.get_neighbor_edges(dual)[1:]) == 4

J = next(
    W for W, m in V.get_neighbors().items()
    if not W.surface.is_product and m > 1
)
U, phi, _ = next(edge for edge in J.get_neighbor_edges() if edge[0] == V)
dual = J.dual_kernel(phi)
edges = U.get_neighbor_edges(dual)
assert edges[0][0] == J and U.kernel_set(*edges[0][2]) == dual
W, psi, K = next(edge for edge in edges[1:] if edge[0] == J)
assert U.kernel_set(*K) != dual

# Return to J along the non-dual edge, then continue excluding each dual isogeny.
dual = U.dual_kernel(psi)
for step in range(6):
    choices = W.get_neighbor_edges(dual)[1:]
    assert len(choices) == 4
    U, phi, K = choices[step % 4]
    assert W.kernel_set(*K) != dual
    if not W.surface.is_product and not U.surface.is_product:
        assert all(
            U.tau()(phi(P)) == phi(W.tau()(P)) for P in W.four_torsion_basis()
        )
    dual = W.dual_kernel(phi)
    W = U

print('Passed: BFS distances and edge multiplicities; walks exclude dual isogenies and preserve RM.')

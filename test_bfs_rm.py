from collections import deque

from sage.all import Graph, Integers, Matrix, is_prime, kronecker

from helpers import get_symplectic_basis
from richelot_rm.richelot_product_isogenies import get_arbitrary_square_example
from richelot_rm.richelot_vertex_RM import RMVertex


def golden_ratio_action_on_symplectic_torsion(ell=2, e=1):
    r"""
    Return the action of the golden ratio on ell^e-torsion in a fixed basis.

    INPUT:

    - ``ell`` -- prime (default: ``2``)
    - ``e`` -- positive integer exponent (default: ``1``)

    OUTPUT: a 4x4 matrix over `\ZZ / ell^e\ZZ`
    """
    Zle = Integers(ell**e)
    return Matrix(Zle, [[0, 1, 0, 0], [1, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 1]])


def get_initial_vertex(p, e):
    r"""
    Return an initial RM vertex for the BFS at level $2^e$ and prime ``p``.

    INPUT:

    - ``p`` -- prime for the base field
    - ``e`` -- positive integer exponent

    OUTPUT: an ``RMVertex``
    """
    r = e + 2
    M = 2**r
    assert (p + 1) % M == 0, "p + 1 must be divisible by 2^(e + 2)."
    square = get_arbitrary_square_example(p)
    E1, E2 = square.E1, square.E2
    torsion_generators = get_symplectic_basis(E1, E2, M)
    rm_action = golden_ratio_action_on_symplectic_torsion(2, r)
    return RMVertex(square, r, torsion_generators, rm_action)


def gen_rm_hash_prime(e, d):
    """Generate a prime p = 4 * 2^e * f - 1 such that p ≡ 3 (mod 4) and d is a square mod p."""
    f = 1
    p = 4 * 2**e * f - 1
    while p % 4 != 3 or kronecker(d, p) != 1 or not is_prime(p):
        f += 1
        p = 4 * 2**e * f - 1
    return p, f


def bfs_rm(initial_vertex, max_depth):
    r"""
    BFS from ``initial_vertex`` in the RM (2,2)-isogeny graph.

    INPUT:

    - ``initial_vertex`` -- starting ``RMVertex``
    - ``max_depth`` -- BFS depth cap; each step consumes one factor of 2 from
      the RM level, so vertices at depth ``d`` have ``r = initial_r - d``

    OUTPUT: pair ``(visited, adjacency)`` where ``visited`` is the set of all
    reached vertices and ``adjacency`` maps each expanded vertex to its
    neighbor list (as returned by ``get_neighbors()``)
    """
    visited = {initial_vertex}
    queue = deque([(initial_vertex, 0)])
    adjacency = {}

    while queue:
        vertex, depth = queue.popleft()
        if depth >= max_depth:
            continue

        neighbors = vertex.get_neighbors()
        adjacency[vertex] = neighbors
        print(f"Depth {depth}: {vertex} -> {len(neighbors)} neighbor(s)")

        for neighbor in neighbors:
            if neighbor not in visited:
                visited.add(neighbor)
                queue.append((neighbor, depth + 1))

    return visited, adjacency


def test_bfs_rm():
    """Run BFS from one starting vertex in the RM graph and save a colored picture."""
    e = 4
    d = 5
    p, f = gen_rm_hash_prime(e, d)
    initial_vertex = get_initial_vertex(p, e)

    # Initial r = e + 2; stop when frontier vertices reach r = 3.
    max_depth = e

    visited, adjacency = bfs_rm(initial_vertex, max_depth)

    G = Graph(adjacency)

    labels = {v: v.get_type() for v in G.vertices()}
    non_start = [v for v in G.vertices() if v != initial_vertex]
    plot = G.plot(
        layout="spring",
        iterations=500,  # Increase iterations for a better spring layout settle
        vertex_size=200, # Adjust vertex size to fit labels better
        figsize=(10, 10), # Larger figure size distributes the nodes more widely
        vertex_labels=labels,
        vertex_colors={
            "#99ffa8": [initial_vertex],
            "#9dc3ff": non_start,
        },
    )
    filename = f"bfs_rm_d={d}_e={e}_f={f}.png"
    plot.save(filename)
    print(f"Saved {filename} with {len(visited)} vertices and {G.num_edges()} edges.")


test_bfs_rm()

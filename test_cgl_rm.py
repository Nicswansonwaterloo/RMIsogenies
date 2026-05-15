from sage.all import Graph, Integers, Matrix, is_prime, kronecker, randint

from helpers import get_symplectic_basis
from richelot_rm.richelot_product_isogenies import get_arbitrary_square_example
from richelot_rm.richelot_vertex_RM import RMVertex


def golden_ratio_action_on_symplectic_torsion(ell=2, e=1):
    r"""
    Return the action of the golden ratio on ell^e-torsion in a fixed basis.

    INPUT:

    - ``ell`` -- prime (default: ``2``)
    - ``e`` -- positive integer exponent (default: ``1``)

    OUTPUT: a 4x4 matrix over ``\ZZ / ell^e\ZZ``
    """
    Zle = Integers(ell**e)
    return Matrix(Zle, [[0, 1, 0, 0], [1, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 1]])


def get_initial_vertex(p, e):
    r"""
    Return an initial RM vertex for the CGL walk at level $2^e$ and prime ``p``.

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


def take_random_walk(initial_vertex, e, allow_backtrack=False):
    r"""
    Take a random walk in the RM graph and return the walk and graph.

    INPUT:

    - ``initial_vertex`` -- starting ``RMVertex``
    - ``e`` -- number of 2,2-isogeny steps to take
    - ``allow_backtrack`` -- boolean (default: ``False``)

    OUTPUT: pair ``(walk, G)`` where ``walk`` is a list of vertices and ``G`` is a graph
    """
    current_vertex = initial_vertex
    # We keep track of neighbors we do not visit for testing.
    graph_dict = {}
    walk = []
    for step in range(e):
        neighbors = current_vertex.get_neighbors()
        # neighbors are returned in a deterministic order with the backtracking neighbor first.
        random_choice = (
            randint(0, len(neighbors) - 1)
            if allow_backtrack
            else randint(1, len(neighbors) - 1)
        )

        next_vertex = neighbors[random_choice]
        graph_dict[current_vertex] = neighbors
        current_vertex = next_vertex
        walk.append(current_vertex)

    return walk, Graph(graph_dict)


def gen_rm_hash_prime(e, d):
    """Generate a prime of the form p = 4 * 2^e * f - 1 such that p is congruent to 3 mod 4 and d is a square mod p."""
    f = 1
    p = 4 * 2**e * f - 1
    while p % 4 != 3 or kronecker(d, p) != 1 or not is_prime(p):
        f += 1
        p = 4 * 2**e * f - 1
    return p, f


def test_non_backtracking_random_walk_cgl():
    """Draw a non-backtracking random walk picture for the fixed test case."""
    # form parameters for hash function.
    d = 5
    e = 32
    p, f = gen_rm_hash_prime(e, d)
    initial_vertex = get_initial_vertex(p, e)

    # take random walk.
    walk, G = take_random_walk(initial_vertex, e, allow_backtrack=False)

    # Nicely print the walk.
    labels = {v: v.get_type() for v in G.vertices()}
    non_walk = [v for v in G.vertices() if v not in walk and v != initial_vertex]
    plot = G.plot(
        vertex_labels=labels,
        vertex_colors={
            "#9dc3ff": walk,
            "#fadb87": non_walk,
            "#99ffa8": [initial_vertex],
        },
    )
    plot.save(f"random_non_backtracking_walk_d={d}_e={e}_f={f}.png")


test_non_backtracking_random_walk_cgl()
